# Plan: Honor init_k for Subclustering with init_clusters

## Problem Statement

Currently when `init_clusters` is provided, subclustering is determined by cell count:
```r
k_sub <- max(1, ceiling(n_cells / cells_per_subcluster))
```

This creates too many subclusters for large datasets. For example, with 10,000 cells and `cells_per_subcluster=100`, you get ~100 subclusters regardless of whether that subdivision is meaningful.

## Goal

Honor the `init_k` parameter (increase default to 20) to ensure exactly `init_k` total subclusters across all init_clusters, selecting the most biologically meaningful subdivisions.

## Proposed Approach: Global Iterative Splitting

### Core Idea

Instead of allocating subclusters per init_cluster and then cutting, we use **iterative splitting** that competes globally:

1. **Build dendrograms**: Run hclust for each init_cluster
2. **Start unified**: Each init_cluster begins as a single subcluster
3. **Iteratively split**: At each step, find the **globally most distinct** split across all current subclusters (where both children have >= `cells_per_subcluster` cells)
4. **Stop at init_k**: Continue until we have exactly `init_k` subclusters

This ensures:
- Large init_clusters don't automatically get more subclusters
- Only the most distinct splits (globally) are selected
- All subclusters guaranteed >= `cells_per_subcluster` cells

### Distinctness Metric (Size-Independent)

Use the **merge height** from Ward's D2 hierarchical clustering:
- The height at which two clusters merge represents their dissimilarity
- Ward's method inherently normalizes for cluster sizes
- Higher merge height = more distinct split

This is size-independent because Ward's D2 measures the increase in total within-cluster variance when merging, normalized by the harmonic mean of cluster sizes.

### Algorithm

```
Input:
  - init_clusters: cell -> init_cluster_id mapping
  - init_k: target total subclusters (hard requirement)
  - cells_per_subcluster: minimum cells per subcluster

Data structures:
  - hclust_list: hclust object for each init_cluster
  - current_nodes: list of (node_id, init_cluster_id, cells)
    Each entry represents a current subcluster as a node in its dendrogram

Initialization:
  1. For each init_cluster:
     a. Run exp_hclust() to get dendrogram
     b. Add root node (all cells) to current_nodes
  2. Total subclusters = n_init_clusters

Iteration (while total < init_k):
  1. For each node in current_nodes:
     a. Check if it can be split:
        - Get its two children in the dendrogram
        - Both children must have >= cells_per_subcluster cells
     b. If splittable, record the split's merge height (distinctness)

  2. If no valid splits exist:
     Error: "Cannot create init_k subclusters while maintaining
            cells_per_subcluster minimum"

  3. Select the split with HIGHEST merge height (most distinct)

  4. Apply the split:
     - Remove parent node from current_nodes
     - Add both children to current_nodes
     - Total subclusters += 1

Output: current_nodes (exactly init_k subclusters)
```

### Why This Works

1. **Guaranteed minimum cells**: We only allow splits where BOTH children have >= `cells_per_subcluster`. This is checked before each split, not inferred from k.

2. **Global competition**: A split in init_cluster A competes with splits in init_cluster B. The most distinct split wins, regardless of which init_cluster it's in.

3. **Size-independent**: Ward's merge height is the distinctness metric. A small init_cluster with a very distinct internal division will be split before a large init_cluster with a less distinct division.

4. **No over-splitting large clusters**: Large init_clusters only get more subclusters if they genuinely have more distinct internal structure.

### Helper Functions Needed

```r
#' Get children of a dendrogram node with their cell counts
#'
#' @param hc hclust object
#' @param node_id Node identifier (negative = leaf, positive = internal)
#' @return List with left/right children info, or NULL if leaf
get_node_children <- function(hc, node_id) { ... }

#' Check if a node can be validly split
#'
#' @param hc hclust object
#' @param node_id Current node
#' @param min_size Minimum cells per child
#' @return List(can_split, height, left_child, right_child) or NULL
check_split <- function(hc, node_id, min_size) { ... }

#' Get all cells belonging to a dendrogram node
#'
#' @param hc hclust object
#' @param node_id Node identifier
#' @return Vector of cell names
get_node_cells <- function(hc, node_id) { ... }
```

## Implementation Plan

### Step 1: Modify parameter defaults (main.R)

**Line 73**: Change `init_k = 3` to `init_k = 20`

**Line 67**: Update documentation for `cells_per_subcluster` to note it's only used when `init_clusters = NULL`

Add documentation for how `init_k` is used with `init_clusters`.

### Step 2: Create helper functions for dendrogram traversal

New functions in `main.R` for working with hclust dendrograms:

```r
#' Get the number of cells (leaves) under a dendrogram node
#'
#' In hclust, merge matrix uses negative indices for leaves and positive
#' indices for internal nodes (previous merges).
#'
#' @param hc hclust object
#' @param node_id Node identifier: negative = leaf (single cell),
#'                positive = internal node (merge index)
#' @return Number of cells under this node
get_node_size <- function(hc, node_id) {
    if (node_id < 0) {
        # Leaf node = single cell
        return(1)
    } else {
        # Internal node: sum of both children
        left <- hc$merge[node_id, 1]
        right <- hc$merge[node_id, 2]
        return(get_node_size(hc, left) + get_node_size(hc, right))
    }
}

#' Get all cell indices under a dendrogram node
#'
#' @param hc hclust object
#' @param node_id Node identifier
#' @return Vector of cell indices (positions in hc$labels)
get_node_members <- function(hc, node_id) {
    if (node_id < 0) {
        # Leaf: the cell index is -node_id
        return(-node_id)
    } else {
        # Internal node: combine children
        left <- hc$merge[node_id, 1]
        right <- hc$merge[node_id, 2]
        return(c(get_node_members(hc, left), get_node_members(hc, right)))
    }
}

#' Check if a node can be split while respecting minimum size
#'
#' @param hc hclust object
#' @param node_id Node to potentially split
#' @param min_size Minimum cells per resulting subcluster
#' @return List(can_split, height, left_id, right_id, left_size, right_size)
#'         or NULL if node is a leaf
check_valid_split <- function(hc, node_id, min_size) {
    if (node_id < 0) {
        # Leaf nodes cannot be split
        return(list(can_split = FALSE))
    }

    left_id <- hc$merge[node_id, 1]
    right_id <- hc$merge[node_id, 2]
    left_size <- get_node_size(hc, left_id)
    right_size <- get_node_size(hc, right_id)
    height <- hc$height[node_id]

    can_split <- (left_size >= min_size) && (right_size >= min_size)

    return(list(
        can_split = can_split,
        height = height,
        left_id = left_id,
        right_id = right_id,
        left_size = left_size,
        right_size = right_size
    ))
}

#' Iteratively split subclusters to reach target k
#'
#' Globally selects the most distinct (highest merge height) valid split
#' at each step until we have init_k subclusters.
#'
#' @param hclust_list Named list of hclust objects per init_cluster
#' @param init_k Target total number of subclusters
#' @param cells_per_subcluster Minimum cells per subcluster
#' @param verbose Print progress
#' @return List of subclusters, each with: init_cluster, node_id, cells
select_subclusters_iterative <- function(hclust_list, init_k,
                                          cells_per_subcluster, verbose = TRUE) {

    n_init_clusters <- length(hclust_list)

    # Validate
    if (init_k < n_init_clusters) {
        stop(glue("init_k ({init_k}) must be >= number of init_clusters ({n_init_clusters})"))
    }

    # Initialize: each init_cluster starts as one subcluster (root node)
    # For hclust with n leaves, root is node n-1 (last merge)
    current_subclusters <- list()
    for (cl in names(hclust_list)) {
        hc <- hclust_list[[cl]]
        if (is.null(hc)) {
            # No hclust = too few cells, treat as single subcluster
            current_subclusters[[length(current_subclusters) + 1]] <- list(
                init_cluster = cl,
                node_id = NA,  # Special: entire cluster, no dendrogram
                hc = NULL
            )
        } else {
            n_leaves <- length(hc$labels)
            root_id <- n_leaves - 1  # Root is the last merge
            current_subclusters[[length(current_subclusters) + 1]] <- list(
                init_cluster = cl,
                node_id = root_id,
                hc = hc
            )
        }
    }

    # Iteratively split until we reach init_k
    while (length(current_subclusters) < init_k) {

        # Find the best valid split across all current subclusters
        best_split <- NULL
        best_height <- -Inf
        best_idx <- NA

        for (i in seq_along(current_subclusters)) {
            sc <- current_subclusters[[i]]

            if (is.null(sc$hc) || is.na(sc$node_id)) {
                # Can't split: no dendrogram or already a leaf
                next
            }

            split_info <- check_valid_split(sc$hc, sc$node_id, cells_per_subcluster)

            if (split_info$can_split && split_info$height > best_height) {
                best_split <- split_info
                best_height <- split_info$height
                best_idx <- i
                best_hc <- sc$hc
                best_init_cluster <- sc$init_cluster
            }
        }

        if (is.null(best_split)) {
            stop(glue(
                "Cannot create {init_k} subclusters: only {length(current_subclusters)} ",
                "possible while maintaining cells_per_subcluster={cells_per_subcluster}. ",
                "Reduce init_k or reduce cells_per_subcluster."
            ))
        }

        # Apply the split: remove parent, add both children
        current_subclusters[[best_idx]] <- NULL  # Remove parent

        # Add left child
        current_subclusters[[length(current_subclusters) + 1]] <- list(
            init_cluster = best_init_cluster,
            node_id = best_split$left_id,
            hc = best_hc
        )

        # Add right child
        current_subclusters[[length(current_subclusters) + 1]] <- list(
            init_cluster = best_init_cluster,
            node_id = best_split$right_id,
            hc = best_hc
        )

        if (verbose) {
            log_message(glue(
                "  Split in {best_init_cluster}: height={round(best_height, 3)}, ",
                "sizes={best_split$left_size}+{best_split$right_size}, ",
                "total subclusters={length(current_subclusters)}"
            ), verbose = verbose)
        }
    }

    # Convert node_ids to actual cell names
    result <- lapply(current_subclusters, function(sc) {
        if (is.null(sc$hc)) {
            # No dendrogram - cells should be retrieved from init_clusters
            return(list(
                init_cluster = sc$init_cluster,
                cells = NULL  # Will be filled in by caller
            ))
        }

        if (sc$node_id < 0) {
            # Leaf node
            cell_idx <- -sc$node_id
            cells <- sc$hc$labels[cell_idx]
        } else {
            # Internal node
            cell_indices <- get_node_members(sc$hc, sc$node_id)
            cells <- sc$hc$labels[cell_indices]
        }

        return(list(
            init_cluster = sc$init_cluster,
            cells = cells
        ))
    })

    return(result)
}
```

### Step 3: Refactor subclustering logic (main.R, lines 238-335)

Replace the current cell-count-based approach with global iterative splitting:

```r
} else if (!is.null(init_clusters)) {

    log_message('Initializing from user-defined clusters with global iterative splitting..', verbose = verbose)
    log_message(glue('Target: {init_k} total subclusters, min {cells_per_subcluster} cells each'), verbose = verbose)

    # Validate init_clusters (keep existing validation)
    if (!is.vector(init_clusters) || is.null(names(init_clusters))) {
        stop('init_clusters must be a named vector with cell IDs as names')
    }

    missing_cells <- setdiff(names(init_clusters), colnames(count_mat))
    if (length(missing_cells) > 0) {
        stop(glue('init_clusters contains {length(missing_cells)} cell IDs not in count_mat'))
    }

    cluster_ids <- unique(init_clusters)
    n_init_clusters <- length(cluster_ids)
    log_message(glue('Processing {n_init_clusters} external clusters..'), verbose = verbose)

    # Phase 1: Run hclust for each init_cluster
    hclust_list <- list()
    cell_ids_list <- list()
    gexp_roll_list <- list()

    for (cl in cluster_ids) {
        cell_ids <- names(init_clusters)[init_clusters == cl]
        n_cells <- length(cell_ids)
        cell_ids_list[[as.character(cl)]] <- cell_ids

        log_message(glue('  Cluster {cl}: {n_cells} cells'), verbose = verbose)

        if (n_cells < cells_per_subcluster) {
            # Too small to subcluster - will be kept as single subtree
            hclust_list[[as.character(cl)]] <- NULL
            log_message(glue('    -> Below minimum size, keeping as single subtree'), verbose = verbose)
            next
        }

        count_mat_sub <- count_mat[, cell_ids, drop = FALSE]
        sc_refs_sub <- sc_refs[cell_ids]

        clust_sub <- exp_hclust(
            count_mat = count_mat_sub,
            lambdas_ref = lambdas_ref,
            gtf = gtf,
            sc_refs = sc_refs_sub,
            ncores = ncores,
            verbose = FALSE
        )

        hclust_list[[as.character(cl)]] <- clust_sub$hc
        gexp_roll_list[[as.character(cl)]] <- clust_sub$gexp_roll_wide
    }

    # Phase 2: Iteratively select most distinct splits globally
    # This competes across all init_clusters - large clusters only get more
    # subclusters if they have the most distinct splits globally
    log_message('Selecting subclusters by global distinctness (merge height)..', verbose = verbose)

    selected_subclusters <- select_subclusters_iterative(
        hclust_list = hclust_list,
        init_k = init_k,
        cells_per_subcluster = cells_per_subcluster,
        verbose = verbose
    )

    # Phase 3: Convert selected subclusters to subtrees
    all_subtrees <- list()
    subtree_id <- 1

    for (sc in selected_subclusters) {
        # Handle init_clusters that had no hclust (too small)
        if (is.null(sc$cells)) {
            sc$cells <- cell_ids_list[[sc$init_cluster]]
        }

        all_subtrees[[subtree_id]] <- list(
            cells = sc$cells,
            size = length(sc$cells),
            sample = as.character(subtree_id),
            members = subtree_id
        )
        subtree_id <- subtree_id + 1
    }

    subtrees <- all_subtrees

    # Log final allocation per init_cluster
    allocation_summary <- table(sapply(selected_subclusters, `[[`, "init_cluster"))
    log_message(glue('Final allocation: {paste(names(allocation_summary), allocation_summary, sep=":", collapse=", ")}'), verbose = verbose)

    # Combine and save smoothed expression from all clusters (keep existing code)
    if (length(gexp_roll_list) > 0) {
        common_genes <- Reduce(intersect, lapply(gexp_roll_list, colnames))
        if (length(common_genes) > 0) {
            all_gexp_roll_common <- lapply(gexp_roll_list, function(x) x[, common_genes, drop = FALSE])
            gexp_roll_wide <- do.call(rbind, all_gexp_roll_common)
            fwrite(
                as.data.frame(gexp_roll_wide) %>% tibble::rownames_to_column('cell'),
                glue('{out_dir}/gexp_roll_wide.tsv.gz'),
                sep = '\t',
                nThread = min(4, ncores)
            )
        }
    }

    log_message(glue('Created {length(subtrees)} subtrees from {n_init_clusters} external clusters'), verbose = verbose)
```

### Step 4: Update documentation

Update roxygen comments for:
- `init_k`: "Number of total subclusters when using init_clusters, or clusters when using default expression-based initialization"
- `cells_per_subcluster`: Add note "(only used when init_clusters = NULL)"

### Step 5: Add logging

Log the heterogeneity scores and allocation decisions for transparency.

## Files to Modify

1. **`R/main.R`**:
   - Line 73: Change default `init_k = 3` to `init_k = 20`
   - Line 34: Update documentation for `init_k` to explain it controls total subclusters
   - Line 67: Update documentation for `cells_per_subcluster` to clarify it's the minimum cells per subcluster
   - Lines 238-335: Replace cell-count-based subclustering with global iterative splitting
   - Add new helper functions (~lines 800-900):
     - `get_node_size()` - count cells under a dendrogram node
     - `get_node_members()` - get cell indices under a node
     - `check_valid_split()` - check if a split respects minimum size
     - `select_subclusters_iterative()` - main algorithm for global selection

## Testing Considerations

1. **Error cases**: Test that appropriate errors are thrown when:
   - `init_k < n_init_clusters` (too few subclusters)
   - `init_k > max_possible` given `cells_per_subcluster` constraint

2. **Minimum size guarantee**: Verify ALL resulting subclusters have >= `cells_per_subcluster` cells
   - This is the critical guarantee - check with various dendrogram structures

3. **Global competition**: Verify that a small init_cluster with very distinct internal structure
   gets split before a large init_cluster with less distinct structure
   - Create synthetic test case where small cluster has higher merge heights

4. **Size independence**: Verify subclusters of different sizes are equally likely to be selected
   based on merge height, not cluster size

5. **Edge cases**:
   - `init_k = n_init_clusters` (no splits needed, each cluster stays whole)
   - `init_k = 1` with multiple init_clusters (should error)
   - Small clusters with `n_cells < cells_per_subcluster` (kept as single subtree, counts toward init_k)
   - Unbalanced dendrograms where one side has < `cells_per_subcluster`

6. **Dendrogram traversal correctness**:
   - Verify `get_node_size()` returns correct counts
   - Verify `get_node_members()` returns correct cell indices
   - Verify cell names are correctly mapped from indices

## Design Decisions (Confirmed)

1. **`init_k` is a hard requirement**: Error if impossible to achieve
2. **No cell-count-based option**: Remove the old scaling behavior entirely
3. **Minimum cells per subcluster**: Each subcluster must have at least `cells_per_subcluster` cells - guaranteed by checking BEFORE each split
4. **Distinctness metric**: Ward's D2 merge height - size-independent measure of how distinct a split is
5. **Global competition**: Splits compete across ALL init_clusters, not allocated per-cluster first

This approach guarantees:
- Exactly `init_k` total subclusters (or error if impossible)
- Every subcluster has >= `cells_per_subcluster` cells (checked at split time)
- Most distinct splits are selected globally (large clusters don't automatically win)
- Size-independent selection (small distinct clusters can "beat" large less-distinct ones)
