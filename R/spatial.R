#' @title Spatial Analysis Functions
#'
#' @description Functions for spatial coordinate handling, neighbor graph
#'   construction, distance calculations, and neighborhood composition analysis.
#'
#' @name spatial
NULL

#' Auto-detect or verify coordinate columns
#'
#' Resolves spatial coordinate column names from metadata. If explicit column
#' names are provided, verifies they exist. Otherwise, auto-detects by trying
#' common spatial coordinate column pairs in order.
#'
#' @param meta A \code{data.table} or \code{data.frame} with metadata columns.
#' @param x_col Character or NULL. Explicit X coordinate column name to verify.
#' @param y_col Character or NULL. Explicit Y coordinate column name to verify.
#'
#' @return A character vector of length 2: \code{c(x_col, y_col)}.
#'
#' @details Auto-detection tries the following pairs in order:
#' \enumerate{
#'   \item \code{spatial_x}, \code{spatial_y}
#'   \item \code{x}, \code{y}
#'   \item \code{x_centroid}, \code{y_centroid}
#' }
#'
#' @examples
#' \dontrun{
#' meta <- data.table(spatial_x = runif(10), spatial_y = runif(10))
#' coords <- get_coord_columns(meta)
#' # Returns c("spatial_x", "spatial_y")
#'
#' # Or specify explicitly:
#' coords <- get_coord_columns(meta, x_col = "spatial_x", y_col = "spatial_y")
#' }
#'
#' @export
get_coord_columns <- function(meta, x_col = NULL, y_col = NULL) {
  col_names <- names(meta)


  # Priority 1: user-specified via arguments
  if (!is.null(x_col) && !is.null(y_col) &&
      nzchar(x_col) && nzchar(y_col)) {
    if (!x_col %in% col_names)
      stop(sprintf("X column '%s' not found in metadata. Available: %s",
                   x_col, paste(head(col_names, 20), collapse = ", ")))
    if (!y_col %in% col_names)
      stop(sprintf("Y column '%s' not found in metadata. Available: %s",
                   y_col, paste(head(col_names, 20), collapse = ", ")))
    return(c(x_col, y_col))
  }

  # Priority 2: auto-detect common pairs
  candidates <- list(
    c("spatial_x", "spatial_y"),
    c("x", "y"),
    c("x_centroid", "y_centroid")
  )
  for (pair in candidates) {
    if (all(pair %in% col_names)) {
      message(sprintf("  Auto-detected coordinate columns: %s, %s", pair[1], pair[2]))
      return(pair)
    }
  }

  stop("Could not find spatial coordinate columns in metadata.\n",
       "  Tried: spatial_x/spatial_y, x/y, x_centroid/y_centroid\n",
       "  Provide x_col and y_col arguments to specify custom column names.\n",
       "  Available columns: ", paste(head(col_names, 30), collapse = ", "))
}


#' Build k-nearest neighbor graph
#'
#' Constructs a k-nearest neighbor graph from spatial coordinates using
#' RANN for fast kNN queries.
#'
#' @param coords Numeric matrix of spatial coordinates (n x 2).
#' @param k Integer. Number of neighbors (default: 20).
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{indices}}{Integer matrix (n x k) of neighbor indices.}
#'   \item{\code{distances}}{Numeric matrix (n x k) of distances to neighbors.}
#' }
#'
#' @details Uses \code{RANN::nn2} for fast kNN queries. The self-neighbor
#'   (distance 0) is automatically excluded.
#'
#' @examples
#' \dontrun{
#' coords <- matrix(runif(200), ncol = 2)
#' knn <- build_knn_graph(coords, k = 10)
#' }
#'
#' @importFrom RANN nn2
#' @export
build_knn_graph <- function(coords, k = 20) {
  message(sprintf("Building %d-nearest neighbor graph for %d cells...", k, nrow(coords)))

  # RANN::nn2 is very fast for kNN queries
  nn_result <- RANN::nn2(coords, coords, k = k + 1)  # +1 because cell is its own neighbor

  # Remove self (first column)
  list(
    indices = nn_result$nn.idx[, -1, drop = FALSE],
    distances = nn_result$nn.dists[, -1, drop = FALSE]
  )
}


#' Build radius-based neighbor graph
#'
#' Constructs a neighbor graph where cells are connected if they are within
#' a specified radius.
#'
#' @param coords Numeric matrix of spatial coordinates (n x 2).
#' @param radius Numeric. Search radius in coordinate units (typically microns).
#'
#' @return A list of length n, where each element is an integer vector of
#'   neighbor indices within the radius.
#'
#' @details Uses \code{RANN::nn2} for an initial broad kNN search, then filters
#'   to only neighbors within the specified radius. The estimated k for the
#'   initial search is based on average cell density.
#'
#' @examples
#' \dontrun{
#' coords <- matrix(runif(200), ncol = 2)
#' neighbors <- build_radius_graph(coords, radius = 50)
#' }
#'
#' @importFrom RANN nn2
#' @export
build_radius_graph <- function(coords, radius) {
  message(sprintf("Building radius neighbor graph (r=%s)...", radius))

  n <- nrow(coords)
  neighbors <- vector("list", n)

  # Use RANN for initial broad search, then filter by radius
  # Estimate k based on typical density
  avg_density <- n / (diff(range(coords[, 1])) * diff(range(coords[, 2])))
  k_est <- min(ceiling(avg_density * pi * radius^2 * 2), n - 1)

  nn_result <- RANN::nn2(coords, coords, k = k_est + 1)

  for (i in seq_len(n)) {
    within_radius <- nn_result$nn.dists[i, ] <= radius & nn_result$nn.dists[i, ] > 0
    neighbors[[i]] <- nn_result$nn.idx[i, within_radius]
  }

  return(neighbors)
}


#' Calculate distance to nearest cell of a given type
#'
#' For each cell, computes the Euclidean distance to the nearest cell of a
#' specified target type using kNN search.
#'
#' @param coords Numeric matrix of spatial coordinates for all cells (n x 2).
#' @param cell_types Character vector of cell type labels (length n).
#' @param target_type Character. The cell type to measure distance to.
#'
#' @return A numeric vector of length n with distances. Returns \code{NA} for
#'   all cells if no target cells are found.
#'
#' @examples
#' \dontrun{
#' coords <- matrix(runif(200), ncol = 2)
#' types <- sample(c("A", "B", "C"), 100, replace = TRUE)
#' dists <- calculate_distance_to_type(coords, types, "A")
#' }
#'
#' @importFrom RANN nn2
#' @export
calculate_distance_to_type <- function(coords, cell_types, target_type) {
  target_mask <- cell_types == target_type
  target_coords <- coords[target_mask, , drop = FALSE]

  if (nrow(target_coords) == 0) {
    warning(sprintf("No cells of type '%s' found", target_type))
    return(rep(NA_real_, nrow(coords)))
  }

  # Find nearest target cell for each cell
  nn_result <- RANN::nn2(target_coords, coords, k = 1)

  return(nn_result$nn.dists[, 1])
}


#' Get neighbor cell types for a set of query cells
#'
#' Retrieves the cell type labels of all neighbors for specified query cells,
#' based on a precomputed kNN graph.
#'
#' @param cell_types Character vector of cell type labels for all cells.
#' @param query_indices Integer vector of indices for query cells.
#' @param knn_result List. Result from \code{\link{build_knn_graph}}, containing
#'   \code{indices} and \code{distances} matrices.
#'
#' @return A \code{data.table} with columns \code{query_cell}, \code{neighbor_cell},
#'   and \code{neighbor_type}.
#'
#' @importFrom data.table rbindlist data.table
#' @export
get_neighbor_cell_types <- function(cell_types, query_indices, knn_result) {
  results <- data.table::rbindlist(lapply(query_indices, function(i) {
    neighbor_idx <- knn_result$indices[i, ]
    data.table::data.table(
      query_cell = i,
      neighbor_cell = neighbor_idx,
      neighbor_type = cell_types[neighbor_idx]
    )
  }))
  return(results)
}


#' Calculate neighborhood composition for a cell type
#'
#' Computes the proportions of different cell types among the neighbors of
#' all cells of a specified query type.
#'
#' @param cell_types Character vector of cell type labels for all cells.
#' @param query_cell_type Character. Cell type whose neighborhoods to analyze.
#' @param knn_result List. Result from \code{\link{build_knn_graph}}.
#'
#' @return A \code{data.table} with columns \code{query_cell_type},
#'   \code{neighbor_type}, \code{count}, and \code{proportion}.
#'
#' @examples
#' \dontrun{
#' coords <- matrix(runif(200), ncol = 2)
#' types <- sample(c("A", "B", "C"), 100, replace = TRUE)
#' knn <- build_knn_graph(coords, k = 10)
#' comp <- calculate_neighbor_composition(types, "A", knn)
#' }
#'
#' @importFrom data.table data.table setcolorder
#' @export
calculate_neighbor_composition <- function(cell_types, query_cell_type, knn_result) {
  query_idx <- which(cell_types == query_cell_type)

  if (length(query_idx) == 0) {
    warning(sprintf("No cells of type '%s' found", query_cell_type))
    return(data.table::data.table())
  }

  # Get all neighbor types
  neighbor_types <- as.vector(knn_result$indices[query_idx, ])
  neighbor_types <- cell_types[neighbor_types]

  # Count and calculate proportions
  counts <- table(neighbor_types)
  dt <- data.table::data.table(
    neighbor_type = names(counts),
    count = as.integer(counts),
    proportion = as.numeric(counts) / sum(counts)
  )

  dt[, query_cell_type := query_cell_type]
  data.table::setcolorder(dt, c("query_cell_type", "neighbor_type", "count", "proportion"))

  return(dt)
}
