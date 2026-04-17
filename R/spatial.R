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
    if (!x_col %in% col_names) {
      stop(sprintf(
        "X column '%s' not found in metadata. Available: %s",
        x_col, paste(head(col_names, 20), collapse = ", ")
      ))
    }
    if (!y_col %in% col_names) {
      stop(sprintf(
        "Y column '%s' not found in metadata. Available: %s",
        y_col, paste(head(col_names, 20), collapse = ", ")
      ))
    }
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

  stop(
    "Could not find spatial coordinate columns in metadata.\n",
    "  Tried: spatial_x/spatial_y, x/y, x_centroid/y_centroid\n",
    "  Provide x_col and y_col arguments to specify custom column names.\n",
    "  Available columns: ", paste(head(col_names, 30), collapse = ", ")
  )
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
  nn_result <- RANN::nn2(coords, coords, k = k + 1) # +1 because cell is its own neighbor

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


#' Check spatial autocorrelation in RIPPLE model residuals
#'
#' Computes Moran's I on Poisson GLM residuals for one or more genes of
#' interest within a specific target cell type. High Moran's I indicates
#' that nearby cells have correlated residuals, which means the GLM's
#' independence assumption is violated and per-sample p-values may be
#' anti-conservative.
#'
#' This is a diagnostic tool — it does not fix the autocorrelation, but
#' lets you assess its severity for specific genes. RIPPLE's per-sample
#' fitting + Fisher's method already mitigates pseudoreplication at the
#' sample level (N = biological replicates, not cells).
#'
#' @param input A Seurat, SCE, or SpatialExperiment object (or path to
#'   an \code{.rds} file containing one).
#' @param genes Character vector of gene names to check.
#' @param celltype_column Cell type column name.
#' @param target_celltype Which cell type to assess (the target, not the query).
#' @param query_celltype Query cell type (needed for distance calculation).
#' @param sample_column Sample/replicate column name (default: \code{"sample_id"}).
#' @param k Number of nearest neighbors for the spatial weights matrix
#'   (default: 20).
#' @param max_distance_um Maximum distance to consider (default: 200).
#' @param x_column X coordinate column (default: NULL, auto-detect).
#' @param y_column Y coordinate column (default: NULL, auto-detect).
#' @param verbose Print progress (default: TRUE).
#'
#' @return A \code{data.table} with columns:
#' \describe{
#'   \item{gene}{Gene name.}
#'   \item{sample_id}{Sample identifier.}
#'   \item{morans_i}{Observed Moran's I statistic.}
#'   \item{morans_expected}{Expected Moran's I under no autocorrelation.}
#'   \item{morans_pvalue}{P-value from \code{spdep::moran.test()}.}
#'   \item{interpretation}{One of "none", "weak", "moderate", "strong".}
#'   \item{n_cells}{Number of cells used.}
#' }
#'
#' @details
#' For each gene and sample, the function:
#' \enumerate{
#'   \item Subsets to target cells within \code{max_distance_um} of the
#'     nearest query cell.
#'   \item Fits the same Poisson GLM as \code{\link{fit_poisson}}.
#'   \item Extracts deviance residuals.
#'   \item Builds a k-nearest-neighbor spatial weights matrix.
#'   \item Computes Moran's I via \code{spdep::moran.test()}.
#' }
#'
#' Interpretation thresholds:
#' \itemize{
#'   \item \code{|I| < 0.05}: "none" — independence assumption holds.
#'   \item \code{0.05 <= |I| < 0.15}: "weak" — minor concern.
#'   \item \code{0.15 <= |I| < 0.30}: "moderate" — p-values may be inflated.
#'   \item \code{|I| >= 0.30}: "strong" — interpret with caution.
#' }
#'
#' @examples
#' \dontrun{
#' autocor <- check_spatial_autocorrelation(
#'   input           = my_spe,
#'   genes           = c("MIF", "CD74", "PDCD1"),
#'   celltype_column = "cell_type",
#'   target_celltype = "CD8_T",
#'   query_celltype  = "tumor",
#'   sample_column   = "patient"
#' )
#' autocor[order(-morans_i)]
#' }
#'
#' @importFrom data.table data.table rbindlist
#' @export
check_spatial_autocorrelation <- function(input,
                                          genes,
                                          celltype_column,
                                          target_celltype,
                                          query_celltype,
                                          sample_column = "sample_id",
                                          k = 20,
                                          max_distance_um = 200,
                                          x_column = NULL,
                                          y_column = NULL,
                                          verbose = TRUE) {
  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("Package 'spdep' is required for spatial autocorrelation testing.\n",
      "Install with: install.packages('spdep')",
      call. = FALSE
    )
  }

  .msg <- function(...) if (isTRUE(verbose)) message(...)

  # Load data
  data <- .resolve_input(input, require_expr = FALSE, verbose = verbose)
  count_matrix <- data$counts
  cell_data <- data$meta
  rm(data)

  # Resolve coordinates

  coord_cols <- get_coord_columns(cell_data, x_col = x_column, y_col = y_column)
  coords <- as.matrix(cell_data[, ..coord_cols])

  # Validate columns
  if (!celltype_column %in% names(cell_data)) {
    stop("celltype_column '", celltype_column, "' not found in metadata.",
      call. = FALSE
    )
  }
  if (!sample_column %in% names(cell_data)) {
    stop("sample_column '", sample_column, "' not found in metadata.",
      call. = FALSE
    )
  }

  # Calculate distances to query
  query_mask <- cell_data[[celltype_column]] == query_celltype
  if (sum(query_mask) < 1) {
    stop("No query cells found for '", query_celltype, "'.", call. = FALSE)
  }
  query_coords <- coords[query_mask, , drop = FALSE]
  nn <- RANN::nn2(query_coords, coords, k = 1)
  cell_data[, dist_to_query := as.vector(nn$nn.dists)]

  # Subset to target cells within max_distance
  target_mask <- cell_data[[celltype_column]] == target_celltype &
    cell_data$dist_to_query <= max_distance_um
  target_data <- cell_data[target_mask]
  target_barcodes <- target_data$barcode
  target_coords <- coords[target_mask, , drop = FALSE]

  .msg(
    "Target cells (", target_celltype, " within ", max_distance_um,
    " um): ", nrow(target_data)
  )

  # Total counts for offset
  target_counts <- count_matrix[, target_barcodes, drop = FALSE]
  total_counts <- Matrix::colSums(target_counts)

  # Validate genes exist
  available_genes <- intersect(genes, rownames(target_counts))
  missing_genes <- setdiff(genes, rownames(target_counts))
  if (length(missing_genes) > 0) {
    .msg("Genes not found in data: ", paste(missing_genes, collapse = ", "))
  }
  if (length(available_genes) == 0) {
    stop("None of the specified genes found in the count matrix.", call. = FALSE)
  }

  samples <- unique(target_data[[sample_column]])
  .msg("Samples: ", length(samples))
  .msg("Genes to check: ", paste(available_genes, collapse = ", "))

  # Per gene x sample: fit GLM, extract residuals, compute Moran's I
  results <- data.table::rbindlist(lapply(available_genes, function(g) {
    .msg("  Gene: ", g)

    data.table::rbindlist(lapply(samples, function(s) {
      samp_idx <- which(target_data[[sample_column]] == s)
      if (length(samp_idx) < max(k + 1, 30)) {
        return(data.table::data.table(
          gene = g, sample_id = s,
          morans_i = NA_real_, morans_expected = NA_real_,
          morans_pvalue = NA_real_, interpretation = "insufficient_cells",
          n_cells = length(samp_idx)
        ))
      }

      samp_counts <- as.numeric(target_counts[g, target_barcodes[samp_idx]])
      samp_dist <- target_data$dist_to_query[samp_idx]
      samp_total <- total_counts[target_barcodes[samp_idx]]
      samp_coords <- target_coords[samp_idx, , drop = FALSE]

      # Fit Poisson GLM (same as fit_poisson but we need the full model)
      valid <- !is.na(samp_counts) & !is.na(samp_dist) &
        samp_total > 0 & is.finite(samp_dist)
      if (sum(valid) < max(k + 1, 30) || sum(samp_counts[valid] > 0) < 5) {
        return(data.table::data.table(
          gene = g, sample_id = s,
          morans_i = NA_real_, morans_expected = NA_real_,
          morans_pvalue = NA_real_, interpretation = "insufficient_cells",
          n_cells = sum(valid)
        ))
      }

      y <- samp_counts[valid]
      d <- samp_dist[valid]
      log_total <- log(samp_total[valid])
      xy <- samp_coords[valid, , drop = FALSE]

      fit <- tryCatch(
        suppressWarnings(stats::glm(y ~ d + offset(log_total),
          family = stats::poisson()
        )),
        error = function(e) NULL
      )

      if (is.null(fit) || !fit$converged) {
        return(data.table::data.table(
          gene = g, sample_id = s,
          morans_i = NA_real_, morans_expected = NA_real_,
          morans_pvalue = NA_real_, interpretation = "glm_failed",
          n_cells = length(y)
        ))
      }

      # Deviance residuals
      resid <- stats::residuals(fit, type = "deviance")

      # Build spatial weights (k-NN)
      effective_k <- min(k, nrow(xy) - 1)
      knn <- spdep::knearneigh(xy, k = effective_k)
      nb <- spdep::knn2nb(knn)
      lw <- spdep::nb2listw(nb, style = "W")

      # Moran's I test
      mt <- tryCatch(
        spdep::moran.test(resid, lw, alternative = "two.sided"),
        error = function(e) NULL
      )

      if (is.null(mt)) {
        return(data.table::data.table(
          gene = g, sample_id = s,
          morans_i = NA_real_, morans_expected = NA_real_,
          morans_pvalue = NA_real_, interpretation = "moran_failed",
          n_cells = length(y)
        ))
      }

      mi <- mt$estimate["Moran I statistic"]
      me <- mt$estimate["Expectation"]
      mp <- mt$p.value

      interp <- if (abs(mi) < 0.05) {
        "none"
      } else if (abs(mi) < 0.15) {
        "weak"
      } else if (abs(mi) < 0.30) {
        "moderate"
      } else {
        "strong"
      }

      data.table::data.table(
        gene = g, sample_id = s,
        morans_i = unname(mi), morans_expected = unname(me),
        morans_pvalue = mp, interpretation = interp,
        n_cells = length(y)
      )
    }))
  }))

  .msg("\nSummary:")
  .msg("  Total tests: ", nrow(results[!is.na(morans_i)]))
  if (nrow(results[!is.na(morans_i)]) > 0) {
    .msg("  Median Moran's I: ", round(stats::median(results$morans_i, na.rm = TRUE), 4))
    interp_tab <- table(results$interpretation)
    .msg("  Interpretation: ", paste(names(interp_tab), interp_tab,
      sep = "=", collapse = ", "
    ))
  }

  results
}
