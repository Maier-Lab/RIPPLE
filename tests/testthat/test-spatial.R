test_that("get_coord_columns auto-detects spatial_x/spatial_y", {
  meta <- data.frame(spatial_x = 1:10, spatial_y = 1:10, other = 1:10)
  result <- get_coord_columns(meta)
  expect_equal(result, c("spatial_x", "spatial_y"))
})

test_that("get_coord_columns auto-detects x/y", {
  meta <- data.frame(x = 1:10, y = 1:10, other = 1:10)
  result <- get_coord_columns(meta)
  expect_equal(result, c("x", "y"))
})

test_that("get_coord_columns uses explicit columns", {
  meta <- data.frame(cx = 1:10, cy = 1:10)
  result <- get_coord_columns(meta, x_col = "cx", y_col = "cy")
  expect_equal(result, c("cx", "cy"))
})

test_that("get_coord_columns errors on missing columns", {
  meta <- data.frame(a = 1:10, b = 1:10)
  expect_error(get_coord_columns(meta))
})

test_that("build_knn_graph returns correct structure", {
  coords <- matrix(runif(200), ncol = 2)
  result <- build_knn_graph(coords, k = 5)
  expect_equal(ncol(result$indices), 5)
  expect_equal(nrow(result$indices), 100)
})

test_that("calculate_distance_to_type computes distances", {
  coords <- matrix(c(0,0, 1,0, 2,0, 10,10), ncol = 2, byrow = TRUE)
  cell_types <- c("A", "B", "A", "B")

  dists <- calculate_distance_to_type(coords, cell_types, "B")

  # Cell 1 (A at 0,0): nearest B is at (1,0), dist = 1
  expect_equal(dists[1], 1.0)
  # Cell 2 (B at 1,0): distance to itself should be 0 or to other B
  expect_true(dists[2] >= 0)
})

test_that("check_spatial_autocorrelation works on synthetic data", {
  skip_if_not_installed("spdep")
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("S4Vectors")

  set.seed(42)
  n <- 200

  # Build a synthetic SPE with query + target cells
  x <- runif(n, 0, 500)
  y <- runif(n, 0, 500)
  ct <- c(rep("query", 30), rep("target", n - 30))

  # Gene with spatially correlated expression near query
  dists_to_query <- sapply(seq_len(n), function(i) {
    min(sqrt((x[i] - x[1:30])^2 + (y[i] - y[1:30])^2))
  })
  gene_counts <- rpois(n, lambda = exp(-2 - 0.005 * dists_to_query))
  bg_counts <- rpois(n, lambda = 2)

  counts <- rbind(gene_counts, bg_counts)
  rownames(counts) <- c("GRADIENT_GENE", "FLAT_GENE")
  colnames(counts) <- paste0("cell_", seq_len(n))
  counts <- methods::as(counts, "CsparseMatrix")

  meta <- data.frame(
    cell_type = ct,
    sample_id = "s1",
    row.names = colnames(counts)
  )
  coords_mat <- cbind(x, y)
  rownames(coords_mat) <- colnames(counts)

  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(meta),
    spatialCoords = coords_mat
  )

  result <- check_spatial_autocorrelation(
    input           = spe,
    genes           = c("GRADIENT_GENE", "FLAT_GENE"),
    celltype_column = "cell_type",
    target_celltype = "target",
    query_celltype  = "query",
    sample_column   = "sample_id",
    k               = 10,
    verbose         = FALSE
  )

  expect_s3_class(result, "data.table")
  expect_true(all(c("gene", "sample_id", "morans_i", "morans_pvalue",
                    "interpretation") %in% names(result)))
  expect_equal(nrow(result), 2)  # 2 genes x 1 sample
  expect_true(all(!is.na(result$morans_i)))
})
