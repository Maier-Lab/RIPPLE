make_synthetic <- function(n_genes = 10, n_cells = 20, seed = 42) {
  set.seed(seed)
  counts <- matrix(
    rpois(n_genes * n_cells, lambda = 2),
    nrow = n_genes, ncol = n_cells
  )
  rownames(counts) <- paste0("gene", seq_len(n_genes))
  colnames(counts) <- paste0("cell", seq_len(n_cells))

  meta <- data.frame(
    cell_type = rep(c("A", "B"), length.out = n_cells),
    sample_id = rep(c("s1", "s2"), each = n_cells / 2),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- colnames(counts)

  coords <- matrix(runif(n_cells * 2, 0, 100), ncol = 2)
  rownames(coords) <- colnames(counts)
  colnames(coords) <- c("x", "y")

  list(counts = counts, meta = meta, coords = coords)
}

test_that("make_ripple_input builds a SpatialExperiment", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("S4Vectors")

  syn <- make_synthetic()
  spe <- make_ripple_input(
    counts       = syn$counts,
    metadata     = syn$meta,
    coords       = syn$coords,
    output_class = "SpatialExperiment"
  )

  expect_s4_class(spe, "SpatialExperiment")
  expect_equal(ncol(spe), 20)
  expect_equal(nrow(spe), 10)
  expect_true(all(c("cell_type", "sample_id") %in% names(SummarizedExperiment::colData(spe))))
})

test_that("make_ripple_input builds a Seurat object", {
  skip_if_not_installed("Seurat")

  syn <- make_synthetic()
  seu <- make_ripple_input(
    counts       = syn$counts,
    metadata     = syn$meta,
    coords       = syn$coords,
    output_class = "Seurat"
  )

  expect_s4_class(seu, "Seurat")
  expect_equal(ncol(seu), 20)
  expect_true(all(c("cell_type", "sample_id", "x", "y") %in% colnames(seu@meta.data)))
})

test_that(".resolve_input handles Seurat input", {
  skip_if_not_installed("Seurat")

  syn <- make_synthetic()
  seu <- make_ripple_input(
    counts       = syn$counts,
    metadata     = syn$meta,
    coords       = syn$coords,
    output_class = "Seurat"
  )

  res <- ripple:::.resolve_input(seu, verbose = FALSE)

  expect_true(is.list(res))
  expect_true(all(c("counts", "meta") %in% names(res)))
  expect_equal(ncol(res$counts), 20)
  expect_true("barcode" %in% names(res$meta))
  expect_true(all(c("cell_type", "sample_id") %in% names(res$meta)))
})

test_that(".resolve_input handles SpatialExperiment input", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("S4Vectors")

  syn <- make_synthetic()
  spe <- make_ripple_input(
    counts       = syn$counts,
    metadata     = syn$meta,
    coords       = syn$coords,
    output_class = "SpatialExperiment"
  )

  res <- ripple:::.resolve_input(spe, verbose = FALSE)

  expect_true(is.list(res))
  expect_equal(ncol(res$counts), 20)
  expect_true("barcode" %in% names(res$meta))
  # Spatial coordinates should have been folded into metadata
  expect_true(any(c("x", "y") %in% names(res$meta)) ||
              any(grepl("^x|^y", names(res$meta), ignore.case = TRUE)))
})

test_that(".resolve_input returns normalized expression when requested", {
  skip_if_not_installed("Seurat")

  syn <- make_synthetic()
  seu <- make_ripple_input(
    counts       = syn$counts,
    metadata     = syn$meta,
    coords       = syn$coords,
    output_class = "Seurat"
  )

  res <- ripple:::.resolve_input(seu, require_expr = TRUE, verbose = FALSE)
  expect_true("expr" %in% names(res))
  expect_equal(dim(res$expr), dim(res$counts))
  expect_true(all(res$expr >= 0))
})

test_that("read_ripple_csv round-trips from disk", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("S4Vectors")

  syn <- make_synthetic()

  tmp_dir <- tempfile("ripple_csv_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE))

  # Write counts.csv: gene column + cells
  counts_df <- data.frame(
    gene = rownames(syn$counts),
    as.data.frame(syn$counts)
  )
  data.table::fwrite(counts_df, file.path(tmp_dir, "counts.csv"))

  # Write metadata.csv: barcode column + metadata
  meta_df <- data.frame(
    barcode = rownames(syn$meta),
    syn$meta
  )
  data.table::fwrite(meta_df, file.path(tmp_dir, "metadata.csv"))

  # Write coords.csv: barcode column + x, y
  coords_df <- data.frame(
    barcode = rownames(syn$coords),
    x = syn$coords[, 1],
    y = syn$coords[, 2]
  )
  data.table::fwrite(coords_df, file.path(tmp_dir, "coords.csv"))

  spe <- read_ripple_csv(tmp_dir, output_class = "SpatialExperiment")
  expect_s4_class(spe, "SpatialExperiment")
  expect_equal(ncol(spe), 20)
})

test_that(".resolve_input errors on unsupported input", {
  expect_error(
    ripple:::.resolve_input(42, verbose = FALSE),
    "Unsupported input type"
  )
})

test_that("make_ripple_input aligns cells between counts and metadata", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("S4Vectors")

  syn <- make_synthetic(n_cells = 20)
  # Drop 5 cells from metadata
  subset_meta <- syn$meta[1:15, ]
  subset_coords <- syn$coords[1:15, ]

  expect_message(
    spe <- make_ripple_input(
      counts       = syn$counts,
      metadata     = subset_meta,
      coords       = subset_coords,
      output_class = "SpatialExperiment"
    ),
    "Subsetting counts"
  )

  expect_equal(ncol(spe), 15)
})
