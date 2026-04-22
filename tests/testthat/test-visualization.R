test_that("plot_gene_category_dotplot builds a ggplot from curated panel", {
  skip_if_not_installed("ggplot2")

  results <- data.table::data.table(
    gene = c("G1", "G2", "G3", "G4", "G5", "G6"),
    median_coef = c(-0.005, -0.001, 0.003, 0.002, -0.004, 0.001),
    fisher_fdr = c(1e-10, 0.2, 1e-5, 0.8, 1e-20, 0.03)
  )

  panel <- list(
    "Category A" = c("G1", "G2"),
    "Category B" = c("G3", "G4"),
    "Category C" = c("G5", "G6")
  )

  p <- plot_gene_category_dotplot(results, panel, query_label = "Test")

  expect_s3_class(p, "ggplot")
  expect_true("category" %in% names(p$data))
  expect_setequal(
    as.character(levels(p$data$category)),
    c("Category A", "Category B", "Category C")
  )
  expect_equal(nrow(p$data), 6)
  expect_true(all(p$data$is_sig == (p$data$fisher_fdr < 0.05)))
})

test_that("plot_gene_category_dotplot drops missing genes and warns when empty", {
  results <- data.table::data.table(
    gene = c("G1", "G2"),
    median_coef = c(-0.001, 0.001),
    fisher_fdr = c(0.01, 0.5)
  )

  # Genes not in results — should error
  expect_error(
    plot_gene_category_dotplot(
      results,
      list("Missing" = c("Z1", "Z2"))
    ),
    "None of the genes"
  )

  # Partial match — keeps overlap only
  p <- plot_gene_category_dotplot(
    results,
    list("A" = c("G1", "Z_absent"), "B" = c("G2"))
  )
  expect_equal(nrow(p$data), 2)
})

test_that("plot_gene_category_dotplot validates inputs", {
  results <- data.table::data.table(
    gene = "G1", median_coef = 0, fisher_fdr = 0.5
  )

  expect_error(
    plot_gene_category_dotplot(results, c("G1", "G2")),
    "named list"
  )
  expect_error(
    plot_gene_category_dotplot(results, list(c("G1"))),
    "named list"
  )
  expect_error(
    plot_gene_category_dotplot(
      data.table::data.table(gene = "G1"),
      list("A" = "G1")
    ),
    "Required column"
  )
})
