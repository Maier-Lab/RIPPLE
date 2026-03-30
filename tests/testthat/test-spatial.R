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
