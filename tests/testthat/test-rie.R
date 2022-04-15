test_that("rie errors on non-binary predictions", {
  expect_error(dplyr::filter(hpc_cv, Resample == "Fold01") %>%
    rie(obs, VF), "A multiclass `truth` input was provided, but only `binary` is supported.")
})

test_that("rie works on a test data.frame", {
  expect_identical(two_class_example %>%
    rie(truth, Class1) %>% dplyr::mutate(.estimate = base::round(.estimate, 3)), tibble::tibble(
    .metric = "rie",
    .estimator = "binary",
    .estimate = 1.931
  ))
})

test_that("rie works on a test vector", {
  expect_identical(base::round(rie_vec(two_class_example$truth, two_class_example$Class1), 3), 1.931)
})
