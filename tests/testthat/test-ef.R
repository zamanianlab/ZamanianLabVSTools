data("hpc_cv")

test_that('ef errors on non-binary predictions', {
  expect_error(dplyr::filter(hpc_cv, Resample == "Fold01") %>%
                 ef(obs, VF), "A multiclass `truth` input was provided, but only `binary` is supported.")
})

data("two_class_example")

test_that('ef works on a test data.frame', {
  expect_identical(two_class_example %>%
                     ef(truth, Class1) %>% mutate(.estimate = round(.estimate, 3)), tibble::tibble(.metric = 'ef',
                                                                                                    .estimator = 'binary',
                                                                                                    .estimate = 1.938))
})

test_that('ef works on a test vector', {
  expect_identical(round(ef_vec(two_class_example$truth, two_class_example$Class1), 3),  1.938)
})
