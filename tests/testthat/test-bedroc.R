data("hpc_cv")

test_that('bedroc errors on non-binary predictions', {
  expect_error(dplyr::filter(hpc_cv, Resample == "Fold01") %>%
                 bedroc(obs, VF), "A multiclass `truth` input was provided, but only `binary` is supported.")
})

data("two_class_example")

test_that('bedroc works on a test data.frame', {
  expect_identical(two_class_example %>%
                     bedroc(truth, Class1) %>% mutate(.estimate = round(.estimate, 3)), tibble::tibble(.metric = 'bedroc',
                                                                                                   .estimator = 'binary',
                                                                                                   .estimate = 0.997))
})

test_that('bedroc works on a test vector', {
  expect_identical(round(bedroc_vec(two_class_example$truth, two_class_example$Class1), 3),  0.997)
})
