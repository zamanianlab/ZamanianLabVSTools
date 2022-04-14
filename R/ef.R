# control what type of data can be used
finalize_estimator_internal.ef <- function(metric_dispatcher, x, estimator) {
  yardstick::validate_estimator(estimator, estimator_override = "binary")

  if (!is.null(estimator)) {
    return(estimator)
  }

  lvls <- levels(x)

  if (length(lvls) > 2) {
    stop("A multiclass `truth` input was provided, but only `binary` is supported.")
  }

  "binary"
}

#' Implementation of enrichment factor
#' @description Calculates a metric often used in virtual screening - enrichment factor.
#' @param data A data.frame containing the truth and estimate columns.
#' @param ... A set of unquoted column names or one or more dplyr selector functions to choose which variables contain the class probabilities. If truth is binary, only 1 column should be selected. Otherwise, there should be as many columns as factor levels of truth.
#' @param truth The column identifier for the true class results (that is a factor). This should be an unquoted column name although this argument is passed by expression and supports quasiquotation (you can unquote column names). For _vec() functions, a factor vector.
#' @param estimator One of "binary", "hand_till", "macro", or "macro_weighted" to specify the type of averaging to be done. "binary" is only relevant for the two class case. The others are general methods for calculating multiclass metrics. The default will automatically choose "binary" or "hand_till" based on truth.
#' @param event_level A single string. Either "first" or "second" to specify which level of truth to consider as the "event". This argument is only applicable when estimator = "binary". The default uses an internal helper that generally defaults to "first", however, if the deprecated global option yardstick.event_first is set, that will be used instead with a warning.
#' @param na_rm	A logical value indicating whether NA values should be stripped before the computation proceeds.
#' @param cutoff The percent value for enrichment (defaults to 0.01 or 1%)
#'
#' @return
#'
#' A tibble with columns .metric, .estimator, and .estimate and 1 row of values.
#'
#' For grouped data frames, the number of rows returned will be the same as the number of groups.
#'
#' For ef_vec(), a single numeric value (or NA).
#'
#' @seealso
#'
#' [bedroc()] for computing the Boltzmann-enhanced discrimination of ROC and [rie()] for computing the Robust Initial Enhancement.
#'
#' @export
ef <- function(data, ...) {
  UseMethod("ef")
}

ef <- yardstick::new_prob_metric(ef, direction = "maximize")

#' @rdname ef
#' @export
ef.data.frame <- function(data,
                          truth,
                          ...,
                          cutoff = 0.01,
                          estimator = NULL,
                          na_rm = TRUE,
                          event_level = "first") {
  estimate <- yardstick::dots_to_estimate(data, !!!enquos(...))

  yardstick::metric_summarizer(
    metric_nm = "ef",
    metric_fn = ef_vec,
    data = data,
    truth = !!enquo(truth),
    estimate = !!estimate,
    estimator = estimator,
    na_rm = na_rm,
    event_level = event_level,
    metric_fn_options = list(cutoff = cutoff)
  )
}

#' @export
#' @rdname ef
ef_vec <- function(truth,
                   estimate,
                   cutoff = 0.01,
                   estimator = NULL,
                   event_level = "first",
                   na_rm = TRUE) {
  estimator <- yardstick::finalize_estimator(truth, estimator, metric_class = "ef")

  ef_impl <- function(truth, estimate, cutoff = 0.01) {
    event <- truth_event_level(truth, event_level)

    N <- length(truth)

    A <- sum(truth == event)

    df <- dplyr::bind_cols(truth = truth, estimate = estimate) %>%
      dplyr::arrange(-estimate) %>%
      dplyr::slice_head(prop = cutoff) %>%
      dplyr::group_by(truth) %>%
      dplyr::tally()

    a <- dplyr::filter(df, truth == event)$n
    n <- sum(df$n)

    (a / n) / (A / N)
  }

  yardstick::metric_vec_template(
    metric_impl = ef_impl,
    truth = truth,
    estimate = estimate,
    estimator = estimator,
    na_rm = na_rm,
    cls = c("factor", "numeric"),
    cutoff = cutoff
  )
}


