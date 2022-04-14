# control what type of data can be used
finalize_estimator_internal.bedroc <- function(metric_dispatcher, x, estimator) {
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

#' Implementation of Boltzmann-enhanced discrimination of receiver operating characteristic (BEDROC)
#' @description Calculates a metric often used in virtual screening - BEDROC
#' @param data A data.frame containing the truth and estimate columns.
#' @param ... A set of unquoted column names or one or more dplyr selector functions to choose which variables contain the class probabilities. If truth is binary, only 1 column should be selected. Otherwise, there should be as many columns as factor levels of truth.
#' @param truth The column identifier for the true class results (that is a factor). This should be an unquoted column name although this argument is passed by expression and supports quasiquotation (you can unquote column names). For _vec() functions, a factor vector.
#' @param estimator One of "binary", "hand_till", "macro", or "macro_weighted" to specify the type of averaging to be done. "binary" is only relevant for the two class case. The others are general methods for calculating multiclass metrics. The default will automatically choose "binary" or "hand_till" based on truth.
#' @param event_level A single string. Either "first" or "second" to specify which level of truth to consider as the "event". This argument is only applicable when estimator = "binary". The default uses an internal helper that generally defaults to "first", however, if the deprecated global option yardstick.event_first is set, that will be used instead with a warning.
#' @param na_rm	A logical value indicating whether NA values should be stripped before the computation proceeds.
#' @param alpha The early recognition parameter (defaults to 20 or as in enrichvs::bedroc())
#'
#' @references
#'
#' Truchon J-F, Bayly CI. Evaluating virtual screening methods: good and bad metrics for the “early recognition” problem. J Chem Inf Model. 2007 Mar;47(2):488–508.
#'
#' @return
#'
#' A tibble with columns .metric, .estimator, and .estimate and 1 row of values.
#'
#' For grouped data frames, the number of rows returned will be the same as the number of groups.
#'
#' For bedroc_vec(), a single numeric value (or NA).
#'
#' #' @seealso
#'
#' [ef()] for computing the enrichment factor.
#'
#' @export
bedroc <- function(data, ...) {
  UseMethod("bedroc")
}

bedroc <- yardstick::new_prob_metric(bedroc, direction = "maximize")

#' @rdname bedroc
#' @export
bedroc.data.frame <- function(data,
                              truth,
                              ...,
                              alpha = 20,
                              estimator = NULL,
                              na_rm = TRUE,
                              event_level = "first") {
  estimate <- yardstick::dots_to_estimate(data, !!!enquos(...))

  yardstick::metric_summarizer(
    metric_nm = "bedroc",
    metric_fn = bedroc_vec,
    data = data,
    truth = !!enquo(truth),
    estimate = !!estimate,
    estimator = estimator,
    na_rm = na_rm,
    event_level = event_level,
    metric_fn_options = list(alpha = alpha)
  )
}

#' @export
#' @rdname bedroc
bedroc_vec <- function(truth,
                       estimate,
                       alpha = 20,
                       estimator = NULL,
                       event_level = "first",
                       na_rm = TRUE) {
  estimator <- yardstick::finalize_estimator(truth, estimator, metric_class = "bedroc")

  bedroc_impl <- function(truth, estimate, alpha = 20) {
    event <- truth_event_level(truth, event_level)

    N <- length(truth)

    n <- length(which(truth == event))
    ord <- order(estimate, decreasing = TRUE)
    m_rank <- which(truth[ord] == event)
    s <- sum(exp(-alpha * m_rank / N))
    ra <- n / N
    ri <- (N - n) / N
    random_sum <- ra * exp(-alpha / N) * (1.0 - exp(-alpha)) / (1.0 - exp(-alpha / N))
    fac <- ra * sinh(alpha / 2.0) / (cosh(alpha / 2.0) - cosh(alpha / 2.0 - alpha * ra))
    cte <- 1.0 / (1 - exp(alpha * ri))
    s / random_sum * fac + cte
  }

  yardstick::metric_vec_template(
    metric_impl = bedroc_impl,
    truth = truth,
    estimate = estimate,
    estimator = estimator,
    na_rm = na_rm,
    cls = c("factor", "numeric"),
    alpha = alpha
  )
}
