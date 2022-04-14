#' Get the correct event level
#'
#' @param truth A factor of the truth data
#' @param event_level The level of the factor that contains the event (the class of interest)
#'
#' @return String
#' @export
#'
#' @examples
#' x <- factor(c('TRUE', 'FALSE', 'TRUE', 'TRUE'))
#' truth_event_level(x, 'first')
truth_event_level <- function(truth, event_level) {
  if (identical(event_level, "first")) {
    levels(truth)[[1]]
  } else {
    levels(truth)[[2]]
  }
}
