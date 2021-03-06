#' Active/decoy predictions
#'
#' A dataset containing predictions on actives/decoys in a docking
#' run (~1200 compounds, ~3000 targets) generated by a random forest
#'
#' @format A data frame with 685800 rows and 9 variables:
#' \describe{
#'   \item{.config}{the preprocessor ID + the model ID}
#'   \item{.row}{the row from the training data}
#'   \item{truth}{whether a compound/target pairing is a known active}
#'   \item{.pred_TRUE}{probability that it's a true compound/target pairing}
#'   \item{.pred_FALSE}{probability that it's a false compound/target pairing}
#'   \item{.pred_class}{class prediction}
#'   ...
#' }
#' @source Internal docking run with GNINA and model fitting with tidymodels
"vs_rf_predictions"
