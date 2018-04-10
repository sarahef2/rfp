#' @title print.survForest
#' @description Print a survForest object
#' @param x A fitted survForest object
#' @param ... ...
#' @examples
#' x = matrix(rnorm(1000), ncol = 10)
#' y = exp(rowMeans(x))
#' c = rbinom(100, 1, 0.9)
#' fit = survForest(x, y, c, ntrees = 10)
#' fit

print.survForest<- function(x, ...)
{
  if (class(x)[2] == "fit")
    #.Call("survForestPrint", x$parameters)
    survForestPrint(x$parameters)

  if (class(x)[2] == "predict")
  {
    cat(paste("survForest prediction for", ncol(x$surv), "subjects on", nrow(x$surv), "unique time points: \n"))
    cat(paste("Minimum time point:", round(min(x$timepoints), 6), "\n"))
    cat(paste("Maximum time point:", round(max(x$timepoints), 6), "\n"))
  }
}

