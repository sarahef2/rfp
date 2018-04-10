#' @title predict.survForest
#' @description The prediction function for personalized survival forest
#' @param object A fitted RLT object
#' @param testx Testing data
#' @param use.cores number of cores (by default use number of threads -1)
#' @param ... ...
#' @return The predicted values.
#' @examples
#' x = matrix(rnorm(1000), ncol = 10)
#' y = exp(rowMeans(x))
#' c = rbinom(100, 1, 0.9)
#' fit = survForest(x, y, c, ntrees = 10)
#' predict(fit, x)

predict.survForest <- function(object, testx, use.cores = 0, ...)
{
  # check test data
  if (missing(testx)) stop("testx is missing")
  if (!is.data.frame(testx) & !is.matrix(testx)) stop("x must be a matrix or dataframe")

  if (ncol(testx) != object$parameters$p)
    stop("dimension of testx is different from training data")

  testx = data.matrix(testx)
  storage.mode(testx) <- "double"

  pred = list()

  # pred[['surv']] = .Call("survForestPredict",
  #                       testx,
  #                       object$FittedForest,
  #                       object$y.point,
  #                       object$censor,
  #                       object$ncat,
  #                       object$subject.weight,
  #                       object$ObsTrack,
  #                       object$NodeRegi,
  #                       object$parameters,
  #                       as.integer(use.cores))[-1, ]
  pred[['surv']] = survForestPredict(
                         testx,
                         object$FittedForest,
                         object$y.point,
                         object$censor,
                         object$ncat,
                         object$subject.weight,
                         object$ObsTrack,
                         object$NodeRegi,
                         object$parameters,
                         as.integer(use.cores))[-1, ]

  pred[['timepoints']] = c(object$timepoints)

  class(pred) <- c("survForest", "predict")

  return(pred)
}

