#' @title survForest
#' @description Main function to fit survival forests
#' @param x A matrix or data.frame for features
#' @param y Response variable, a numeric/factor vector or a Surv object
#' @param censor The censoring indicator if survival model is used
#' @param ntrees Number of trees default is \code{ntrees = 500}
#' @param mtry Number of variables used at each internal node
#' @param split.gen How the cutting points are generated
#' @param split.rule How to compare the splits
#' @param nsplit Number of random cutting points to compare for each variable at an internal node
#' @param nmin Minimum number of observations reqired in an internal node to perform a split. Set this to twice of the desired terminal node size.
#' @param alpha Minimum number of observations required for each child node as a portion of the parent node. Must be within \code{(0, 0.5]}.
#' @param replacement Whether the in-bag samples are sampled with replacement
#' @param resample.prob Proportion of in-bag samples
#' @param subject.weight Subject weights
#' @param variable.weight Variable weights when randomly sample \code{mtry} to select the splitting rule
#' @param importance Should importance measures be calculated
#' @param use.cores Number of cores for parallel computing
#' @param verbose Printing additional information
#' @return A \code{survForest} object;
#' @examples
#' x = matrix(rnorm(1000), ncol = 10)
#' y = exp(rowMeans(x))
#' c = rbinom(100, 1, 0.9)
#' fit = survForest(x, y, c, ntrees = 10)

survForest <- function(x, y, censor,
        ntrees = 500,
        mtry = max(1, as.integer(ncol(x)/3)),
        split.gen = "random",
        split.rule = "logrank",
        nsplit = 1,
        nmin = max(1, as.integer(log(nrow(x)))),
        alpha = 0,
        replacement = FALSE,
        resample.prob = 0.632,
        subject.weight = NULL,
        variable.weight = NULL,
        importance = TRUE,
        nimpute = 1,
        use.cores = 0,
        verbose = FALSE)
{
  # check inputs

  if (missing(x)) stop("x is missing")
  if (missing(y)) stop("y is missing")
  if (missing(censor)) stop("censor is missing")

  if (!is.data.frame(x) & !is.matrix(x)) stop("x must be a matrix or dataframe")
  if (!is.vector(y)) stop("y must be a vector")
  if (!is.vector(censor)) stop("censor must be a vector")

  split.gen.mode = c("random", "rank", "best")
  match.arg(split.gen, split.gen.mode)

  split.rule.mode = c("logrank", "suplogrank")
  match.arg(split.rule, split.rule.mode)

  xnames <- colnames(x)
  n = nrow(x)
  p = ncol(x)

  if (any(is.na(x))) stop("NA not permitted in x")
  if (any(is.na(y))) stop("NA not permitted in y")
  if (any(is.na(censor))) stop("NA not permitted in censor")

  #if (nrow(x) != length(y)) stop("number of observations in x and y does not match")
  if (nrow(x) != length(y)) stop("number of observations in x and y do not match")
  #if (nrow(x) != length(censor)) stop("number of observations in x and censor does not match")
  if (nrow(x) != length(censor)) stop("number of observations in x and censor do not match")

  # prepare x, continuous and categorical
  if (is.data.frame(x))
  {
    # data.frame, check for categorical variables
    xlevels <- lapply(x, function(x) if (is.factor(x)) levels(x) else 0)
    ncat <- sapply(xlevels, length)
    x <- data.matrix(x)
  }else{
  #numerical matrix for x, all continuous
        ncat <- rep(1, p)
        xlevels <- as.list(rep(0, p))
  }

  storage.mode(ncat) <- "integer"
  storage.mode(x) <- "double"

  if (max(ncat) > 53)
        stop("Cannot handle categorical predictors with more than 53 categories")

  # prepare y, convert y to integers that represent the unique time failure time points
  if (!is.numeric(y)) stop("y should be numerical")
  if (!is.numeric(y)) stop("y should be numerical")

  if (any(y <= 0)) stop("y should be positive")

  timepoints = sort(unique(y[censor == 1]))
  y.point = rep(NA, length(y))
  for (i in 1:length(y))
  {
    if (censor[i] == 1)
      y.point[i] = match(y[i], timepoints)
    else
      y.point[i] = sum(y[i] >= timepoints)
  }

  y.point <- data.matrix(y.point)
  storage.mode(y.point) <- "integer"

  # print(sort(y.point))

  censor <- as.integer(censor)
  storage.mode(censor) <- "integer"

  interval = timepoints - c(0, timepoints[-length(timepoints)])
  interval <- data.matrix(c(0, interval))
  storage.mode(interval) <- "double"

  # set up and check parameters
  use.cores = max(1, use.cores)
  nmin = max(1, floor(nmin))
  resample.prob = max(0, min(resample.prob, 1))
  alpha = max(0, min(alpha, 0.5))

  if (split.gen == "random" && alpha > 0)
  {
    warning("Cannot use alpha > 0 when split.gen = random. The child node size requirment in this mode is not exact. alpha reset to zero.")
    alpha = 0
  }

  # subject weights for calculating the score of each split
  use.sub.weight = !is.null(subject.weight)
  if (is.null(subject.weight)) subject.weight = rep(1/n, n) else {subject.weight = subject.weight/sum(subject.weight)}
  if (length(subject.weight) != n) {warning("Subject weights length must be n, reset to equal weights"); subject.weight = rep(1/n, n); use.sub.weight = FALSE}
  if (any(subject.weight<=0)) {warning("Subject weights cannot be 0 or negative, reset to equal weights"); subject.weight = rep(1/n, n); use.sub.weight = FALSE}
  storage.mode(subject.weight) <- "double"

  # variable weights for mtry
  use.var.weight = !is.null(variable.weight)
  if (is.null(variable.weight)) variable.weight = rep(1/p, p) else {variable.weight = variable.weight/sum(variable.weight)}
  if (length(variable.weight) != p) {warning("Variable weights length must by p, reset to equal weights"); subject.weight = rep(1/p, p); use.var.weight = FALSE}
  if (any(variable.weight<0)) {warning("Variable weights cannot be negative, reset to equal weights"); subject.weight = rep(1/p, p); use.var.weight = FALSE}
  storage.mode(variable.weight) <- "double"

  if (use.var.weight && split.gen.mode == "best")
    stop("Cannot use the best splitting rule if variable weights is used. Switch split.gen to random or rank")

  # variable importantce

  resample.prob = max(0, min(1, resample.prob))

  if(importance & !replacement & resample.prob == 1)
    stop("Cannot perform importance if all observations are sampled without replacement")

  nimpute = max(1, nimpute)

  parameters = list("n" = as.integer(n),
                    "p" = as.integer(p),
                    "nfail" = as.integer(length(timepoints)),
                    "ntrees" = as.integer(ntrees),
                    "mtry" = as.integer(mtry),
                    "split.gen" = as.integer(match(split.gen, split.gen.mode)),
                    "split.rule" = as.integer(match(split.rule, split.rule.mode)),
                    "nsplit" = as.integer(nsplit),
                    "nmin" = as.integer(nmin),
                    "alpha" = as.double(alpha),
                    "replacement" = as.integer(replacement),
                    "resample.prob" = as.double(resample.prob),
                    "use.sub.weight" = as.integer(use.sub.weight),
                    "use.var.weight" = as.integer(use.var.weight),
                    "importance" = as.integer(importance),
                    "nimpute" = as.integer(nimpute),
                    "verbose" = as.integer(verbose))

  # fit model

  # fit = .Call("survForestFit",
  #             x,
  #             y.point,
  #             censor,
  #             ncat,
  #             interval,
  #             subject.weight,
  #             variable.weight,
  #             parameters,
  #             as.integer(use.cores))
  fit=survForestFit(x,y.point,censor,ncat,interval,subject.weight,variable.weight,parameters,as.integer(use.cores))

  fit[["ncat"]] = ncat
  fit[["subject.weight"]] = subject.weight

  if(use.var.weight)
    fit[["variable.weight"]] = variable.weight
  else
    fit[["variable.weight"]] = NULL

  fit[["parameters"]] = parameters
  fit[["y"]] = y
  fit[["censor"]] = censor
  fit[["timepoints"]] = timepoints
  fit[["y.point"]] = y.point

  class(fit) <- c("survForest", "fit")
  return(fit)
}
