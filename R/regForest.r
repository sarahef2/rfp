#' @title regForest
#' @description Main function to fit regression forests
#' @param x A matrix or data.frame for features
#' @param y Response variable, a numeric/factor vector or a Surv object
#' @param ntrees Number of trees default is \code{ntrees = 500}
#' @param mtry Number of variables used at each internal node
#' @param split.gen How the cutting points are generated
#' @param split.rule How to compare the splits
#' @param nsplit Number of random cutting points to compare for each variable at an internal node
#' @param nmin Minimum number of observations required in an internal node to perform a split. Will split an internal node if it has 2*\code{nmin} observations. (CHECK)
#' @param nmin.control Should the terminal node size be forced to be at least nmin?  Default FALSE will split any node with 2*\code{nmin} observations without regard to nmin.
#' @param alpha Minimum number of observations required for each child node as a portion of the parent node. Must be within \code{(0, 0.5]}.
#' @param replacement Whether the in-bag samples are sampled with replacement
#' @param resample.prob Proportion of in-bag samples
#' @param subject.weight Subject weights
#' @param variable.weight Variable weights when randomly sample \code{mtry} to select the splitting rule
#' @param importance Should importance measures be calculated
#' @param nimpute Number of imputations for calculating the variable importance
#' @param use.cores Number of cores for parallel computing
#' @param verbose Printing additional information
#' @return A \code{regForest} object;
#' @examples
#' x = matrix(rnorm(1000), ncol = 10)
#' y = exp(rowMeans(x))
#' fit = regForest(x, y, c, ntrees = 10)

regForest <- function(x, y,
                      ntrees = 500,
                      mtry = max(1, as.integer(ncol(x)/3)),
                      split.gen = "random",
                      split.rule = "var",
                      nsplit = 1,
                      nmin = max(1, as.integer(log(nrow(x)))),
                      nmin.control = FALSE,
                      alpha = 0,
                      replacement = TRUE,
                      resample.prob = 1,
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

  if (!is.data.frame(x) & !is.matrix(x)) stop("x must be a matrix or dataframe")
  if (!is.vector(y)) stop("y must be a vector")

  split.gen.mode = c("random", "rank", "best")
  match.arg(split.gen, split.gen.mode)

  split.rule.mode = c("var")
  match.arg(split.rule, split.rule.mode)

  xnames <- colnames(x)
  n = nrow(x)
  p = ncol(x)

  if (any(is.na(x))) stop("NA not permitted in x")
  if (any(is.na(y))) stop("NA not permitted in y")

  if (nrow(x) != length(y)) stop("number of observations in x and y do not match")

  # prepare x, continuous and categorical
  if (is.data.frame(x))
  {
    # data.frame, check for categorical variables
    xlevels <- lapply(x, function(x) if (is.factor(x)) levels(x) else 0)
    ncat <- sapply(xlevels, length)
    x <- data.matrix(x)
  }else{
    # numerical matrix for x, all continuous
    ncat <- rep(1, p)
    xlevels <- as.list(rep(0, p))
  }

  storage.mode(ncat) <- "integer"
  storage.mode(x) <- "double"

  if (max(ncat) > 53)
        stop("Cannot handle categorical predictors with more than 53 categories")

  # prepare y
  if (!is.numeric(y)) stop("y should be numerical")

  storage.mode(y) <- "double"

  # set up and check parameters
  use.cores = as.integer(max(1, use.cores))
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
  if(any(variable.weight<0)) {
    warning("Variable weights must be greater than 0.  Setting all variable weights less than 0 as 0.")
    variable.weight = ifelse(variable.weight < 0, 0, variable.weight)
  }
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
                    "ntrees" = as.integer(ntrees),
                    "mtry" = as.integer(mtry),
                    "split.gen" = as.integer(match(split.gen, split.gen.mode)),
                    "split.rule" = as.integer(match(split.rule, split.rule.mode)),
                    "nsplit" = as.integer(nsplit),
                    "nmin" = as.integer(nmin),
                    "nmin.control" = as.integer(nmin.control),
                    "alpha" = as.double(alpha),
                    "replacement" = as.integer(replacement),
                    "resample.prob" = as.double(resample.prob),
                    "use.sub.weight" = as.integer(use.sub.weight),
                    "use.var.weight" = as.integer(use.var.weight),
                    "importance" = as.integer(importance),
                    "nimpute" = as.integer(nimpute),
                    "verbose" = as.integer(verbose))

  # fit model

  fit = regForestFit(x, y, ncat,
                     subject.weight,
                     variable.weight,
                     parameters,
                     use.cores)

  fit[["ncat"]] = ncat
  fit[["subject.weight"]] = subject.weight

  if(use.var.weight)
    fit[["variable.weight"]] = variable.weight
  else
    fit[["variable.weight"]] = NULL

  fit[["parameters"]] = parameters
  fit[["y"]] = y

  class(fit) <- c("RLT", "fit", "regForest")
  return(fit)
}
