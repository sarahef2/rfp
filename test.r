library(RLT)
library(MASS)
#library(survival)
library(randomForestSRC)

#use.cores = 6
use.cores = 1

P = 100
trainN = 10000
testN = 1000

rho = 0.25
V <- rho^abs(outer(1:P, 1:P, "-"))

set.seed(5)
X = as.matrix(mvrnorm(trainN + testN, mu=rep(0,P), Sigma=V))
Link <- function(x) x[,1] + x[,20] +x[,100] + x[,15]

set.seed(5)
T = exp(0.5*rnorm(trainN + testN, mean = Link(X)))
set.seed(5)
C = pmin(3, exp(0.5+0.5*rnorm(trainN + testN, mean = X[,5])))
Y = pmin(T, C)
Censor = (T <= C)

#plot(survfit(Surv(Y, Censor)~1))
mean(Censor)

# training and testing data
trainX = X[1:trainN, ]
trainY = Y[1:trainN]
trainCensor = Censor[1:trainN]

testX = X[-(1:trainN), ]
testY = Y[-(1:trainN)]
testCensor = Censor[-(1:trainN)]


# test the package
pre = Sys.time()
fit = survForest(trainX, trainY, trainCensor, subject.weight = runif(nrow(trainX), 1, 2), ntrees = 1, nmin = 15, split.gen = "best", split.rule = "logrank", nsplit = 10,
                 use.cores = use.cores, replacement = TRUE,resample.prob = 1)
Sys.time() - pre
pre = Sys.time()
fit = survForest_old(trainX, trainY, trainCensor, subject.weight = runif(nrow(trainX), 1, 2), ntrees = 1, nmin = 15, split.gen = "best", split.rule = "logrank", nsplit = 10,
                 use.cores = use.cores, replacement = TRUE,resample.prob = 1)
Sys.time() - pre
#Time difference of 20.11461 secs

pre = Sys.time()
fit = survForest(trainX, trainY, trainCensor, subject.weight = runif(nrow(trainX), 1, 2), ntrees = 2, nmin = 15, split.gen = "best", split.rule = "logrank", nsplit = 10,
                 use.cores = use.cores, replacement = TRUE,resample.prob = 1)
Sys.time() - pre
pre = Sys.time()
fit = survForest_old(trainX, trainY, trainCensor, subject.weight = runif(nrow(trainX), 1, 2), ntrees = 2, nmin = 15, split.gen = "best", split.rule = "logrank", nsplit = 10,
                 use.cores = use.cores, replacement = TRUE,resample.prob = 1)
Sys.time() - pre
#Time difference of 37.95218 secs

# test the package
pre = Sys.time()
fit = survForest(trainX, trainY, trainCensor, subject.weight = runif(nrow(trainX), 1, 2), ntrees = 1, nmin = 25, split.gen = "best", split.rule = "logrank", nsplit = 10,
                 use.cores = use.cores, replacement = TRUE,resample.prob = 1)
Sys.time() - pre
pre = Sys.time()
fit = survForest_old(trainX, trainY, trainCensor, subject.weight = runif(nrow(trainX), 1, 2), ntrees = 1, nmin = 25, split.gen = "best", split.rule = "logrank", nsplit = 10,
                 use.cores = use.cores, replacement = TRUE,resample.prob = 1)
Sys.time() - pre
#Time difference of 16.88461 secs

pre = Sys.time()
fit = survForest(trainX, trainY, trainCensor, subject.weight = runif(nrow(trainX), 1, 2), ntrees = 1, nmin = 15, split.gen = "random", split.rule = "logrank", nsplit = 1,
                 use.cores = use.cores, replacement = TRUE,resample.prob = 1)
Sys.time() - pre
pre = Sys.time()
fit = survForest_old(trainX, trainY, trainCensor, subject.weight = runif(nrow(trainX), 1, 2), ntrees = 1, nmin = 15, split.gen = "random", split.rule = "logrank", nsplit = 1,
                 use.cores = use.cores, replacement = TRUE,resample.prob = 1)
Sys.time() - pre
#Time difference of 5.707115 secs

pre = Sys.time()
fit = survForest(trainX, trainY, trainCensor, subject.weight = runif(nrow(trainX), 1, 2), ntrees = 1, nmin = 15, split.gen = "random", split.rule = "logrank", nsplit = 100,
                 use.cores = use.cores, replacement = TRUE,resample.prob = 1)
Sys.time() - pre
pre = Sys.time()
fit = survForest_old(trainX, trainY, trainCensor, subject.weight = runif(nrow(trainX), 1, 2), ntrees = 1, nmin = 15, split.gen = "random", split.rule = "logrank", nsplit = 100,
                 use.cores = use.cores, replacement = TRUE,resample.prob = 1)
Sys.time() - pre
#Time difference of 10.39944 secs


fit$FittedForest
table(fit$NodeRegi)

pre = Sys.time()
pred = predict(fit, testX, use.cores = use.cores)
Sys.time() - pre

plot(Link(X[-(1:trainN), ]), pred$surv[100,])

# categorical variables

trainX = data.frame(trainX)
testX = data.frame(testX)
for (j in 1:P)
{
  trainX[,j] = as.factor(trainX[,j]>0)
  testX[,j] = as.factor(testX[,j]>0)
}

pre = Sys.time()
fit = survForest(trainX, trainY, trainCensor, subject.weight = runif(nrow(trainX), 1, 2), ntrees = 10, nmin = 5, split.gen = "best", split.rule = "logrank", nsplit = 10,
                 verbose = TRUE, use.cores = use.cores, replacement = TRUE)
Sys.time() - pre

pre = Sys.time()
pred = predict(fit, testX, use.cores = use.cores)
Sys.time() - pre

plot(Link(X[-(1:trainN), ]), pred$surv[100,])







pre = Sys.time()
fit2 = rfsrc(Surv(Y, C)~., data.frame(trainX, "Y" = trainY, "C" = trainCensor), ntree = 10, nodesize = 10, importance = FALSE, mtry = P/3)
Sys.time() - pre


