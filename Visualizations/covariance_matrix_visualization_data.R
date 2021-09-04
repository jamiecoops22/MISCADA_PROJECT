rm(list=ls())
library(spBayes)
library(fields)
library(MBA)
library(tidyverse)

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

set.seed(321)

n <- 5*5*5
cs1 <- c(1:5)
times <- c()
for (i in cs1){
  times <- c(times, i, i, i, i, i)
}
coords <- cbind(runif(n, 0, 1), runif(n, 0, 1), times)

ord <- order(coords[,3])
coords <- coords[ord,]

X <- cbind(rnorm(n, 0, 1), rnorm(n, 0, 1))

b0 <- 1
b1 <- 5
B <- c(b0, b1)

sigma_squ <- 1
a <- 50
c <- 25
kap <- 0.75
tau_squ <- 0.1


C <- corrmat(coords, coords, sigma_squ, a, c, kap)

w <- rmvn(1, rep(0,n), C)
eps <- rnorm(n, 0, tau_squ)

y <- c()
for (i in 1:n){
  yi <- t(X[i,])%*%B + w[i] + eps[i]
  y <- c(y, yi)
}

synth_full <- data.frame(coords = coords,
                         x = X,
                         w = w,
                         y = y)

synth_full <- synth_full %>%
  rename(time = coords.times)

View(synth_full)