rm(list=ls())
library(spBayes)
library(fields)
library(MBA)
library(StempCens)
library(dplyr)

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

set.seed(321)

n <- 10*10*10
times <- as.matrix(seq(1,10))


cs2 <- seq(105,150,by=5)
cs3 <- matrix(0, nrow = 1, ncol = 2)
for (i in 1:10){
  for (j in 1:10){
    cs3 <- rbind(cs3, c(cs2[i], cs2[j]))
  }
}
cs3 <- cs3[-1,]
s_coords <- cs3

sigma2 <- 1
tau2 <- 0.1
phi <- 5
rho <- 1

Ms <- as.matrix(dist(s_coords))     # Spatial distances
Mt <- as.matrix(dist(times))      # Temporal distances
Cov <- CovarianceM(phi,rho,tau2,sigma2,distSpa=Ms,disTemp=Mt,kappa=0,type.S="exponential")
Cov

cs4 <- matrix(0, nrow = 1, ncol = 3)
for (i in 1:10){
  for (j in 1:10){
    for (k in 1:10){
      cs4 <- rbind(cs4, c(cs2[i], cs2[j], k))
    }
  }
}
cs4 <- cs4[-1,]
coords <- cs4

#x1 <- (runif(n, 1,5))
incpt <- as.matrix(rep(1, n))
#x <- cbind(incpt, x1)
x <- incpt
B <- as.matrix(1)
w <- rmvn(1, rep(0,n), Cov)
y <- rnorm(n, x%*%B + w, sqrt(tau2))

exp_synth_full <- data.frame(coords = coords,
                             x = x,
                             w = w,
                             y = y)

exp_synth_full <- exp_synth_full %>%
  rename(
    Day = coords.3,
    Longitude = coords.1,
    Latitude = coords.2
  )



set.seed(1)
ho <- sort(sample(1:nrow(coords), 100))

exp_synth_holdout <- data.frame(coords = coords[ho,],
                                x = x[ho,],
                                w = w[ho],
                                y = y[ho])

exp_synth_holdout <- exp_synth_holdout %>%
  rename(
    Day = coords.3
  )

exp_synth_model <- data.frame(coords = coords[-ho,],
                              x = x[-ho,],
                              w = w[-ho],
                              y = y[-ho])

exp_synth_model <- exp_synth_model %>%
  rename(
    Day = coords.3
  )

exp_synth_model <- exp_synth_model %>%
  rename(
    Longitude = coords.1,
    Latitude = coords.2
  )

exp_synth_holdout <- exp_synth_holdout %>%
  rename(
    Longitude = coords.1,
    Latitude = coords.2
  )

View(exp_synth_model)

write.table(coords[ho,], "coords.ho", sep="\t", row.names=F, col.names=F)
write.table(y[ho], "y.ho", sep="\t", row.names=F, col.names=F)
write.table(x[ho,], "x.ho", sep="\t", row.names=F, col.names=F)
write.table(w[ho], "w.ho", sep="\t", row.names=F, col.names=F)

write.table(coords[-ho,], "coords.mod", sep="\t", row.names=F, col.names=F)
write.table(y[-ho], "y.mod", sep="\t", row.names=F, col.names=F)
write.table(x[-ho,], "x.mod", sep="\t", row.names=F, col.names=F)
write.table(w[-ho], "w.mod", sep="\t", row.names=F, col.names=F)