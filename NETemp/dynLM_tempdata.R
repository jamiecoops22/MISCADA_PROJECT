library(fields)
library(geoR)
library(MBA)
library(spBayes)


ne.temp <- NETemp.dat
ne.temp <- ne.temp[reduced,]
y.t <- ne.temp[, 4:27]

N.t <- ncol(y.t)
n <- nrow(y.t)


miss <- sample(1:N.t, 10)
holdout.station.ids <- sample(1:n, 24)
y.t.holdout <- c()
for (i in 1:24){
  holdout.station.id <- holdout.station.ids[i]
  y.t.holdout_part <- y.t[holdout.station.id, miss]
  y.t[holdout.station.id, miss] <- NA
  y.t.holdout <- rbind(y.t.holdout, y.t.holdout_part)
}

coords <- as.matrix(ne.temp[, c("UTMX", "UTMY")]/1000)
max.d <- max(iDist(coords))
plot(coords, xlab = "Easting (km)", ylab = "Northin (km)")

p <- 2
starting <- list(beta = rep(0, N.t * p), phi = rep(3/(0.5 * max.d), N.t), 
                 sigma.sq = rep(2.1, N.t), tau.sq = rep(0.9, N.t), 
                 sigma.eta = diag(rep(0.01, p)))
tuning <- list(phi = rep(0.75, N.t))
priors <- list(beta.0.Norm = list(rep(0, p), diag(1000, p)), phi.Unif = list(rep(3/(0.9 * max.d), N.t), rep(3/(0.05 * max.d), N.t)), 
               sigma.sq.IG = list(rep(2, N.t), rep(1, N.t)), tau.sq.IG = list(rep(2, N.t), rep(1, N.t)),
               sigma.eta.IW = list(2, diag(0.001, p)))

mods <- lapply(paste(colnames(y.t), "elev", sep = "~"), as.formula)
n.samples <- 1000
dlm <-  spDynLM(mods, data = cbind(y.t, ne.temp[, "elev", drop = FALSE]),
                coords = coords, starting = starting, tuning = tuning, priors = priors,
                get.fitted = TRUE, cov.model = "exponential", n.samples = n.samples,
                n.report = 100)

### post-sampling

burn.in <- floor(0.75 * n.samples)
quant <- function(x) {
  quantile(x, prob = c(0.5, 0.025, 0.975))
}

beta <- apply(dlm$p.beta.samples[burn.in:n.samples, ], 2, quant)
beta.0 <- beta[, grep("Intercept", colnames(beta))]
beta.1 <- beta[, grep("elev", colnames(beta))]
par(mfrow = c(3, 2))

plot(1:N.t, beta.0[1, ], pch = 19, cex = 0.5, xlab = "months", ylab = "beta.0",
     ylim = range(beta.0))
arrows(1:N.t, beta.0[1, ], 1:N.t, beta.0[3, ], length = 0.02, angle = 90)
arrows(1:N.t, beta.0[1, ], 1:N.t, beta.0[2, ], length = 0.02, angle = 90)

plot(1:N.t, beta.1[1, ], pch = 19, cex = 0.5, xlab = "months", ylab = "beta.1",
     ylim = range(beta.1))
arrows(1:N.t, beta.1[1, ], 1:N.t, beta.1[3, ], length = 0.02, angle = 90)
arrows(1:N.t, beta.1[1, ], 1:N.t, beta.1[2, ], length = 0.02, angle = 90)

theta <- apply(dlm$p.theta.samples[burn.in:n.samples, ], 2, quant)
sigma.sq <- theta[, grep("sigma.sq", colnames(theta))]
tau.sq <- theta[, grep("tau.sq", colnames(theta))]
phi <- theta[, grep("phi", colnames(theta))]
par(mfrow = c(3, 1))

plot(1:N.t, sigma.sq[1, ], pch = 19, cex = 0.5, xlab = "months",
     ylab = "sigma.sq", ylim = range(sigma.sq))
arrows(1:N.t, sigma.sq[1, ], 1:N.t, sigma.sq[3, ], length = 0.02,
         angle = 90)
arrows(1:N.t, sigma.sq[1, ], 1:N.t, sigma.sq[2, ], length = 0.02,
         angle = 90)
plot(1:N.t, tau.sq[1, ], pch = 19, cex = 0.5, xlab = "months", ylab = "tau.sq",
       ylim = range(tau.sq))
arrows(1:N.t, tau.sq[1, ], 1:N.t, tau.sq[3, ], length = 0.02, angle = 90)
arrows(1:N.t, tau.sq[1, ], 1:N.t, tau.sq[2, ], length = 0.02, angle = 90)
plot(1:N.t, 3/phi[1, ], pch = 19, cex = 0.5, xlab = "months", ylab = "eff. range (km)",
       ylim = range(3/phi))
arrows(1:N.t, 3/phi[1, ], 1:N.t, 3/phi[3, ], length = 0.02, angle = 90)
arrows(1:N.t, 3/phi[1, ], 1:N.t, 3/phi[2, ], length = 0.02, angle = 90)


y.hat <- apply(dlm$p.y.samples[, burn.in:n.samples], 1, quant)
y.hat.med <- matrix(y.hat[1, ], ncol = N.t)
y.hat.up <- matrix(y.hat[3, ], ncol = N.t)
y.hat.low <- matrix(y.hat[2, ], ncol = N.t)


y.obs <- as.vector(as.matrix(y.t[-holdout.station.ids, -miss]))



y.obs.hat.med <- as.vector(y.hat.med[-holdout.station.ids, -miss])
y.obs.hat.up <- as.vector(y.hat.up[-holdout.station.ids, -miss])
y.obs.hat.low <- as.vector(y.hat.low[-holdout.station.ids, -miss])

y.ho <- as.matrix(y.t.holdout)

y.ho.hat.med <- as.vector(y.hat.med[holdout.station.ids, miss])
y.ho.hat.up <- as.vector(y.hat.up[holdout.station.ids, miss])
y.ho.hat.low <- as.vector(y.hat.low[holdout.station.ids, miss])

par(mfrow = c(1, 2))
plot(y.obs, y.obs.hat.med, pch = 19, cex = 0.5, xlab = "observed",
       ylab = "fitted", main = "Observed vs. fitted")
arrows(y.obs, y.obs.hat.med, y.obs, y.obs.hat.up, length = 0.02,
         angle = 90)
arrows(y.obs, y.obs.hat.med, y.obs, y.obs.hat.low, length = 0.02,
         angle = 90)
lines(-50:50, -50:50, col = "blue")
plot(y.ho, y.ho.hat.med, pch = 19, cex = 0.5, xlab = "observed",
       ylab = "predicted", , main = "Observed vs. predicted")
arrows(y.ho, y.ho.hat.med, y.ho, y.ho.hat.up, length = 0.02, angle = 90)
arrows(y.ho, y.ho.hat.med, y.ho, y.ho.hat.low, length = 0.02, angle = 90)
lines(-50:50, -50:50, col = "blue")


y.ho.hat.med
as.vector(y.ho)

dlm_ypred <- y.ho.hat.med
dlm_ytrue <- as.vector(y.ho)
dlm_ytrue
dlm_ypred

dlm_RMSE <- RMSE(dlm_ypred, dlm_ytrue)
dlm_MAE <- MAE(dlm_ypred, dlm_ytrue)
dlm_R2 <- R2_Score(dlm_ypred, dlm_ytrue)

dlm_RMSE
dlm_MAE
dlm_R2

dlmsig2_avg <- mean(sigma.sq[1,])
dlmtau2_avg <- mean(tau.sq[1,])
dlmphi_avg <- 3/mean(phi[1,])
dlmbeta0_avg <- mean(beta.0[1,])
dlmbeta1_avg <- mean(beta.1[1,])

dlm_tempresults <- data.frame(sigma2_estimate = c(dlmsig2_avg), 
                            tau2_estimate = c(dlmtau2_avg), 
                            phi_estimate = c(dlmphi_avg), 
                            beta0_estimate = c(dlmbeta0_avg),
                            beta1_estimate = c(dlmbeta1_avg),
                            RMSE = c(dlm_RMSE),
                            MAE = c(dlm_MAE),
                            R2 = c(dlm_R2)
)

View(dlm_tempresults)

write.csv(dlm_tempresults,"dlm_tempdata_results.csv", row.names = TRUE)
