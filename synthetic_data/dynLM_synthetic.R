library(fields)
library(geoR)
library(MBA)
library(spBayes)
library(MLmetrics)


View(exp_synth_full)

exp_full_mat <- as.matrix(exp_synth_full)
synth_mat <- matrix(0, nrow = 100, ncol = 13)
for (i in 1:100){
  synth_mat[i,1] <- 1
  synth_mat[i,2] <- exp_full_mat[i*10, 1]
  synth_mat[i,3] <- exp_full_mat[i*10, 2]
  for (j in 1:10){
    synth_mat[i,j+3] <- exp_full_mat[(i*10)-(10-j), 6]
  }
}
View(synth_mat)

exp_full_dlm <- data.frame(x = synth_mat[,1],
                           long = synth_mat[,2],
                           lat = synth_mat[,3],
                           y.1 = synth_mat[,4],
                           y.2 = synth_mat[,5],
                           y.3 = synth_mat[,6],
                           y.4 = synth_mat[,7],
                           y.5 = synth_mat[,8],
                           y.6 = synth_mat[,9],
                           y.7 = synth_mat[,10],
                           y.8 = synth_mat[,11],
                           y.9 = synth_mat[,12],
                           y.10 = synth_mat[,13])

View(exp_full_dlm)

y.t <- exp_full_dlm[, 4:13]

N.t <- ncol(y.t)
n <- nrow(y.t)



miss <- sample(1:N.t, 5)
holdout.station.ids <- sample(1:n, 20)
y.t.holdout <- c()
for (i in 1:20){
  holdout.station.id <- holdout.station.ids[i]
  y.t.holdout_part <- y.t[holdout.station.id, miss]
  y.t[holdout.station.id, miss] <- NA
  y.t.holdout <- rbind(y.t.holdout, y.t.holdout_part)
}
####

coords <- as.matrix(exp_full_dlm[, c("long", "lat")])
max.d <- max(iDist(coords))
plot(coords, xlab = "Easting (km)", ylab = "Northin (km)")

3/(0.9 * max.d)
3/(0.05 * max.d)

p <- 1
starting <- list(beta = rep(0, N.t * p), phi = rep(3/(0.5 * max.d), N.t), 
                 sigma.sq = rep(1, N.t), tau.sq = rep(0.1, N.t), 
                 sigma.eta = as.matrix(rep(0.01, p)))
tuning <- list(phi = rep(0.75, N.t))
priors <- list(beta.0.Norm = list(rep(0, p), diag(1000, p)), phi.Unif = list(rep(3/(0.9 * max.d), N.t), rep(3/(0.05 * max.d), N.t)), 
               sigma.sq.IG = list(rep(2, N.t), rep(1, N.t)), tau.sq.IG = list(rep(2, N.t), rep(1, N.t)),
               sigma.eta.IW = list(2, diag(0.001, p)))

mods <- lapply(paste(colnames(y.t), 1, sep = "~"), as.formula)
#mods <- ((colnames(y.t)) as.formula)
n.samples <- 1000
dlm <-  spDynLM(mods, data = cbind(y.t, exp_full_dlm[, "x", drop = FALSE]),
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
par(mfrow = c(2, 1))

plot(1:N.t, beta.0[1, ], pch = 19, cex = 0.5, xlab = "months", ylab = "beta.0",
     ylim = range(beta.0))
arrows(1:N.t, beta.0[1, ], 1:N.t, beta.0[3, ], length = 0.02, angle = 90)
arrows(1:N.t, beta.0[1, ], 1:N.t, beta.0[2, ], length = 0.02, angle = 90)



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

View(y.t)


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
plot(as.vector(y.ho), y.ho.hat.med, pch = 19, cex = 0.5, xlab = "observed",
     ylab = "predicted", , main = "Observed vs. predicted")
arrows(y.ho, y.ho.hat.med, y.ho, y.ho.hat.up, length = 0.02, angle = 90)
arrows(y.ho, y.ho.hat.med, y.ho, y.ho.hat.low, length = 0.02, angle = 90)
lines(-50:50, -50:50, col = "blue")


y.ho.hat.med
as.vector(y.ho)

synth_dlm_ypred <- y.ho.hat.med
synth_dlm_ytrue <- as.vector(y.ho)
synth_dlm_ytrue
synth_dlm_ypred

synth_dlm_RMSE <- RMSE(synth_dlm_ypred, synth_dlm_ytrue)
synth_dlm_MAE <- MAE(synth_dlm_ypred, synth_dlm_ytrue)
synth_dlm_R2 <- R2_Score(synth_dlm_ypred, synth_dlm_ytrue)

synth_dlm_RMSE
synth_dlm_MAE
synth_dlm_R2

synth_dlmsig2_avg <- mean(sigma.sq[1,])
synth_dlmtau2_avg <- mean(tau.sq[1,])
synth_dlmphi_avg <- 3/mean(phi[1,])
synth_dlmbeta0_avg <- mean(beta.0[1,])


synth_dlm_tempresults <- data.frame(sigma2_estimate = c(synth_dlmsig2_avg), 
                              tau2_estimate = c(synth_dlmtau2_avg), 
                              phi_estimate = c(synth_dlmphi_avg), 
                              beta_estimate = c(synth_dlmbeta0_avg),
                              RMSE = c(synth_dlm_RMSE),
                              MAE = c(synth_dlm_MAE),
                              R2 = c(synth_dlm_R2)
)

View(synth_dlm_tempresults)

write.csv(synth_dlm_tempresults,"new_dlm_synthdata_results.csv", row.names = TRUE)
