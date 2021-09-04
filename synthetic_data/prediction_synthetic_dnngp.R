library(SpatialTools)
library(pracma)
library(MLmetrics)

rho <- 1

pred_mat <- matrix(0, nrow = 100, ncol = 3)
for (i in 1:100){
  for (j in 1:3){
    pred_mat[i,j] <- exp_synth_holdout[i,j]
  }
}
pred_mat

ref_set <- matrix(0, nrow = 900, ncol = 3)
for (i in 1:900){
  for (j in 1:3){
    ref_set[i,j] <- exp_synth_model[i,j]
  }
}
ref_set

ys_true <- matrix(0, nrow = 100, ncol = 1)
for (i in 1:100){
  ys_true[i] <- exp_synth_holdout$y[i]
}
ys_true

y_N0 <- function(nns, ref_set, model_set){
  inds <- c()
  for (i in 1:nrow(nns)){
    ind <- vecMatch(ref_set, nns[i,])
    inds <- c(inds, ind)
  }
  return(model_set$y[inds])
}

predict_s0 <- function(s0, ref_set, m, params, model_data){
  sig2_avg <- params[1]
  tau2_avg <- params[2]
  phi_avg <- params[3]
  beta_avg <- params[4]
  
  nns <- nnsfunc(m, ref_set, s0)
  
  ccs <- rbind(s0, nns)
  ccsS <- ccs[,1:2]
  ccsT <- ccs[,3]
  
  fcm <- cov.st(coords = ccsS, time = ccsT, sp.type = "exponential",
         sp.par = c(sig2_avg, phi_avg), error.var = tau2_avg, t.type = "ar1", t.par = (1-1e-6))
  Eps_sn <- fcm$V[,1]
  Eps_sn <- Eps_sn[-1]
  Eps_ss <- fcm$V[1,1]
  Eps_nn <- cov.st(coords = nns[,1:2], time = nns[,3], sp.type = "exponential",
                   sp.par = c(sig2_avg, phi_avg), error.var = tau2_avg, t.type = "ar1", t.par = (1-1e-6))
  Eps_nn <- Eps_nn$V
  
  y_n <- y_N0(nns, ref_set, model_data)
  x_N0 <- as.matrix(rep(1, nrow(nns)))
  B <- as.matrix(beta_avg)
  
  c <- Eps_sn
  
  m <- dot(Eps_sn, solve(Eps_nn, y_n - x_N0%*%B))
  V <- Eps_ss - c%*%(solve(Eps_nn, c))
  y_new <- beta_avg + m + sqrt(V)*rnorm(1)
  return(y_new)
}

set.seed(11)
predict_s0(pred_mat[1,], ref_set, 15, c(sig2_avg, tau2_avg, phi_avg, beta_avg), exp_synth_model)

predict_all <- function(pred_set, ref_set, m, params, model_data){
  set.seed(11)
  ys_new <- matrix(0, nrow = nrow(pred_set), ncol = 1)
  for (i in 1:nrow(pred_set)){
    ys_new[i] <- predict_s0(pred_set[i,], ref_set, m, params, model_data)
  }
  return(ys_new)
}
ys_pred <- predict_all(pred_mat, ref_set, 15, c(sig2_avg, tau2_avg, phi_avg, beta_avg), exp_synth_model)


dnngp_RMSE <- RMSE(ys_pred, ys_true)
dnngp_MAE <- MAE(ys_pred, ys_true)
dnngp_R2 <- R2_Score(ys_pred, ys_true)

dnngp_synth <- data.frame(sigma2_estimate = c(sig2_avg), 
                          tau2_estimate = c(tau2_avg), 
                          phi_estimate = c(phi_avg), 
                          beta_estimate = c(beta_avg),
                          RMSE = c(dnngp_RMSE),
                          MAE = c(dnngp_MAE),
                          R2 = c(dnngp_R2)
                          )

View(dnngp_synth)

write.csv(dnngp_synth,"dnngp_synthetic_results.csv", row.names = TRUE)
