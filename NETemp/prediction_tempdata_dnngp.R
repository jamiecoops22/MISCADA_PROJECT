library(SpatialTools)
library(pracma)
library(MLmetrics)

rho <- 1

rpred_mat <- matrix(0, nrow = 240, ncol = 3)
for (i in 1:240){
  for (j in 1:3){
    rpred_mat[i,j] <- temp_data_holdout[i,j]
  }
}
rpred_mat

rref_set <- matrix(0, nrow = 2160, ncol = 3)
for (i in 1:2160){
  for (j in 1:3){
    rref_set[i,j] <- temp_data_model[i,j]
  }
}
rref_set

rys_true <- matrix(0, nrow = 240, ncol = 1)
for (i in 1:240){
  rys_true[i] <- temp_data_holdout$y[i]
}
rys_true

y_N0 <- function(nns, ref_set, model_set){
  inds <- c()
  for (i in 1:nrow(nns)){
    ind <- vecMatch(ref_set, nns[i,])
    inds <- c(inds, ind)
  }
  return(model_set$y[inds])
}

x_N0 <- function(nns, ref_set, model_set){
  inds <- c()
  for (i in 1:nrow(nns)){
    ind <- vecMatch(ref_set, nns[i,])
    inds <- c(inds, ind)
  }
  return(model_set$x.2[inds])
}

x_s0 <- function(s0, pred_set, holdout_set){
  ind <- vecMatch(pred_set, s0)
  return(holdout_set$x.2[ind])
}

rpredict_s0 <- function(s0, pred_set, ref_set, m, params, model_data, holdout_data){
  sig2_avg <- params[1]
  tau2_avg <- params[2]
  phi_avg <- params[3]
  beta0_avg <- params[4]
  beta1_avg <- params[5]
  
  nns <- nnsfunc(m, ref_set, s0)
  
  
  ccs <- rbind(s0, nns)
  ccsS <- ccs[,1:2]
  ccsT <- ccs[,3]
  
  fcm <- cov.st(coords = ccsS, time = ccsT, sp.type = "exponential",
                sp.par = c(sig2_avg, phi_avg), error.var = tau2_avg, t.type = "ar1", t.par = .5)
  Eps_sn <- fcm$V[,1]
  Eps_sn <- Eps_sn[-1]
  Eps_ss <- fcm$V[1,1]
  Eps_nn <- cov.st(coords = nns[,1:2], time = nns[,3], sp.type = "exponential",
                   sp.par = c(sig2_avg, phi_avg), error.var = tau2_avg, t.type = "ar1", t.par = .5)
  Eps_nn <- Eps_nn$V
  
  y_n <- y_N0(nns, ref_set, model_data)
  x_n <- cbind((rep(1, nrow(nns))), x_N0(nns, ref_set, model_data))
  B <- as.matrix(c(beta0_avg, beta1_avg))
  x_s <- cbind(1, x_s0(s0, pred_set, holdout_data))
  
  c <- Eps_sn
  
  m <- dot(Eps_sn, solve(Eps_nn, y_n - x_n%*%B))
  V <- Eps_ss - c%*%(solve(Eps_nn, c))
  y_new <- x_s%*%B + m + sqrt(V)*rnorm(1)
  return(y_new)
}

set.seed(1)
rpredict_s0(rpred_mat[3,], rpred_mat, rref_set, 15, c(rsig2_avg, rtau2_avg, rphi_avg, rbeta0_avg, rbeta1_avg), temp_data_model, temp_data_holdout)

rpredict_all <- function(pred_set, ref_set, m, params, model_data, holdout_data){
  set.seed(1)
  ys_new <- matrix(0, nrow = nrow(pred_set), ncol = 1)
  for (i in 1:nrow(pred_set)){
    ys_new[i] <- rpredict_s0(pred_set[i,], pred_set, ref_set, m, params, model_data, holdout_data)
  }
  return(ys_new)
}
rys_pred <- rpredict_all(rpred_mat, rref_set, 15, c(rsig2_avg, rtau2_avg, rphi_avg, rbeta0_avg, rbeta1_avg), temp_data_model, temp_data_holdout)


a2rdnngp_RMSE <- RMSE(rys_pred, rys_true)
a2rdnngp_MAE <- MAE(rys_pred, rys_true)
a2rdnngp_R2 <- R2_Score(rys_pred, rys_true)

a2rdnngp_RMSE
a2rdnngp_MAE
a2rdnngp_R2

a2rdnngp_temp <- data.frame(sigma2_estimate = c(rsig2_avg), 
                          tau2_estimate = c(rtau2_avg), 
                          phi_estimate = c(rphi_avg), 
                          beta0_estimate = c(rbeta0_avg),
                          beta1_estimate = c(rbeta1_avg),
                          RMSE = c(a2rdnngp_RMSE),
                          MAE = c(a2rdnngp_MAE),
                          R2 = c(a2rdnngp_R2)
)

View(a2rdnngp_temp)

write.csv(a2rdnngp_temp,"dnngp_tempdata_results.csv", row.names = TRUE)
