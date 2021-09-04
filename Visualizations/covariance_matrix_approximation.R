library(spBayes)
library(GpGp)
library(matlib)

# useful functions
eucldist <- function(c1, c2){
  d1 <- (c1[1] - c2[1])
  d2 <- (c1[2] - c2[2])
  d  <- sqrt(d1**2 + d2**2)
  return(d)
}

vecMatch <- function(x, want) {
  if (length(x) == length(want)){
    if (x[1]==want[1]&x[2]==want[2]&x[3]==want[3]){
      return(1)
    }
    else {
      return(NULL)
    }
  }
  tfv <- apply(x, 1, function(x, want) isTRUE(all.equal(x, want)), want)
  which(tfv == TRUE)
}


hist_set <- function(m, ref_set, start_point){
  times <- ref_set[,3]
  r <- length(times)
  sqrm <- sqrt(m)
  coords <- ref_set[,1:2]
  sind <- vecMatch(ref_set, start_point)
  #print(sind)
  start_index <- sind
  #print(start_index)
  if (length(start_index) == 0){
    # start_point not in ref set
    spirs <- FALSE
  }
  else {
    spirs <- TRUE
  }
  #print(spirs)
  start_time <- start_point[3]
  history_set <- matrix(c(0,0,0),nrow=1, ncol=3)
  hist_times <- which(times < start_time)
  #print(hist_times)
  history_set <- rbind(history_set, ref_set[hist_times,])
  current_times <- which(times == start_time)
  #print(current_times)
  if (spirs){
    is <- c()
    for (i in current_times){
      if (i < start_index){
        is <- c(is, i)
      }
    }
    #print(is)
    if (length(is) != 0){
      history_set <- rbind(history_set, ref_set[is,])
    }
  }
  else {
    history_set <- rbind(history_set, ref_set[current_times,])
  }
  #history_set <- rbind(history_set, ref_set[current_times,])
  #print(history_set)
  return((history_set <- history_set[-1,]))
}


nn_from_hist <- function(m, history_set, start_point){
  ### DEFINE TIME GAP
  #print(history_set)
  time_gap <- 1
  if (length(history_set) == 3){
    return(history_set)
  }
  h <- length(history_set[,3])
  times <- history_set[,3]
  coords <- history_set[,1:2]
  rm <- sqrt(m)
  if (h <= m){
    #print('here')
    return(unique.matrix(history_set))
  }
  start_time <- start_point[3]
  nns <- matrix(c(0,0,0),nrow=1, ncol=3)
  # previous times
  for (i in (start_time-time_gap):(start_time-(rm+1)*time_gap)){
    eds <- rep(1e6, h)
    for (j in seq(1, h)){
      #print('times_j')
      #print(times[j])
      #print('i')
      #print(i)
      if (format(round(times[j], 5), nsmall = 5) == format(round(i, 5), nsmall = 5)){
        #print('YES')
        ed <- eucldist(start_point[1:2], coords[j,])
        eds[j] <- ed
      }
    }
    cor_inds <- which(eds < 1e6)
    #print(cor_inds)
    if (length(cor_inds) != 0){
      space_order_i <- order(eds)[1:min(length(cor_inds),rm)]
      #print(space_order_i)
      nearest_spatial_i <- history_set[space_order_i,]
      nns <- rbind(nns, (nearest_spatial_i))
    }
  }
  nns <- nns[-1,]
  # current time
  current_times <- which(times == start_time)
  #print('current times')
  #print(current_times)
  if (length(current_times) != 0){
    edst <- rep(1e6, length(current_times))
    for (i in current_times){
      ed <- eucldist(start_point[1:2], coords[i,])
      edst[i] <- ed
    }
    tinds <- which(edst < 1e6)
    time_order <- order(edst)[1:min(length(tinds),rm)]
    nearest_times <- history_set[time_order,]
    nns <- rbind(nns, nearest_times)
  }
  return(nns)
}

# overall nearest neighbor function

nnsfunc <- function(m, ref_set, start_point){
  history_set <- hist_set(m, ref_set, start_point)
  nns <- nn_from_hist(m, history_set, start_point)
  #print(nns)
  if (length(nns) <= 3){
    return(nns)
  }
  else {
    return(unique.matrix(nns))
  }
}

# correlation function as of DNNGP article
corrf <- function(c1, c2, sigma_squ, a, c, kap){
  s1 <- c1[1:2]
  s2 <- c2[1:2]
  t1 <- c1[3]
  t2 <- c2[3]
  h <- eucldist(s1, s2)
  u <- abs(t1-t2)
  f <- sigma_squ/((a*(u**2) + 1)**kap)
  s <- exp((-c*h)/((a*(u**2)+1)**(kap/2)))
  return(f*s)
}

# covariance matrix for all location possibilities
corrmat <- function(locs1, locs2, sigma_squ, a, c, kap){
  locs1_test <- rbind(locs1, c(0, 0, 0))
  locs2_test <- rbind(locs2, c(0, 0, 0))
  if (length(locs1_test[,1]) == 2){
    nrows <- 1
  }
  else {
    nrows <- length(locs1[,1])
  }
  if (length(locs2_test[,1]) == 2){
    ncols <- 1
  }
  else {
    ncols <- length(locs2[,1])
  }
  cMat <- matrix(0, nrow = nrows, ncol = ncols)
  for (i in 1:nrows){
    for (j in 1:ncols){
      if (nrows == 1 & ncols == 1){
        cMat[i,j] <- corrf(locs1, locs2, sigma_squ, a, c, kap)
      }
      else if (nrows == 1 & ncols != 1){
        cMat[i,j] <- corrf(locs1, locs2[j,], sigma_squ, a, c, kap)
      }
      else if (nrows != 1 & ncols == 1){
        cMat[i,j] <- corrf(locs1[i,], locs2, sigma_squ, a, c, kap)
      }
      else {
        cMat[i,j] <- corrf(locs1[i,], locs2[j,], sigma_squ, a, c, kap)
      }
    }
  }
  return(cMat)
}


# function to compute a_N(li)

a_vector <- function(l_i, nns, sigma_squ, a, c, kap){
  C_NN <- corrmat(nns, nns, sigma_squ, a, c, kap)
  C_Nl <- corrmat(nns, l_i, sigma_squ, a, c, kap)
  if (length(C_NN) == 1){
    iC <- 1/(C_NN)
  }
  else {
    # print('nns')
    # print(nns)
    # print('C_NN')
    # print(C_NN)
    iC <- inv(C_NN)
  }
  return(iC %*% C_Nl)
}


# function to compute diagonal element of matrix F i.e. f_li

f_li <- function(l_i, nns, sigma_squ, a, c, kap){
  corr <- corrf(l_i, l_i, sigma_squ, a, c, kap)
  C_lN <- corrmat(l_i, nns, sigma_squ, a, c, kap)
  C_NN <- corrmat(nns, nns, sigma_squ, a, c, kap)
  C_Nl <- corrmat(nns, l_i, sigma_squ, a, c, kap)
  #print((C_NN))
  if (length(C_NN) == 1){
    iC <- 1/(C_NN)
  }
  else {
    iC <- inv(C_NN)
  }
  return(corr - ((C_lN %*% iC) %*% C_Nl))
}

# function to compute F matrix

F_mat <- function(ref_set, m, sigma_squ, a, c, kap){
  r <- length(ref_set[,1])
  Fmatrix <- matrix(0, nrow = r, ncol = r)
  Fmatrix[1,1] <- corrf(ref_set[1,], ref_set[1,], sigma_squ, a, c, kap)
  for (i in 2:r){
    l_i <- ref_set[i,]
    nns <- nnsfunc(m, ref_set, l_i)
    #print(nns)
    f <- f_li(l_i, nns, sigma_squ, a, c, kap)
    #print(f)
    Fmatrix[i,i] <- f
  }
  return(Fmatrix)
}

# function to compute V matrix

V_mat <- function(ref_set, m, sigma_squ, a, c, kap){
  r <- length(ref_set[,1])
  V <- matrix(1e6, nrow = r, ncol = r)
  V[1,1] <- 1
  for (i in 2:r){
    V[i,1] <- 0
  }
  for (j in 2:r){
    nns_j <- nnsfunc(m, ref_set, ref_set[j,])
    inds <- c()
    for (i in 1:r){
      if (i == j){
        V[i,j] <- 1
      }
      else {
        ind_match <- vecMatch(nns_j, ref_set[i,])
        if (length(ind_match) == 0){
          # li not in N_lj
          V[i,j] <- 0
        }
        else {
          inds <- c(inds, i)
        }
      }
    }
    a_Nlj <- a_vector(ref_set[j,], nns_j, sigma_squ, a, c, kap)
    for (i in 1:length(inds)){
      V[inds[i],j] <- -a_Nlj[i]
    }
  }
  return(V)
}

### K FUNCTION

K_func <- function(ref_set, m, sigma_squ, a, c, kap){
  V <- V_mat(ref_set, m, sigma_squ, a, c, kap)
  #view(V)
  F_matrix <- F_mat(ref_set, m, sigma_squ, a, c, kap)
  #view(F_matrix)
  iF <- inv(F_matrix)
  Vt <- t(V)
  iK <- (Vt %*% (iF %*% V))
  #view(iK)
  return(iK)
}


K_func_test <- function(V, F_matrix){
  iF <- inv(F_matrix)
  Vt <- t(V)
  iK <- (Vt %*% iF) %*% V
  return(iK)
}


C_func_K_known <- function(l_i, l_j, ref_set, K, m, sigma_squ, a, c, kap){
  p <- vecMatch(ref_set, l_i)
  if (length(p) != 0){
    i_in_rs <- TRUE
  }
  else {
    i_in_rs <- FALSE
  }
  q <- vecMatch(ref_set, l_j)
  if (length(q) != 0){
    j_in_rs <- TRUE
  }
  else {
    j_in_rs <- FALSE
  }
  if (i_in_rs){
    if (j_in_rs){
      return(K[p,q])
    }
    else {
      nns_lj <- nnsfunc(m, ref_set, l_j)
      a_lj <- a_vector(l_j, nns_lj, sigma_squ, a, c, kap)
      return((t(a_lj))%*%K[,p])
    }
  }
  else {
    if (j_in_rs){
      nns_li <- nnsfunc(m, ref_set, l_i)
      a_li <- a_vector(l_i, nns_li, sigma_squ, a, c, kap)
      return((t(a_li))%*%K[,q])
    }
    else {
      nns_li <- nnsfunc(m, ref_set, l_i)
      a_li <- a_vector(l_i, nns_li, sigma_squ, a, c, kap)
      nns_lj <- nnsfunc(m, ref_set, l_j)
      a_lj <- a_vector(l_j, nns_lj, sigma_squ, a, c, kap)
      if (vecMatch(l_i, l_j) == 1){
        extr <- f_li(l_i, nns_li, sigma_squ, a, c, kap)
      }
      else {
        extr <- 0
      }
      main <- ((t(a_li))%*%K)%*%a_lj
      return(main + extr)
    }
  }
}

C_func <- function(l_i, l_j, ref_set, m, sigma_squ, a, c, kap){
  K <- K_func(ref_set, m, sigma_squ, a, c, kap)
  C <- C_func_K_known(l_i, l_j, ref_set, K, m, sigma_squ, a, c, kap)
  return(C)
}

C_mat <- function(ref_set, m, sigma_squ, a, c, kap){
  r <- length(ref_set[,1])
  K <- K_func(ref_set, m, sigma_squ, a, c, kap)
  print('K found')
  C <- matrix(0, nrow = r, ncol = r)
  count <- 0
  for (i in 1:r){
    for (j in 1:r){
      C[i,j] <- C_func_K_known(ref_set[i,], ref_set[j,], ref_set, K, m, sigma_squ, a, c, kap)
      count <- count + 1
      print(count)
    }
  }
  return(C)
}

### visualization of matrices

sigma_squ <- 1
a <- 50
c <- 25
kap <- 0.75
tau_squ <- 0.1

# full cov mat is C

cs <- cbind(synth_full$coords.V1, synth_full$coords.V2, synth_full$time)

V_test <- V_mat(cs, 4, sigma_squ, a, c, kap)
F_test <- F_mat(cs, 4, sigma_squ, a, c, kap)
K_test <- K_func_test(V_test, F_test)
iK <- K_test
iC <- inv(C)

C_plot_x <- c()
C_plot_y <- c()
K_plot_x <- c()
K_plot_y <- c()
for (i in 1:125){
  for (j in 1:125){
    if (iC[i,j] != 0){
      C_plot_x <- c(C_plot_x, i)
      C_plot_y <- c(C_plot_y, -j)
    }
    if (iK[i,j] != 0){
      K_plot_x <- c(K_plot_x, i)
      K_plot_y <- c(K_plot_y, -j)
    }
  }
}

par(mfrow = c(1,2))
plot(C_plot_x, C_plot_y, xaxt='n', yaxt='n', xlab="", ylab = "", main = expression(C^{-1}))
plot(K_plot_x, K_plot_y,xaxt='n', yaxt='n', xlab="", ylab = "", main = expression(tilde(C)^{-1}))