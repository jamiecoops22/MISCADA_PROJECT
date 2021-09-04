library(coda)

### guess sig2 2 (1), tau2 1 (1), phi 0.003 (0.001, 0.02), beta (10,0), tuning 0.00003, seed 1

rtchain1 <- read.table("C:/project-files/TEMP_DATA/chain1/temp-expchain-theta", sep = "" , header = F , nrows = 3,
                      na.strings ="", stringsAsFactors= F)
rbchain1 <- read.table("C:/project-files/TEMP_data/chain1/temp-expchain-beta", sep = "" , header = F , nrows = 2,
                      na.strings ="", stringsAsFactors= F)
View(rbchain1)

nc <- ncol(rtchain1)

rchain1_mat <- matrix(0, nrow = nc, ncol = 3)
for (i in 1:3){
  for (j in 1:nc){
    rchain1_mat[j,i] <- rtchain1[i,j]
  }
}
rchain1_mat[,3] <- 3/rchain1_mat[,3]
rchain1_mat <- cbind(rchain1_mat, t(rbchain1))
View(rchain1_mat)

rchain1 <- mcmc(data = rchain1_mat)
varnames(rchain1) <- c('sigma2', 'tau2', 'phi', 'beta0', 'beta1')


plot(rchain1, density = FALSE)

# remove burn in 

rchain1_b <- mcmc(data = rchain1_mat[2000:5000,])
summary(rchain1_b)

rsig2_c1 <- summary(rchain1_b)$statistics[1]
rtau2_c1 <- summary(rchain1_b)$statistics[2]
rphi_c1 <- summary(rchain1_b)$statistics[3]
rbeta0_c1 <- summary(rchain1_b)$statistics[4]
rbeta1_c1 <- summary(rchain1_b)$statistics[5]


### guess sig2 1.5, tau2 1.5, phi 0.005, beta (5,1), tuning 0.00003, seed 2

# rtchain2 <- read.table("C:/project-files/TEMP_DATA/chain2/temp-expchain2-theta", sep = "" , header = F , nrows = 3,
#                        na.strings ="", stringsAsFactors= F)
# rbchain2 <- read.table("C:/project-files/TEMP_data/chain2/temp-expchain2-beta", sep = "" , header = F , nrows = 2,
#                        na.strings ="", stringsAsFactors= F)

rtchain2 <- read.table("C:/project-files/TEMP_DATA/attempt2/chain1/rda2-chain1-theta", sep = "" , header = F , nrows = 3,
                       na.strings ="", stringsAsFactors= F)
rbchain2 <- read.table("C:/project-files/TEMP_data/attempt2/chain1/rda2-chain1-beta", sep = "" , header = F , nrows = 2,
                       na.strings ="", stringsAsFactors= F)

View(rbchain2)

nc <- ncol(rtchain2)

rchain2_mat <- matrix(0, nrow = nc, ncol = 3)
for (i in 1:3){
  for (j in 1:nc){
    rchain2_mat[j,i] <- rtchain2[i,j]
  }
}
rchain2_mat[,3] <- 3/rchain2_mat[,3]
rchain2_mat <- cbind(rchain2_mat, t(rbchain2))

rchain2 <- mcmc(data = rchain2_mat)
varnames(rchain2) <- c('sigma2', 'tau2', 'phi', 'beta0', 'beta1')

par(mfrow = c(1,5))
plot(rchain2, density = FALSE)

# remove burn in 

rchain2_b <- mcmc(data = rchain2_mat[2000:5000,])
summary(rchain2_b)

rsig2_c2 <- summary(rchain2_b)$statistics[1]
rtau2_c2 <- summary(rchain2_b)$statistics[2]
rphi_c2 <- summary(rchain2_b)$statistics[3]
rbeta0_c2 <- summary(rchain2_b)$statistics[4]
rbeta1_c2 <- summary(rchain2_b)$statistics[5]


### guess sig2 2 (1), tau2 1 (1), phi 0.003 (0.001, 0.02), beta (10,0), tuning 0.00003, seed 1

rtchain3 <- read.table("C:/project-files/TEMP_DATA/chain3/temp-expchain3-theta", sep = "" , header = F , nrows = 3,
                     na.strings ="", stringsAsFactors= F)
rbchain3 <- read.table("C:/project-files/TEMP_data/chain3/temp-expchain3-beta", sep = "" , header = F , nrows = 2,
                       na.strings ="", stringsAsFactors= F)
View(rbchain3)

nc <- ncol(rtchain3)

rchain3_mat <- matrix(0, nrow = nc, ncol = 3)
for (i in 1:3){
  for (j in 1:nc){
    rchain3_mat[j,i] <- rtchain3[i,j]
  }
}
rchain3_mat[,3] <- 3/rchain3_mat[,3]
rchain3_mat <- cbind(rchain3_mat, t(rbchain3))

rchain3 <- mcmc(data = rchain3_mat)
varnames(rchain3) <- c('sigma2', 'tau2', 'phi', 'beta0', 'beta1')

plot(rchain3, density = FALSE)

# remove burn in 

rchain3_b <- mcmc(data = rchain3_mat[2000:5000,])
summary(rchain3_b)

rsig2_c3 <- summary(rchain3_b)$statistics[1]
rtau2_c3 <- summary(rchain3_b)$statistics[2]
rphi_c3 <- summary(rchain3_b)$statistics[3]
rbeta0_c3 <- summary(rchain3_b)$statistics[4]
rbeta1_c3 <- summary(rchain3_b)$statistics[5]


### AVERAGE PARAMETER ESTIMATES

rsig2_avg <- mean(c(rsig2_c1, rsig2_c2, rsig2_c3))
rtau2_avg <- mean(c(rtau2_c1, rtau2_c2, rtau2_c3))
rphi_avg  <- mean(c(rphi_c1, rphi_c2, rphi_c3))
rbeta0_avg <- mean(c(rbeta0_c1, rbeta0_c2, rbeta0_c3))
rbeta1_avg <- mean(c(rbeta1_c1, rbeta1_c2, rbeta1_c3))

c(rsig2_c1, rsig2_c2, rsig2_c3)
c(rtau2_c1, rtau2_c2, rtau2_c3)
c(rphi_c1, rphi_c2, rphi_c3)
c(rbeta0_c1, rbeta0_c2, rbeta0_c3)
c(rbeta1_c1, rbeta1_c2, rbeta1_c3)

rsig2_avg
rtau2_avg
rphi_avg
rbeta0_avg
rbeta1_avg
