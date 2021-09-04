library(coda)

### guess sig2 1, tau2 0.1, phi 5, beta 1, tuning 0.0015, seed 1

tchain1 <- read.table("C:/project-files/SYNTH_DATA/exporesults/chain1/synth-chain1-theta", sep = "" , header = F , nrows = 3,
                         na.strings ="", stringsAsFactors= F)
bchain1 <- read.table("C:/project-files/SYNTH_DATA/exporesults/chain1/synth-chain1-beta", sep = "" , header = F , nrows = 1,
                      na.strings ="", stringsAsFactors= F)
View(chain1)

nc <- ncol(tchain1)

chain1_mat <- matrix(0, nrow = nc, ncol = 3)
for (i in 1:3){
  for (j in 1:nc){
    chain1_mat[j,i] <- tchain1[i,j]
  }
}

chain1_mat <- cbind(chain1_mat, t(bchain1))

chain1 <- mcmc(data = chain1_mat)
varnames(chain1) <- c('sigma2', 'tau2', 'phi', 'beta')

plot(chain1, density = FALSE)

# remove burn in 

chain1_b <- mcmc(data = chain1_mat[2000:5000,])
summary(chain1_b)

sig2_c1 <- summary(chain1_b)$statistics[1]
tau2_c1 <- summary(chain1_b)$statistics[2]
phi_c1 <- summary(chain1_b)$statistics[3]
beta_c1 <- summary(chain1_b)$statistics[4]

### chain 2 guess sig2 1.1, tau2 0.15, phi 6, beta 1.5, tuning 0.0015, seed 3

tchain2 <- read.table("C:/project-files/SYNTH_DATA/exporesults/chain2/synth-chain2-theta", sep = "" , header = F , nrows = 3,
                     na.strings ="", stringsAsFactors= F)
bchain2 <- read.table("C:/project-files/SYNTH_DATA/exporesults/chain2/synth-chain2-beta", sep = "" , header = F , nrows = 1,
                      na.strings ="", stringsAsFactors= F)
View(chain1)

nc <- ncol(tchain2)

chain2_mat <- matrix(0, nrow = nc, ncol = 3)
for (i in 1:3){
  for (j in 1:nc){
    chain2_mat[j,i] <- tchain2[i,j]
  }
}

chain2_mat <- cbind(chain2_mat, t(bchain2))

chain2 <- mcmc(data = chain2_mat)
varnames(chain2) <- c('sigma2', 'tau2', 'phi', 'beta')

plot(chain2, density = FALSE)

# remove burn in 

chain2_b <- mcmc(data = chain2_mat[2000:5000,])
summary(chain2_b)

sig2_c2 <- summary(chain2_b)$statistics[1]
tau2_c2 <- summary(chain2_b)$statistics[2]
phi_c2 <- summary(chain2_b)$statistics[3]
beta_c2 <- summary(chain2_b)$statistics[4]



### chain 3 guess sig2 0.95, tau2 0.095, phi 5.5, tuning 0.0015, beta 0.9

tchain3 <- read.table("C:/project-files/SYNTH_DATA/exporesults/chain3/synth-chain3-theta", sep = "" , header = F , nrows = 3,
                     na.strings ="", stringsAsFactors= F)
bchain3 <- read.table("C:/project-files/SYNTH_DATA/exporesults/chain3/synth-chain3-beta", sep = "" , header = F , nrows = 1,
                      na.strings ="", stringsAsFactors= F)

nc <- ncol(tchain3)

chain3_mat <- matrix(0, nrow = nc, ncol = 3)
for (i in 1:3){
  for (j in 1:nc){
    chain3_mat[j,i] <- tchain3[i,j]
  }
}

chain3_mat <- cbind(chain3_mat, t(bchain3))

chain3 <- mcmc(data = chain3_mat)
varnames(chain3) <- c('sigma2', 'tau2', 'phi', 'beta')

plot(chain3, density = FALSE)

# remove burn in 

chain3_b <- mcmc(data = chain3_mat[2000:5000,])
summary(chain3_b)

sig2_c3 <- summary(chain3_b)$statistics[1]
tau2_c3 <- summary(chain3_b)$statistics[2]
phi_c3 <- summary(chain3_b)$statistics[3]
beta_c3 <- summary(chain3_b)$statistics[4]

### AVERAGE PARAMETER ESTIMATES

sig2_avg <- mean(c(sig2_c1, sig2_c2, sig2_c3))
tau2_avg <- mean(c(tau2_c1, tau2_c2, tau2_c3))
phi_avg  <- mean(c(phi_c1, phi_c2, phi_c3))
beta_avg <- mean(c(beta_c1, beta_c2, beta_c3))

sig2_avg
tau2_avg
phi_avg
beta_avg
