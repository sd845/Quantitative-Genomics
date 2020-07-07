pheno_data <- read.csv("pheno+pop.csv", header = FALSE)
n <- length(pheno_data[,2])
cat("Sample size (n): ", n)

hist(pheno_data[,2], main="histogram of HDL measurements",
     xlab="HDL measurements", ylab = "frequency")

geno_input <- read.csv("genotypes_V2.csv", header = FALSE)

#number of SNPs
N <- ncol(geno_input)
cat("Number of SNPs (N): ", N)

#importing MASS
library(MASS)

#converting to matrix
geno_data <- as.matrix(geno_input)

Xa <- matrix(NA, nrow= nrow(geno_data), ncol=ncol(geno_data))
Xd <- matrix(NA, nrow= nrow(geno_data), ncol=ncol(geno_data))

#given genotype coding: 0 = A1A1; 1 = A1A2; 2 = A2A2

#coding to Xa and Xd
for (i in 1:nrow(geno_data)){
  for (j in 1:ncol(geno_data)){
    if (geno_data[i,j]==0){
      Xa[i,j] <- -1
      Xd[i,j] <- -1
    }
    if (geno_data[i,j]==1){
      Xa[i,j] <- 0
      Xd[i,j] <- 1
    }
    if (geno_data[i,j]==2){
      Xa[i,j] <- 1
      Xd[i,j] <- -1
    }
  }
}

#calculating pvalues
pval_calculator <- function(pheno_input, xa_input, xd_input){
  
  n_samples <- length(xa_input)
  X_mx <- cbind(1,xa_input,xd_input)
  
  MLE_beta <- ginv(t(X_mx) %*% X_mx) %*% t(X_mx) %*% pheno_input
  y_hat <- X_mx %*% MLE_beta
  
  SSM <- sum((y_hat - mean(pheno_input))^2)
  SSE <- sum((pheno_input - y_hat)^2)
  
  df_M <- 2
  df_E <- n_samples - 3 
  
  MSM <- SSM / df_M
  MSE <- SSE / df_E
  
  Fstatistic <- MSM / MSE
  
  pval <- pf(Fstatistic, df_M, df_E,lower.tail = FALSE)
  return(pval)
}

pval_mx <- c()
for(i in 1:ncol(Xa)){
  pval_mx <- append(pval_mx, pval_calculator(pheno_data$V2, Xa[,i], Xd[,i]))
}

library(ggplot2)
plot_df <- data.frame(genotypes = 1:length(pval_mx), pval = pval_mx)
ggplot(plot_df, aes(genotypes, -log10(pval_mx))) + geom_point() + ggtitle("Manhattan Plot") + xlab("Position") + ylab("-log(p-value)")

N_tests <- length(pval_mx)
p_bf <- 0.05/N_tests
N_sig <- 0
for (i in 1:length(pval_mx)){
  if (pval_mx[i] < p_bf){
    N_sig <- N_sig + 1
  }
}
cat("Number of significant SNPs:", N_sig)

n1 <- 0
n2 <- 0
for (i in 1:nrow(pheno_data)){
  if (substr(pheno_data[i,1],1,2)=="HG"){
    n1 <- n1 + 1
  }
  if (substr(pheno_data[i,1],1,2)=="NA"){
    n2 <- n2 + 1
  }
}
cat("n1 =", n1, "\n")
cat("n2 =", n2, "\n")

Xz <- rep(0, nrow(pheno_data))

for (i in 1:nrow(pheno_data)){
  if (substr(pheno_data[i,1],1,2) =="HG"){
    Xz[i] <- -1}
  else if (substr(pheno_data[i,1],1,2) =="NA") {
    Xz[i] <- 1}
}
pval_calculator_w_covars <- function(pheno_input, xa_input, xd_input, xz_input){
  
  n_samples <-  length(pheno_input)
  
  X_mx <- cbind(rep(1,length(xa_input)), xa_input, xd_input, xz_input)
  MLE_beta <- ginv(t(X_mx) %*% X_mx) %*% t(X_mx) %*% pheno_input
  
  x_h0 <-   cbind(rep(1,length(xa_input)), xz_input)
  MLE_h0 <-  ginv(t(x_h0) %*% x_h0) %*% t(x_h0) %*% pheno_input
  
  y_hat_0 <-  x_h0 %*% MLE_h0
  y_hat_1 <- X_mx %*% MLE_beta
  
  SSE_theta_0 <-  sum((pheno_input-y_hat_0)^2) 
  SSE_theta_1 <- sum((pheno_input-y_hat_1)^2)
  
  
  df_M <- 2
  df_E <- n_samples - 6
  
  numerator <- (SSE_theta_0-SSE_theta_1) / df_M
  denom <- SSE_theta_1 / df_E
  Fstatistic <-numerator / denom
  
  pval  <- pf(Fstatistic, df_M, df_E,lower.tail = FALSE)
  return(pval)
}

pval_mx_covar <- rep(0,ncol(Xa))

for(i in 1:ncol(Xa)){
  pval_mx_covar[i] <- pval_calculator_w_covars(pheno_data$V2, Xa[,i], Xd[,i], Xz)
}

#make your Manhattan plot
plot_df <- data.frame(index = 1:length(pval_mx_covar), pval = pval_mx_covar)
ggplot(plot_df, aes(index, -log10(pval))) + geom_point() + ggtitle("Manhattan Plot with a Covariate") +xlab("Position")+ylab("-log(p-value)")

N_tests <- length(pval_mx_covar)
p_bf <- 0.05/N_tests
N_sig_covar <- 0
for (i in 1:length(pval_mx_covar)){
  if (pval_mx_covar[i] < p_bf){
    N_sig_covar <- N_sig_covar + 1
  }
}
cat("Number of significant SNPs (with covariate):", N_sig_covar)