#reading in Cholesterol levels measured for each of the n individuals
cholestrol_levels = read.csv("2020QG_cholesterol.csv", header = FALSE)

#sample size n
n = dim(cholestrol_levels)[1]
cat("The sample size n is", n)

#histogram of cholestrol levels
hist(cholestrol_levels[,1], main=("Cholesterol Levels"),xlab="Phenotype Measurements", ylab = "Frequency", col="grey", breaks = 30)


#importing genotype data
#SNP genotype data which contains data on N SNPs measured for each of the n individuals in the sample
#genotype of each SNP for an individual is coded as follows: `0' = homozygote, `1' = heterozygote, `2' = homozygote
geno_data = read.csv("2020QG_genotypes.csv", header = FALSE)
#number of genotypes N
N = dim(geno_data)[2]

#performing PCA
geno.pca <- prcomp(geno_data)

#potting pricnipal components
plot(geno.pca$x[,1], geno.pca$x[,2], main = "Genotype PC projections", xlab = "PC1", ylab = "PC2")


library(MASS)
library(ggplot2)

#calculating p-values

Xa = matrix(NA, nrow= n, ncol=N)
Xd = matrix(NA, nrow= n, ncol=N)

#given genotype coding: 0 = A1A1; 1 = A1A2; 2 = A2A2

#coding to Xa and Xd
Xa = geno_data - 1
Xd = 1 - 2*abs(Xa)

#calculating p-values
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

pval_mx = rep(0,N)

for (i in 1:N){
  pval_mx[i] <- pval_calculator(cholestrol_levels[,1], Xa[,i], Xd[,i])
}

N_sig = 0
for (i in 1:N){
  if (pval_mx[i]<0.05/N){
    N_sig = N_sig + 1
  }
}
cat("Number of Bonferroni corrected significant SNPs:",N_sig)

#making a Manhattan Plot for the p-values
#red line indicates Bonferroni corrected p-value cut-off
df <- data.frame(genotypes = 1:N, pval_mx)
g1 <- ggplot(df, aes(genotypes, -log10(pval_mx))) + geom_point() + ggtitle("Manhattan Plot") + xlab("Position") + ylab("-log(p-value)") + geom_abline(intercept = -log10(0.05/N), slope = 0, color="red") + theme_bw()
g1

#making a QQ plot
df_QQ <- data.frame(x = sort(-log10(seq(0, 1, length.out = N))), y = sort(-log10(pval_mx)))
q1 <- ggplot(df_QQ)+geom_point(aes(x,y)) + geom_abline(intercept = 0, slope = 1, color="red") + ggtitle("QQ Plot") + xlab("Theoretical Quantiles") + ylab("Observed Quantiles") + theme_bw()
q1

#accounting for covariates
library(MASS)
pval_calculator_w_covars <- function(pheno_input, xa_input, xd_input, xz_input){
  
  n_samples <-  length(pheno_input)
  
  X_mx <- cbind(rep(1,length(xa_input)), xa_input, xd_input, xz_input)
  MLE_beta <- ginv(t(X_mx) %*% X_mx) %*% t(X_mx) %*% pheno_input
  
  x_h0 =   cbind(rep(1,length(xa_input)), xz_input)
  MLE_h0 =  ginv(t(x_h0) %*% x_h0) %*% t(x_h0) %*% pheno_input
  
  y_hat_0 =  x_h0 %*% MLE_h0
  y_hat_1 = X_mx %*% MLE_beta
  
  SSE_theta_0 =  sum((pheno_input-y_hat_0)^2) 
  SSE_theta_1 = sum((pheno_input-y_hat_1)^2)
  
  
  df_M <- 2
  df_E <- n_samples - 5
  
  numerator <- (SSE_theta_0-SSE_theta_1) / df_M
  denom <- SSE_theta_1 / df_E
  Fstatistic <-numerator / denom
  
  pval  = pf(Fstatistic, df_M, df_E,lower.tail = FALSE)
  return(pval)
}


Xz <- cbind(geno.pca$x[,1], geno.pca$x[,2])

pval_mx_covar = rep(0,N)
for(i in 1:N){
  pval_mx_covar[i] <- pval_calculator_w_covars(cholestrol_levels[,1], Xa[,i], Xd[,i], Xz)
}

N_sig = 0
for (i in 1:N){
  if (pval_mx_covar[i]<0.05/N){
    N_sig = N_sig + 1
  }
}
cat("Number of Bonferroni corrected significant SNPs:",N_sig)

#making a Manhattan Plot for the p-values
df_covar <- data.frame(genotypes = 1:N, pval_mx_covar)
g2 <- ggplot(df_covar, aes(genotypes, -log10(pval_mx_covar))) + geom_point() + ggtitle("Manhattan Plot") + xlab("Position") + ylab("-log(p-value)") + geom_abline(intercept = -log10(0.05/N), slope = 0, color="red") + theme_bw()
g2

#making a QQ plot
df_QQ_covar <- data.frame(x = sort(-log10(seq(0, 1, length.out = N))), y = sort(-log10(pval_mx_covar)))
q2 <- ggplot(df_QQ_covar)+geom_point(aes(x,y)) + geom_abline(intercept = 0, slope = 1, color="red") + ggtitle("QQ Plot") + xlab("Theoretical Quantiles") + ylab("Observed Quantiles") + theme_bw()
q2

p_bf = 0.05/N
cat("p-value cutoff using Bonferroni correction:", p_bf, "\n")

N_sig_bf = 0
for (i in 1:N){
  if (pval_mx_covar[i]<p_bf){
    N_sig_bf = N_sig_bf + 1
  }
}
cat("Number of significant SNPs with Bonferroni correction:",N_sig_bf, "\n")

pval1 = pval_mx_covar[5000:15000]
pval2 = pval_mx_covar[15000:25000]
most_sig_pval1 = pval1[which(pval1==min(pval1))]
most_sig_pval2 = pval2[which(pval2==min(pval2))]
cat("most significant pvalues for peak 1:", most_sig_pval1[1], "at position", which(pval1==min(pval1))[1] + 5000, "\n")
cat("most significant pvalues for peak 2:", most_sig_pval2[1], "at position", which(pval2==min(pval2))[1] + 15000, "\n")

#reading in heart disease phenotypes
heart_disease = read.csv("2020QG_heartdisease.csv", header = FALSE)

sample_size = dim(heart_disease)[1]
cat("sample size, n:", sample_size)

#histogram of heart disease phenotype data
hist(heart_disease[,1], main=("Heart Disease Phenotype"),xlab="Phenotype Measurements", ylab = "Frequency", col="grey", breaks = 30)

library(MASS)
library(ggplot2)
library(ggthemes)

gamma_inv_calc <- function(X_mx, beta_t){
  #initialize gamma
  # K is the part which goes into the exponent
  K <- X_mx %*% beta_t
  gamma_inv <- exp(K)/(1+exp(K))
  return(gamma_inv)
}

W_calc <- function(gamma_inv){
  W <- diag(as.vector(gamma_inv * (1- gamma_inv)))
  return(W)
}

beta_update <- function(X_mx, W, Y, gamma_inv, beta){
  beta_up <- beta + ginv(t(X_mx)%*%W%*%X_mx)%*%t(X_mx)%*%(Y-gamma_inv)
  return(beta_up)
}

dev_calc <- function(Y, gamma_inv){
  deviance <- 2*( sum(Y[Y==1]*log(Y[Y==1]/gamma_inv[Y==1])) + sum((1-Y[Y==0])*log((1-Y[Y==0])/(1-gamma_inv[Y==0]))) )  
  return(deviance)
}

log_likelihood_calc <- function(Y, gamma_inv){
  log_likelihood <- sum(Y[Y==1]*log(gamma_inv[Y==1])) + sum((1-Y[Y==0])*log(1-gamma_inv[Y==0]))
  return (log_likelihood)
}


logistic_IRLS_keep_track_betas<- function(Xa,Xd,Y = Y, beta.initial.vec = c(0,0,0), d.stop.th = 1e-6, it.max = 100) {
  
  #Create the X matrix
  X_mx <- cbind(rep(1,length(Xa)), Xa, Xd)
  X_h0 <- cbind(rep(1,length(Xa)))
  
  #initializing beta vector at t=0
  beta_t <- c(0,0,0)
  beta_h0 <- c(0)
  
  # initialize deviance at d[t]
  dt <- 0
  dt_h0 <- 0
  
  beta_history = as.matrix(t(beta.initial.vec))
  
  #initialize gamma
  gamma_inv <- gamma_inv_calc(X_mx, beta_t)
  gamma_inv_h0 <- gamma_inv_calc(X_h0, beta_h0)
  
  for(i in 1:it.max) {
    dpt1 <- dt #store previous deviance
    dpt1_h0 <- dt_h0
    
    # create empty matrix W
    W <- W_calc(gamma_inv)
    W_h0 <- W_calc(gamma_inv_h0)
    
    beta_t <- beta_update(X_mx, W, Y, gamma_inv, beta_t)
    beta_h0 <- beta_update(X_h0, W_h0, Y, gamma_inv_h0, beta_h0)
    
    beta_history = rbind(beta_history, t(beta_t))
    
    #update gamma since it's a function of beta
    gamma_inv <- gamma_inv_calc(X_mx, beta_t)
    gamma_inv_h0 <- gamma_inv_calc(X_h0, beta_h0)
    
    #calculate new deviance
    dt <- dev_calc(Y, gamma_inv)
    dt_h0 <- dev_calc(Y, gamma_inv_h0)
    
    absD <- abs(dt - dpt1)
    absD_h0 <- abs(dt_h0 - dpt1_h0)
    
    if(absD < d.stop.th & absD_h0< d.stop.th) {
      #cat("Convergence at iteration:", i, "at threshold:", d.stop.th, "\n")
      logl = c(log_likelihood_calc(Y, gamma_inv_h0), log_likelihood_calc(Y, gamma_inv))
      return(list(beta_t,logl,beta_history))
    }	
  }
  cat("Convergence not reached after iteration:", i, "at threshold:", d.stop.th, "\n")
  return(list(beta_t= c(NA,NA,NA),logl=c(NA, NA), beta_history))
}

pval_mx_2 = rep(0,N)

for(j in seq(1, ncol(Xa))){
  logl_alt = logistic_IRLS_keep_track_betas(Xa[,j], Xd[,j], heart_disease[,1])[[2]][2]
  logl_null = logistic_IRLS_keep_track_betas(Xa[,j], Xd[,j], heart_disease[,1])[[2]][1]
  LRT = 2*logl_alt - 2*logl_null
  pval_mx_2[j] = pchisq(LRT, df=2)
}

N_sig_bf_1 = 0
for (i in 1:N){
  if (pval_mx_2[i]<p_bf){
    N_sig_bf_1 = N_sig_bf_1 + 1
  }
}
cat("Number of significant SNPs with Bonferroni correction:",N_sig_bf_1, "\n")

#making a Manhattan Plot for the p-values
df_2 <- data.frame(genotypes = 1:N, pval_mx_2)
g3 <- ggplot(df_2, aes(genotypes, -log10(pval_mx_2))) + geom_point() + ggtitle("Manhattan Plot") + xlab("Position") + ylab("-log(p-value)") + geom_abline(intercept = -log10(0.05/N), slope = 0, color="red")
g3

#making a QQ plot
df_QQ_2 <- data.frame(x = sort(-log10(seq(0, 1, length.out = N))), y = sort(-log10(pval_mx_2)))
q3 <- ggplot(df_QQ_2)+geom_point(aes(x,y)) + geom_abline(intercept = 0, slope = 1, color="red") + ggtitle("QQ Plot") + xlab("Theoretical Quantiles") + ylab("Observed Quantiles") + theme_bw()
q3

logistic_IRLS_keep_track_betas_w_covars<- function(Xa,Xd, Xz, Y = Y, beta.initial.vec = c(0,0,0), d.stop.th = 1e-6, it.max = 100) {
  
  #Create the X matrix
  X_mx <- cbind(rep(1,length(Xa)), Xa, Xd, Xz)
  X_h0 <- cbind(rep(1,length(Xa)), Xz)
  
  #initializing beta vector at t=0
  beta_t <- rep(0,ncol(Xz)+3)
  beta_h0 <- rep(0,ncol(Xz)+1)
  
  # initialize deviance at d[t]
  dt <- 0
  dt_h0 <- 0
  
  beta_history = as.matrix(t(rep(0,ncol(Xz)+3)))
  
  #initialize gamma
  gamma_inv <- gamma_inv_calc(X_mx, beta_t)
  gamma_inv_h0 <- gamma_inv_calc(X_h0, beta_h0)
  
  for(i in 1:it.max) {
    dpt1 <- dt #store previous deviance
    dpt1_h0 <- dt_h0
    
    # create empty matrix W
    W <- W_calc(gamma_inv)
    W_h0 <- W_calc(gamma_inv_h0)
    
    beta_t <- beta_update(X_mx, W, Y, gamma_inv, beta_t)
    beta_h0 <- beta_update(X_h0, W_h0, Y, gamma_inv_h0, beta_h0)
    
    beta_history = rbind(beta_history, t(beta_t))
    
    #update gamma since it's a function of beta
    gamma_inv <- gamma_inv_calc(X_mx, beta_t)
    gamma_inv_h0 <- gamma_inv_calc(X_h0, beta_h0)
    
    #calculate new deviance
    dt <- dev_calc(Y, gamma_inv)
    dt_h0 <- dev_calc(Y, gamma_inv_h0)
    
    absD <- abs(dt - dpt1)
    absD_h0 <- abs(dt_h0 - dpt1_h0)
    
    if(absD < d.stop.th & absD_h0< d.stop.th) {
      #cat("Convergence at iteration:", i, "at threshold:", d.stop.th, "\n")
      logl = c(log_likelihood_calc(Y, gamma_inv_h0), log_likelihood_calc(Y, gamma_inv))
      return(list(beta_t,logl,beta_history))
    }	
  }
  cat("Convergence not reached after iteration:", i, "at threshold:", d.stop.th, "\n")
  return(list(beta_t= c(NA,NA,NA),logl=c(NA, NA), beta_history))
}

pval_mx_covar_2 = rep(0,N)

for(j in seq(1, ncol(Xa))){
  logl_alt = logistic_IRLS_keep_track_betas_w_covars(Xa[,j], Xd[,j], Xz, heart_disease[,1])[[2]][2]
  logl_null = logistic_IRLS_keep_track_betas_w_covars(Xa[,j], Xd[,j], Xz, heart_disease[,1])[[2]][1]
  LRT = 2*logl_alt - 2*logl_null
  pval_mx_covar_2[j] = pchisq(LRT, df=2)
}

#making a Manhattan Plot for the p-values
df_covar_2 <- data.frame(genotypes = 1:N, pval_mx_covar_2)
g4 <- ggplot(df_covar_2, aes(genotypes, -log10(pval_mx_covar_2))) + geom_point() + ggtitle("Manhattan Plot") + xlab("Position") + ylab("-log(p-value)") + geom_abline(intercept = -log10(0.05/N), slope = 0, color="red")
g4

#making a QQ plot
df_QQ_covar_2 <- data.frame(x = sort(-log10(seq(0, 1, length.out = N))), y = sort(-log10(pval_mx_covar_2)))
q4 <- ggplot(df_QQ_covar_2)+geom_point(aes(x,y)) + geom_abline(intercept = 0, slope = 1, color="red") + ggtitle("QQ Plot") + xlab("Theoretical Quantiles") + ylab("Observed Quantiles") + theme_bw()
q4

p_bf = 0.05/N
cat("p-value cutoff using Bonferroni correction:", p_bf, "\n")

N_sig_bf_2 = 0
for (i in 1:N){
  if (pval_mx_covar_2[i]<p_bf){
    N_sig_bf_2 = N_sig_bf_2 + 1
  }
}
cat("Number of significant SNPs with Bonferroni correction:",N_sig_bf_2, "\n")

pval1 = pval_mx_covar_2[0:10000]
pval2 = pval_mx_covar_2[30000:length(pval_mx_covar_2)]
most_sig_pval1 = pval1[which(pval1==min(pval1))]
most_sig_pval2 = pval2[which(pval2==min(pval2))]
cat("most significant pvalues for peak 1:", most_sig_pval1[1], "at position", which(pval1==min(pval1))[1] + 0, "\n")
cat("most significant pvalues for peak 2:", most_sig_pval2[1], "at position", which(pval2==min(pval2))[1] + 30000, "\n")