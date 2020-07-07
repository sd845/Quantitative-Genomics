## reading in the data files
#reading in phenotype data or measurements of 5 gene expression levels
pheno_import = read.csv("./data_files/phenotypes.csv", header = TRUE)
#reading in genotype data file
geno_import = read.csv("./data_files/genotypes.csv", header = TRUE)
#reading in covariates to account for population structure and sex
covar_data = read.csv("./data_files/covars.csv", header = TRUE)
#reading in chromosome number, start and end positions
gene_info = read.csv("./data_files/gene_info.csv", header = TRUE)
#reading in SNP info with chromosome number and rsID
SNP_info = read.csv("./data_files/SNP_info.csv")

## examining phenotype and genotype data
#Number of samples
n = nrow(pheno_import)
cat("Sample Size: ", n, "\n")

#converting genotype and phenotype data to matrix
geno_data = as.matrix(geno_import[1:n,2:ncol(geno_import)])
pheno_data = as.matrix(pheno_import[1:n,2:ncol(pheno_import)])

#Number of phenotypes
n_pheno = ncol(pheno_data)
cat("Number of Phenotypes: ", n_pheno, "\n")

#number of SNPs
N = ncol(geno_data)
cat("Number of SNPs: ", N, "\n")

par(mfrow=c(2,3))

hist(pheno_data[,1], main=("ERAP2"),xlab="Gene Expression Measurements", ylab = "Frequency", col="grey", breaks = 30)
hist(pheno_data[,2], main=("PEX6"),xlab="Gene Expression Measurements", ylab = "Frequency", col="grey", breaks = 30)
hist(pheno_data[,3], main=("FAHD1"),xlab="Gene Expression Measurements", ylab = "Frequency", col="grey", breaks = 30)
hist(pheno_data[,4], main=("GFM1"),xlab="Gene Expression Measurements", ylab = "Frequency", col="grey", breaks = 30)
hist(pheno_data[,5], main=("MARCH7"),xlab="Gene Expression Measurements", ylab = "Frequency", col="grey", breaks = 30)

From the histograms, we find that each phenotype follows a normal distribution so we can proceed with using linear regression for our eQTL analysis. Linear regression assumes that the residuals follow a normal distribution.

## Building our Xa and Xd matrices

Xa = matrix(NA, nrow= n, ncol=N)
Xd = matrix(NA, nrow= n, ncol=N)

#given genotype coding: 0 = A1A1; 1 = A1A2; 2 = A2A2

#coding to Xa and Xd using Brute force
for (i in 1:n){
  for (j in 1:N){
    if (geno_data[i,j]==0){
      Xa[i,j] = -1
      Xd[i,j] = -1
    }
    if (geno_data[i,j]==1){
      Xa[i,j] = 0
      Xd[i,j] = 1
    }
    if (geno_data[i,j]==2){
      Xa[i,j] = 1
      Xd[i,j] = -1
    }
  }
}

## Principal component analysis
library(ggplot2)
#pca.result <- prcomp(Xa%*%t(Xa))
pca.result <- prcomp(Xa)

pca.result$sdev
(pca.result$sdev / sum(pca.result$sdev))*100
#summary(pca.result)

#constructing a pca_df
pcaDf <- data.frame(pc1=pca.result$x[,1], pc2=pca.result$x[,2])
ggplot(pcaDf,aes(pc1,pc2))+geom_point(aes(pc1,pc2))+xlab("Principal Component 1")+ylab("Principal Component 2") + theme_bw()
ggsave("PCA.png")


## calculating pvalues without a covariate
#importing MASS
library(MASS)

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

pval_mx = matrix(NA, nrow = N, ncol=n_pheno)

for (i in 1:n_pheno){
  for(j in 1:N){
    pval_mx[j,i] <- pval_calculator(pheno_data[,i], Xa[,j], Xd[,j])
  }
}

## creating Manhattan plots without covariate
library(ggplot2)

df_no_covar <- data.frame(genotypes = 1:N, pval_mx)

#using bonferroni correction
p_bf = 0.05/N

g1 <- ggplot(df_no_covar, aes(genotypes, -log10(df_no_covar[,2]))) + geom_point() + ggtitle("Manhattan Plot for Phenotype 1") + xlab("Position") + ylab("-log(p-value)") + geom_abline(intercept = -log10(p_bf), slope = 0, color="red") + ylim(0,75)
g2 <- ggplot(df_no_covar, aes(genotypes, -log10(df_no_covar[,3]))) + geom_point() + ggtitle("Manhattan Plot for Phenotype 2") + xlab("Position") + ylab("-log(p-value)") + geom_abline(intercept = -log10(p_bf), slope = 0, color="red")+ ylim(0,75)
g3 <- ggplot(df_no_covar, aes(genotypes, -log10(df_no_covar[,4]))) + geom_point() + ggtitle("Manhattan Plot for Phenotype 3") + xlab("Position") + ylab("-log(p-value)") + geom_abline(intercept = -log10(p_bf), slope = 0, color="red")+ ylim(0,75)
g4 <- ggplot(df_no_covar, aes(genotypes, -log10(df_no_covar[,5]))) + geom_point() + ggtitle("Manhattan Plot for Phenotype 4") + xlab("Position") + ylab("-log(p-value)") + geom_abline(intercept = -log10(p_bf), slope = 0, color="red")+ ylim(0,75)
g5 <- ggplot(df_no_covar, aes(genotypes, -log10(df_no_covar[,6]))) + geom_point() + ggtitle("Manhattan Plot for Phenotype 5") + xlab("Position") + ylab("-log(p-value)") + geom_abline(intercept = -log10(p_bf), slope = 0, color="red")+ ylim(0,75)
print(g1)
print(g2)
print(g3)
print(g4)
print(g5)

There is one peak in the Manhattan plots for the first three phenotypes which could be indicative of a cu]ausal polymorphism. For phenotypes 4 and 5, the Manhattan plots show no significant SNPs.

## To calculate the number of significant SNPs for each phenotype

N_sig_1 = 0
N_sig_2 = 0
N_sig_3 = 0
for (i in 1:N){
  if (pval_mx[i,1] < p_bf){
    N_sig_1 = N_sig_1 + 1
  }
  else if (pval_mx[i,2] < p_bf){
    N_sig_2 = N_sig_2 + 1
  }
  else if (pval_mx[i,3] < p_bf){
    N_sig_3 = N_sig_3 + 1
  }
}
cat("Number of Significant SNPs for Phenotype 1:", N_sig_1, "\n")
cat("Number of Significant SNPs for Phenotype 2:", N_sig_2, "\n")
cat("Number of Significant SNPs for Phenotype 3:", N_sig_3, "\n")

## Including covariates

We include three covariates for our analysis which account for population structure and sex. We follow the coding method given below:
  
  For covariate 1 (Xz1)
HG = 1
NA = 2

For covariate 2 (Xz2)
CEU (Utah residents with European ancestry) = 1
FIN (Finns) = 2
GBR (British) = 3
TSI (Toscani) = 4

For covariate 3 (Xz3)
MALE = 1
FEMALE = 2

Xz1 = rep(0, n)
Xz2 = rep(0, n)
Xz3 = rep(0, n)

for (i in 1:n){
  if (substr(pheno_data[i,1],1,2) =="HG"){
    Xz1[i] = -1}
  else if (substr(pheno_data[i,1],1,2) =="NA") {
    Xz1[i] = 1}
}

for (i in 1:n){
  if (covar_data[i,2] =="GBR"){
    Xz2[i] = 1}
  else if (covar_data[i,2] =="FIN") {
    Xz2[i] = 2}
  else if (covar_data[i,2] =="CEU") {
    Xz2[i] = 3}
  else if (covar_data[i,2] =="TSI") {
    Xz2[i] = 4}
}

for (i in 1:n){
  if (covar_data[i,3] =="MALE"){
    Xz3[i] = 1}
  else if (covar_data[i,3] =="FEMALE") {
    Xz3[i] = 2}
}

#calculating pvalue accounting for covariates

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
  df_E <- n_samples - 6
  
  numerator <- (SSE_theta_0-SSE_theta_1) / df_M
  denom <- SSE_theta_1 / df_E
  Fstatistic <-numerator / denom
  
  pval  = pf(Fstatistic, df_M, df_E,lower.tail = FALSE)
  return(pval)
}


Xz <- cbind(Xz1,Xz2,Xz3)
pval_mx_covar = matrix(NA, nrow = N, ncol=n_pheno)

for (i in 1:n_pheno){
  for(j in 1:N){
    pval_mx_covar[j,i] <- pval_calculator_w_covars(pheno_data[,i], Xa[,j], Xd[,j], Xz)
  }
}

## To calculate the number of significant SNPs for each phenotype- with covariate

N_sig_1 = 0
N_sig_2 = 0
N_sig_3 = 0
for (i in 1:N){
  if (pval_mx_covar[i,1] < p_bf){
    N_sig_1 = N_sig_1 + 1
  }
  else if (pval_mx_covar[i,2] < p_bf){
    N_sig_2 = N_sig_2 + 1
  }
  else if (pval_mx_covar[i,3] < p_bf){
    N_sig_3 = N_sig_3 + 1
  }
}
cat("Number of Significant SNPs for Phenotype 1:", N_sig_1, "\n")
cat("Number of Significant SNPs for Phenotype 2:", N_sig_2, "\n")
cat("Number of Significant SNPs for Phenotype 3:", N_sig_3, "\n")

## creating Manhattan Plots

chr = rep(NA,N)

for (i in 1:N){
  rs_id = colnames(geno_data)[i] 
  chr[i] = SNP_info[SNP_info$id == rs_id,]$chromosome
}

df_w_covar <- data.frame(genotypes = 1:N, pval_mx_covar,chr)

#cbp1 <- c("#FF0000", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#cbp1 <- c("coral", "blue", "pink", "aquamarine", "darkolivegreen", "orange", "maroon", "khaki", "red", "green", "coral", "blue", "pink", "aquamarine", "darkolivegreen", "orange", "maroon", "khaki", "red", "green", "black", "grey")
cbp1 = rep(c("grey", "black"), 11)


g1 <- ggplot(df_w_covar, aes(genotypes, -log10(df_w_covar[,2]), color = factor(chr))) + geom_point(show.legend = FALSE) + ggtitle("ERAP2: Manhattan Plot") + xlab("Position") + ylab("-log(p-value)") + geom_abline(intercept = -log10(p_bf), slope = 0, color="black") + ylim(0,75) + scale_colour_manual(values = cbp1) + theme_bw()

g2 <- ggplot(df_w_covar, aes(genotypes, -log10(df_w_covar[,3]), color = factor(chr))) + geom_point(show.legend = FALSE) + ggtitle("PEX6: Manhattan Plot") + xlab("Position") + ylab("-log(p-value)") + geom_abline(intercept = -log10(p_bf), slope = 0, color="black") + ylim(0,75) + scale_colour_manual(values = cbp1) + theme_bw()

g3 <- ggplot(df_w_covar, aes(genotypes, -log10(df_w_covar[,4]), color = factor(chr))) + geom_point(show.legend = FALSE) + ggtitle("FAHD1: Manhattan Plot") + xlab("Position") + ylab("-log(p-value)") + geom_abline(intercept = -log10(p_bf), slope = 0, color="black") + ylim(0,75) + scale_colour_manual(values = cbp1) + theme_bw()

g4 <- ggplot(df_w_covar, aes(genotypes, -log10(df_w_covar[,5]), color = factor(chr))) + geom_point(show.legend = FALSE) + ggtitle("GFM1: Manhattan Plot") + xlab("Position") + ylab("-log(p-value)") + geom_abline(intercept = -log10(p_bf), slope = 0, color="black") + ylim(0,75) + scale_colour_manual(values = cbp1) + theme_bw()

g5 <- ggplot(df_w_covar, aes(genotypes, -log10(df_w_covar[,6]), color = factor(chr))) + geom_point(show.legend = FALSE) + ggtitle("MARCH7: Manhattan Plot") + xlab("Position") + ylab("-log(p-value)") + geom_abline(intercept = -log10(p_bf), slope = 0, color="black") + ylim(0,75) + scale_color_manual(values = cbp1) + theme_bw()

print(g1)
print(g2)
print(g3)
print(g4)
print(g5)


## Multiple testing: QQ Plots
df_w_covar1<- data.frame(normalQuantiles = -log10(seq(0, 1, length.out = N)), -log10(pval_mx_covar))

q1 <- ggplot(df_w_covar1)+geom_point(aes(sort(normalQuantiles), sort(df_w_covar1[,2]))) + geom_abline(intercept = 0, slope = 1, color="red") + ggtitle("ERAP2: QQ Plot") + xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + theme_bw()

q2 <- ggplot(df_w_covar1)+geom_point(aes(sort(normalQuantiles), sort(df_w_covar1[,3]))) + geom_abline(intercept = 0, slope = 1, color="red") + ggtitle("PEX6: QQ Plot") + xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + theme_bw()

q3 <- ggplot(df_w_covar1)+geom_point(aes(sort(normalQuantiles), sort(df_w_covar1[,4]))) + geom_abline(intercept = 0, slope = 1, color="red") + ggtitle("FAHD1: QQ Plot") + xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + theme_bw()

q4 <- ggplot(df_w_covar1)+geom_point(aes(sort(normalQuantiles), sort(df_w_covar1[,5]))) + geom_abline(intercept = 0, slope = 1, color="red") + ggtitle("GFM1: QQ Plot") + xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + theme_bw()

q5 <- ggplot(df_w_covar1)+geom_point(aes(sort(normalQuantiles), sort(df_w_covar1[,6]))) + geom_abline(intercept = 0, slope = 1, color="red") + ggtitle("MARCH7: QQ Plot") + xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + theme_bw()

print(q1)
print(q2)
print(q3)
print(q4)
print(q5)

library(ggpubr)
f1 <- ggarrange(g1, q1, g2, q2, ncol = 2, nrow = 2)
ggsave("f1.png")
f2 <- ggarrange(g3, q3, g4, q4, ncol = 2, nrow = 2)
ggsave("f2.png")
f3 <- ggarrange(g5, q5, g5, q5, ncol = 2, nrow = 2)
ggsave("f3.png")

## Correcting p-values
N_tests = length(pval_mx_covar)
p_bf = 0.05/N_tests
N_sig_covar = 0
for (i in 1:length(pval_mx_covar)){
  if (pval_mx_covar[i] < p_bf){
    N_sig_covar = N_sig_covar + 1
  }
}
cat("Number of significant SNPs (with covariate):", N_sig_covar)

## Getting chromosome ID and position 
# Phenotype 1
i1 <- which(pval_mx_covar[,1] == min(pval_mx_covar[,1]))
rs_id1 <- colnames(geno_data)[i1] 
pos1 <- SNP_info[SNP_info$id %in% rs_id1,]
gene1 <- gene_info[gene_info$chromosome==pos1[1,1],]$symbol
gene1 <- cbind(pos1, symbol=gene1)

idx1 <- c(min(which(pval_mx_covar[,1] < p_bf)), max(which(pval_mx_covar[,1] < p_bf)))
rs1 <- colnames(geno_data)[idx1] 
position1 <- SNP_info[SNP_info$id %in% rs1,]$position
chr1 <- SNP_info[SNP_info$id %in% rs1,]$chromosome

cat("For Phenotype 1 \n")
cat("location start:", position1[1],"on chromosome", chr1[1], "\n")
cat("location end:", position1[2],"on chromosome", chr1[2], "\n\n")


# Phenotype 2
i2 <- which(pval_mx_covar[,2] == min(pval_mx_covar[,2]))
rs_id2 <- colnames(geno_data)[i2] 
pos2 <- SNP_info[SNP_info$id %in% rs_id2,]
gene2 <- gene_info[gene_info$chromosome==pos2[1,1],]$symbol
gene2 <- cbind(pos2, symbol=gene2)

idx2 <- c(min(which(pval_mx_covar[,2] < p_bf)), max(which(pval_mx_covar[,2] < p_bf)))
rs2 <- colnames(geno_data)[idx2] 
position2 <- SNP_info[SNP_info$id %in% rs2,]$position
chr2 <- SNP_info[SNP_info$id %in% rs2,]$chromosome

cat("For Phenotype 2 \n")
cat("location start:", position2[1],"on chromosome", chr2[1], "\n")
cat("location end:", position2[2],"on chromosome", chr2[2], "\n\n")

# Phenotype 3
i3 <- which(pval_mx_covar[,3] == min(pval_mx_covar[,3]))
rs_id3 <- colnames(geno_data)[i3] 
pos3 <- SNP_info[SNP_info$id %in% rs_id3,]
gene3 <- gene_info[gene_info$chromosome==pos3[1,1],]$symbol
gene3 <- cbind(pos3, symbol=gene3)

idx3 <- c(min(which(pval_mx_covar[,3] < p_bf)), max(which(pval_mx_covar[,3] < p_bf)))
rs3 <- colnames(geno_data)[idx3] 
position3 <- SNP_info[SNP_info$id %in% rs3,]$position
chr3 <- SNP_info[SNP_info$id %in% rs3,]$chromosome

cat("For Phenotype 3 \n")
cat("location start:", position3[1],"on chromosome", chr3[1], "\n")
cat("location end:", position3[2],"on chromosome", chr3[2], "\n\n")
