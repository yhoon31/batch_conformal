#BH procedure
BH <- function(p_values, alpha) 
{
  ##inputs
  #p_values : vector of p values
  #alpha : target level
  
  ##output
  #rejection threshold
  
  p_sorted <- sort(p_values)
  m <- length(p_values)
  if(min(p_sorted*m/(1:m)) <= alpha)
  {
    threshold <- p_sorted[max(which(p_sorted <= (1:m) / m * alpha))]
  }
  else
  {
    threshold <- 0
  }
  return(threshold)
  
}

##group size generation
group_size_unif <- function(K,m_l,m_u)
{##inputs
  #K : number of groups
  #m_l, m_u : lower and upper bound of group sizes
  return(sample(m_l:m_u, K,replace=TRUE))
}

##batch conformal p-value
batch_conformal <- function(S_cal, S_test, eta)
{
  n <- length(S_cal)
  m <- length(S_test)
  probs <- sapply(1:(n+1), function(i) choose(i+eta-2,eta-1)*choose(n+m-i-eta+1,m-eta)/choose(n+m,m))
  t <- sum(S_cal < sort(S_test)[eta])
  return(sum(probs[(t+1):(n+1)]))
}

##two sample z test
two_sample <- function(v, w, sigma_v, sigma_w) {
  n_v <- length(v)
  n_w <- length(w)
  mean_v <- mean(v)
  mean_w <- mean(w)
  z <- (mean_v - mean_w) / sqrt((sigma_v^2 / n_v) + (sigma_w^2 / n_w))
  p_value <- pnorm(z)
  return(p_value)
}

##two sample t test
two_sample_t <- function(v, w)
{
  t_result <- t.test(v, w, alternative = "less")
  return(t_result$p.value)
}

##---------------------------------------------------------------------------------------------


#simulation for one-dimensional setting
sim_1 <- function(mu_0,mu_1,m_l,m_u,null,K,alpha_set)
{
  m_vec <- group_size_unif(K,m_l,m_u)
  X_cal <- rnorm(n,mu_0,sigma)
  
  X_test_list <- lapply(1:K, function(k) null[k]*rnorm(m_vec[k],mu_0,sigma) + (1-null[k])*rnorm(m_vec[k],mu_1,sigma))
  
  #batch conformal p value
  p_values <- sapply(X_test_list, function(v) batch_conformal(X_cal,v,max(round(0.5*length(v)),1)))
  p_threshold <- sapply(alpha_set,function(alpha) BH(p_values,alpha))
  rej <- lapply(p_threshold, function(t) p_values <= t)
  V <- sapply(rej, function(v) sum(null==1 & v==TRUE))
  W <- sapply(rej, function(v) sum(null==0 & v==TRUE))
  R <- sapply(rej, function(v) sum(v==TRUE))
  fdp <- V/max(R,1)
  power <- W/max(R,1)

  return(rbind(fdp,power))
}



n <- 100 #calibration size

mu_0 <- 0 #mean under H0
sigma <- 3 #variance
alpha_set <- seq(0.05,0.3,by=0.05) #target levels

m_l <- 30
m_u <- 50
#m_l <- 1
#m_u <- 1

p_null <- 0.5 #proportion of true nulls


mu_1 <- 3 #mean under H1 / signal strength
K <- 50 #number of groups
num_null <- round(K*p_null) #number of nulls
null <- c(rep(0,num_null), rep(1,K-num_null)) #true null indicators

tn <- 500
fdp <- matrix(0,nrow=tn,ncol=length(alpha_set))
tpp <- matrix(0,nrow=tn,ncol=length(alpha_set))
fdp_z <- matrix(0,nrow=tn,ncol=length(alpha_set))
tpp_z <- matrix(0,nrow=tn,ncol=length(alpha_set))


pb <- txtProgressBar(min = 0, max = tn, style = 3)
for(i in 1:tn)
{
  m_vec <- group_size_unif(K,m_l,m_u)
  X_cal <- rnorm(n,mu_0,sigma)
  
  X_test_list <- lapply(1:K, function(k) null[k]*rnorm(m_vec[k],mu_0,sigma) + (1-null[k])*rnorm(m_vec[k],mu_1,sigma))
  
  #batch conformal p value
  p_values <- sapply(X_test_list, function(v) batch_conformal(X_cal,v,max(round(0.5*length(v)),1)))
  p_threshold <- sapply(alpha_set,function(alpha) BH(p_values,alpha))
  rej <- lapply(p_threshold, function(t) p_values <= t)
  V <- sapply(rej, function(v) sum(null==1 & v==TRUE))
  W <- sapply(rej, function(v) sum(null==0 & v==TRUE))
  R <- sapply(rej, function(v) sum(v==TRUE))
  fdp[i,] <- V/max(R,1)
  tpp[i,] <- W/max(R,1)
  
  #two sample z test p value
  p_values_z <- sapply(X_test_list, function(v) two_sample(X_cal,v,sigma,sigma))
  p_threshold_z <- sapply(alpha_set,function(alpha) BH(p_values_z,alpha))
  rej_z <- lapply(p_threshold_z, function(t) p_values_z <= t)
  V_z <- sapply(rej_z, function(v) sum(null==1 & v==TRUE))
  W_z <- sapply(rej_z, function(v) sum(null==0 & v==TRUE))
  R_z <- sapply(rej_z, function(v) sum(v==TRUE))
  fdp_z[i,] <- V_z/max(R_z,1)
  tpp_z[i,] <- W_z/max(R_z,1)

  setTxtProgressBar(pb, i) 
}


fdr <- colMeans(fdp)
fdr_se <- apply(fdp,2,sd)/sqrt(tn)
power <- colMeans(tpp)
power_se <- apply(tpp,2,sd)/sqrt(tn)
fdr_z <- colMeans(fdp_z)
fdr_z_se <- apply(fdp_z,2,sd)/sqrt(tn)
power_z <- colMeans(tpp_z)
power_z_se <- apply(tpp_z,2,sd)/sqrt(tn)

res <- cbind(fdr,fdr_se,power,power_se,fdr_z,fdr_z_se,power_z,power_z_se)


colz <- c("#1976D2","#D32F2F")
par(mfrow=c(1,2), cex.lab = 2, cex.main=2, cex.axis=2, mar=c(2,2,2,3), xpd=NA, oma = c(3, 4, 1, 0.5),font.main=1)

matplot(alpha_set,res[,c(1,5)],ylim=c(0,0.3),type="l",lty=1,lwd=2,xlab=expression(alpha),ylab="FDR",main="",col=colz)
arrows(alpha_set, res[,1]-res[,2], alpha_set, res[,1]+res[,2], length=0.01, angle=90, code=3,lwd=2,col=colz[1])
arrows(alpha_set, res[,5]-res[,6], alpha_set, res[,5]+res[,6], length=0.01, angle=90, code=3,lwd=2,col=colz[2])
segments(0.05,0.05*p_null,0.3,0.3*p_null,lty=2,lwd=2)

matplot(alpha_set,res[,c(3,7)],ylim=c(0,1),type="l",lty=1,lwd=2,xlab=expression(alpha),ylab="Power",main="",col=colz)
arrows(alpha_set, res[,3]-res[,4], alpha_set, res[,3]+res[,4], length=0.01, angle=90, code=3,lwd=2,col=colz[1])
arrows(alpha_set, res[,7]-res[,8], alpha_set, res[,7]+res[,8], length=0.01, angle=90, code=3,lwd=2,col=colz[2])


##--------------------------------------------------------------------------------------------------------

##simulations under different settings
mu_1_set <- c(1,2,3)
K_set <- c(20,50,200)
p_null <- 0.7


fdr <- matrix(0,nrow=9,ncol=length(alpha_set))
power <- matrix(0,nrow=9,ncol=length(alpha_set))
fdr_se <- matrix(0,nrow=9,ncol=length(alpha_set))
power_se <- matrix(0,nrow=9,ncol=length(alpha_set))

tn <- 50

t <- 1
for(K in K_set)
{
  for(mu_1 in mu_1_set)
  {
    num_null <- round(K*p_null)
    null <- c(rep(0,num_null), rep(1,K-num_null))
    fdp <- matrix(0,nrow=tn,ncol=length(alpha_set))
    tpp <- matrix(0,nrow=tn,ncol=length(alpha_set))
    
    
    pb <- txtProgressBar(min = 0, max = tn, style = 3)
    for(i in 1:tn)
    {
      res <- sim_1(mu_0,mu_1,m_l,m_u,null,K,alpha_set)
      fdp[i,] <- res[1,]
      tpp[i,] <- res[2,]
      setTxtProgressBar(pb, i) 
    }

    fdr[t,] <- colMeans(fdp)
    fdr_se[t,] <- apply(fdp,2,sd)/sqrt(tn)
    power[t,] <- colMeans(tpp)
    power_se[t,] <- apply(tpp,2,sd)/sqrt(tn)
    t <- t+1
  }
}

col1=c("#A0C8F0", "#5A8BBA", "#2C5572")

par(mfrow=c(2,3), cex.lab = 2, cex.main=2, cex.axis=2, mar=c(3,3,3,3), xpd=NA, oma = c(5.5, 3.5, 1, 0.5),font.main=1)

matplot(alpha_set,t(fdr[1:3,]),ylim=c(0,0.3),type="l",lty=1,lwd=2,xlab="",ylab="FDR",main="K = 20",col=col1)
arrows(alpha_set, fdr[1,]-fdr_se[1,], alpha_set, fdr[1,]+fdr_se[1,], length=0.02, angle=90, code=3,lwd=2,col=col1[1])
arrows(alpha_set, fdr[2,]-fdr_se[2,], alpha_set, fdr[2,]+fdr_se[2,], length=0.02, angle=90, code=3,lwd=2,col=col1[2])
arrows(alpha_set, fdr[3,]-fdr_se[3,], alpha_set, fdr[3,]+fdr_se[3,], length=0.02, angle=90, code=3,lwd=2,col=col1[3])
segments(0.05,0.05*p_null,0.3,0.3*p_null,lty=2,lwd=2)

matplot(alpha_set,t(fdr[4:6,]),ylim=c(0,0.3),type="l",lty=1,lwd=2,xlab="",ylab="",main="K = 50",col=col1)
arrows(alpha_set, fdr[4,]-fdr_se[4,], alpha_set, fdr[4,]+fdr_se[4,], length=0.02, angle=90, code=3,lwd=2,col=col1[1])
arrows(alpha_set, fdr[5,]-fdr_se[5,], alpha_set, fdr[5,]+fdr_se[5,], length=0.02, angle=90, code=3,lwd=2,col=col1[2])
arrows(alpha_set, fdr[6,]-fdr_se[6,], alpha_set, fdr[6,]+fdr_se[6,], length=0.02, angle=90, code=3,lwd=2,col=col1[3])
segments(0.05,0.05*p_null,0.3,0.3*p_null,lty=2,lwd=2)

matplot(alpha_set,t(fdr[7:9,]),ylim=c(0,0.3),type="l",lty=1,lwd=2,xlab="",ylab="",main="K = 200",col=col1)
arrows(alpha_set, fdr[7,]-fdr_se[7,], alpha_set, fdr[7,]+fdr_se[7,], length=0.02, angle=90, code=3,lwd=2,col=col1[1])
arrows(alpha_set, fdr[8,]-fdr_se[8,], alpha_set, fdr[8,]+fdr_se[8,], length=0.02, angle=90, code=3,lwd=2,col=col1[2])
arrows(alpha_set, fdr[9,]-fdr_se[9,], alpha_set, fdr[9,]+fdr_se[9,], length=0.02, angle=90, code=3,lwd=2,col=col1[3])
segments(0.05,0.05*p_null,0.3,0.3*p_null,lty=2,lwd=2)

matplot(alpha_set,t(power[1:3,]),ylim=c(0,1),type="l",lty=1,lwd=2,xlab=expression(alpha),ylab="Power",main="",col=col1)
arrows(alpha_set, power[1,]-power_se[1,], alpha_set, power[1,]+power_se[1,], length=0.02, angle=90, code=3,lwd=2,col=col1[1])
arrows(alpha_set, power[2,]-power_se[2,], alpha_set, power[2,]+power_se[2,], length=0.02, angle=90, code=3,lwd=2,col=col1[2])
arrows(alpha_set, power[3,]-power_se[3,], alpha_set, power[3,]+power_se[3,], length=0.02, angle=90, code=3,lwd=2,col=col1[3])

matplot(alpha_set,t(power[4:6,]),ylim=c(0,1),type="l",lty=1,lwd=2,xlab=expression(alpha),ylab="",main="",col=col1)
arrows(alpha_set, power[4,]-power_se[4,], alpha_set, power[4,]+power_se[4,], length=0.02, angle=90, code=3,lwd=2,col=col1[1])
arrows(alpha_set, power[5,]-power_se[5,], alpha_set, power[5,]+power_se[5,], length=0.02, angle=90, code=3,lwd=2,col=col1[2])
arrows(alpha_set, power[6,]-power_se[6,], alpha_set, power[6,]+power_se[6,], length=0.02, angle=90, code=3,lwd=2,col=col1[3])

matplot(alpha_set,t(power[7:9,]),ylim=c(0,1),type="l",lty=1,lwd=2,xlab=expression(alpha),ylab="",main="",col=col1)
arrows(alpha_set, power[7,]-power_se[7,], alpha_set, power[7,]+power_se[7,], length=0.02, angle=90, code=3,lwd=2,col=col1[1])
arrows(alpha_set, power[8,]-power_se[8,], alpha_set, power[8,]+power_se[8,], length=0.02, angle=90, code=3,lwd=2,col=col1[2])
arrows(alpha_set, power[9,]-power_se[9,], alpha_set, power[9,]+power_se[9,], length=0.02, angle=90, code=3,lwd=2,col=col1[3])

legend1=expression(paste(delta," = 1     "))
legend2=expression(paste(delta," = 2     "))
legend3=expression(paste(delta," = 3     "))
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(legend1, legend2, legend3), col = col1,lty=c(1,1,1),lwd=2,box.lty=0,cex=1.7,bty="n",horiz = T,xpd = TRUE)




##---------------------------------------------------------------------------------------------

#simulation under non-normal distribution

library(sn)
library(VGAM)

sim_1 <- function(mu_0,mu_1,m_l,m_u,null,K,alpha_set)
{
  m_vec <- group_size_unif(K,m_l,m_u)
  X_cal <- rnorm(n,mu_0,sigma)+runif(n,-1,1)
  
  X_test_list <- lapply(1:K, function(k) null[k]*rnorm(m_vec[k],mu_0,sigma) + (1-null[k])*rnorm(m_vec[k],mu_1,sigma) +runif(m_vec[k],-1,1))
  
  #batch conformal p value
  p_values <- sapply(X_test_list, function(v) batch_conformal(X_cal,v,round(0.5*length(v))))
  p_threshold <- sapply(alpha_set,function(alpha) BH(p_values,alpha))
  rej <- lapply(p_threshold, function(t) p_values <= t)
  V <- sapply(rej, function(v) sum(null==1 & v==TRUE))
  W <- sapply(rej, function(v) sum(null==0 & v==TRUE))
  R <- sapply(rej, function(v) sum(v==TRUE))
  fdp <- V/max(R,1)
  power <- W/max(R,1)
  
  return(rbind(fdp,power))
}


set.seed(1234)
n <- 100 #calibration size

mu_0 <- 0 #mean under H0
sigma <- 3 #variance
alpha_set <- seq(0.05,0.3,by=0.05) #target levels

m_l <- 30
m_u <- 50

p_null <- 0.7 #proportion of true nulls


mu_1 <- 1 #mean under H1 / signal strength
K <- 50 #number of groups
num_null <- round(K*p_null) #number of nulls
null <- c(rep(1,num_null), rep(0,K-num_null)) #true null indicators

tn <- 500
fdp <- matrix(0,nrow=tn,ncol=length(alpha_set))
tpp <- matrix(0,nrow=tn,ncol=length(alpha_set))
fdp_z <- matrix(0,nrow=tn,ncol=length(alpha_set))
tpp_z <- matrix(0,nrow=tn,ncol=length(alpha_set))


pb <- txtProgressBar(min = 0, max = tn, style = 3)
for(i in 1:tn)
{
  m_vec <- group_size_unif(K,m_l,m_u)
  X_cal <- mu_0+ifelse(runif(n) > 0.5, rcauchy(n, location = 0, scale = 1), runif(n, min = -1, max = 1))
  X_test_list <- lapply(1:K, function(k) null[k]*(mu_0+ifelse(runif(m_vec[k]) > 0.5, rcauchy(m_vec[k], location = 0, scale = 1), runif(m_vec[k], min = -1, max = 1))) + (1-null[k])*(mu_1+ifelse(runif(m_vec[k]) > 0.5, rcauchy(m_vec[k], location = 0, scale = 1), runif(m_vec[k], min = -1, max = 1))))
  
  #batch conformal p value
  p_values <- sapply(X_test_list, function(v) batch_conformal(X_cal,v,round(0.5*length(v))))
  p_threshold <- sapply(alpha_set,function(alpha) BH(p_values,alpha))
  rej <- lapply(p_threshold, function(t) p_values <= t)
  V <- sapply(rej, function(v) sum(null==1 & v==TRUE))
  W <- sapply(rej, function(v) sum(null==0 & v==TRUE))
  R <- sapply(rej, function(v) sum(v==TRUE))
  fdp[i,] <- V/max(R,1)
  tpp[i,] <- W/max(R,1)
  
  #two sample t test p value
  p_values_z <- sapply(X_test_list, function(v) two_sample_t(X_cal,v))
  p_threshold_z <- sapply(alpha_set,function(alpha) BH(p_values_z,alpha))
  rej_z <- lapply(p_threshold_z, function(t) p_values_z <= t)
  V_z <- sapply(rej_z, function(v) sum(null==1 & v==TRUE))
  W_z <- sapply(rej_z, function(v) sum(null==0 & v==TRUE))
  R_z <- sapply(rej_z, function(v) sum(v==TRUE))
  fdp_z[i,] <- V_z/max(R_z,1)
  tpp_z[i,] <- W_z/max(R_z,1)
  
  setTxtProgressBar(pb, i) 
}

fdr <- colMeans(fdp)
fdr_se <- apply(fdp,2,sd)/sqrt(tn)
power <- colMeans(tpp)
power_se <- apply(tpp,2,sd)/sqrt(tn)
fdr_z <- colMeans(fdp_z)
fdr_z_se <- apply(fdp_z,2,sd)/sqrt(tn)
power_z <- colMeans(tpp_z)
power_z_se <- apply(tpp_z,2,sd)/sqrt(tn)

res <- cbind(fdr,fdr_se,power,power_se,fdr_z,fdr_z_se,power_z,power_z_se)

colz <- c("#1976D2","#D32F2F")
par(mfrow=c(1,2), cex.lab = 2, cex.main=2, cex.axis=2, mar=c(2,2,2,3), xpd=NA, oma = c(3, 4, 1, 0.5),font.main=1)

matplot(alpha_set,res[,c(1,5)],ylim=c(0,0.3),type="l",lty=1,lwd=2,xlab=expression(alpha),ylab="FDR",main="",col=colz)
arrows(alpha_set, res[,1]-res[,2], alpha_set, res[,1]+res[,2], length=0.01, angle=90, code=3,lwd=2,col=colz[1])
arrows(alpha_set, res[,5]-res[,6], alpha_set, res[,5]+res[,6], length=0.01, angle=90, code=3,lwd=2,col=colz[2])
segments(0.05,0.05*p_null,0.3,0.3*p_null,lty=2,lwd=2)

matplot(alpha_set,res[,c(3,7)],ylim=c(0,1),type="l",lty=1,lwd=2,xlab=expression(alpha),ylab="Power",main="",col=colz)
arrows(alpha_set, res[,3]-res[,4], alpha_set, res[,3]+res[,4], length=0.01, angle=90, code=3,lwd=2,col=colz[1])
arrows(alpha_set, res[,7]-res[,8], alpha_set, res[,7]+res[,8], length=0.01, angle=90, code=3,lwd=2,col=colz[2])


##----------------------------------------------------------------------


##simulation with multivariate outcome
library(MASS)
set.seed(1234)
n_train <- 100
n_cal <- 50
p <- 10

#parameters
set.seed(123)
beta_1 <- matrix(runif(p*3),nrow=p,ncol=3)
beta_2 <- matrix(0.5*runif(p*3),nrow=p,ncol=3)
beta_3 <- matrix(0.1*runif(p),nrow=p,ncol=1)
Sig <- matrix(c(2,1,0,1,2,1,0,1,2),nrow=3,ncol=3)

n_group <- 50
group_num <- c(25,10,5,5,5) #true null proportion = 0.5
tr <- rep(c(0,1,2,3,4),group_num)
group_sizes <- 5+rpois(50,20)
null <- as.numeric(tr==0)

#training data
X_train <- matrix(runif(n_train*p),nrow=n_train,ncol=p)
Y_123_train <- t(sapply(1:n_train,function(j) mvrnorm(1,mu=(X_train%*%beta_1)[j,],Sigma=Sig)))
Y_4_train <- sapply(1:n_train,function(j) rchisq(1,sum(Y_123_train[j,1:2]^2)))
Y_5_train <- sapply(1:n_train,function(j) rbeta(1,abs(Y_123_train[j,1]),abs(Y_123_train[j,1])+0.5*abs(Y_123_train[j,1])))

Y_train <- cbind(Y_123_train,Y_4_train,Y_5_train)
colnames(Y_train) <- c("","","","","")

#constructing scores
library(randomForest)

##score 1
rf_1 <- randomForest(X_train,Y_train[,1])
rf_2 <- randomForest(X_train,Y_train[,2])
rf_3 <- randomForest(X_train,Y_train[,3])
rf_4 <- randomForest(X_train,Y_train[,4])
rf_5 <- randomForest(X_train,Y_train[,5])

muhat_train <- matrix(0,nrow=n_train,ncol=5)
muhat_train[,1] <- predict(rf_1,X_train)
muhat_train[,2] <- predict(rf_2,X_train)
muhat_train[,3] <- predict(rf_3,X_train)
muhat_train[,4] <- predict(rf_4,X_train)
muhat_train[,5] <- predict(rf_5,X_train)

sigmahat <- cov(muhat_train)

tn <- 500
fdp <- matrix(0,nrow=tn,ncol=length(alpha_set))
tpp <- matrix(0,nrow=tn,ncol=length(alpha_set))
alpha_set <- seq(0.05,0.3,by=0.05)

pb <- txtProgressBar(min = 0, max = tn, style = 3)
for(i in 1:tn)
{
  
  #calibration data
  X_cal <- matrix(runif(n_cal*p),nrow=n_cal,ncol=p)
  Y_123_cal <- t(sapply(1:n_cal,function(j) mvrnorm(1,mu=(X_cal%*%beta_1)[j,],Sigma=Sig)))
  Y_4_cal <- sapply(1:n_cal,function(j) rchisq(1,0.1*sum(Y_123_cal[j,1:2]^2)))
  Y_5_cal <- sapply(1:n_cal,function(j) rbeta(1,abs(Y_123_cal[j,1]),abs(Y_123_cal[j,1])+0.5*abs(Y_123_cal[j,1])))
  
  Y_cal <- cbind(Y_123_cal,Y_4_cal,Y_5_cal)
  
  #test data
  X_test <- lapply(1:n_group, function(k) matrix(runif(group_sizes[k]*p),nrow=group_sizes[k],ncol=p))
  Y_test <- list()
  for(k in 1:n_group)
  {
    Y_123_k <- t(sapply(1:group_sizes[k],function(j) mvrnorm(1,mu=(X_test[[k]]%*%beta_1)[j,]+tr[k]*(X_test[[k]]%*%beta_2)[j,]^2,Sigma=Sig)))
    Y_4_k <- sapply(1:group_sizes[k],function(j) rchisq(1,3*tr[k]*abs(sum(X_test[[k]][j,]*beta_3))+0.1*sum(Y_123_k[j,1:2]^2)))
    Y_5_k <- sapply(1:group_sizes[k],function(j) rbeta(1,abs(Y_123_k[j,1]),abs(Y_123_k[j,1])+0.5*abs(Y_123_k[j,1])))
    
    Y_test[[k]] <- cbind(Y_123_k,Y_4_k,Y_5_k)
  }
  
  ##main procedure
  muhat_cal <- matrix(0,nrow=n_cal,ncol=5)
  muhat_cal[,1] <- predict(rf_1,X_cal)
  muhat_cal[,2] <- predict(rf_2,X_cal)
  muhat_cal[,3] <- predict(rf_3,X_cal)
  muhat_cal[,4] <- predict(rf_4,X_cal)
  muhat_cal[,5] <- predict(rf_5,X_cal)
  
  score_cal <- diag((Y_cal-muhat_cal)%*%solve(sigmahat,t(Y_cal-muhat_cal)))
  
  
  p_val <- rep(0,n_group)
  for(k in 1:n_group)
  {
    muhat_k <- matrix(0,nrow=group_sizes[k],ncol=5)
    muhat_k[,1] <- predict(rf_1,X_test[[k]])
    muhat_k[,2] <- predict(rf_2,X_test[[k]])
    muhat_k[,3] <- predict(rf_3,X_test[[k]])
    muhat_k[,4] <- predict(rf_4,X_test[[k]])
    muhat_k[,5] <- predict(rf_5,X_test[[k]])
    
    score_k <- diag((Y_test[[k]]-muhat_k)%*%solve(sigmahat,t(Y_test[[k]]-muhat_k)))
    eta_k <- max(round(group_sizes[k]/2),1)
    
    p_val[k] <- batch_conformal(score_cal,score_k,eta_k)
    
  }

  p_threshold <- sapply(alpha_set, function(alpha) BH(p_val,alpha))
  rej <- lapply(p_threshold, function(t) as.numeric(p_val <= t))
  V <- sapply(rej, function(v) sum(null==1 & v==1))
  R <- sapply(rej, function(v) sum(v==1))
  W <- sapply(rej, function(v) sum(null==0 & v==1))
  fdp[i,] <- V/max(R,1)
  tpp[i,] <- W/max(R,1)
  setTxtProgressBar(pb, i)
}


fdr <- colMeans(fdp)
fdr_se <- apply(fdp,2,sd)/sqrt(tn)
power <- colMeans(tpp)
power_se <- apply(tpp,2,sd)/sqrt(tn)

dev.off()
plot.new()
plot(alpha_set,colMeans(fdp),type="l",ylim=c(0,0.2),xlab=expression(alpha),ylab="FDR")
segments(0.05,0.05,0.3,0.3*0.5,lty=2)





##score 2
rf_1 <- randomForest(as.data.frame(X_train),Y_train[,1])
rf_2 <- randomForest(as.data.frame(cbind(X_train,Y_train[,1])),Y_train[,2])
rf_3 <- randomForest(as.data.frame(cbind(X_train,Y_train[,1:2])),Y_train[,3])
rf_4 <- randomForest(as.data.frame(cbind(X_train,Y_train[,1:3])),Y_train[,4])
rf_5 <- randomForest(as.data.frame(cbind(X_train,Y_train[,1:4])),Y_train[,5])

muhat_train <- matrix(0,nrow=n_train,ncol=5)
muhat_train[,1] <- predict(rf_1,X_train)
muhat_train[,2] <- predict(rf_2,as.data.frame(cbind(X_train,Y_train[,1])))
muhat_train[,3] <- predict(rf_3,as.data.frame(cbind(X_train,Y_train[,1:2])))
muhat_train[,4] <- predict(rf_4,as.data.frame(cbind(X_train,Y_train[,1:3])))
muhat_train[,5] <- predict(rf_5,as.data.frame(cbind(X_train,Y_train[,1:4])))

sigmahat <- cov(muhat_train)


tn <- 500
fdp <- matrix(0,nrow=tn,ncol=length(alpha_set))
tpp <- matrix(0,nrow=tn,ncol=length(alpha_set))
alpha_set <- seq(0.05,0.3,by=0.05)

pb <- txtProgressBar(min = 0, max = tn, style = 3)
for(i in 1:tn)
{
  
  #calibration data
  X_cal <- matrix(runif(n_cal*p),nrow=n_cal,ncol=p)
  Y_123_cal <- t(sapply(1:n_cal,function(j) mvrnorm(1,mu=(X_cal%*%beta_1)[j,],Sigma=Sig)))
  Y_4_cal <- sapply(1:n_cal,function(j) rchisq(1,0.1*sum(Y_123_cal[j,1:2]^2)))
  Y_5_cal <- sapply(1:n_cal,function(j) rbeta(1,abs(Y_123_cal[j,1]),abs(Y_123_cal[j,1])+0.5*abs(Y_123_cal[j,1])))
  
  Y_cal <- cbind(Y_123_cal,Y_4_cal,Y_5_cal)
  colnames(Y_cal) <- c("","","","","")
  
  #test data
  
  X_test <- lapply(1:n_group, function(k) matrix(runif(group_sizes[k]*p),nrow=group_sizes[k],ncol=p))
  Y_test <- list()
  for(k in 1:n_group)
  {
    Y_123_k <- t(sapply(1:group_sizes[k],function(j) mvrnorm(1,mu=(X_test[[k]]%*%beta_1)[j,]+tr[k]*(X_test[[k]]%*%beta_2)[j,]^2,Sigma=Sig)))
    Y_4_k <- sapply(1:group_sizes[k],function(j) rchisq(1,3*tr[k]*abs(sum(X_test[[k]][j,]*beta_3))+0.1*sum(Y_123_k[j,1:2]^2)))
    Y_5_k <- sapply(1:group_sizes[k],function(j) rbeta(1,abs(Y_123_k[j,1]),abs(Y_123_k[j,1])+0.5*abs(Y_123_k[j,1])))
    
    Y_test[[k]] <- cbind(Y_123_k,Y_4_k,Y_5_k)
    colnames(Y_test[[k]]) <- c("","","","","")
  }
  
  ##main procedure
  muhat_cal <- matrix(0,nrow=n_cal,ncol=5)
  muhat_cal[,1] <- predict(rf_1,as.data.frame(X_cal))
  muhat_cal[,2] <- predict(rf_2,as.data.frame(cbind(X_cal,Y_cal[,1])))
  muhat_cal[,3] <- predict(rf_3,as.data.frame(cbind(X_cal,Y_cal[,1:2])))
  muhat_cal[,4] <- predict(rf_4,as.data.frame(cbind(X_cal,Y_cal[,1:3])))
  muhat_cal[,5] <- predict(rf_5,as.data.frame(cbind(X_cal,Y_cal[,1:4])))
  
  score_cal <- diag((Y_cal-muhat_cal)%*%solve(sigmahat,t(Y_cal-muhat_cal)))
  
  
  p_val <- rep(0,n_group)
  for(k in 1:n_group)
  {
    muhat_k <- matrix(0,nrow=group_sizes[k],ncol=5)
    muhat_k[,1] <- predict(rf_1,as.data.frame(X_test[[k]]))
    muhat_k[,2] <- predict(rf_2,as.data.frame(cbind(X_test[[k]],Y_test[[k]][,1])))
    muhat_k[,3] <- predict(rf_3,as.data.frame(cbind(X_test[[k]],Y_test[[k]][,1:2])))
    muhat_k[,4] <- predict(rf_4,as.data.frame(cbind(X_test[[k]],Y_test[[k]][,1:3])))
    muhat_k[,5] <- predict(rf_5,as.data.frame(cbind(X_test[[k]],Y_test[[k]][,1:4])))
    
    score_k <- diag((Y_test[[k]]-muhat_k)%*%solve(sigmahat,t(Y_test[[k]]-muhat_k)))
    
    eta_k <- round(group_sizes[k]/2)
    
    p_val[k] <- batch_conformal(score_cal,score_k,eta_k)
    
  }
  
  p_threshold <- sapply(alpha_set, function(alpha) BH(p_val,alpha))
  rej <- lapply(p_threshold, function(t) as.numeric(p_val <= t))
  V <- sapply(rej, function(v) sum(null==1 & v==1))
  R <- sapply(rej, function(v) sum(v==1))
  W <- sapply(rej, function(v) sum(null==0 & v==1))
  fdp[i,] <- V/max(R,1)
  tpp[i,] <- W/max(R,1)
  setTxtProgressBar(pb, i)
}

fdr <- colMeans(fdp)
fdr_se <- apply(fdp,2,sd)/sqrt(tn)
power <- colMeans(tpp)
power_se <- apply(tpp,2,sd)/sqrt(tn)

plot(alpha_set,colMeans(fdp),type="l",ylim=c(0,0.2),xlab=expression(alpha),ylab="FDR")
segments(0.05,0.05,0.3,0.3*0.5,lty=2)


##score 3
rf_1 <- randomForest(as.data.frame(X_train),Y_train[,5])
rf_2 <- randomForest(as.data.frame(cbind(X_train,Y_train[,5])),Y_train[,4])
rf_3 <- randomForest(as.data.frame(cbind(X_train,Y_train[,4:5])),Y_train[,3])
rf_4 <- randomForest(as.data.frame(cbind(X_train,Y_train[,3:5])),Y_train[,2])
rf_5 <- randomForest(as.data.frame(cbind(X_train,Y_train[,2:5])),Y_train[,1])

muhat_train <- matrix(0,nrow=n_train,ncol=5)
muhat_train[,5] <- predict(rf_1,X_train)
muhat_train[,4] <- predict(rf_2,as.data.frame(cbind(X_train,Y_train[,5])))
muhat_train[,3] <- predict(rf_3,as.data.frame(cbind(X_train,Y_train[,4:5])))
muhat_train[,2] <- predict(rf_4,as.data.frame(cbind(X_train,Y_train[,3:5])))
muhat_train[,1] <- predict(rf_5,as.data.frame(cbind(X_train,Y_train[,2:5])))

sigmahat <- cov(muhat_train)


tn <- 500
fdp <- matrix(0,nrow=tn,ncol=length(alpha_set))
tpp <- matrix(0,nrow=tn,ncol=length(alpha_set))
alpha_set <- seq(0.05,0.3,by=0.05)

pb <- txtProgressBar(min = 0, max = tn, style = 3)
for(i in 1:tn)
{
  
  #calibration data
  X_cal <- matrix(runif(n_cal*p),nrow=n_cal,ncol=p)
  Y_123_cal <- t(sapply(1:n_cal,function(j) mvrnorm(1,mu=(X_cal%*%beta_1)[j,],Sigma=Sig)))
  Y_4_cal <- sapply(1:n_cal,function(j) rchisq(1,0.1*sum(Y_123_cal[j,1:2]^2)))
  Y_5_cal <- sapply(1:n_cal,function(j) rbeta(1,abs(Y_123_cal[j,1]),abs(Y_123_cal[j,1])+0.5*abs(Y_123_cal[j,1])))
  
  Y_cal <- cbind(Y_123_cal,Y_4_cal,Y_5_cal)
  colnames(Y_cal) <- c("","","","","")
  
  #test data
  
  X_test <- lapply(1:n_group, function(k) matrix(runif(group_sizes[k]*p),nrow=group_sizes[k],ncol=p))
  Y_test <- list()
  for(k in 1:n_group)
  {
    Y_123_k <- t(sapply(1:group_sizes[k],function(j) mvrnorm(1,mu=(X_test[[k]]%*%beta_1)[j,]+tr[k]*(X_test[[k]]%*%beta_2)[j,]^2,Sigma=Sig)))
    Y_4_k <- sapply(1:group_sizes[k],function(j) rchisq(1,3*tr[k]*abs(sum(X_test[[k]][j,]*beta_3))+0.1*sum(Y_123_k[j,1:2]^2)))
    Y_5_k <- sapply(1:group_sizes[k],function(j) rbeta(1,abs(Y_123_k[j,1]),abs(Y_123_k[j,1])+0.5*abs(Y_123_k[j,1])))
    
    Y_test[[k]] <- cbind(Y_123_k,Y_4_k,Y_5_k)
    colnames(Y_test[[k]]) <- c("","","","","")
  }
  
  ##main procedure
  muhat_cal <- matrix(0,nrow=n_cal,ncol=5)
  muhat_cal[,5] <- predict(rf_1,as.data.frame(X_cal))
  muhat_cal[,4] <- predict(rf_2,as.data.frame(cbind(X_cal,Y_cal[,5])))
  muhat_cal[,3] <- predict(rf_3,as.data.frame(cbind(X_cal,Y_cal[,4:5])))
  muhat_cal[,2] <- predict(rf_4,as.data.frame(cbind(X_cal,Y_cal[,3:5])))
  muhat_cal[,1] <- predict(rf_5,as.data.frame(cbind(X_cal,Y_cal[,2:5])))
  
  score_cal <- diag((Y_cal-muhat_cal)%*%solve(sigmahat,t(Y_cal-muhat_cal)))
  
  
  p_val <- rep(0,n_group)
  for(k in 1:n_group)
  {
    muhat_k <- matrix(0,nrow=group_sizes[k],ncol=5)
    muhat_k[,5] <- predict(rf_1,as.data.frame(X_test[[k]]))
    muhat_k[,4] <- predict(rf_2,as.data.frame(cbind(X_test[[k]],Y_test[[k]][,5])))
    muhat_k[,3] <- predict(rf_3,as.data.frame(cbind(X_test[[k]],Y_test[[k]][,4:5])))
    muhat_k[,2] <- predict(rf_4,as.data.frame(cbind(X_test[[k]],Y_test[[k]][,3:5])))
    muhat_k[,1] <- predict(rf_5,as.data.frame(cbind(X_test[[k]],Y_test[[k]][,2:5])))
    
    score_k <- diag((Y_test[[k]]-muhat_k)%*%solve(sigmahat,t(Y_test[[k]]-muhat_k)))
    
    eta_k <- round(group_sizes[k]/2)
    
    p_val[k] <- batch_conformal(score_cal,score_k,eta_k)
    
  }
  
  p_threshold <- sapply(alpha_set, function(alpha) BH(p_val,alpha))
  rej <- lapply(p_threshold, function(t) as.numeric(p_val <= t))
  V <- sapply(rej, function(v) sum(null==1 & v==1))
  R <- sapply(rej, function(v) sum(v==1))
  W <- sapply(rej, function(v) sum(null==0 & v==1))
  fdp[i,] <- V/max(R,1)
  tpp[i,] <- W/max(R,1)
  setTxtProgressBar(pb, i)
}

fdr <- colMeans(fdp)
fdr_se <- apply(fdp,2,sd)/sqrt(tn)
power <- colMeans(tpp)
power_se <- apply(tpp,2,sd)/sqrt(tn)

plot(alpha_set,colMeans(fdp),type="l",ylim=c(0,0.2),xlab=expression(alpha),ylab="FDR")
segments(0.05,0.05,0.3,0.3*0.5,lty=2)



