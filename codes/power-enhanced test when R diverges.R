library("MASS")
library("expm")
library("sphunif")
library("mvtnorm")
library("energy")
library("HHG")
library("dHSIC")
library("dcov")
library("semidist")
library("extraDistr")
###Basic Function###
vm <- function(x, y, R, py, y_ca) {
  n <- length(x)
  num_geq_x <- n - rank(x, ties.method = "min") + 1
  vm_rec <- numeric(R)
  for (r in 1:R) {
    if(py[r]==0){vm_rec[r]=0}
    if(py[r]!=0){
      sum_count_r <- sum(num_geq_x[y == y_ca[r]])
      sum_vmm <- sum_count_r / (n * py[r])
      vmm2 <- ((sum_vmm / (n + 1)) - 0.5)^2 / ((1 / (12 * (n + 1))) * ((1 / py[r]) - 1))
      vm_rec[r] <- vmm2
    }
  }
  sum(py * vm_rec)
}

max_vm <- function(x, y, R, py, y_ca) {
  n <- length(x)
  num_geq_x <- n - rank(x, ties.method = "min") + 1
  vm_rec <- numeric(R)
  for (r in 1:R) {
    if(py[r]==0){vm_rec[r]=0}
    if(py[r]!=0){
      sum_count_r <- sum(num_geq_x[y == y_ca[r]])
      sum_vmm <- sum_count_r / (n * py[r])
      vmm2 <- ((sum_vmm / (n + 1)) - 0.5)^2 / ((1 / (12 * (n + 1))) * ((1 / py[r]) - 1))
      vm_rec[r] <- vmm2
    }
  }
  max(vm_rec)
}

sdcov_fast<- function(X, y) {
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  } else if (!is.matrix(X)) {
    X <- as.matrix(X, ncol = 1)
  }
  n <- nrow(X)
  group_idx <- split(seq_len(n), y)
  group_idx <- group_idx[lengths(group_idx) > 1]
  
  if (length(group_idx) == 0) {
    return(0)
  }
  sum_x <- 2 * sum(dist(X))
  sum_y_term <- 0
  for (idx in group_idx) {
    n_r <- length(idx)
    sum_grp <- 2 * sum(dist(X[idx, , drop = FALSE]))
    sum_y_term <- sum_y_term + sum_grp / (n_r - 1)
  }
  sdcov <- sum_x / (n * (n - 1)) - sum_y_term / n
  return(sdcov)
}

sim=1000
pow1=matrix(NA,nrow=sim,ncol=6)
pow2=matrix(NA,nrow=sim,ncol=6)
pow3=matrix(NA,nrow=sim,ncol=6)
pow4=matrix(NA,nrow=sim,ncol=6)
aa=c(0,0.4,0.8,1.2,1.6,2.0)
n=1000
R=100
p=1
for (m in 1:sim) {
  #print(m)
  for (o in 1:length(aa)) {
    prob_set=rep(1/R,R)
    #prob_set=sapply(1:R, function(i) (2*(1+(i-1)/(R-1)))/(3*R))
    y=sample(rep(1:R), n, replace = TRUE, prob = prob_set)
    x=matrix(NA,nrow = n,ncol = p)
    for (i in 1:n) {
      if(y[i]==100){
        x[i,]=rnorm(1,mean = aa[o],sd=1)
      }
      if(y[i]!=100){
        x[i,]=rnorm(1,mean = aa[o]*0.01,sd=1)
      }
    }
    
    ###G-sum Test###
    y_ca=rep(1:R)
    p_y=sapply(1:length(y_ca),function(i) length(which(y==y_ca[i]))/(length(y)))
    sta1=vm(x,y,R,p_y,y_ca)
    sta_use1=(R*(sta1-1)^2)/2
    pval1 <-1-pchisq(sta_use1,df=1)
    pow1[m,o]=ifelse(pval1<=0.05,1,0)
    
    ###Power-enhanced Test###
    sta2=max_vm(x,y,R,p_y,y_ca)
    rec=sta2-2*log(R)+log(log(R))
    p_max=1-exp((-1/sqrt(pi))*exp(-rec/2))
    T_final=-2*log(pval1)-2*log(p_max)
    pval2 <-1-pchisq(T_final,df=4)
    pow2[m,o]=ifelse(pval2<=0.05,1,0)
    
    ###SD Test with Asymptotic Distribution###
    sdp=((n*sdcov_fast(as.matrix(x),as.factor(y)))/(sqrt(R-1)*dcov2d(x,x)))/sqrt(2)
    pval3=1-pchisq(sdp^2,df=1)
    pow3[m,o] <- ifelse(pval3 <= 0.05, 1, 0)
    
    ###Permutation Test for SD###
    pval4=sd_test(as.matrix(x), as.factor(y), test_type = "perm", num_perm = 200)$pvalue
    pow4[m,o] <- ifelse(pval4 <= 0.05, 1, 0)
    
  }
}
colMeans(pow1)
colMeans(pow2)
colMeans(pow3)
colMeans(pow4)