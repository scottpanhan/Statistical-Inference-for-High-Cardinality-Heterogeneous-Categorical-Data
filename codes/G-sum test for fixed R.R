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
###G-Sum Test###
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
###VM Test###
vm2 <- function(x, y, R, py, y_ca) {
  n <- length(x)
  num_geq_x <- n - rank(x, ties.method = "min") + 1
  vm_rec <- numeric(R)
  for (r in 1:R) {
    if(py[r]==0){vm_rec[r]=0}
    if(py[r]!=0){
      sum_count_r <- sum(num_geq_x[y == y_ca[r]])
      sum_vmm <- sum_count_r / (n * py[r])
      vmm2 <- ((sum_vmm / (n + 1)) - 0.5)^2
      vm_rec[r] <- vmm2
    }
  }
  sum(py * vm_rec)
}

sim=1000
pow1=matrix(NA,nrow=sim,ncol=4)
pow2=matrix(NA,nrow=sim,ncol=4)
pow3=matrix(NA,nrow=sim,ncol=4)
pow4=matrix(NA,nrow=sim,ncol=4)
pow5=matrix(NA,nrow=sim,ncol=4)
pow6=matrix(NA,nrow=sim,ncol=4)
pow7=matrix(NA,nrow=sim,ncol=4)
pow8=matrix(NA,nrow=sim,ncol=4)
pow9=matrix(NA,nrow=sim,ncol=4)
pow10=matrix(NA,nrow=sim,ncol=4)
aa=c(0,1,2,3)
n=100
R=2
p=1
for (m in 1:sim) {
  print(m)
  for (o in 1:length(aa)) {
    prob_set=sapply(1:R, function(i) (2*(1+(i-1)/(R-1)))/(3*R))
    y=sample(rep(1:R), n, replace = TRUE, prob = prob_set)
    mu_x=aa[o]*c(-1/3,0)
    x=matrix(NA,nrow = n,ncol = p)
    for (i in 1:n) {
      x[i,]=mu_x[y[i]]+rnorm(1)
    }
    y_ca=rep(1:R)
    p_y=sapply(1:length(y_ca),function(i) length(which(y==y_ca[i]))/(length(y)))
    ###G-sum Test###
    sta1=vm(x,y,R,p_y,y_ca)
    sig_Z=diag(1,R)
    for (i in 1:(R-1)) {
      for (j in (i+1):R) {
        sig_Z[i,j]=sig_Z[j,i]=-sqrt(p_y[i]*p_y[j])/sqrt((1-p_y[i])*(1-p_y[j]))
      }
    }
    AA=diag(p_y)
    sig_sqr=sqrtm(AA)
    BB=sig_sqr%*%sig_Z%*%sig_sqr
    value_BB=eigen(BB)$values
    pval1 <-1-p_wschisq(sta1,weights = value_BB,dfs = rep(1,R))
    pow1[m,o]=ifelse(pval1<=0.05,1,0)
    
    ###VM Test###
    sta2=n*vm2(x,y,R,p_y,y_ca)
    pval2 <-1-p_wschisq(sta2,weights = rep(1/12,R-1),dfs = rep(1,R-1))
    pow2[m,o]=ifelse(pval2<=0.05,1,0)
    
    ###Permutation tests for G-sum and VM methods###
    pr_sta1 <- c()
    pr_sta2 <- c()
    for (k in 1:200) {
      pr_sim <- sample(1:n, n, replace = FALSE)
      y_pe <- y[pr_sim]
      pr_sta1[k] <- vm(x,y_pe,R,p_y,y_ca)
      pr_sta2[k] <- n*vm2(x,y_pe,R,p_y,y_ca)
    }
    p_e1 <- mean(pr_sta1 >= sta1)
    pow3[m,o] <- ifelse(p_e1 <= 0.05, 1, 0)
    p_e2 <- mean(pr_sta2 >= sta2)
    pow4[m,o] <- ifelse(p_e2 <= 0.05, 1, 0)
    
    ###DC Test###
    y_dum=switch_cat_repr(as.factor(y))
    dtp=dcov.test(x, y_dum, R = 200)$p.value
    pow5[m,o]=ifelse(dtp<=0.05,1,0)
    
    ###HHG Test###
    dist_y <- as.matrix(dist(y_dum))
    dist_x <- as.matrix(dist(x))
    hhgp=hhg.test(dist_x, dist_y, nr.perm = 200)$perm.pval.hhg.sc
    pow6[m,o]=ifelse(hhgp<=0.05,1,0)
    
    ###SD Test###
    sdp=sd_test(as.matrix(x), as.factor(y), test_type = "perm", num_perm = 200)$pvalue
    pow7[m,o]=ifelse(sdp<=0.05,1,0)
    
    ###MINT###
    try({mintp=MINTsemiauto(x, as.factor(y), kmax=8, B1 = 50, B2 =200)$pvalue},silent = TRUE)
    pow8[m,o]=ifelse(mintp<=0.05,1,0)
    
    ###HSIC###
    hsicp=dhsic.test(dist_x, dist_y,B=200)$p.value
    pow9[m,o]=ifelse(hsicp<=0.05,1,0)
    
    ###MV Test###
    mvp=mv_test(x, as.factor(y),test_type = "perm", num_perm = 200)$pvalue
    pow10[m,o]=ifelse(mvp<=0.05,1,0)
  }
}
colMeans(pow1)
colMeans(pow2)
colMeans(pow3)
colMeans(pow4)
colMeans(pow5)
colMeans(pow6)
colMeans(pow7)
colMeans(pow8)
colMeans(pow9)
colMeans(pow10)