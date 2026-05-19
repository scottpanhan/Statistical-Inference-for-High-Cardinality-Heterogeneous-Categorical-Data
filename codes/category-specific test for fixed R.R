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
  vm_rec
}

sim=1000
aa=c(0,1,2,3)
n=50
R=100
p=1
pow=array(NA, dim = c(sim, R, length(aa)))
for (m in 1:sim) {
  print(m)
  pow_rec=matrix(NA,nrow = R,length(aa))
  for (o in 1:length(aa)) {
    prob_set=rep(1/R,R)
    #prob_set=sapply(1:R, function(i) (2*(1+(i-1)/(R-1)))/(3*R))
    y=sample(rep(1:R), n, replace = TRUE, prob = prob_set)
    mu_x=aa[o]*c(0,0,0,0,0,0,0,0,-0.5,0.5)
    x=matrix(NA,nrow = n,ncol = p)
    for (i in 1:n) {
      x[i,]=mu_x[y[i]]+rnorm(1)
    }
    y_ca=rep(1:R)
    p_y=sapply(1:length(y_ca),function(i) length(which(y==y_ca[i]))/(length(y)))
    sta1=vm(x,y,R,p_y,y_ca)
    
    
    for (rr in 1:R) {
      cate=sta1[rr]^2
      pval=1-pchisq(cate,df=1)
      pow_rec[rr,o]=ifelse(pval<=0.05,1,0)
    }
    
  }
  pow[m,,]=pow_rec
}
average_power <- apply(pow, c(2,3), mean)
average_power