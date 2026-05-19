library("knockoff")
library("splines")
library("IHW")
library("adaptMT")
library("MASS")
library("expm")
library("sphunif")
library("mvtnorm")
library("corpcor")  
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

CAFDR <- function(T_n, tilT_n, pTn, ptilTn, sigma, alpha, tau = 0.5) {
  
  p <- length(T_n)
  weighted_kernel_est <- function(t_val, sigma_val, T_obs, T_cal, sigma_obs, h_t, h_s) {
    
    kernel_weights <- dnorm(sigma_val - sigma_obs, mean = 0, sd = h_s)
    numerator <- sum(kernel_weights * (dnorm(t_val - T_obs, mean = 0, sd = h_t) + dnorm(t_val - T_cal, mean = 0, sd = h_t)))
    denominator <- 2 * sum(kernel_weights)
    
    return(numerator / denominator)
  }

  estimate_pi_sigma <- function(sigma_val, p_T_obs, p_T_cal, sigma_obs, h_s, tau_val) {
    
    kernel_weights <- dnorm(sigma_val - sigma_obs, mean = 0, sd = h_s)
    
    numerator <- sum(kernel_weights * ( (p_T_obs > tau_val) + (p_T_cal > tau_val) ))
    denominator <- 2 * (1 - tau_val) * sum(kernel_weights)
    
    pi_hat <- 1 - (numerator / denominator)
    return(pi_hat)
  }
  
  f_hat_sigma <- numeric(p)
  f_hat_sigma_tilde <- numeric(p)
  pi_hat_sigma <- numeric(p)
  
  h_sigma <- bw.SJ(sigma)
  h_T <- bw.SJ(c(T_n,tilT_n))
  
  for (j in 1:p) {
    f_hat_sigma[j] <- weighted_kernel_est(T_n[j], sigma[j], T_n, tilT_n, sigma, h_T, h_sigma)
    f_hat_sigma_tilde[j] <- weighted_kernel_est(tilT_n[j], sigma[j], T_n, tilT_n, sigma, h_T, h_sigma)
    pi_hat_sigma[j] <- estimate_pi_sigma(sigma[j], pTn, ptilTn, sigma, h_sigma, tau)
  }

  eta <- 0.001
  pi_hat_star <- pi_hat_sigma
  pi_hat_star[pi_hat_sigma <= 0] <- eta
  pi_hat_star[pi_hat_sigma > 0.5] <- 0.5 - eta

  # Null density f0(t)
  f0_T <- d_wschisq(T_n,weights = value_BB,dfs = rep(1,R))
  f0_T_tilde <- d_wschisq(tilT_n,weights = value_BB,dfs = rep(1,R))

  m <- 0.999
  Clfdr_hat_T_mod <- pmin((1 - pi_hat_star) * f0_T / f_hat_sigma, m)
  Clfdr_hat_T_tilde_mod <- pmin((1 - pi_hat_star) * f0_T_tilde / f_hat_sigma_tilde, m)
  
  term <- (1 - Clfdr_hat_T_mod) / Clfdr_hat_T_mod
  s <- 1 + (((1 - pi_hat_star) / (0.5 - pi_hat_star)) * term)
  s <- 1 / s 
  
  term_tilde <- (1 - Clfdr_hat_T_tilde_mod) / Clfdr_hat_T_tilde_mod
  s_tilde <- 1 + (((1 - pi_hat_star) / (0.5 - pi_hat_star)) * term_tilde)
  s_tilde <- 1 / s_tilde 

  Q1 <- which(s < s_tilde)
  Q2 <- which(s_tilde < s)

  t_values <- sort(unique(c(s, s_tilde)))
  
  M_t <- sapply(t_values, function(t) {
    numerator <- 1 + sum(s_tilde[Q2] <= t)
    denominator <- max(1, sum(s[Q1] <= t))
    return(numerator / denominator)
  })
  
  valid_t <- t_values[M_t <= alpha]
  gamma <- if (length(valid_t) > 0) max(valid_t) else -Inf

  R <- Q1[s[Q1] <= gamma]

  return(sort(R))
}

group.func = function(data1, data2, f0){
  m = length(data1)
  mix.dens = density(c(data1,data2),from=min(c(data1,data2))-10, 
                     to=max(c(data1,data2))+10, n=1000)
  drx = dnorm(data1,f0[1],f0[2])/ lin.itp(data1, mix.dens$x, mix.dens$y)
  drm = dnorm(data2,f0[1],f0[2])/ lin.itp(data2, mix.dens$x, mix.dens$y)
  
  pvs = 2*(1-pnorm(abs(c(data1,data2)),f0[1],f0[2]))
  pp0 = sum(pvs>0.5)/(0.5*2*m)
  if (pp0==Inf | is.na(pp0)){
    pp0 = 1
  }
  
  s = pp0 * drx
  st = pp0 * drm
  
  s = pmin(s,0.999)
  st = pmin(st,0.999)
  pp = 1-pp0
  
  s = (0.5-pp)/pp0*s/(1-s)
  st = (0.5-pp)/pp0*st/(1-st)
  s = s/(1+s)
  st = st/(1+st)
  return(list(s=s,st=st))
}

Adadetect = function(x,y,al){
  m = length(x)
  t = sort(unique(c(x, y)))
  fr = NA
  for (j in 1:length(t)){
    fr[j] = (sum(y<=t[j])+1)/(m+1) /(sum(x<=t[j])/m)
  }
  t.hat = max(na.omit(t[fr<=al]))
  decision = which(x<=t.hat)
  return(decision)
}

lin.itp<-function(x, X, Y){
  x.N<-length(x)
  X.N<-length(X)
  y<-rep(0, x.N)
  for (k in 1:x.N){
    i<-max(which((x[k]-X)>=0))
    if (i<X.N)
      y[k]<-Y[i]+(Y[i+1]-Y[i])/(X[i+1]-X[i])*(x[k]-X[i])
    else 
      y[k]<-Y[i]
  }
  return(y)
}

generate_corr_matrix<-function(p,rho){
  corr_matrix<-matrix(NA,nrow = p,ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      if(i==j){
        corr_matrix[i,j]<-1
      }else{
        corr_matrix[i,j]<-rho^abs(i-j)
      }
    }
  }
  return(corr_matrix)
}

compute_fdr_power <- function(rejection_set, true_signal) {

	n_discoveries <- length(rejection_set)
	
	if (n_discoveries == 0) {
		fdr <- 0
		power <- 0
	} else {
		n_true_discoveries <- sum(true_signal[rejection_set])
		n_false_discoveries <- n_discoveries - n_true_discoveries
		
		fdr <- n_false_discoveries / n_discoveries

		n_true_signals <- sum(true_signal)
		if (n_true_signals > 0) {
			power <- n_true_discoveries / n_true_signals
		} else {
			power <- NA
		}
	}
	
	return(list(fdr = fdr, power = power))
}

generate_data <- function(c_val, distribution = "normal", n = 100, p = 1000, l = 50, R = 5, seed = 123) {
	set.seed(seed)
	Y <- sample(1:R, n, replace = TRUE, prob = rep(1/R, R))
	S <- c(1, 5)
	Sigma <- generate_corr_matrix(p,0.5)
	l_half <- floor(l/2)
	D1 <- diag(runif(l_half, 1, 2))
	Sigma11 <- generate_corr_matrix(l_half,0.8)
	Sigma1 <- D1 %*% Sigma11 %*% D1
	D2 <- diag(runif(l_half, 5, 6))
	Sigma22 <- generate_corr_matrix(l_half,0.2)
	Sigma2 <- D2 %*% Sigma22 %*% D2
	Sigma <- make.positive.definite(Sigma)
	Sigma1 <- make.positive.definite(Sigma1)
	Sigma2 <- make.positive.definite(Sigma2)
	X <- matrix(0, n, p)
	for (i in 1:n) {
		if (distribution == "normal") {
			z_i <- mvrnorm(1, mu = rep(0, p), Sigma = Sigma)
		} else if (distribution == "t") {
			z_i <- rmvt(1, sigma = Sigma, df = 3)
		}
		
		if (Y[i] %in% S) {
			u1 <- runif(l_half, 0.5, 2)
			u2 <- runif(l_half, 0.5, 2)
			if (distribution == "normal") {
				eta1 <- mvrnorm(1, mu = rep(0, l_half), Sigma = Sigma1)
				eta2 <- mvrnorm(1, mu = rep(0, l_half), Sigma = Sigma2)
			} else if (distribution == "t") {
				eta1 <- rmvt(1, sigma = Sigma1, df = 3)
				eta2 <- rmvt(1, sigma = Sigma2, df = 3)
			}
			mu_i <- numeric(l)
			mu_i[1:l_half] <- c_val * u1 + eta1
			mu_i[(l_half+1):l] <- c_val * u2 + eta2
			X[i, 1:l] <- mu_i + z_i[1:l]
			X[i, (l+1):p] <- z_i[(l+1):p]
			
		} else {
			X[i, ] <- z_i
		}
	}
	true_signal <- rep(0, p)
	true_signal[1:l] <- 1
	return(list(
		X = X,
		Y = Y,
		true_signal = true_signal,
		c = c_val,
		distribution = distribution,
		Sigma = Sigma,
		Sigma1 = Sigma1,
		Sigma2 = Sigma2
	))
}

###Simulation###
n = 100
p = 1000
l = 50  
R = 5
c = 1.5
seed=123
data <- generate_data(c,"normal",n,p,l,R,seed)
x <- data$X
y <- data$Y
y_ca=rep(1:5)
p_y=sapply(1:R,function(i) length(which(y==y_ca[i]))/n)
T_n=c()
for (j in 1:p) {
  T_n[j]=vm(x[,j],y,R,p_y,y_ca)
}
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
tilT_n=r_wschisq(p,weights = value_BB,dfs = rep(1,R))
pp_Tn=c()
pp_tilTn=c()
for (j in 1:p) {
  pp_Tn[j]=1-p_wschisq(T_n[j],weights = value_BB,dfs = rep(1,R))
  pp_tilTn[j]=1-p_wschisq(tilT_n[j],weights = value_BB,dfs = rep(1,R))
}
Sigma_xx=c()
for (j in 1:p) {
  Sigma_xx[j]=sd(x[,j])
}
###CASH###
Rejection_CA=CAFDR(T_n, tilT_n, pp_Tn, pp_tilTn, Sigma_xx, alpha=0.15, tau = 0.5)
Rejection_CA
###Knockoff###
alphakn_threshold <- function(W, fdr, offset) {
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}
x_kn=create.second_order(x,method = "equi") 
T_n2=c()
for (j in 1:p) {
  T_n2[i]=vm(x_kn[,j],y,R,p_y,y_ca)
}
W_sta=abs(T_n)-abs(T_n2)
thr=alphakn_threshold(W_sta, 0.15, 1)
Rejection_kn<-which(W_sta>=thr)
Rejection_kn
###BH###
sorted_pvalues <- sort(pp_Tn, decreasing = F, index.return = T)
cutoff <- max(which(sorted_pvalues$x <= (1:p)*0.15 / p))
Rejection_BH <- sorted_pvalues$ix[1:cutoff]
Rejection_BH
###IHW###
ihw_fdr <- ihw(pp_Tn, Sigma_xx, 0.15)
Rejection_IHW=which(rejected_hypotheses(ihw_fdr)==TRUE)
Rejection_IHW
###AdaDetect###
f0 <- c(0,1)
score.pool = group.func(T_n,tilT_n,f0)
pools = score.pool$s
poolst = score.pool$st
Rejection_ada= Adadetect(pools,poolst,0.15)
Rejection_ada

#############Simulation for R=5###########
sim=1000
fdr1 <- c()
fdr2 <- c()
fdr3 <- c()
fdr4 <- c()
fdr5 <- c()
pow1 <- c()
pow2 <- c()
pow3 <- c()
pow4 <- c()
pow5 <- c()
for (m in 1:sim) {
		print(m)
		n = 100
		p = 1000
		l = 50  
		R = 5
		c = 0.25 #0.5,0.75,1,1.25,1.5
		seed=123
		data <- generate_data(c,"normal",n,p,l,R,seed+m)
		x <- data$X
		y <- data$Y
		truesignal <- data$true_signal
		#############CASH###############
		y_ca=rep(1:5)
		p_y=sapply(1:R,function(i) length(which(y==y_ca[i]))/n)
		T_n=c()
		for (j in 1:p) {
			T_n[j]=vm(x[,j],y,R,p_y,y_ca)
		}
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
		tilT_n=r_wschisq(p,weights = value_BB,dfs = rep(1,R))
		pp_Tn=c()
		pp_tilTn=c()
		for (j in 1:p) {
			pp_Tn[j]=1-p_wschisq(T_n[j],weights = value_BB,dfs = rep(1,R))
			pp_tilTn[j]=1-p_wschisq(tilT_n[j],weights = value_BB,dfs = rep(1,R))
		}
		Sigma_xx=c()
		for (j in 1:p) {
			Sigma_xx[j]=sd(x[,j])
		}
		rjCAFDR <- CAFDR(T_n, tilT_n, pp_Tn, pp_tilTn, Sigma_xx, alpha=0.15, tau = 0.5)
		metrics_cafdr <- compute_fdr_power(rjCAFDR, truesignal)
		fdr1[m] <- metrics_cafdr$fdr
		pow1[m] <- metrics_cafdr$power
		###Knockoff###
		alphakn_threshold <- function(W, fdr, offset) {
			ts = sort(c(0, abs(W)))
			ratio = sapply(ts, function(t)
				(offset + sum(W <= -t)) / max(1, sum(W >= t)))
			ok = which(ratio <= fdr)
			ifelse(length(ok) > 0, ts[ok[1]], Inf)
		}
		x_kn=create.second_order(x,method = "equi") 
		T_n2=c()
		for (j in 1:p) {
			T_n2[i]=vm(x_kn[,j],y,R,p_y,y_ca)
		}
		W_sta=abs(T_n)-abs(T_n2)
		thr=alphakn_threshold(W_sta, 0.15, 1)
		Rejection_kn<-which(W_sta>=thr)
		metrics_knockoff <- compute_fdr_power(Rejection_kn, truesignal)
		fdr2[m] <- metrics_knockoff$fdr
		pow2[m] <- metrics_knockoff$power
		###BH###
		sorted_pvalues <- sort(pp_Tn, decreasing = F, index.return = T)
		cutoff <- max(which(sorted_pvalues$x <= (1:p)*0.15 / p))
		Rejection_BH <- sorted_pvalues$ix[1:cutoff]
		metrics_bh <- compute_fdr_power(Rejection_BH, truesignal)
		fdr3[m] <- metrics_bh$fdr
		pow3[m] <- metrics_bh$power
		###IHW###
		ihw_fdr <- ihw(pp_Tn, Sigma_xx, 0.15)
		Rejection_IHW=which(rejected_hypotheses(ihw_fdr)==TRUE)
		metrics_ihw <- compute_fdr_power(Rejection_IHW, truesignal)
		fdr4[m] <- metrics_ihw$fdr
		pow4[m] <- metrics_ihw$power
		###AdaPT###
		f0 <- c(0,1)
		score.pool = group.func(T_n,tilT_n,f0)
		pools = score.pool$s
		poolst = score.pool$st
		Rejection_ada= Adadetect(pools,poolst,0.15)
		metrics_ada <- compute_fdr_power(Rejection_ada, truesignal)
		fdr5[m] <- metrics_ada$fdr
		pow5[m] <- metrics_ada$power
}
avg_fdr1 <- mean(fdr1, na.rm = TRUE)
avg_fdr2 <- mean(fdr2, na.rm = TRUE)
avg_fdr3 <- mean(fdr3, na.rm = TRUE)
avg_fdr4 <- mean(fdr4, na.rm = TRUE)
avg_fdr5 <- mean(fdr5, na.rm = TRUE)
avg_pow1 <- mean(pow1, na.rm = TRUE)
avg_pow2 <- mean(pow2, na.rm = TRUE)
avg_pow3 <- mean(pow3, na.rm = TRUE)
avg_pow4 <- mean(pow4, na.rm = TRUE)
avg_pow5 <- mean(pow5, na.rm = TRUE)