#rm(list=ls())
setwd("/Users/wangl29/Dropbox/Research/Work with Li/Life gained from intervention/Codes")
library(survival)
library(survminer)
library(survey)
library(matrixStats)
library(MASS)
library(data.table)
source("subfunctions.R")
seed = read.table("seed.txt", header = T)
seed1 = seed[,1]
seed2 = seed[,2]



################### # Generate a finite population with size N=200,000 ####################  
#n_simu=500
#pop_beta1 = matrix(1, n_simu, 3)
#pop_beta2 = matrix(1, n_simu, 3)
#rt = rep(0, n_simu)
#seed=runif(n_simu, min=0, max=10000)
#simu=1
#for(simu in simu:n_simu){
#set.seed(seed[simu])
#path = "/Volumes/wangl29/Life_Intervention/rlts/t5_gamma0/"
#set.seed(4837.442); N=5e5; t1=5; gamma2 = 0
#path = "/Volumes/wangl29/Life_Intervention/rlts/t12_gamma0/"
#set.seed(3834.606); N=5e5; t1=12; gamma2 = 0
#path = "/Volumes/wangl29/Life_Intervention/rlts/t46_gamma0/"
#set.seed(2728.588); N=5e5; t1=46; gamma2 = 0
#path = "/Volumes/wangl29/Life_Intervention/rlts/t5_gamma_0057/"
#set.seed(9838.613); N=5e5; t1=5; gamma2=-0.0825 
path = "/Volumes/wangl29/Life_Intervention/rlts/t5_gamma0057/"
set.seed(8265.993); N=5e5; t1=5; gamma2=0.0825  


# Finite population data generation
c1 = 7  # end of follow-up for survey
c2 = 15 # end of follow-up for trial

# Log-HR's for the two events
# event 1
beta1 = c(log(-log(.95)/c2),-0.2, 0.3); gamma1 = -0.8 # intervention effect on event 1 of interest
#event 2
beta2 = c(log(-log(.7)/c2), -0.1, 0.2); 

# Three covariates 
x1 = rnorm(N, sd=2)
x2=rchisq(N, df=1.3)
# Prevention service variable
y = as.numeric(runif(N)>0.5)

# HR's 
theta1_y0 = exp(cbind(1, x1, x2)%*%c(beta1)) #event 1, y=0
theta2_y0 = exp(cbind(1, x1, x2)%*%c(beta2)) #event 2, y=0
theta1_y = theta1_y0*exp(gamma1*y) #event 1, actual y
theta2_y = theta2_y0*exp(gamma2*y) #event 2, actual y

theta_y = cbind(theta1_y, theta2_y)
theta_y0 = cbind(theta1_y0, theta2_y0)

# Generate time to event (with time-varying effect of the intervention)
alpha = 1.2 # time varying parameter 
# Using formula (11)
u = runif(N)
# before intervention effect disappear
t_d.l = (-log(u)/rowSums(theta_y))^(1/alpha)
# after intervention effect disappear
t_d.h = ((-log(u)-t1^alpha*rowSums(theta_y-theta_y0))/rowSums(theta_y0))^(1/alpha)
# threshold
thred = t1^alpha*rowSums(theta_y)
t_d = t_d.h
t_d[(-log(u)-thred)<=0]=t_d.l[(-log(u)-thred)<=0]
quantile(t_d, probs = seq(0,1,0.1))
# combined baseline hazard
lambda = alpha*t_d^(alpha-1)*theta_y
lambda[t_d>t1] = (alpha*t_d^(alpha-1)*theta_y0)[t_d>t1,]
rnd.d = runif(N)
# indicator of event 1
d1 = as.numeric(rnd.d<=lambda[,1]/rowSums(lambda))
d1.gamma=d1; d1.gamma[t_d>t1]=0
# indicator of event 2
d2 = 1-d1
d2.gamma=d2; d2.gamma[t_d>t1]=0
rm(lambda)
######################### event time with no intervention ######################### 
t_d0 = (-log(u)/rowSums(theta_y0))^(1/alpha)
lambda_y0 = alpha*t_d0^(alpha-1)*theta_y0
d1_y0 = as.numeric(rnd.d<=lambda_y0[,1]/rowSums(lambda_y0))
d2_y0 = 1-d1_y0
rm(lambda_y0)
################################################################
# Censoring
# Survey sample
c_s = rweibull(N, shape=1, scale=-c2/log(0.8))
c_s[c_s>=c2]=c2
t_s  = apply(cbind(t_d, c_s),  1, FUN=min)
t_s0 = apply(cbind(t_d0, c_s), 1, FUN=min)
d1_s  = d1; d2_s= d2
d1_s[t_s<t_d]=0; d2_s[t_s<t_d]=0; 
d1_s0 = d1_y0; d2_s0= d2_y0
d1_s0[t_s0<t_d0]=0; d2_s0[t_s0<t_d0]=0; 
# Trial
c_t = rweibull(N, shape=1, scale=-c1/log(0.8))
c_t[c_t>=c1]=c1
t_t = apply(cbind(t_d, c_t), 1, FUN=min)
d1.gamma_t = d1.gamma
d1.gamma_t[t_t<t_d]=0
d2.gamma_t = d2.gamma
d2.gamma_t[t_t<t_d]=0
pop = data.frame(x1, x2,  y, t_d, t_d0, d1, d2, d1_y0, d2_y0, d1.gamma,
                 t_s, d1_s, d2_s,t_s0, d1_s0, d2_s0, 
                 t_t, d1.gamma_t,  d2.gamma_t, w=1)
# Check distributions of time variables
quantile(t_d,        probs = seq(0.1,0.9,0.1))
quantile(t_d[y==0],  probs = seq(0.1,0.9,0.1))
quantile(t_d[y==1],  probs = seq(0.1,0.9,0.1))
quantile(t_d0,       probs = seq(0.1,0.9,0.1))
quantile(t_d0[y==1], probs = seq(0.1,0.9,0.1))
round(quantile(x1, probs = seq(0.1,0.9,0.1)),4)
round(quantile(x1[y==0], probs = seq(0.1,0.9,0.1)),4)

# Cox regression model
# Full models
fm_fit.cox1       = as.formula("Surv(t_d, d1) ~ x1+x2+y")
fm_fit.cox1.beta  = as.formula("Surv(t_d, d1) ~ x1+x2")
fm_fit.cox2.beta  = as.formula("Surv(t_d, d2) ~ x1+x2")
fm_fit.cox10.beta = as.formula("Surv(t_d0, d1_y0) ~ x1+x2")
fm_fit.cox20.beta = as.formula("Surv(t_d0, d2_y0) ~ x1+x2")
fm_fit.cox1.gamma = as.formula("Surv(t_d, d1.gamma) ~ y")
fm_fit.cox2.gamma = as.formula("Surv(t_d, d2.gamma) ~ y")


# Finite population parameters
# no censoring
cox_pop1          = coxph(fm_fit.cox1,       data = pop)
cox_pop1.beta     = coxph(fm_fit.cox1.beta,  data = pop)
cox_pop2.beta     = coxph(fm_fit.cox2.beta,  data = pop)
cox_pop1_d0.beta  = coxph(fm_fit.cox10.beta, data = pop)
cox_pop2_d0.beta  = coxph(fm_fit.cox20.beta, data = pop)
cox_pop1_y0.beta  = coxph(fm_fit.cox1.beta,  data = pop[y==0,])
cox_pop2_y0.beta  = coxph(fm_fit.cox2.beta,  data = pop[y==0,])
cox_pop1.gamma  = coxph(fm_fit.cox1.gamma,   data = pop)
cox_pop1$coefficients
rbind(cox_pop1.beta$coefficients,
      cox_pop2.beta$coefficients)
rbind(cox_pop1_d0.beta$coefficients,
      cox_pop2_d0.beta$coefficients)
rbind(cox_pop1_y0.beta$coefficients,
      cox_pop2_y0.beta$coefficients)
pop_beta1 = c(cox_pop1.beta$coefficients, cox_pop1.gamma$coefficients);pop_beta1

if(gamma2!=0){
  cox_pop2.gamma = coxph(fm_fit.cox2.gamma, data = pop)
  pop_beta2 = c(cox_pop2.beta$coefficients, cox_pop2.gamma$coefficients);pop_beta2
}else{
  pop_beta2 = c(cox_pop2.beta$coefficients, 0);pop_beta2
}

# Unique event times in the population
t_star = quantile(t_d, probs=seq(0.1,0.9,0.1))
# cumulative baseline hazards
cbind(cum_bsln_hzd(surv.fit=cox_pop1,         t_star=t_star),
      cum_bsln_hzd(surv.fit=cox_pop1.beta,    t_star=t_star, 
                   rel_hzd = exp(model.matrix(cox_pop1)%*%c(cox_pop1$coefficients))),
      cum_bsln_hzd(surv.fit=cox_pop1.beta,    t_star=t_star),
      cum_bsln_hzd(surv.fit=cox_pop1_d0.beta, t_star=t_star),
      cum_bsln_hzd(surv.fit=cox_pop1_y0.beta, t_star=t_star))


######################################################################################
# Population parameters with a cutoff time
# A set of variables of x1 and x2 for absolute risk prediction
x0=cbind(quantile(x1, probs=c(.10, .25, .5, .75, .90)), quantile(x2, probs=c(.10, .25, .5, .75, .90)))
# cutoff time 
t_tau=c(7, 15)
t_star.c = sort(pop$t_s0[((pop$d1_s0==1)|(pop$d2_s0==1))&(pop$y==0)])

######################################################################################
# Model fitted to the trial
fm_fit.cox.t1 = as.formula("Surv(t_t, d1.gamma_t) ~ y")
fm_fit.cox.t2 = as.formula("Surv(t_t, d2.gamma_t) ~ y")
# Model fitted to the survey sample
fm_fit.cox.s1 = as.formula("Surv(t_s0, d1_s0) ~ x1+x2")
fm_fit.cox.s2 = as.formula("Surv(t_s0, d2_s0) ~ x1+x2")

cox_pop.s1= coxph(fm_fit.cox.s1,  data = pop)
cox_pop.s2= coxph(fm_fit.cox.s2,  data = pop)
cox_pop.t1= coxph(fm_fit.cox.t1 , data = pop)
pop_beta1.c = c(cox_pop.s1$coeff, cox_pop.t1$coeff)
if(gamma2!=0){
  cox_pop.t2 = coxph(fm_fit.cox.t2, data = pop)
  pop_beta2.c = c(cox_pop.s2$coeff, cox_pop.t2$coeff)
}else{
  pop_beta2.c = c(cox_pop.s2$coeff, 0)
}



date()
bh=basehaz(cox_pop.s1, centered = FALSE)          #baseline cumulative hazard function for all covariates set to 0 (centered = false
bh_y0=basehaz(cox_pop1_y0.beta, centered = FALSE) #baseline cumulative hazard function for all covariates set to 0 (centered = false
# baseline hazards among the untreated population and the whole population
plot(bh[bh[,2]%in%bh_y0[,2],1],bh_y0[bh_y0[,2]%in%bh[,2],1],
     xlab="Whole population", 
     ylab="Untreated popualtion y=0",
     main="Baseline hazard in the whole Pop vs. y=0 Pop")
abline(0,1,col="red")
  # Baseline cumulative hazard at the time to event 1
  popLambda0_s1=cum_bsln_hzd(surv.fit=cox_pop.s1, t_star=t_star.c)
  popLambda0_s2=cum_bsln_hzd(surv.fit=cox_pop.s2, t_star=t_star.c)
  popLambda0.s = cbind(popLambda0_s1, popLambda0_s2)
  rm(popLambda0_s1, popLambda0_s2)

  Delta_mu_x0.s = Delta_mu_inf(Lambda0=popLambda0.s, t=t_star.c, t1, x =x0, 
                               betas = cbind(pop_beta1.c, pop_beta2.c))
  Delta_mu_tau1.s = colSums(Delta_mu_x0.s$Delta_mu[t_star.c<=t_tau[1],])
  Delta_mu_tau2.s = colSums(Delta_mu_x0.s$Delta_mu[t_star.c<=t_tau[2],])
  rm(Delta_mu_x0.s, popLambda0.s)
  pop_Delta_mu = rbind(Delta_mu_tau1.s,Delta_mu_tau2.s)
  #write.table(pop_Delta_mu, paste0(path, "pop_Delta_mu.txt"), sep = ",", row.names = F, col.names = F)
  pop_param = list(NSIMU_tot = NSIMU_tot,
                   NSIMU = NSIMU,
                   x0 = x0,
                   betas = rbind(beta1, beta2),
                   pop_beta = rbind(pop_beta1, pop_beta2),
                   pop_Delta_mu = pop_Delta_mu,
                   t_tau = t_tau
  )  #/rbind(Delta_mu_tau1,Delta_mu_tau2)
  date()
  
  rm(cox_pop1, cox_pop.s1, cox_pop.s2, cox_pop.t1, cox_pop1_d0.beta, cox_pop1_y0.beta,
     cox_pop1.beta, cox_pop2.beta, cox_pop2_d0.beta, cox_pop2_y0.beta, cox_pop1.gamma)
  save(pop_param, file=paste0(path, 'pop_param_2.RData'))
  
####################   Selection probabilities for trial and survey sample ####################  
  n_t =  2000 # sample size of the trial
  n_s = 8000 # sample size of the survey
  fs = n_s/N;fs*100 # sample fraction of the trial
  ft = n_t/N;ft*100 # sample fraction of the survey
    
  # Selection probabilities for the survey sample (PPS, with size "odds.s"), depending on x1 and x2
  gamma_s = c(0.1, -0.07)
  odds.s = exp(as.matrix(pop[,c("x1", "x2")])%*%c(gamma_s));
  pi.s = odds.s*N*fs/sum(odds.s)
  max(abs(pi.s))/min(abs(pi.s)) #poisson sampling with identified intercept
  quantile(pi.s)
  
  #Selection probabilities for the trial (PPS, with size "odds.c"), depending on x1 and x2
  gamma_t = c(-0.1, 0.07)
  odds.t = exp(as.matrix(pop[,c("x1", "x2")])%*%c(gamma_t));
  pi.t = odds.t*N*ft/sum(odds.t)
  max(abs(pi.t))/min(abs(pi.t)) #poisson sampling with identified intercept
  quantile(pi.t)
  
  ######################################## Simulations ######################################## 
  NSIMU = 10 # number of simulation runs 
  # Simulation parameters
  skip_simu = NULL
  # Number of simulation runs
  # Matrices for inference of Delta_mu and relative hazards storage
  Delta_mu_est = array(0, c(NSIMU, length(t_tau),nrow(x0)))
  Delta_mu_var = array(0, c(NSIMU, length(t_tau),nrow(x0)))
  if (gamma2==0){
    beta_est = matrix(0, NSIMU, 5)
    beta_var = matrix(0, NSIMU, 5)
  }else{
    beta_est = matrix(0, NSIMU, 6)
    beta_var = matrix(0, NSIMU, 6)
  }
  
  # Start of the simulation
  simu=1
  start.time=date()
  for(simu in c(simu:NSIMU)){
    # Selecting a trial
    set.seed(seed1[simu])
    samp.t = sam.pps(popul=pop, Msize = odds.t, n=n_t)
    # Selecting a survey sample    
    set.seed(seed2[simu])
    samp.s = sam.pps(popul=pop, Msize = odds.s, n=n_s)
    # Fit Cox regression models to the samples for events 1 and 2
    # Trial
    cox_t1 = coxph(fm_fit.cox.t1, data = samp.t, robust=T, ties="breslow")#;cox_c.d1
    # Survey sample
    ds.s = svydesign(ids=~1, data=samp.s, weights=~wt) # Survey design
    cox_s.d1 = svycoxph(fm_fit.cox.s1, ds.s)#;cox_s.d1
    cox_s.d2 = svycoxph(fm_fit.cox.s2, ds.s)#;cox_s.d2
    # Log-relative hazards
    if(gamma2==0){
      beta_est[simu,] = c(cox_s.d1$coefficients, cox_t1$coefficients, 
                          cox_s.d2$coefficients)
      beta.simu = cbind(c(cox_s.d1$coefficients, cox_t1$coefficients), 
                        c(cox_s.d2$coefficients, 0))
      beta_w = array(0, c(n_s+n_t, 3, 2))
      beta_w.out.s1 = beta_pw.cox(x.mtrx = model.matrix(cox_s.d1), rel_hzd = exp(cox_s.d1$linear.predictors), 
                                  dat = samp.s, pw = "wt", t="t_s0", d="d1_s0")
      beta_w.out.s2 = beta_pw.cox(x.mtrx = model.matrix(cox_s.d2), rel_hzd = exp(cox_s.d2$linear.predictors), 
                                  dat = samp.s, pw = "wt", t="t_s0", d="d2_s0")
      beta_w.out.t1 = beta_pw.cox(x.mtrx = model.matrix(cox_t1), rel_hzd = exp(cox_t1$linear.predictors), 
                                  dat = samp.t, pw = "wt", t="t_t", d="d1.gamma_t")
      
      beta_var[simu,]=c(diag(var(samp.s$wt*beta_w.out.s1$beta_pw)*(n_s-1)),
                        diag(var(samp.t$wt*beta_w.out.t1$beta_pw)*(n_t-1)),
                        diag(var(samp.s$wt*beta_w.out.s2$beta_pw)*(n_s-1)))
      beta_w[1:n_s,1:2,1] = beta_w.out.s1$beta_pw
      beta_w[1:n_s,1:2,2] = beta_w.out.s2$beta_pw
      beta_w[(1:n_t+n_s),3,1] = beta_w.out.t1$beta_pw
      
      
    }else{
      cox_t2 = coxph(fm_fit.cox.t2, data = samp.t, robust=T, ties="breslow")#;cox_c.d1
      beta_est[simu,] = c(cox_s.d1$coefficients, cox_t1$coefficients, 
                          cox_s.d2$coefficients, cox_t2$coefficients)
      beta.simu = cbind(c(cox_s.d1$coefficients, cox_t1$coefficients), 
                        c(cox_s.d2$coefficients, cox_t2$coefficients))
      beta_w = array(0, c(n_s+n_t, 3, 2))
      beta_w.out.s1 = beta_pw.cox(x.mtrx = model.matrix(cox_s.d1), rel_hzd = exp(cox_s.d1$linear.predictors), 
                                  dat = samp.s, pw = "wt", t="t_s0", d="d1_s0")
      beta_w.out.s2 = beta_pw.cox(x.mtrx = model.matrix(cox_s.d2), rel_hzd = exp(cox_s.d2$linear.predictors), 
                                  dat = samp.s, pw = "wt", t="t_s0", d="d2_s0")
      beta_w.out.t1 = beta_pw.cox(x.mtrx = model.matrix(cox_t1), rel_hzd = exp(cox_t1$linear.predictors), 
                                  dat = samp.t, pw = "wt", t="t_t", d="d1.gamma_t")
      beta_w.out.t2 = beta_pw.cox(x.mtrx = model.matrix(cox_t2), rel_hzd = exp(cox_t2$linear.predictors), 
                                  dat = samp.t, pw = "wt", t="t_t", d="d2.gamma_t")
      
      beta_var[simu,]=c(diag(var(samp.s$wt*beta_w.out.s1$beta_pw)*(n_s-1)),
                        diag(var(samp.t$wt*beta_w.out.t1$beta_pw)*(n_t-1)),
                        diag(var(samp.s$wt*beta_w.out.s2$beta_pw)*(n_s-1)),
                        diag(var(samp.t$wt*beta_w.out.t2$beta_pw)*(n_t-1)))
      beta_w[1:n_s,1:2,1] = beta_w.out.s1$beta_pw
      beta_w[1:n_s,1:2,2] = beta_w.out.s2$beta_pw
      beta_w[(1:n_t+n_s),3,1] = beta_w.out.t1$beta_pw
      beta_w[(1:n_t+n_s),3,2] = beta_w.out.t2$beta_pw
    }
    # Unique event times in the survey sample
    t_star.s = sort(samp.s$t_s0[(samp.s$d1_s0+samp.s$d2_s0)==1])
    
    # Baseline hazard at the event times for events 1 and 2, and their derivative  w.r.t. weights 
    lambda.s1_w = lambda_w(surv.fit=cox_s.d1, beta_wt=beta_w.out.s1$beta_pw)
    lambda.s2_w = lambda_w(surv.fit=cox_s.d2, beta_wt=beta_w.out.s2$beta_pw)
    # Cumulative baseline hazard at the event times for events 1 and 2, and their derivative  w.r.t. weights 
    Lambda0.s1_w = Lambda_w(lambda_out=lambda.s1_w, t_star.s)  
    Lambda0.s2_w = Lambda_w(lambda_out=lambda.s2_w, t_star.s)
    rm(lambda.s1_w, lambda.s2_w)
    
    S_Delta_mu = Delta_mu_inf(Lambda0=cbind(Lambda0.s1_w$Lambda, Lambda0.s2_w$Lambda), 
                              t=t_star.s, t1, x = x0, 
                              betas = beta.simu,
                              Lambda0_w = array(c(Lambda0.s1_w$Lambda_wt, Lambda0.s2_w$Lambda_wt),
                                                c(n_s, length(t_star.s), 2)), 
                              beta_w=beta_w, deviate=T)
    
    rm(Lambda0.s1_w, Lambda0.s2_w)
    # Estimate of LYG
    Delta_mu_est[simu,1,] = colSums(S_Delta_mu$Delta_mu[t_star.s<=t_tau[1],]) # at t=7
    Delta_mu_est[simu,2,] = colSums(S_Delta_mu$Delta_mu[t_star.s<=t_tau[2],]) # at t=15
    Delta_mu_w1 = sapply(1:nrow(x0), function(i)rowSums(S_Delta_mu$Delta_mu_w[,t_star.s<=t_tau[1],i]))
    Delta_mu_w2 = sapply(1:nrow(x0), function(i)rowSums(S_Delta_mu$Delta_mu_w[,t_star.s<=t_tau[2],i]))
    # Variance of LYG
    Delta_mu_var[simu,1,] = diag(var(c(samp.s$wt, samp.t$wt)*Delta_mu_w1)*(n_s+n_t-1))  # at t=7
    Delta_mu_var[simu,2,] = diag(var(c(samp.s$wt, samp.t$wt)*Delta_mu_w2)*(n_s+n_t-1))  # at t=15
    rm(S_Delta_mu)
    print(simu)
  }   # End of the simualtion
  end.time = date()
  start.time; end.time
  
  

  #########################################  Simulation Inference ##############################     
  if(simu<NSIMU) simu=simu-1
  # Simulation mean of the Cox regression coefficients
  apply(beta_est[1:simu,], 2,mean)
  apply(beta_var[1:simu,], 2, mean)/apply(beta_est[1:simu,], 2, var)
  # Simulation mean of Delta_mu
  Delta_mu_est.mxt=Delta_mu_est
  dim(Delta_mu_est.mxt)=c(NSIMU, length(t_tau)*nrow(x0))
  matrix(apply(Delta_mu_est.mxt[1:simu,], 2, mean), 2,)/rbind(Delta_mu_tau1.s,Delta_mu_tau2.s)
  # Variances of Delta_mu
  apply(Delta_mu_est[1:simu,1,], 2, var)  # Empirical variances (simulation variance)
  apply(Delta_mu_var[1:simu,1,], 2, mean) # Analytical vairanes (simulation mean of the TL variance estimates)
  apply(Delta_mu_est[1:simu,2,], 2, var)  # Empirical variances (simulation variance)
  apply(Delta_mu_var[1:simu,2,], 2, mean) # Analytical vairanes (simulation mean of the TL variance estimates)
  
  apply(Delta_mu_var[1:simu,1,], 2, mean)/apply(Delta_mu_est[1:simu,1,], 2, var) # Variance ratio
  apply(Delta_mu_var[1:simu,2,], 2, mean)/apply(Delta_mu_est[1:simu,2,], 2, var) # Variance ratio
  




