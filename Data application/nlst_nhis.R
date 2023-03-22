#rm(list=ls())

#Use NLST with long-term follow-up to quantify the effect of screening on mortality
setwd("/Users/wangl29/Dropbox/Research/Codes/Data application")
library(survival)
library(dynsurv)
library(survey)
library(survminer)
library(matrixStats)
library(MASS)
library(data.table)

source("subfunctions2.R")
# Read nlst data
nlst <- read.csv("prsn.csv")
nlst <- subset(nlst,mortality_exitdays_2015>0,select=c(pid,rndgroup,mortality_exitdays_2015,finaldeathlc_2015,gender))
nlst$lcd <- ifelse(is.na(nlst$finaldeathlc_2015)==0 & nlst$finaldeathlc_2015==1,1,0)
nlst$othcd <- ifelse(is.na(nlst$finaldeathlc_2015)==0 & nlst$finaldeathlc_2015==0,1,0)
nlst$ct <- ifelse(nlst$rndgroup==1,1,0)

nlst$mtlty_yr = ceiling(nlst$mortality_exitdays_2015/365.25)
nlst$lcd_y6 = nlst$lcd
nlst$lcd_y6[nlst$mtlty_yr>6]=0
fit_gamma <- coxph(Surv(mtlty_yr,lcd_y6)~ct,data=nlst, robust=T, ties="breslow")
nlst$wt=1
beta_w.out.t = beta_pw.cox(x.mtrx  = model.matrix(fit_gamma), 
                           rel_hzd = exp(model.matrix(fit_gamma)%*%fit_gamma$coefficients), 
                           dat     = nlst, pw = "wt", t="mtlty_yr", d="lcd_y6")

n_t = nrow(beta_w.out.t$beta_pw)

fit_gamma2 <- coxph(Surv(mtlty_yr,othcd)~ct,data=nlst, robust=T, ties="breslow")


load(file="nhis97_01.imputed.1.RData")

#cause-specific death data to the end of 2006 (nhis$deathage2, nhis$lung.cancer.death)
#overall mortality data to 2015
nhis$other.death <- ifelse(nhis$died==1 & nhis$deathyear2<2007 & nhis$lung.cancer.death==0,1,0)
nhis$died10 <- ifelse(nhis$died==1 & nhis$deathyear2<2007,1,0)
nhis$flwup <- nhis$deathyear-nhis$year
nhis$flwup_yr = ceiling(nhis$flwup)
nhis <- subset(nhis,is.na(wt_mort5)==0)
master  <- svydesign(id=~psu, strata=~strata, weights=~wt_mort5, data=nhis, nest=TRUE)
#ever-smokers with known history, 40-84 years old, no prior lung cancers
design <- subset(master, analypop==1 & age>=40&age<=84 & lung.cancer.before==0)
design2 <- subset(master, analypop==1 & age>=40&age<=84 & lung.cancer.before==0 & is.na(deathage2)==0) #removed missing cause-specific mortality
elig_nhis = design2$variables
#design2=update(design2, age.c=cut(age,c(39, seq(45,85,10))))
#elig_nhis = design2$variables

#LCDRAT model

# try different functions of age to minimize the AIC
modlcd <- svycoxph(Surv(flwup_yr, lung.cancer.death) ~ 
                     age+female+race+edu6+emp+I(bmi <= 18.5)+I(cpd > 20)+
                     I(pkyr.cat>=30 & pkyr.cat<40)+I(pkyr.cat>=40 & pkyr.cat<50)+I(pkyr.cat>=50)+ 
                     I(log(bmi))+I(log(qtyears+1))+I(log(smkyears)), 
                   design=design2, model = TRUE, x=TRUE, y = TRUE) 
#AIC(modlcd)
rr=exp(model.matrix(modlcd)%*%modlcd$coefficients)
design2=update(design2, rr = rr)
svyquantile(~rr, design2, c(.1, .5, .9))
wt_p = sum(elig_nhis$wt_mort5)*c(.1, .5, .9)
q.indx = apply(abs(outer(cumsum(elig_nhis$wt_mort5[order(rr)]),wt_p, FUN="-")), 2, which.min)
sort(rr)[q.indx]

betas.modlcd = cbind(modlcd$coefficients, as.matrix(summary(modlcd)$conf.int[,-2]))
beta_w.out.s1 = beta_pw.cox(x.mtrx = model.matrix(modlcd), 
                            rel_hzd = exp(model.matrix(modlcd)%*%modlcd$coefficients), 
                            dat = elig_nhis, pw = "wt_mort5", t="flwup_yr", d="lung.cancer.death")
n_s = nrow(beta_w.out.s1$beta_pw)
n_beta1 = nrow(betas.modlcd)

#overall mortality model fit to other cause mortality
modoth <- svycoxph(Surv(flwup_yr, other.death) ~ 
                     age+female + race + edu6 + 
                     I(bmi <= 18.5)+I(bmi > 18.5 & bmi <= 20)+I(bmi > 25 & bmi <= 30)+I(bmi > 30 & bmi <= 35)+I(bmi > 35)+ 
                     emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                     I(log(qtyears + 1)) + I(log(cpd)) + I(sqrt(pkyr.cat)),
                   design=design2, model = TRUE, x=TRUE, y = TRUE)
betas.modoth = cbind(modoth$coefficients, as.matrix(summary(modoth)$conf.int[,-2]))
n_beta2 = nrow(betas.modoth)

#x0=rbind(55, 1, 0, 0, 0, 1, 1, 0, )
beta_w.out.s2 = beta_pw.cox(x.mtrx = model.matrix(modoth), 
                            rel_hzd = exp(model.matrix(modoth)%*%modoth$coefficients), 
                            dat = elig_nhis, pw = "wt_mort5", t="flwup_yr", d="other.death")

lambda.s1_w = lambda_w(surv.fit=modlcd, beta_wt=beta_w.out.s1$beta_pw)
lambda.s2_w = lambda_w(surv.fit=modoth, beta_wt=beta_w.out.s2$beta_pw)
# Cumulative baseline hazard at the event times for events 1 and 2, and their derivative  w.r.t. weights 
t_star.s=sort(unique(elig_nhis$flwup_yr))
Lambda0.s1_w = Lambda_w(lambda_out=lambda.s1_w, t_star.s)  
Lambda0.s2_w = Lambda_w(lambda_out=lambda.s2_w, t_star.s)

betas.modlcd_all = betas.modlcd[,1]
betas.modoth_all = betas.modoth[,1]

betas.modlcd_all[setdiff(row.names(betas.modoth), row.names(betas.modlcd))] <- 0
betas.modoth_all[setdiff(row.names(betas.modlcd), row.names(betas.modoth))] <- 0

betas = cbind(c(betas.modlcd_all, fit_gamma$coefficients),
              c(betas.modoth_all, 0))
n_beta = nrow(betas)


beta_w = array(0, c(n_s+n_t, n_beta, 2))
beta_w[1:n_s,1:(n_beta-1),1] = cbind(beta_w.out.s1$beta_pw, matrix(0,n_s, (n_beta-1-n_beta1)))
beta_w[1:n_s,1:(n_beta-1),2] = cbind(beta_w.out.s2$beta_pw, matrix(0,n_s, (n_beta-1-n_beta2)))[,mtch.order]
beta_w[(1:n_t+n_s),n_beta,1] = beta_w.out.t$beta_pw

nhis_ds = elig_nhis[,c("psu", "strata", "wt_mort5")]

#overall mortality model fit to other cause mortality

modall <- svycoxph(Surv(flwup_yr, died) ~ 
                     age+female + race + edu6 + 
                     I(bmi <= 18.5)+I(bmi > 18.5 & bmi <= 20)+I(bmi > 25 & bmi <= 30)+I(bmi > 30 & bmi <= 35)+I(bmi > 35)+ 
                     emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                     I(log(qtyears + 1)) + I(log(cpd)) + I(sqrt(pkyr.cat)),
                   design=design2, model = TRUE, x=TRUE, y = TRUE)
betas.modall = cbind(modall$coefficients, summary(modall)$conf.int[,-2])
#round(cbind(betas.modoth, betas.modall10, betas.modall), 3)
beta_w.out.s = beta_pw.cox(x.mtrx = model.matrix(modall), 
                           rel_hzd = exp(model.matrix(modall)%*%modall$coefficients), 
                           dat = elig_nhis, pw = "wt_mort5", t="flwup_yr", d="died")
lambda.s_w = lambda_w(surv.fit=modall, beta_wt=beta_w.out.s$beta_pw)
Lambda0.s_w = Lambda_w(lambda_out=lambda.s_w, t_star.s)

betas.modall_all = betas.modall[,1]
betas.modall_all[setdiff(row.names(betas.modlcd), row.names(betas.modall))] <- 0
betas = cbind(betas,
              c(betas.modall_all, 0))

mtch.order = match(names(betas.modlcd_all), names(betas.modoth_all))
hr_ci = data.frame(beta_lcd=betas.modlcd_all, hr_lcd=0, hr_lcd.lw=0, hr_lcd.up=0,
                   beta_oth=betas.modoth_all[mtch.order], hr_oth=0, hr_oth.lw=0, hr_oth.up=0,
                   beta_all=betas.modall_all[mtch.order], hr_all=0, hr_all.lw=0, hr_all.up=0
)

mtch.order = match(row.names(hr_ci), row.names(betas.modoth))
mtch.order = mtch.order[!is.na(mtch.order)]
hr_ci[which(hr_ci$beta_lcd!=0),2:4]=betas.modlcd[,-1]
hr_ci[which(hr_ci$beta_oth!=0),6:8]=betas.modoth[mtch.order,-1]
hr_ci[which(hr_ci$beta_all!=0),10:12]=betas.modall[mtch.order,-1]

hr_ci=rbind(hr_ci,ct = c(fit_gamma$coefficients, summary(fit_gamma)$conf.in[-2],rep(0,8)))
write.csv(hr_ci,"beta_out.csv")

rm(beta_w.out.s)
rm(nlst, elig_nhis, design, design2, master, nhis)
rm(beta_w.out.t, beta_w.out.s1, beta_w.out.s2)

load("nhis_15_18_imputed_1.RData")
nhis_1518 = nhis; rm(nhis)
c("analypop", "age", "lung.cancer.before")%in%names(nhis_1518)
elig.indx = (nhis_1518$analypop==1 & nhis_1518$age>=40&nhis_1518$age<=84 & nhis_1518$lung.cancer.before==0)
elig_nhis_1518 = nhis_1518[elig.indx,]
elig_nhis_1518$combdty_n = rowSums(elig_nhis_1518[,c("emp", "hypertension", "chd", "angina", "heartattack", 
                                                     "heartdisease", "stroke","diab", "bron","kidney", "liver", 
                                                     "speceq", "prior.cancer")])
x0.1 = model.matrix(pid~age+female+race+edu6+emp+I(bmi <= 18.5)+I(cpd > 20)+
                      I(pkyr.cat>=30 & pkyr.cat<40)+I(pkyr.cat>=40 & pkyr.cat<50)+I(pkyr.cat>=50)+ 
                      I(log(bmi))+I(log(qtyears+1))+I(log(smkyears)), data = elig_nhis_1518)[,-1]
x0.2 = model.matrix(pid~age+female + race + edu6 + 
                      I(bmi <= 18.5)+I(bmi > 18.5 & bmi <= 20)+I(bmi > 25 & bmi <= 30)+I(bmi > 30 & bmi <= 35)+I(bmi > 35)+ 
                      emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                      I(log(qtyears + 1)) + I(log(cpd)) + I(sqrt(pkyr.cat)), data = elig_nhis_1518)[,-1]

# PURE risk
rr.1518 = exp(x0.1%*%modlcd$coefficients)
pureR1518=absR_w(beta_est=modlcd$coefficients, Lambda=Lambda0.s1_w$Lambda[19], x0=x0.1)$absR

# ABSOLUTE risk 
last_u=function(u,t,i){
  dif=u-t[i]
  dif[dif>0]=-10^6
  which.max(dif)
}
indx1=sapply(1:length(t_star.s), last_u, u=lambda.s1_w$u, t=t_star.s)
lambda0_1 = lambda.s1_w$lambda[indx1]
rr.1518_oth = exp(x0.2%*%modoth$coefficients)

absR1518 = colSums(outer(c(lambda0_1), c(rr.1518))*(exp(-outer(c(Lambda0.s1_w$Lambda), c(rr.1518))-outer(c(Lambda0.s2_w$Lambda), c(rr.1518_oth)))))
rm(lambda.s1_w, lambda.s2_w)


n_nhis = length(absR1518)
Delta_out=NULL
n_break = ifelse(j==28, 29, 30)
lw = 1200*(j-1)+ 0:(n_break-1)*40+1
up = 1:30*40+1200*(j-1)
if(j==28) up=c(1:(n_break-1)*40+1200*(j-1), n_nhis)
for(i in 1:n_break){
  x0.indx = order(absR1518)[lw[i]:up[i]]
  x0 = cbind(x0.1[x0.indx,],
             x0.2[x0.indx,!(names(betas.modoth)%in%names(betas.modlcd))])
  n_x0 = length(x0.indx)
  S_Delta_mu = Delta_mu_inf(Lambda0=cbind(Lambda0.s1_w$Lambda, Lambda0.s2_w$Lambda), 
                            t=t_star.s, t1=6, x = x0, 
                            betas = betas[,1:2],
                            Lambda0_w = array(c(Lambda0.s1_w$Lambda_wt, Lambda0.s2_w$Lambda_wt),
                                              c(n_s, length(t_star.s), 2)), 
                            beta_w=beta_w, deviate=T)
  
  Delta_mu.est = colSums(S_Delta_mu$Delta_mu)
  Delta_mu_w = sapply(1:n_x0, function(i)rowSums(S_Delta_mu$Delta_mu_w[,,i]))
  
  nhis_Delta_mu = cbind(nhis_ds, Delta_mu_w[1:n_s,])
  names(nhis_Delta_mu)[1:n_x0+3] = paste("Delta_mu_w.", c(1:n_x0), sep="")
  Delta.dsgn  <- svydesign(id=~psu, strata=~strata, weights=~wt_mort5, data=nhis_Delta_mu, nest=TRUE)
  Delta_mu.var = diag(vcov(svytotal(as.formula(paste0("~",paste(paste0("Delta_mu_w.", c(1:n_x0)), collapse = "+"))), 
                                    Delta.dsgn))+
                        var(Delta_mu_w[1:n_t+n_s,])*(n_t-1))
  S_Delta_mu.all = Delta_mu_3m(Lambda0=cbind(Lambda0.s1_w$Lambda, Lambda0.s2_w$Lambda)[1:6,], 
                           Lambda0.all = Lambda0.s_w$Lambda, t=t_star.s, t1=6, x=x0, betas=betas)
  Delta_mu.est1 = colSums(S_Delta_mu.all$Delta_mu)
  
  #Delta_mu.var=0
  Delta_out = rbind(Delta_out, cbind(Delta_mu.est, Delta_mu.var, 
                                     cbind(c(rr.1518), c(absR1518), c(pureR1518),
                                           elig_nhis_1518[,c("pid", "age", "wt", "combdty_n")])[x0.indx,],
                                     Delta_mu.est1)
  )
  print(i)
}
colnames(Delta_out)[c(3:5)]=c("rr", "absR", "pureR")


write.table(Delta_out, paste0("Delta1518_out_", j, ".txt"))

