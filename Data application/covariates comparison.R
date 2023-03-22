#rm(list=ls())
library(survey)
setwd("/Users/wangl29/Dropbox/Research/Codes/Data application")
nlst <- read.csv("prsn.csv")
nrow(nlst)
range(nlst$age)
nlst$age_c = cut(nlst$age, breaks = c(39, seq(54,69,5), 79, 85), include.lowest = T)
table(nlst$age_c)
nlst$race6=nlst$race
nlst$race=nlst$race6
nlst$race[nlst$ethnic==1]=3
nlst$race[nlst$race%in%c(4:7)]=4
nlst$race[nlst$race>7]=1
table(nlst$race)
nlst$female=nlst$gender-1
table(nlst$female)
nlst$edu6=nlst$educat
nlst$edu6[nlst$educat==1]=2
nlst$edu6[nlst$educat>7]=3
nlst$edu6=nlst$edu6-1
table(nlst$edu6)
nlst$weight.kg = nlst$weight/2.20462
nlst$weight.kg[is.na(nlst$weight.kg)]=mean(nlst$weight.kg[!is.na(nlst$weight.kg)])
nlst$height.m = nlst$height/39.3701
nlst$height.m[is.na(nlst$height.m)]=mean(nlst$height.m[!is.na(nlst$height.m)])
nlst$bmi = nlst$weight.kg/(nlst$height.m)^2
nlst$bmi_c=cut(nlst$bmi, c(0, 18.5, 24.9, 29.9, 200))
nlst$bmi_c = as.factor(nlst$bmi_c)
levels(nlst$bmi_c) = c("Under", "Normal", "Over", "Obese")
table(nlst$bmi_c)
nlst$pkyr_c = cut(nlst$pkyr, breaks = c(0,30, 50, 70, 600))
table(nlst$pkyr_c)
nlst$smkyr_c = cut(nlst$smokeyr, breaks=c(0,20,40,50,80))
table(nlst$smkyr_c)
#hist((nlst$age-nlst$age_quit)[!is.na(nlst$age_quit)], probs=seq(0,1,.1))
nlst$qtyr = nlst$age-nlst$age_quit
nlst$qtyr[nlst$qtyr<0]=0
nlst$qtyr_c = cut(nlst$qtyr, breaks = c(0,5,10,70), include.lowest = T)
table(nlst$qtyr_c)
table(nlst$smokeday, useNA = "ifany")
nlst$cpd_c = cut(nlst$smokeday, breaks =c(0,10,30,50,300))
# diabete, hypertension, heart disease/attach, bronch, stroke
table(nlst$diagdiab, useNA = "ifany")[2]/nrow(nlst)
table(nlst$diaghype, useNA = "ifany")[2]/nrow(nlst)
table(nlst$diaghear, useNA = "ifany")[2]/nrow(nlst)
table(nlst$diagbron, useNA = "ifany")[2]/nrow(nlst)
table(nlst$diagstro, useNA = "ifany")[2]/nrow(nlst)
diagcopd=nlst$diagcopd; diagcopd[is.na(nlst$diagcopd)]=0
diagchro=nlst$diagchro; diagchro[is.na(nlst$diagchro)]=0
nlst$bron = as.numeric(diagcopd+diagchro>0)
table(nlst$bron)
nlst_pcan = nlst[,which(startsWith(names(nlst), "canc"))[1:15]]
nlst_pcan[is.na(nlst_pcan)]=0
nlst$prior.cancer = as.numeric(rowSums(nlst_pcan)>0)
table(nlst$prior.cancer)
nlst$chd=NA; nlst$liver=NA;#nlst$canckidn
# revert na to 0
naT0 = function(vars){
  vars[is.na(vars)]=0
  vars
}
nlst$combdt8 = naT0(nlst$diagdiab)+naT0(nlst$diagemph)+naT0(nlst$diaghype)+naT0(nlst$diaghear)+naT0(nlst$bron)+
  naT0(nlst$diagstro)+naT0(nlst$prior.cancer)

################################## NHIS 1997-2001 ########################################################
load(file="data_example/nhis97_01.imputed.1.RData")
nhis_9701 = nhis; rm(nhis)
elig.indx = (nhis_9701$analypop==1 & nhis_9701$age>=40&nhis_9701$age<=84 & nhis_9701$lung.cancer.before==0)
elig_nhis_9701 = nhis_9701[elig.indx,]
elig_nhis_9701$age_c=cut(elig_nhis_9701$age, breaks = c(39, seq(54,69,5), 79, 85), include.lowest = T)
table(elig_nhis_9701$age_c)
table(elig_nhis_9701$race)
elig_nhis_9701$edu6=elig_nhis_9701$edu
elig_nhis_9701$edu6[elig_nhis_9701$edu6=="Missing"]="3"
elig_nhis_9701$edu6[elig_nhis_9701$edu6=="1"]="2"
elig_nhis_9701$edu6=as.numeric(elig_nhis_9701$edu6)-1
table(elig_nhis_9701$edu6)
elig_nhis_9701$bmi_c = as.factor(cut(elig_nhis_9701$bmi, c(0, 18.5, 24.9, 29.9, 200)))
levels(elig_nhis_9701$bmi_c) = c("Under", "Normal", "Over", "Obese")
table(elig_nhis_9701$bmi_c)
elig_nhis_9701$pkyr_c=cut(elig_nhis_9701$packyears,breaks = c(0,30, 50, 70, 600))
table(elig_nhis_9701$pkyr_c)
elig_nhis_9701$smkyr_c = cut(elig_nhis_9701$smkyears, breaks=c(0,20,40,50,80))
table(elig_nhis_9701$smkyr_c)
elig_nhis_9701$qtyr_c = cut(elig_nhis_9701$qtyears, breaks = c(0,5,10,70), include.lowest = T)
table(elig_nhis_9701$qtyr_c)
elig_nhis_9701$cpd_c = cut(elig_nhis_9701$cpd, breaks =c(0,10,30,50,300))
elig_nhis_9701$diaghear=as.numeric((elig_nhis_9701$heartattack + elig_nhis_9701$heartdisease)>0)
elig_nhis_9701$combdt8 = naT0(elig_nhis_9701$diab)+naT0(elig_nhis_9701$emp)+naT0(elig_nhis_9701$hypertension)+
  naT0(elig_nhis_9701$diaghear)+naT0(elig_nhis_9701$bron)+naT0(elig_nhis_9701$stroke)+naT0(elig_nhis_9701$prior.cancer)
elig_nhis_9701$combdt12 = elig_nhis_9701$combdt8+naT0(elig_nhis_9701$chd)+naT0(elig_nhis_9701$angina)+naT0(elig_nhis_9701$liver)+
  naT0(elig_nhis_9701$kidney)+  naT0(elig_nhis_9701$speceq)
ds.nhis9701  <- svydesign(id=~psu, strata=~strata, weights=~wt_mort5, data=elig_nhis_9701, nest=TRUE)
ds.nhis9701.fmr = subset(ds.nhis9701, former==1)
load("nhis_15_18_imputed_1.RData")
nhis_1518 = nhis; rm(nhis)
elig.indx = (nhis_1518$analypop==1 & nhis_1518$age>=40&nhis_1518$age<=84 & nhis_1518$lung.cancer.before==0)
elig_nhis_1518 = nhis_1518[elig.indx,]
elig_nhis_1518$age_c=cut(elig_nhis_1518$age, breaks = c(39, seq(54,69,5), 79, 85), include.lowest = T)
table(elig_nhis_1518$age_c)
table(elig_nhis_1518$race) 
table(elig_nhis_1518$edu6)
elig_nhis_1518$bmi_c = as.factor(cut(elig_nhis_1518$bmi, c(0, 18.5, 24.9, 29.9, 200)))
levels(elig_nhis_1518$bmi_c) = c("Under", "Normal", "Over", "Obese")
table(elig_nhis_1518$bmi_c)
elig_nhis_1518$pkyr_c=cut(elig_nhis_1518$packyears,breaks = c(0,30, 50, 70, 600))
table(elig_nhis_1518$pkyr_c)
elig_nhis_1518$smkyr_c = cut(elig_nhis_1518$smkyears, breaks=c(0,20,40,50,80))
table(elig_nhis_1518$smkyr_c)
elig_nhis_1518$qtyr_c = cut(elig_nhis_1518$qtyears, breaks = c(0,5,10,70), include.lowest = T)
table(elig_nhis_1518$qtyr_c)
elig_nhis_1518$cpd_c = cut(elig_nhis_1518$cpd, breaks =c(0,10,30,50,300))
elig_nhis_1518$diaghear=as.numeric((elig_nhis_1518$heartattack + elig_nhis_1518$heartdisease)>0)
elig_nhis_1518$combdt8 = naT0(elig_nhis_1518$diab)+naT0(elig_nhis_1518$emp)+naT0(elig_nhis_1518$hypertension)+
  naT0(elig_nhis_1518$diaghear)+naT0(elig_nhis_1518$bron)+naT0(elig_nhis_1518$stroke)+naT0(elig_nhis_1518$prior.cancer)
elig_nhis_1518$combdt12 = elig_nhis_1518$combdt8+naT0(elig_nhis_1518$chd)+naT0(elig_nhis_1518$angina)+naT0(elig_nhis_1518$liver)+naT0(elig_nhis_1518$kidney)+naT0(elig_nhis_1518$speceq)
ds.nhis1518  <- svydesign(id=~psu, strata=~strata, weights=~wt, data=elig_nhis_1518, nest=TRUE)
ds.nhis1518.fmr = subset(ds.nhis1518, former==1)

#######################################
dist_fun = function(var.name, ds.fmr=F){
  if(!ds.fmr){
    out=cbind(table(nlst[,var.name])/nrow(nlst),
              (svytable(as.formula(paste0("~", var.name)), ds.nhis9701)/sum(elig_nhis_9701$wt_mort5)),
              (svytable(as.formula(paste0("~", var.name)), ds.nhis1518)/sum(elig_nhis_1518$wt)))
  }else{
    out=cbind(table(nlst[,var.name])/sum(table(nlst[,var.name])),
              svytable(as.formula(paste0("~", var.name)), ds.nhis9701.fmr)/sum(elig_nhis_9701$wt_mort5[elig_nhis_9701$former==1]),
              svytable(as.formula(paste0("~", var.name)), ds.nhis1518.fmr)/sum(elig_nhis_1518$wt[elig_nhis_1518$former==1]))
  }
  rownames(out) = paste0(var.name, "_", rownames(out))
  out
}
prop_fun = function(var.name1, var.name2){
  out = matrix(c(table(nlst[,var.name1], useNA = "ifany")[2]/nrow(nlst),
                 svymean(as.formula(paste0("~",var.name2)), ds.nhis9701), 
                 svymean(as.formula(paste0("~",var.name2)), ds.nhis1518)),1,3)
  rownames(out) = var.name2
  out
}
distr3 = rbind(dist_fun("age_c"), 
               dist_fun("female"),
               dist_fun("race"),
               dist_fun("edu6"),
               dist_fun("bmi_c"),
               dist_fun("pkyr_c"),
               dist_fun("smkyr_c"),
               dist_fun("qtyr_c", ds.fmr=T),
               dist_fun("cpd_c"),
               prop_fun("diagdiab", "diab"),
               prop_fun("diagemph", "emp"),
               prop_fun("diaghype", "hypertension"),
               prop_fun("diaghear", "diaghear"),
               prop_fun("bron", "bron"),
               prop_fun("diagstro", "stroke"),
               prop_fun("prior.cancer", "prior.cancer"),
               c(mean(nlst$combdt8),  svymean(~combdt8,    ds.nhis9701), svymean(~combdt8,    ds.nhis1518)),
               c(NA, svymean(~chd,      ds.nhis9701), svymean(~chd,      ds.nhis1518)),
               c(NA, svymean(~angina,   ds.nhis9701), svymean(~angina,   ds.nhis1518)),
               c(NA, svymean(~liver,    ds.nhis9701), svymean(~liver,    ds.nhis1518)),
               c(NA, svymean(~kidney,   ds.nhis9701), svymean(~kidney,   ds.nhis1518)),
               c(NA, svymean(~speceq,   ds.nhis9701), svymean(~speceq,   ds.nhis1518)),
               c(NA, svymean(~combdt12, ds.nhis9701), svymean(~combdt12, ds.nhis1518))
)

rownames(distr3)[45:51]=c("combdt8", "chd", "angina", "kidney", "liver", "special_eqp", "combdt12")
round(distr3*100, 1)

