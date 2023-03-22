rm(list=ls())
library(ggplot2)
library(survey)
library(RColorBrewer)
########################################################################################################
path = paste0("/Users/wangl29/Dropbox/Research/Codes/Data application/")
Delta_out = NULL
for (j in 1:28){
  Delta_tmp  = read.table(paste0(path, "Delta1518_out_" , j, ".txt"), row.names = NULL)[,-1]
  Delta_out  = rbind(Delta_out, Delta_tmp)
  #k=k+1
  print(j)
}
names(Delta_out)[c(1,10)]=c("Delta_mu.est0", "Delta_mu.est")
Delta_out$age_c = cut(Delta_out$age, breaks = ,c(40, seq(44, 84, 5)), include.lowest = T)

load("nhis_15_18_imputed_1.RData")
##load("nhis_15_18_imputed_2.RData")
nhis_1518 = nhis; rm(nhis)

nhis_1518$uspstf50_80 = 0
nhis_1518$uspstf50_80[(nhis_1518$age>=50)&(nhis_1518$age<=80)&(nhis_1518$packyears>=20)&(nhis_1518$qtyears<=15)]=1
table(nhis_1518$uspstf50_80, nhis_1518$uspstf.eligible)
nhis_1518$uspstf.eligible[is.na(nhis_1518$uspstf.eligible)]=FALSE
ds.nhis = svydesign(id=~psu, strata = ~strata, weights=~wt, nest = T, data = nhis_1518)
svymean(~uspstf.eligible, ds.nhis)
svymean(~uspstf50_80, ds.nhis)
mean(nhis_1518$uspstf50_80)
mean(nhis_1518$uspstf.eligible)

elig.indx = (nhis_1518$analypop==1 & nhis_1518$age>=40&nhis_1518$age<=84 & nhis_1518$lung.cancer.before==0)
elig_nhis_1518 = nhis_1518[elig.indx,]
c("analypop", "age", "lung.cancer.before")%in%names(nhis_1518)
elig_nhis_1518$combdty_n = rowSums(elig_nhis_1518[,c("emp", "hypertension", "chd", "angina", "heartattack", 
                                                     "heartdisease", "stroke","diab", "bron","kidney", "liver", 
                                                     "speceq", "prior.cancer")])
table(as.numeric(elig_nhis_1518$uspstf.eligible), elig_nhis_1518$age)
table(as.numeric(elig_nhis_1518$uspstf.eligible), elig_nhis_1518$packyears>=20)
table(as.numeric(elig_nhis_1518$uspstf.eligible), elig_nhis_1518$qtyears<=15)
table(elig_nhis_1518$uspstf50_80, elig_nhis_1518$age)
ds.nhis = svydesign(id=~psu, strata = ~strata, weights=~wt, nest = T, data = elig_nhis_1518)
# selection threshold
p.sel = svymean(~uspstf50_80, ds.nhis)
mean(elig_nhis_1518$uspstf50_80)
svymean(~uspstf.eligible, ds.nhis)
mean(elig_nhis_1518$uspstf.eligible)
table(elig_nhis_1518$uspstf.eligible, elig_nhis_1518$uspstf50_80, elig_nhis_1518$age<55)
absR_out = cbind(elig_nhis_1518[,c("pid", "age","wt", "combdty_n")], Delta_out$rr, Delta_out$absR)
#write.table(Delta_out, paste0(path, "Delta1518_out_all.txt"), row.names = F)
ds = svydesign(id=~1,weights=~wt, data=Delta_out)
q.indx=100-p.sel[[1]]*100
thred_Delta = as.numeric(svyquantile(~Delta_mu.est, ds, quantile=q.indx/100))
thred_absR  = as.numeric(svyquantile(~absR        , ds, quantile=q.indx/100))
benf1 =as.numeric(Delta_out$Delta_mu.est>=thred_Delta)
absr1 =as.numeric(Delta_out$absR        >=thred_absR)
slct = rep("rsk", nrow(Delta_out))
slct[(benf1*absr1)==1]="both"
slct[(benf1==1)&(absr1==0)]="bnf"
slct[(benf1==0)&(absr1==0)]="none"
Delta_out$slct = slct
freq.tab = svytable(~slct, design=svydesign(id=~1,weights=Delta_out$wt))/sum(Delta_out$wt)
print(freq.tab)
thred_Delta*365.25
range_Delta_mu = cbind(range(Delta_out$Delta_mu.est), sapply(1:4,  function(i) range(Delta_out$Delta_mu.est[slct==unique(slct)[i]])))
range_absR     = cbind(range(Delta_out$absR),         sapply(1:4,  function(i) range(Delta_out$absR[slct==unique(slct)[i]])))


mean(Delta_out$absR<range_absR[2,5])
mean(Delta_out$absR<range_absR[1,5])

##################################Figure 3##########################################
Delta_plot = data.frame(Delta_mu = log(Delta_out$Delta_mu.est*365.25),
                        absR = log(Delta_out$absR),
                        slct = slct,
                        wt   = Delta_out$wt)
set.seed(1)
plot.indx = rep(1:nrow(Delta_plot), round(Delta_plot$wt))[sample(1:sum(Delta_plot$wt), 70000)]
Delta_plot1 = Delta_plot[plot.indx,]
## Use densCols() output to get density at each point
x <- densCols(Delta_plot1$absR,Delta_plot1$Delta_mu, 
              colramp=colorRampPalette(c("black", "white")))
Delta_plot1$dens <- col2rgb(x)[1,] + 1L

## Map densities to colors

cols <-  colorRampPalette(c(brewer.pal(8,"Reds")[3:6],
                            rev(hsv(1, 1, seq(0,1,length.out = 12)))))(256)
Delta_plot1$col <- cols[Delta_plot1$dens]
## Plot it, reordering rows so that densest points are plotted on top
par(mfrow=c(1,1), mar=c(3, 3, 0.6, 0.3))
plot(Delta_mu~absR, data=Delta_plot1[order(Delta_plot1$dens),], pch=20, col=col, cex=0.3,
     xlab="", ylab="", cex.axis=0.8, xaxt="n", yaxt="n", ylim=c(-4,6))
axis(1, at = seq(-10, 0, 2), label = rep("", 6), tck = -0.01)
axis(1, seq(-10, 0, 2), line = -0.9, lwd = 0, cex.axis = 0.8)
axis(2, at = seq(-4,6,2), label = rep("", 6), tck = -0.01)
axis(2, seq(-4,6,2), line = -0.9, lwd = 0, cex.axis = 0.8)
mtext(side=1, text="Absolute risk of lung cancer death \n in 10 years (Log scale)", 
      line=1.8, cex=0.8)
mtext(side=2, text="Life-gained from screening \n in log(days)", 
      line=1, cex=0.8)
abline(v=log(thred_absR), col="blue")
abline(h=log(thred_Delta*365.25), col="blue")
text(-8.5, -3.5, paste0("None\n (",         round(freq.tab[3]*100, 1), "%)"), lwd=1.5, cex=0.7)
text(-2.2, -3.5, paste0("Risk Only\n (",    round(freq.tab[4]*100, 1), "%)"), lwd=1.5, cex=0.7)
text(-8, 5.5,    paste0("Benefit Only\n (", round(freq.tab[1]*100, 1), "%)"), lwd=1.5, cex=0.7)
text(-2, 5.5,    paste0("Both\n (",         round(freq.tab[2]*100, 1), "%)"), lwd=1.5, cex=0.7)
text(log(thred_absR), -2.2, paste0("Risk threshold\n (r=",round(thred_absR, 3),")"), 
                                 lwd=1.5, cex=0.7, col="blue", srt=-90, font=2)
text(-7.5, log(thred_Delta*365.25), paste0("Benefit threshold\n (LYG=",round(thred_Delta*365.25,1)," days)"), 
     lwd=1.5, cex=0.7, col="blue", font=2)

Delta_slct = Delta_out[slct%in%c("bnf", "both"),]
absR_slct  = Delta_out[slct%in%c("rsk", "both"),]

########## Comparing characteristics of people selected by Benifit vs absR ###########
Delta_only = Delta_out[slct%in%c("bnf"),]
absR_only  = Delta_out[slct%in%c("rsk"),]
B_absR = rbind(c(svymean(~absR, svydesign(id=~1,weights=~wt, data=Delta_only)),
                 svymean(~absR, svydesign(id=~1,weights=~wt, data=absR_only ))),
               c(svymean(~Delta_mu.est, svydesign(id=~1,weights=~wt, data=Delta_only))*365.25,
                 svymean(~Delta_mu.est, svydesign(id=~1,weights=~wt, data=absR_only ))*365.25),
               c(svymean(~age, svydesign(id=~1,weights=~wt, data=Delta_only)),
                 svymean(~age, svydesign(id=~1,weights=~wt, data=absR_only ))),
               c(svymean(~combdty_n, svydesign(id=~1,weights=~wt, data=Delta_only)),
                 svymean(~combdty_n, svydesign(id=~1,weights=~wt, data=absR_only ))),
               c(svyquantile(~absR, svydesign(id=~1,weights=~wt, data=Delta_only), quantiles=0.5),
                 svyquantile(~absR, svydesign(id=~1,weights=~wt, data=absR_only ), quantiles=0.5)),
               c(svyquantile(~Delta_mu.est, svydesign(id=~1,weights=~wt, data=Delta_only), quantiles=0.5)*365.25,
                 svyquantile(~Delta_mu.est, svydesign(id=~1,weights=~wt, data=absR_only ), quantiles=0.5)*365.25),
               c(svyquantile(~age, svydesign(id=~1,weights=~wt, data=Delta_only), quantiles=0.5),
                 svyquantile(~age, svydesign(id=~1,weights=~wt, data=absR_only ), quantiles=0.5)),
               c(svyquantile(~combdty_n, svydesign(id=~1,weights=~wt, data=Delta_only), quantiles=0.5),
                 svyquantile(~combdty_n, svydesign(id=~1,weights=~wt, data=absR_only ), quantiles=0.5))
)
row.names(B_absR) = paste0(rep(c("Mean", "Median"), each=4),"_", c("absR", "LYG", "Age", "combdty"))
colnames(B_absR) = c("Benefit", "AbsR")
B_absR

Delta_age = svytable(~age, svydesign(id=~1,weights=~wt, data=Delta_slct))
absR_age  = svytable(~age, svydesign(id=~1,weights=~wt, data=absR_slct ))
Delta_age.c = svytable(~age_c, svydesign(id=~1,weights=~wt, data=Delta_slct))
absR_age.c  = svytable(~age_c, svydesign(id=~1,weights=~wt, data=absR_slct ))

agedist.c = data.frame(age_c=names(absR_age.c), Delta = as.numeric(Delta_age.c), absR = as.numeric(absR_age.c))
sel.prob = data.frame(age_c=names(absR_age.c), 
                      Delta = as.numeric(agedist.c[,2]/svytable(~age_c, svydesign(id=~1,weights=~wt, data=Delta_out))),
                      absR  = as.numeric(agedist.c[,3]/svytable(~age_c, svydesign(id=~1,weights=~wt, data=Delta_out)))
)
sel.prob$age_c = factor(as.factor(sel.prob$age_c), levels=sel.prob$age_c)
#Delta_agedist= as.data.frame(table(Deata_15$age)/nrow(Deata_15)); names(Delta_agedist)[1]="age"
#absR_agedist = as.data.frame(table(absR_15$age)/nrow(absR_15));   names(absR_agedist)[1]="age"
agedist = merge(merge(data.frame(age=40:84), Delta_age, by="age", all=T), absR_age, by="age", all=T)
names(agedist)=c("age", "Delta_mu", "absR")
agedist[is.na(agedist)]=0

agedist[,2]=agedist[,2]/svytable(~age, svydesign(id=~1,weights=~wt, data=Delta_out))
agedist[,3]=agedist[,3]/svytable(~age, svydesign(id=~1,weights=~wt, data=Delta_out))

agedist$age.c=cut(agedist$age,c(40, seq(44, 84, 5)), include.lowest = T)


# Table S5
range_Delta_mu*365.25
range_absR*365.25
Delta_out$id = 1:nrow(Delta_out)
Delta_out = Delta_out[order(Delta_out$Delta_mu.est),]
Delta_out$q.Delta = cumsum(Delta_out$wt)/sum(Delta_out$wt)
Delta_out = Delta_out[order(Delta_out$absR),]
Delta_out$q.absR = cumsum(Delta_out$wt)/sum(Delta_out$wt)
Delta_out = Delta_out[order(Delta_out$id),]
range(Delta_out$q.Delta[Delta_out$slct=="bnf"])
range(Delta_out$q.Delta[Delta_out$slct=="rsk"])
q0.indx=which.min(Delta_out$Delta_mu.est[Delta_out$slct=="rsk"])
q0 = (Delta_out$q.Delta[Delta_out$slct=="rsk"])[q0.indx]; q0*100
q0.absR = (Delta_out$q.absR[Delta_out$slct=="rsk"])[q0.indx]
q1.indx=which.min(abs(Delta_out$q.Delta[Delta_out$slct=="none"]-.25))
q1 = (Delta_out$q.Delta[Delta_out$slct=="none"])[q1.indx]; q1*100
q1.absR = (Delta_out$q.absR[Delta_out$slct=="none"])[q1.indx]
q2.indx=which.min(abs(Delta_out$q.Delta[Delta_out$slct=="rsk"]-.5))
q2 = (Delta_out$q.Delta[Delta_out$slct=="rsk"])[q2.indx]; q2*100
q2.absR = (Delta_out$q.absR[Delta_out$slct=="rsk"])[q2.indx]
q3.indx=which.min(abs(Delta_out$q.Delta[Delta_out$slct=="bnf"]-.77))
q3 = (Delta_out$q.Delta[Delta_out$slct=="bnf"])[q3.indx]; q3*100
q3.absR = (Delta_out$q.absR[Delta_out$slct=="bnf"])[q3.indx]
q4.indx = which.max((Delta_out$q.Delta-Delta_out$q.absR)[Delta_out$slct=="bnf"])
q4 = (Delta_out$q.Delta[Delta_out$slct=="bnf"])[q4.indx]; q4*100
q4.absR = (Delta_out$q.absR[Delta_out$slct=="bnf"])[q4.indx]

q = c(q0, q1, q2, q3, q4)
q.wt = sum(Delta_out$wt)*q
q.wt.indx = apply(abs(outer(cumsum(Delta_out$wt[order(Delta_out$Delta_mu.est)]), 
                            q.wt,FUN="-")), 2, which.min)
q.wt.indx = c(which.min(Delta_out$Delta_mu.est), 
              c(1:nrow(Delta_out))[order(Delta_out$Delta_mu.est)][q.wt.indx],
              which.max(Delta_out$Delta_mu.est))
Delta_out$slct[q.wt.indx]
Delta_out$Delta_mu.est[q.wt.indx]*365.25
tableS5 = data.frame(Q=round(c(0, q, 1)*100,1))
tableS5$Delta_mu.est= round(Delta_out$Delta_mu.est[q.wt.indx]*365.25, 2)
tableS5$Delta_mu.sd = round(sqrt(Delta_out$Delta_mu.var[q.wt.indx])*365.25, 2)
tableS5$absR        = round(Delta_out$absR[q.wt.indx]*100,2)
tableS5$q.absR = round(c(Delta_out$q.absR[which.min(Delta_out$Delta_mu.est)],
                  q0.absR, q1.absR, q2.absR, q3.absR, q4.absR,
                  Delta_out$q.absR[which.max(Delta_out$Delta_mu.est)])*100,1)
svyquantile(~Delta_mu.est, svydesign(id=~1,weights=~wt, data=Delta_out),
            quantile=q)*365.25
tableS5$slct = c("None", "Risk", "None", "Risk", "Benefit", "Benefit", "Both")
x0.1 = model.matrix(pid~age+female+race+edu6+bmi+ qtyears + smkyears + cpd + packyears, data = elig_nhis_1518)[,-1]
x0.1 = as.data.frame(x0.1[match(Delta_out$pid, elig_nhis_1518$pid),])
# check the order
sum(x0.1$age == Delta_out$age)
x0.2 = model.matrix(pid~emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + 
                      bron + kidney + liver + speceq + prior.cancer, data = elig_nhis_1518)[,-1]
x0.2 = x0.2[match(Delta_out$pid, elig_nhis_1518$pid),]

x0.1$race="NH-White"
x0.1$race[x0.1$race1==1]="NH-Black"
x0.1$race[x0.1$race2==1]="Hispanic"
x0.1$race[x0.1$race3==1]="NH-Other"
x0.1$edu6= as.factor(x0.1$edu6)
levels(x0.1$edu6) =  c("<= 11th grade", "High school", "Post-high", 
                       "Assoc. degree/some college", "Bachelor", "Post-graduate")
table(x0.1$edu6)
tableS5.cov = x0.1[q.wt.indx,c(1,2,12,6:11)]
tableS5 = cbind(tableS5, tableS5.cov)
tableS5$n.combdt = rowSums(x0.2[q.wt.indx,c("emp", "hypertension", "chd", "angina", "heartattack", "heartdisease", "stroke",
                                           "diab", "bron", "kidney", "liver", "speceq", "prior.cancer")])
tableS5
write_xlsx(tableS5, "/Users/wangl29/Dropbox/Research/Work with Li/Life gained from intervention/data example/tableS5.xlsx")
for(i in 1:5) print(sum(((Delta_out$absR*100)<tableS5$absR[i])*Delta_out$wt)/sum(Delta_out$wt))

# Figure 4
agedist = merge(merge(data.frame(age=40:84), Delta_age, by="age", all=T), absR_age, by="age", all=T)
names(agedist)=c("age", "Delta_mu", "absR")
agedist[is.na(agedist)]=0

agedist[,2]=agedist[,2]/svytable(~age, svydesign(id=~1,weights=~wt, data=Delta_out))
agedist[,3]=agedist[,3]/svytable(~age, svydesign(id=~1,weights=~wt, data=Delta_out))

agedist$age.c=cut(agedist$age,breaks=c(40, seq(44, 84, 5)), include.lowest = T)

agedist.plot = data.frame(age=rep(agedist$age,2), freq=c(agedist$Delta_mu, agedist$absR),
                          Method=rep(c("Delta_mu", "absR"), each=45))
ggplot( aes(x=age, y=freq, group=Method, color=Method), data=agedist.plot) +
  geom_hline(yintercept = seq(0, .6, .1), color = "white") +
  geom_vline(xintercept = seq(40, 80, 10), color = "white") +
  geom_line(aes(linetype=Method)) +
  scale_y_continuous(name="% Ever-Smokers Selected for Screening", limits=c(0, .6),
                     breaks = seq(0, .6, .1), 
                     labels = seq(0, 60, 10))+
  scale_x_continuous(name="Age (in years)", limits=c(40, 85),
                     breaks = c(seq(40, 80, 5), 84),
                     labels = c(seq(40, 80, 5), 84))+
  scale_colour_manual(labels = c("Risk-based", "Benefit-based"), values = c("red", "darkblue"))+
  scale_linetype_manual(labels = c("Risk-based", "Benefit-based"), values=c("longdash", "solid"))+
  annotate("text", x=71, y=.55, label="Risk-based", col="red")+
  annotate("text", x=55, y=.43, label="Benefit-based", col="darkblue")+
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
  )

age_c.dist.plot = data.frame(age_c=rep(1:9, 2), freq=c(sel.prob$Delta, sel.prob$absR),
                             Method=rep(c("Delta_mu", "absR"), each=9))
#Table S6
TabS6 = round(rbind(age_c.dist.plot$freq[1:9],
                    age_c.dist.plot$freq[10:18])*100,
              1)
