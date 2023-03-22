#################################################################################
#-FUNCTION sam.pps is to select a sample with pps sampling                      #
# INPUT:  popul - the population including response and covariates              #
#         MSize - the size of measure to compute the selection probabilities(>0)#
#         n -     the sample size                                               #
# OUTPUT: sample selected with pps sampling of size n,                          #
#         including response, covariates, and sample weights                    #
#################################################################################
sam.pps<-function(popul,Msize, n){
  N=nrow(popul)
  pps.samID=sample(N,n,replace=F,prob=Msize)   #F but make sure N is large!
  if (dim(popul)[2] == 1){
    sam.pps.data=as.data.frame(popul[pps.samID,])
    names(sam.pps.data) = names(popul)
  }else{sam.pps.data=popul[pps.samID,]}
  sam.pps.data$wt=sum(Msize)/n/Msize[pps.samID]               
  return(sam.pps = sam.pps.data)
}

###############################################################################################
# cum_bsln_hzd is the FUNCTION to estimate cumulative baseline hazard using Breslow's method  #
# INPUT                                                                                       #
#  surv.fit - an object of fitted Cox regression model (output of coxph or svycoxph)          #
#  t_star   - a vector of desired time at which the cumulative baseline hazards are estimated #
# OUPUT                                                                                       #
#  cum_hzd  - a vector of estimates of cumulative baseline hazard at given time points t_star #
#  lambda.u - a vector of baseline hazard at sorted unique event times                        #
#  u        - sorted unique event times                                                       #   
###############################################################################################
cum_bsln_hzd = function(surv.fit, t_star, lambda=F){
  var.names = all.vars(surv.fit$formula)
  t = var.names[1]
  d = var.names[2]
  dat = eval(surv.fit$call$data)
  if(is.null(dat)) dat = eval(surv.fit$survey.design$call$data)
  rel_hzd = exp(model.matrix(surv.fit)%*%c(surv.fit$coefficients))
  dat$rel_hzd = rel_hzd
  dat$wt=1;
  if(!is.null(surv.fit$weights)) dat$wt = surv.fit$weights*surv.fit$n
  dat$w_e = dat$wt*dat$rel_hzd

  dat = dat[order(dat[,t], -dat[,d]),]
  dat$dnom = rev(cumsum(rev(dat$w_e)))
  dat1 = dat[dat[,d]==1,c(t, d, "wt", "dnom")]
  # Find ties in the dataset
  dup = duplicated(dat1[,t])
  last_l0=function(dat,i){
    dif=dat1[,t]-t_star[i]
    dif[dif>0]=-10^6
    which.max(dif)
  }
  if(sum(dup)>0){
    dnom_dat = dat1[!dup, c(t, "dnom")]
    dup_t = unique(dat1[,t][dup])
    dup_dat1 = dat1[dat1[,t]%in%dup_t, c(t, "wt")]
    dup_de = aggregate(dup_dat1$wt, by=list(dup_dat1[,t]), FUN=sum)
    names(dup_de) = c(t, "wt")
    unq_dat1 = dat1[!dat1[,t]%in%dup_t, c(t, "wt")]
    num_dat = rbind(unq_dat1, dup_de)
    num_dat = num_dat[order(num_dat[,t]),]
    #sum(dnom[,t]!=num_dat[,t])
    indx=sapply(1:length(t_star), last_l0, dat=num_dat)
    lambda.u = num_dat$wt/dnom_dat$dnom
  }else{
    indx=sapply(1:length(t_star), last_l0, dat=dat1)
    lambda.u = dat1$wt/dat1$dnom
  }
  cum_hzd = cumsum(lambda.u)[indx]
  if(!lambda){
    return(cum_hzd  = cum_hzd)
  }else{
    return(list(cum_hzd  = cum_hzd,
                lambda.u = lambda.u, 
                u = unique(dat1[,t])))
  }
  
}

#########################################################################################
#Derivative of U and beta (cox regression) w.r.t. pseudo weight
beta_pw.cox = function(x.mtrx, rel_hzd, dat, pw, t="t", d="d", pi.c_est=NULL, post=NULL){
  if(!is.null(post)){
    dat$f = post$f
    f_w = post$f_w
    dat[,pw] = dat[,pw]/dat$f
  }else{dat$f=1}
  n = nrow(dat)
  p = ncol(x.mtrx)
  if(is.null(pi.c_est)) pi.c_est = 1/dat[,pw]
  dat$id = 1:n
  dat$rel_hzd = rel_hzd
  dat$pw_e = dat$f*dat[,pw]*dat$rel_hzd
  
  x.mtrx = as.matrix(x.mtrx[order(dat[,t]),])
  dat = dat[order(dat[,t]),]
  dat$H_dnom = rev(cumsum(rev(dat$pw_e)))# H2
  dat$H_num = as.matrix(cumsum(as.data.frame(c(dat$pw_e)*x.mtrx)[n:1,]))[n:1,] #H1
  
  ties = duplicated(dat[,t])
  if(sum(ties)>0){
    H_uniq = dat[!ties,c(t, "H_dnom", "H_num")]
    H_dat = H_uniq[rep(1:nrow(H_uniq), table(dat[,t])),]
    #h_dat[,t]==dat[,t]
    dat[, c("H_dnom", "H_num")] = H_dat[, c("H_dnom", "H_num")]
    
  }
  
  x.mtrx = as.matrix(x.mtrx[order(dat$id), ])
  dat = dat[order(dat$id), ]
  
  H = as.matrix(dat$H_num/dat$H_dnom)
  
  d_indx = which(dat[,d]==1)
  dat1 = dat[d_indx,]
  # create a two-column matrix sum_wt.t, with the first column being the unique event time, 
  # the second column being the sum of weights for individuals having the event time t
  if(sum(duplicated(dat1[, t]))>0){
    dat1$weight = dat1[,pw]*dat1$f
    
    dat1_tb <- data.table(dat1, key=t)
    names(dat1_tb)[ncol(dat1_tb)]="weight"
    names(dat1_tb)[names(dat1_tb)==t]="t_tb"
    sum_wt.t <- as.data.frame(dat1_tb[, sum(weight), by=list(t_tb)])
    #sum_wt.t = aggregate(dat1[,pw]*dat1$f, by=list(rep(1, nrow(dat1)),dat1[,t]), sum)
    sum_wt.t = sum_wt.t[match(unique(dat1[, t]), sum_wt.t[,1]),]
    unq_t.d_indx = d_indx[!duplicated(dat1[, t])]
  }else{
    sum_wt.t=cbind(dat1[, t], dat1$f*dat1[,pw])
    unq_t.d_indx = d_indx
  }
  
  #U_w_2.k.out =NULL
  if(is.null(post)){
    U_w_2 = 0 # Ai 
    for(j in 1:nrow(sum_wt.t)){
      k=unq_t.d_indx[j]
      U_w_2.k = sum_wt.t[j, 2]*(c((dat[,t]>=dat[,t][k])*dat$rel_hzd)*x.mtrx/dat$H_dnom[k]-
                                  outer(c((dat[,t]>=dat[,t][k])*dat$rel_hzd),as.numeric(as.matrix(dat$H_num)[k,]))/(dat$H_dnom[k])^2)
      U_w_2 = U_w_2.k+U_w_2
    }
    
    #for(k in d_indx){
    #  U_w_2.k = dat$f[k]*dat[k,pw]*(c((dat[,t]>=dat[,t][k])*dat$rel_hzd)*x.mtrx/dat$H_dnom[k]-
    #                         outer(c((dat[,t]>=dat[,t][k])*dat$rel_hzd),as.numeric(as.matrix(dat$H_num)[k,]))/(dat$H_dnom[k])^2)
    #  U_w_2 = U_w_2.k+U_w_2
    #  U_w_2.k.out = cbind(U_w_2.k.out, U_w_2.k[,1])
    #  #print(k)
    #}
    Ui_pw= as.matrix(dat[,d]*(x.mtrx- H)-U_w_2)
  }else{
    if(nrow(f_w)==length(d_indx)){
      fw_indx=d_indx
    }else if(nrow(f_w)==nrow(dat)){
      fw_indx=1:nrow(dat)
    }else{
      stop("The partial derivative matrix of poststratification factor should be given for event cases or the whole sample!")
    }
    U_w_2 = 0 # Ai 
    U_w_4 = 0 # Bi
    #U_w_4.k.out =NULL
    for(j in 1:nrow(sum_wt.t)){
      k=unq_t.d_indx[j]
      U_w_2.k = sum_wt.t[j, 2]*(c(dat$f*(dat[,t]>=dat[,t][k])*dat$rel_hzd)*x.mtrx/dat$H_dnom[k]-
                                  outer(c(dat$f*(dat[,t]>=dat[,t][k])*dat$rel_hzd),as.numeric(as.matrix(dat$H_num)[k,]))/(dat$H_dnom[k])^2)
      U_w_2 = U_w_2.k+U_w_2
      #U_w_2.k.out = cbind(U_w_2.k.out, U_w_2.k[,1])
      U_w_4.k = sum_wt.t[j, 2]*(f_w_mtrx(f_w, (c(dat[,pw]*(dat[,t]>=dat[,t][k])*dat$rel_hzd)*x.mtrx)[fw_indx,])/dat$H_dnom[k]-
                                  outer(c(f_w_mtrx(f_w, c((dat[,pw]*(dat[,t]>=dat[,t][k])*dat$rel_hzd)[fw_indx]))),
                                        as.numeric(as.matrix(dat$H_num)[k,]))/dat$H_dnom[k]/dat$H_dnom[k])
      U_w_4 = U_w_4.k+U_w_4
      #U_w_4.k.out = cbind(U_w_4.k.out,U_w_4.k[,1])
      #print(j)
    }
    #for(k in d_indx){
    #  U_w_2.k = dat$f[k]*dat[k,pw]*(c(dat$f*(dat[,t]>=dat[,t][k])*dat$rel_hzd)*x.mtrx/dat$H_dnom[k]-
    #                                  outer(c(dat$f*(dat[,t]>=dat[,t][k])*dat$rel_hzd),as.numeric(as.matrix(dat$H_num)[k,]))/(dat$H_dnom[k])^2)
    #  U_w_2 = U_w_2.k+U_w_2
    #  U_w_2.k.out=cbind(U_w_2.k.out, U_w_2.k[,1])
    #  U_w_4.k = dat$f[k]*dat[k,pw]*((f_w%*%as.matrix((c(dat[,pw]*(dat[,t]>=dat[,t][k])*dat$rel_hzd)*x.mtrx)[d_indx,]))/dat$H_dnom[k]-
    #                                  outer(c(f_w%*%c((dat[,pw]*(dat[,t]>=dat[,t][k])*dat$rel_hzd)[d_indx])),as.numeric(as.matrix(dat$H_num)[k,]))/(dat$H_dnom[k])^2)
    #  U_w_4 = U_w_4.k+U_w_4
    #  U_w_4.k.out = cbind(U_w_4.k.out,U_w_4.k[,1])
    #}
    if(length(fw_indx)== length(d_indx)){
      id_d = c(dat$id[d_indx], dat$id[dat[,d]==0])
      U_w_3 = rbind(f_w_mtrx(f_w, (dat[,pw]*(x.mtrx- H))[d_indx,]), matrix(0, n-sum(dat[,d]), p))[order(id_d),]
      U_w_4 = rbind(U_w_4, matrix(0, n-sum(dat[,d]), p))[order(id_d),]
    }else{
      U_w_3 = f_w_mtrx(f_w, (c(dat[,pw]*dat[,d])*(x.mtrx- H)))
    }
    Ui_pw= as.matrix(c(dat$f*dat[,d])*(x.mtrx- H)-U_w_2+U_w_3-U_w_4)
    
  }
  U_beta_1 = 0
  for(j in 1:nrow(sum_wt.t)){
    k=unq_t.d_indx[j]
    sum_j = t(as.matrix(x.mtrx*c(dat$pw_e*as.numeric(dat[,t]>=dat[,t][k]))))%*%as.matrix(x.mtrx)
    U_beta_1 =-sum_wt.t[j, 2]*sum_j/dat$H_dnom[k]+U_beta_1
  }
  #for(i in d_indx){
  #  sum_j = t(as.matrix(x.mtrx*c(dat$pw_e*as.numeric(dat[,t]>=dat[,t][i]))))%*%as.matrix(x.mtrx)
  #  U_beta_1 =-dat$f[i]*dat[i,pw]*sum_j/dat$H_dnom[i]+U_beta_1
  #}
  
  U_beta = U_beta_1 +t(c(dat[d_indx,pw])*H[d_indx,])%*%H[d_indx,]
  beta_pw = -Ui_pw%*%solve(U_beta)
  Delta_beta.pw = c(dat[,pw])*beta_pw
  
  #print("beta_pw.cox OK")
  return(list(beta_pw = beta_pw,
              Delta_beta.pw = Delta_beta.pw)
         )
}
###################################################################################################
# lambda_w is the FUNCTION to estimate baseline hazards and their derivatives w.r.t. sample weight#
# INPUT                                                                                           #
#  surv.fit - an object of fitted Cox regression model (output of coxph or svycoxph)              #
#  beta_wt  - derivative of Cox regression coefficients w.r.t. sample weight if different         #
#             from the output of coxph or svycoxph. The default is NULL                           #
# OUPUT                                                                                           #
#  lambda    - a vector of estimates of  baseline hazards at sorted unique event times            #
#  lambda_wt - a matrix of derivative of baseline hazards w.r.t. sample weight                    #
#              at sorted unique event times                                                       #
#  u         - sorted unique event times                                                          #   
###################################################################################################
lambda_w = function(surv.fit, beta_wt=NULL){
  var.names = all.vars(surv.fit$formula)
  t = var.names[1]
  d = var.names[2]
  dat = eval(surv.fit$call$data)
  if(is.null(dat)){
    dat = eval(surv.fit$survey.design$variables)
    dat$wt = 1/surv.fit$survey.design$prob
  } else(dat$wt=1)
  x.mtrx = model.matrix(surv.fit)
  rel_hzd = exp(x.mtrx%*%c(surv.fit$coefficients))
  dat$rel_hzd = rel_hzd
  n = nrow(dat)
  dat$id = 1:n

  Nt=NULL; Zt=NULL; Yt=NULL; t_unq = NULL
  #dat = dat[order(dat$t),]
  dat1 = dat[dat[,d]==1,]
  lambda_dat = data.frame(t = unique(dat1[,t]))
  names(lambda_dat)=t
  dat$Yi_t = outer(dat[,t], lambda_dat[,], FUN=">=")
  colnames(dat$Yi_t)=paste("t=", lambda_dat[,t], sep="")
  row.names(dat$Yi_t) = row.names(dat)
  dat$Ii_t = outer(dat[,t], lambda_dat[,t], FUN="==")
  colnames(dat$Ii_t)=paste("t=", lambda_dat[,t], sep="")
  row.names(dat$Ii_t) = row.names(dat)
  
  for (i in 1:length(lambda_dat[,t])){
    Nt_i = sum((dat$wt*dat[,d])*dat$Ii_t[,i])
    Zt_i = sum((dat$wt*dat$rel_hzd)*dat$Yi_t[,i])
    t_unq = c(t_unq, lambda_dat[,t][i])
    Nt = c(Nt, Nt_i)
    Zt = c(Zt, Zt_i)
  }  

  lambda_dat$Nt = Nt
  lambda_dat$Zt = Zt
  lambda_dat$lambda = Nt/Zt
  if(is.null(beta_wt)) beta_wt = residuals(surv.fit, type="dfbeta")/dat$wt
  Nt_w = c(dat[,d])*dat$Ii_t
  Zt_w = c(dat$rel_hzd)*dat$Yi_t+beta_wt%*%t(t(c(dat$wt*dat$rel_hzd)*dat$Yi_t)%*%x.mtrx)
  lambda_wt = t((lambda_dat$Zt)^-1*(t(Nt_w)-lambda_dat$lambda*t(Zt_w)))
  colnames(lambda_wt) = lambda_dat[,t]
  return(list(lambda = lambda_dat$lambda,
              lambda_wt = lambda_wt,
              u = lambda_dat[,t])
  )
}
##########################################################################################
# Lambda_w is the FUNCTION to estimate cumulative baseline hazard at given time points   #
#                          and its derivative w.r.t. sample weight                       #
# INPUT                                                                                  #
#  lambda_out - an object of output of function lambda_w                                 #
#  t_star     - a vector of time points for cumulative baseline hazard inference         #
# OUPUT                                                                                  #
#  Lambda    - a vector of estimates of cumulative baseline hazards at given time pionts #
#  Lambda_wt - a matrix of derivative of cumulative baseline hazardS w.r.t. sample weight#
#              at given time points                                                      #
##########################################################################################
Lambda_w = function(lambda_out, t_star){
  u = lambda_out$u
  lambda_wt = lambda_out$lambda_wt[,order(u)]
  lambda    = lambda_out$lambda[order(u)]
  u = sort(u)
  t_diff = outer(u, t_star, FUN="-")
  t_diff[t_diff>0]=-10^6
  Lambda = cumsum(lambda)[apply(t_diff, 2, which.max)]
  Lambda_wt = cbind(NULL,rowCumsums(lambda_wt)[,apply(t_diff, 2, which.max)])
  colnames(Lambda_wt) = paste("Lambda(t=", t_star, ")", sep="")
  return(list(Lambda = Lambda,
              Lambda_wt = Lambda_wt))
}


Lambda_t = function(t, alpha, t1, x, betas){
  n_beta = nrow(betas)
  t.id =1:length(t);t.id=t.id[order(t)]
  t = sort(t)
  t.l = t[t<=t1];  t.h = t[t>t1]
  Lambda0 = outer(t^alpha, exp(betas[1,]))
  Lambda = outer(rowSums(exp(x%*%betas)), t.l^alpha)
  rr_y0 = exp(x[,-n_beta]%*%betas[-n_beta,])
  Lambda=cbind(Lambda, outer(rowSums(rr_y0), t.h^alpha)+t1^alpha*rowSums(rr_y0*(exp(outer(x[,n_beta],betas[n_beta,]))-1)))
  Lambda = Lambda[,order(t.id)]
  return(list(Lambda0 = Lambda0, Lambda = t(Lambda)))
}

#Lambda_tmax.0 = Lambda_t(t=t_star, alpha, t1, x=cbind(1, x0, 0), betas=cbind(c(beta1, gamma1), c(beta2, gamma2)))
#Lambda_tmax.1 = Lambda_t(t=t_star, alpha, t1, x=cbind(1, x0, 1), betas=cbind(c(beta1, gamma1), c(beta2, gamma2)))
#
#plot(t_star^alpha*exp(beta1[1]), pop_Lambda0_1$cum_hzd*1.15895, 
#     ylab=expression(paste(hat(Lambda)[0],"Breslow")), xlab=expression(t^alpha*e^beta[0]));abline(0,1,col="red")
##plot(t_star^alpha*exp(beta1[1]), survfit(cox_pop1.beta)$cumhaz)
#plot(pop_Lambda0_1$cum_hzd[-c(1:2)], survfit(cox_pop1.beta)$cumhaz[-1]);abline(0, 1, col="red")
#t_diff = t_star-c(0, t_star[-length(t_star)])
#cum_S = t_diff*(exp(-Lambda_tmax.1)-exp(-Lambda_tmax.0))
#plot(t_star, cumsum(cum_S[,1]), type="l")
#colSums(cum_S[t_star<=15,])
#colSums(cum_S[t_star<=100,])
#colSums(cum_S[t_star<=200,])
#colSums(cum_S[t_star<=400,])

Lambda_inf = function(Lambda0, t, t1, x, y, betas, Lambda0_w = NULL, beta_w=NULL, deviate=F){
  n_s = nrow(Lambda0_w)
  n = nrow(beta_w)
  n_t=n-n_s
  n_beta = nrow(betas)
  x = matrix(x, ncol=n_beta-1)
  n_x = nrow(x)
  K = ncol(betas)
  t1.indx = which.min(abs(t-t1))
  t1=t[t1.indx]
  t.l.indx = which(t<=t1);  t.h.indx = which(t>t1)
  rr = exp(cbind(x, y)%*%betas)
  rr_y0 = exp(x%*%betas[-n_beta,])
  rr_y = exp(outer(y,betas[n_beta,]))
  Lambda=0;Lambda_w=0
  for(k in 1:K){
    if((betas[n_beta,k]==0)|(y==0)){
      Lambda.k = outer(rr[,k], Lambda0[,k])
      if(deviate){
        Lambda_w.k=array(0,c(n, length(t), nrow(x)))
        outer.row = function(i) outer(c(beta_w[1:n_s,-n_beta,k]%*%c(x[i,])), rr[i,k]*Lambda0[,k])
        Lambda_w.k[1:n_s,,] = outer(Lambda0_w[,,k], rr[,k])+
          array(sapply(1:nrow(x), FUN=outer.row), c(n_s,length(t),n_x))
        Lambda_w = Lambda_w+Lambda_w.k
      }
    }else{
      Lambda.k = outer(rr[,k], Lambda0[t.l.indx,k])
      Lambda.k = cbind(Lambda.k, outer(rr_y0[,k], Lambda0[t.h.indx,k])+Lambda0[t1.indx,k]*rr_y0[,k]*(rr_y[,k]-1))
      if(deviate){
        Lambda_w.k=array(0,c(n, length(t), nrow(x)))
        outer.row = function(rr,t.indx,i) outer(c(beta_w[1:n_s,-n_beta,k]%*%c(x[i,])), rr[i,k]*Lambda0[t.indx,k])
        Lambda_w.k[1:n_s,1:length(t.l.indx),]= outer(Lambda0_w[,t.l.indx,k], rr[,k])+
          array(sapply(1:nrow(x), FUN=outer.row, t.indx=t.l.indx, rr=rr), 
                c(n_s,length(t.l.indx),n_x))
        Lambda_w.k[1:n_t+n_s,1:length(t.l.indx),] = outer(outer(beta_w[1:n_t+n_s,n_beta,k], Lambda0[t.l.indx,k]), rr[,k])
        
        if(length(t.l.indx)<length(t)){
          Lambda_w.k.h.t1 = outer(Lambda0_w[,t1.indx,k], (rr_y[,k]-1)*rr_y0[,k])+
            Lambda0[t1.indx,k]*(rr_y[,k]-1)*(beta_w[1:n_s,-n_beta,k]%*%t(rr_y0[,k]*x))
          Lambda_w.k.h.t1 = array(Lambda_w.k.h.t1, c(n_s, nrow(x),1))
          Lambda_w.k.h.t1 = aperm(array(Lambda_w.k.h.t1[,,rep(1, length(t.h.indx))], c(n_s, n_x, length(t.h.indx))), 
                                  c(1,3,2))
          Lambda_w.k.h.tj = outer(Lambda0_w[,t.h.indx,k], rr_y0[,k])+
            array(sapply(1:nrow(x), FUN=outer.row, t.indx=t.h.indx, rr=rr_y0), 
                  c(n_s,length(t.h.indx),n_x))
          Lambda_w.k[1:n_s,1:length(t.h.indx)+length(t.l.indx),] = Lambda_w.k.h.t1+Lambda_w.k.h.tj
          Lambda_w.k[1:n_t+n_s,1:length(t.h.indx)+length(t.l.indx),] = 
            aperm(array(array(Lambda0[t1.indx,k]*rr_y[,k]*outer(beta_w[1:n_t+n_s,n_beta,k], rr_y0[,k]), 
                              c(n_t,n_x,1))[,,rep(1,length(t.h.indx))], c(n_t,n_x,length(t.h.indx))), 
                  c(1, 3, 2))
        }
        Lambda_w = Lambda_w+Lambda_w.k
      }
    } 
    Lambda = Lambda+t(Lambda.k)
    #print(k)
  }
  return(list(Lambda = Lambda, Lambda_w = Lambda_w))
}


Delta_mu_inf = function(Lambda0, t, t1, x, betas, Lambda0_w = NULL, beta_w=NULL, Delta_mu = T, deviate=F){
  n_s = nrow(Lambda0_w)
  n = nrow(beta_w)
  n_t=n-n_s
  n_beta = nrow(betas)
  x = matrix(x, ncol=n_beta-1)
  n_x = nrow(x)
  K = ncol(betas)
  t1.indx = which.min(abs(t-t1))
  t1=t[t1.indx]
  t.l.indx = which(t<=t1);  t.h.indx = which(t>t1)
  Lambda.inf_in = function(y, deviate){
    rr = exp(cbind(x, y)%*%betas)
    rr_y0 = exp(x%*%betas[-n_beta,])
    rr_y = exp(outer(y,betas[n_beta,]))
    Lambda=0;Lambda_w=0
    for(k in 1:K){
      if((betas[n_beta,k]==0)|(y==0)){
        Lambda.k = outer(rr[,k], Lambda0[,k])
        if(deviate){
          Lambda_w.k=array(0,c(n, length(t), nrow(x)))
          outer.row = function(i) outer(c(beta_w[1:n_s,-n_beta,k]%*%c(x[i,])), rr[i,k]*Lambda0[,k])
          Lambda_w.k[1:n_s,,] = outer(Lambda0_w[,,k], rr[,k])+
            array(sapply(1:nrow(x), FUN=outer.row), c(n_s,length(t),n_x))
          Lambda_w = Lambda_w+Lambda_w.k
        }
      }else{
        Lambda.k = outer(rr[,k], Lambda0[t.l.indx,k])
        Lambda.k = cbind(Lambda.k, outer(rr_y0[,k], Lambda0[t.h.indx,k])+Lambda0[t1.indx,k]*rr_y0[,k]*(rr_y[,k]-1))
        if(deviate){
          Lambda_w.k=array(0,c(n, length(t), nrow(x)))
          outer.row = function(rr,t.indx,i) outer(c(beta_w[1:n_s,-n_beta,k]%*%c(x[i,])), rr[i,k]*Lambda0[t.indx,k])
          Lambda_w.k[1:n_s,1:length(t.l.indx),]= outer(Lambda0_w[,t.l.indx,k], rr[,k])+
            array(sapply(1:nrow(x), FUN=outer.row, t.indx=t.l.indx, rr=rr), 
                  c(n_s,length(t.l.indx),n_x))
          Lambda_w.k[1:n_t+n_s,1:length(t.l.indx),] = outer(outer(beta_w[1:n_t+n_s,n_beta,k], Lambda0[t.l.indx,k]), rr[,k])
          if(length(t.l.indx)<length(t)){
            Lambda_w.k.h.t1 = outer(Lambda0_w[,t1.indx,k], (rr_y[,k]-1)*rr_y0[,k])+
              Lambda0[t1.indx,k]*(rr_y[,k]-1)*(beta_w[1:n_s,-n_beta,k]%*%t(rr_y0[,k]*x))
            Lambda_w.k.h.t1 = array(Lambda_w.k.h.t1, c(n_s, nrow(x),1))
            Lambda_w.k.h.t1 = aperm(array(Lambda_w.k.h.t1[,,rep(1, length(t.h.indx))], c(n_s, n_x, length(t.h.indx))), 
                                    c(1,3,2))
            Lambda_w.k.h.tj = outer(Lambda0_w[,t.h.indx,k], rr_y0[,k])+
              array(sapply(1:nrow(x), FUN=outer.row, t.indx=t.h.indx, rr=rr_y0), 
                    c(n_s,length(t.h.indx),n_x))
            Lambda_w.k[1:n_s,1:length(t.h.indx)+length(t.l.indx),] = Lambda_w.k.h.t1+Lambda_w.k.h.tj
            Lambda_w.k[1:n_t+n_s,1:length(t.h.indx)+length(t.l.indx),] = 
              aperm(array(array(Lambda0[t1.indx,k]*rr_y[,k]*outer(beta_w[1:n_t+n_s,n_beta,k], rr_y0[,k]), 
                                c(n_t,n_x,1))[,,rep(1,length(t.h.indx))], c(n_t,n_x,length(t.h.indx))), 
                    c(1, 3, 2))
          }
          Lambda_w = Lambda_w+Lambda_w.k
        }
      } 
      Lambda = Lambda+t(Lambda.k)
      #print(k)
    }
    return(list(Lambda = Lambda, Lambda_w = Lambda_w))
  }
  Lambda_out1=Lambda.inf_in(1, deviate=deviate)
  Lambda_out0=Lambda.inf_in(0, deviate=deviate)
  S1 = exp(-Lambda_out1$Lambda)
  S0 = exp(-Lambda_out0$Lambda)
  if(Delta_mu&deviate){
    t_diff = t-c(0, t[-length(t)])
    S1_w = -array(array(S1,c(1,length(t),n_x))[rep(1,n),,], c(n,length(t),n_x))*Lambda_out1$Lambda_w
    S0_w = -array(array(S0,c(1,length(t),n_x))[rep(1,n),,], c(n,length(t),n_x))*Lambda_out0$Lambda_w
    return(list(Delta_mu = t_diff*(S1-S0),
                Delta_mu_w = array(array(t_diff,c(1,length(t),n_x))[rep(1,n),,], c(n,length(t),n_x))*(S1_w-S0_w),
                S1 = S1, S0 = S0,
                S1_w=S1_w ,S0_w=S0_w,
                #S=cbind(S1[t1.indx,], S0[t1.indx,]),
                Lambda = cbind(Lambda_out1$Lambda,Lambda_out0$Lambda),
                Lambda1_w = Lambda_out1$Lambda_w,
                Lambda0_w = Lambda_out0$Lambda_w
    ))
  }else if(Delta_mu){
    t_diff = t-c(0, t[-length(t)])
    return(list(Delta_mu = t_diff*(S1-S0),
                S1 = S1, S0 = S0,
                #S=cbind(S1[t1.indx,], S0[t1.indx,]),
                Lambda = cbind(Lambda_out1$Lambda,
                               Lambda_out0$Lambda)
    ))
  }else{return(list(S1 = S1, S0 = S0))}
}



Delta_mu_inf.all = function(Lambda0, t, t1, x, betas, Lambda0_w = NULL, beta_w=NULL, Delta_mu = T, deviate=F){
  if(!is.null(Lambda0_w)){
    n_s = nrow(Lambda0_w)
    n = nrow(beta_w)
    n_t=n-n_s}
  x = matrix(x, ncol=n_beta-1)
  n_x = nrow(x)
  K = ncol(betas)-1
  t1.indx = which.min(abs(t-t1))
  t1=t[t1.indx]
  t.l.indx = which(t<=t1);  t.h.indx = which(t>t1)
  rr_y0 = exp(x%*%betas[-n_beta,])
  rr_y1 = cbind(exp(cbind(x,1)%*%betas[,1]), rr_y0[,-1])
  Lambda.y1_tau1 = matrix(rowSums(sapply(1:K, function(k) outer(rr_y1[,k], Lambda0[t.l.indx,k]))),n_x,)
  S1_tau1 = exp(-Lambda.y1_tau1)
  Lambda.y0_tau1 = matrix(rowSums(sapply(1:K, function(k) outer(rr_y0[,k], Lambda0[t.l.indx,k]))),n_x,)
  S0_tau1 = exp(-Lambda.y0_tau1)
  term1 = colSums(t(S1_tau1-S0_tau1)*
                    (t.l.indx-c(0,t.l.indx[-length(t.l.indx)])))
  term2.1 = (exp(-rr_y1[,1]*Lambda0[t1.indx,1])-exp(-rr_y0[,1]*Lambda0[t1.indx,1]))*
            exp(-rr_y0[,2]*Lambda0[t1.indx,2]+rr_y0[,K+1]*Lambda0[t1.indx,K+1])
  term2.2 = colSums(t(exp(-outer(rr_y0[,K+1], Lambda0[t.h.indx,K+1])))*
                      (c(t.h.indx)-c(0,t.h.indx[-length(t.h.indx)])))
  return(list(S1 = S1, S0 = S0))
}
  

absR_w = function(beta_est, Lambda, x0, beta_w=NULL, Lambda_w=NULL){
  n_beta = length(beta_est)
  x0 = matrix(x0, ncol=n_beta)
  beta_est = as.matrix(beta_est)
  rr = exp(x0%*%beta_est)
  absR = 1-exp(-outer(c(Lambda), c(rr)))
  if((!is.null(beta_w))&(!is.null(Lambda_w))){
    if((length(Lambda)>1)&(nrow(x0)==1)){
      n_t = length(Lambda)
      n = nrow(beta_w)
      beta_w_x0 = matrix(rep(beta_w%*%t(x0), n_t), n, n_t)
      absR_w = t(c((1-absR)*c(rr))*(Lambda*t(beta_w_x0)+t(Lambda_w)))
    }
    if((length(Lambda)==1)&(ncol(x0)>=1)){
      absR_w = t(c((1-absR)*rr)*t(Lambda*(beta_w%*%t(x0))+Lambda_w))
    }
    
    return(list(absR = absR,
                absR_w = absR_w))
  }else(return(list(absR = absR)))
}

Delta_mu_3m = function(Lambda0, Lambda0.all, t, t1, x, betas){
  Delta_mu_t1 = Delta_mu_inf(Lambda0=Lambda0, t=t[t<=t1], t1=t1, x=x, betas=betas[,1:2])
  nbeta = nrow(betas)
  rr0 = exp(cbind(x, 0)%*%betas[,1:2]); 
  Lmd12.0 = t(t(rr0)*Lambda0[t1,])
  Lmd1.1 = Lambda0[t1,1]*rr0[,1]*exp(betas[nbeta,1])
  Lmd.0 = outer(c(Lambda0.all[t>t1]), c((exp(cbind(x, 0)%*%betas[,-(1:2)]))))
  Delta_mu_tau = rbind(Delta_mu_t1$Delta_mu,
                       ((exp(-Lmd1.1))-(exp(-Lmd12.0[,1])))*exp(-Lmd12.0[,2]+Lmd.0[1,])*exp(-Lmd.0))
  return(list(Delta_mu=Delta_mu_tau))
}
