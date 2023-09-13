#######
#functions
#######

logbf_to_pp_vector=function(bf,pi=1e-4) {
  library(corrcoverage)
  n=length(bf)
  if(length(pi)==1) { # scalar pi
    if(pi > 1/n)
      pi=1/n
    pi=rep(pi,n)
  }
  if(any(pi == 0)) { # avoid division by zero
    pi[ pi==0 ] = 1e-16
    pi=pi/sum(pi)
  }
  priors=log(pi)
  denom=logsum(bf + priors)
  exp(bf + priors - denom)
}

approx.bf.estimates <- function (z, V, type="quant", suffix=NULL, sdY=1) {
  sd.prior <- if (type == "quant") { 0.15*sdY } else { 0.2 }
  r <- sd.prior^2/(sd.prior^2 + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  ret <- data.frame(V, z, r, lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}

pp_to_logbf_scale_vector <- function(pp, pi=1e-4,abf) {
  
  n <- length(pp)
  
  if (length(pi) == 1) {
    if (pi > 1/n)
      pi <- 1/n
    pi <- rep(pi, n)
  }
  
  if (any(pi == 0)) {
    pi[pi == 0] <- 1e-16
    pi <- pi / sum(pi)
  }
  
  log_prior <- log(pi)
  
  logbf0 <- log(pp) - log_prior
  
  #
  v=logbf0
  ind=which.min(pp)
  #ind=which.max(pp)
  delta=v[ind]-abf[ind]
  #ind=order(v)[1:10]
  #delta=mean(v[ind])-mean(abf[ind])
  v=v-delta
  logbf0=v
  
  return(logbf0)
}

finemap_logBF=function(pip,maxk=10,m){
  
  prior_Pk=function(k,m,maxk=10){
    
    non_normolized_pk=function(k,m){
      p=choose(n = m,k = k) * (1/m)^k * ((m-1)/m)^(m-k)
      p
    }
    
    pki=NULL
    for (i in 1:maxk){
      pki=c(pki,non_normolized_pk(i,m))
    }
    
    pki_norm=pki/sum(pki)
    
    pki_norm[k]
  }
  
  #prior_Pk(k=2,m=1000,maxk=10)
  
  p_gamma_l_1=function(m,maxk=10){
    k=1
    out=0
    for (k in 1:maxk){
      out=out+prior_Pk(k=k,m=m,maxk=maxk)*k/m
    }
    out
  }
  
  p_gamma_l_0=function(m,maxk=10){
    k=1
    out=0
    for (k in 1:maxk){
      out=out+prior_Pk(k=k,m=m,maxk=maxk)*(m-k)/m
    }
    out
  }
  
  prior0=p_gamma_l_1(m=m,maxk=maxk)/p_gamma_l_0(m=m,maxk=maxk)
  logBFs=log((pip/(1-pip))/prior0)
  logBFs
}

##### FULL
library(data.table)
setwd("~/Projects/susie_coloc_data/01_marc_test/")
source("~/Projects/Sanger_OT_susie_coloc/02_CARMA_test/priors.R")

t2=fread("finemapping_cd4_test_v2.txt",data.table=FALSE)

#####
ensg_list=unique((t2$ENSG))
i=1
ensg=ensg_list[i]

for (i in 1:length(ensg_list)){
  ensg=ensg_list[i]
  ind=which(t2$ENSG==ensg)
  st=t2[ind,]
  z=st$beta/st$beta_se
  abf=approx.bf.estimates(z=z,V=st$beta_se^2,sdY=1, type="quant")
  cs=1
  for (cs in 1:5){
    real_bfs=st[,paste0("lbf_cs_",cs)]
    pip=logbf_to_pp_vector(real_bfs,pi=1/length(real_bfs))
    bf_v=pp_to_logbf_scale_vector(pip,pi=1/nrow(st),abf=abf$lABF)
    logBFs=finemap_logBF(pip,maxk=10,m=nrow(st))
    
    t2[ind,paste0("pred_v1_lbf_cs_",cs)]=bf_v
    t2[ind,paste0("pred_v2_lbf_cs_",cs)]=logBFs
    
    if (max(logBFs)>max(bf_v)) t2[ind,paste0("pred_max_lbf_cs_",cs)]=logBFs
    if (max(logBFs)<max(bf_v)) t2[ind,paste0("pred_max_lbf_cs_",cs)]=bf_v
  }
}

cs=1
for (cs in 1:5){
  t2_cs=t2[which(t2$SusieRss_CS==cs),]
  pip=t2[which(t2$SusieRss_CS==cs),"SusieRss_pip"]
  
  #table(is.na(t2_cs[,paste0("pred_v2_lbf_cs_",cs)]))
  plot(t2_cs[,paste0("lbf_cs_",cs)],t2_cs[,paste0("pred_v1_lbf_cs_",cs)])
  #plot(t2_cs[,paste0("lbf_cs_",cs)],t2_cs[,paste0("pred_v2_lbf_cs_",cs)])
  #plot(t2_cs[,paste0("lbf_cs_",cs)],t2_cs[,paste0("pred_max_lbf_cs_",cs)])
  ind=which(is.infinite(t2_cs[,paste0("pred_v2_lbf_cs_",cs)]))
  print("CS:")
  print(cs)
  print(nrow(t2_cs))
  print(cor(t2_cs[,paste0("lbf_cs_",cs)],t2_cs[,paste0("pred_v1_lbf_cs_",cs)]))
  print(mean(t2_cs[,paste0("lbf_cs_",cs)]-t2_cs[,paste0("pred_v1_lbf_cs_",cs)]))
  print(cor(t2_cs[-ind,paste0("lbf_cs_",cs)],t2_cs[-ind,paste0("pred_v2_lbf_cs_",cs)]))
  print(cor(t2_cs[-ind,paste0("lbf_cs_",cs)],t2_cs[-ind,paste0("pred_max_lbf_cs_",cs)]))
}

#####
##### CROPPPPPPPP
library(data.table)
setwd("~/Projects/susie_coloc_data/01_marc_test/")
source("~/Projects/Sanger_OT_susie_coloc/02_CARMA_test/priors.R")

t2=fread("finemapping_cd4_test_v2.txt",data.table=FALSE)

#####
t2=t2[!is.na(t2$SusieRss_CS),]
ensg_list=unique((t2$ENSG))
i=1
ensg=ensg_list[i]

for (i in 1:length(ensg_list)){
  ensg=ensg_list[i]
  ind=which(t2$ENSG==ensg)
  st=t2[ind,]
  z=st$beta/st$beta_se
  abf=approx.bf.estimates(z=z,V=st$beta_se^2,sdY=1, type="quant")
  cs=1
  for (cs in 1:5){
    real_bfs=st[,paste0("lbf_cs_",cs)]
    pip=logbf_to_pp_vector(real_bfs,pi=1/length(real_bfs))
    bf_v=pp_to_logbf_scale_vector(pip,pi=1/nrow(st),abf=abf$lABF)
    logBFs=finemap_logBF(pip,maxk=10,m=nrow(st))
    
    t2[ind,paste0("pred_v1_lbf_cs_",cs)]=bf_v
    t2[ind,paste0("pred_v2_lbf_cs_",cs)]=logBFs
    
    if (sum(is.nan(logBFs))==0){
     if (max(logBFs)>max(bf_v)) t2[ind,paste0("pred_max_lbf_cs_",cs)]=logBFs
     if (max(logBFs)<max(bf_v)) t2[ind,paste0("pred_max_lbf_cs_",cs)]=bf_v
    }
  }
}

cs=1
for (cs in 1:5){
  t2_cs=t2[which(t2$SusieRss_CS==cs),]
  pip=t2[which(t2$SusieRss_CS==cs),"SusieRss_pip"]
  
  #table(is.na(t2_cs[,paste0("pred_v2_lbf_cs_",cs)]))
  plot(t2_cs[,paste0("lbf_cs_",cs)],t2_cs[,paste0("pred_v1_lbf_cs_",cs)])
  #plot(t2_cs[,paste0("lbf_cs_",cs)],t2_cs[,paste0("pred_v2_lbf_cs_",cs)])
  #plot(t2_cs[,paste0("lbf_cs_",cs)],t2_cs[,paste0("pred_max_lbf_cs_",cs)])
  ind=which(is.infinite(t2_cs[,paste0("pred_v2_lbf_cs_",cs)]))
  #ind=which(t2_cs[,paste0("pred_v1_lbf_cs_",cs)]<50)
  print("CS:")
  print(cs)
  print(nrow(t2_cs))
  print(cor(t2_cs[,paste0("lbf_cs_",cs)],t2_cs[,paste0("pred_v1_lbf_cs_",cs)]))
  print(mean(t2_cs[,paste0("lbf_cs_",cs)]-t2_cs[,paste0("pred_v1_lbf_cs_",cs)]))
  print(cor(t2_cs[-ind,paste0("lbf_cs_",cs)],t2_cs[-ind,paste0("pred_v2_lbf_cs_",cs)]))
  print(cor(t2_cs[-ind,paste0("lbf_cs_",cs)],t2_cs[-ind,paste0("pred_max_lbf_cs_",cs)]))
}

#####

ensg_list=unique((t2$ENSG))
i=1
ensg=ensg_list[i]
out_bf_pred=NULL
out_bf_pred_v2=NULL
out_bf_real=NULL
for (i in 1:length(ensg_list)){
  ensg=ensg_list[i]
  st=t2[t2$ENSG==ensg,]
  z=st$beta/st$beta_se
  abf=approx.bf.estimates(z=z,V=st$beta_se^2,sdY=1, type="quant")
  #st$SusieRss_CS[is.na(st$SusieRss_CS)]=1
  n_cs=unique(st$SusieRss_CS)  
  n_cs=n_cs[!is.na(n_cs)]
  cs=n_cs[1]
  for (cs in n_cs){
    ind=which(st$SusieRss_CS==cs)
    st_cs=st[ind,]
    bf_v=pp_to_logbf_scale_vector(st_cs$SusieRss_pip,pi=1/nrow(st_cs),abf=abf$lABF[ind])
    out_bf_pred=c(out_bf_pred,bf_v)
    out_bf_real=c(out_bf_real,st_cs[,paste0("lbf_cs_",cs)])
    
    pip=st_cs$SusieRss_pip
    logBFs=finemap_logBF(pip,maxk=10,m=nrow(st))
    out_bf_pred_v2=c(out_bf_pred_v2,logBFs)
  }
}
plot(out_bf_pred,out_bf_real)
plot(out_bf_pred_v2,out_bf_real)
cor(out_bf_pred,out_bf_real)
cor(out_bf_pred_v2,out_bf_real)


