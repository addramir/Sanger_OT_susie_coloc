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
  #ind=which.min(pp)
  ind=which.max(pp)
  delta=v[ind]-abf[ind]
  #ind=order(v)[1:10]
  #delta=mean(v[ind])-mean(abf[ind])
  v=v-delta
  logbf0=v
  
  return(logbf0)
}

#####
library(data.table)
setwd("Projects/susie_coloc_data/01_marc_test/")
t2=fread("finemapping_cd4_test_v2.txt",data.table=FALSE)

ensg_list=unique((t2$ENSG))
i=1
ensg=ensg_list[i]
out_bf_pred=NULL
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
  }
}


