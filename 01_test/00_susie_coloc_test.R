library(coloc)
data(coloc_test_data)
attach(coloc_test_data)
par(mfrow=c(2,1))
plot_dataset(D3, main="Dataset D3")
plot_dataset(D4, main="Dataset D4")

my.res <- coloc.abf(dataset1=D3, dataset2=D4)
class(my.res)
my.res
sensitivity(my.res,"H4 > 0.9")
check_dataset(D3,req="LD")
check_dataset(D4,req="LD")

S3=runsusie(D3)
summary(S3)
S4=runsusie(D4)
summary(S4)

susie.res=coloc.susie(S3,S4)
print(susie.res$summary)

sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=D3,dataset2=D4)
sensitivity(susie.res,"H4 > 0.9",row=2,dataset1=D3,dataset2=D4)


########
susie.res0=coloc.susie(S3,S4)
susie.res0$summary

names(S3)
#
lst=c("lbf_variable","sets")
S3_test=S3[lst]
class(S3_test)="susie"
S4_test=S4[lst]
class(S4_test)="susie"

res_test=coloc.susie(S3_test,S4_test)
res_test$summary
table(susie.res0$results==res_test$results)


lst=c("pip","sets","alpha")
S3_test=S3[lst]
class(S3_test)="susie"
S4_test=S4[lst]
class(S4_test)="susie"

res_test=coloc.susie(S3_test,S4_test,back_calculate_lbf = TRUE)
res_test$summary


logbf_to_pp=function(bf,pi=1e-4) {
  library(corrcoverage)
  n=ncol(bf)
  if(length(pi)==1) { # scalar pi
    if(pi > 1/n)
      pi=1/n
    pi=rep(pi,n)
  }
  if(any(pi == 0)) { # avoid division by zero
    pi[ pi==0 ] = 1e-16
    pi=pi/sum(pi)
  }
  priors=matrix(log(pi),nrow(bf),ncol(bf),byrow=TRUE)
  denom=matrix(apply(bf + priors,1,logsum), nrow(bf), ncol(bf))
  exp(bf + priors - denom)
}


ppp=logbf_to_pp(S3$lbf_variable,1e-4)

v1=cbind(ppp[1,],S3$alpha[1,])
cor(v1)
plot(log(v1))

v1=cbind(ppp[1,],S3$alpha[1,])
cor(v1)
plot(log(v1))



 
pp_to_logbf <- function(pp, pi=1e-4) {
  n <- ncol(pp)
  
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
  
  error <-log_prior + logbf0- apply(logbf0 + log_prior, 1, logsum)-log(pp)
  
  if (sum(error)<1e-5){
    out=logbf0
  } else{
    out=logbf0
    print("Error!")
  }
  
  out=t(apply(out, 1, FUN=function(x){x-min(x)}))
  return(out)
}

bf_v=pp_to_logbf(ppp,pi=0.002)

plot(bf_v[1,],S3$lbf_variable[1,])

lm(S3$lbf_variable[1,]~bf_v[1,])



approx.bf.estimates <- function (z, V, type="quant", suffix=NULL, sdY=1) {
  sd.prior <- if (type == "quant") { 0.15*sdY } else { 0.2 }
  r <- sd.prior^2/(sd.prior^2 + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  ret <- data.frame(V, z, r, lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}

z=(D3$beta/sqrt(D3$varbeta))
abf=approx.bf.estimates(z=z,V=D3$varbeta,sdY=D3$sdY)
head(abf)
cor(abf$lABF,t(S3$lbf_variable))

plot(abf$lABF,S3$lbf_variable[2,])

ind=which(abf$lABF>=10)
cor(abf$lABF[ind],t(S3$lbf_variable[,ind]))
cor(abf$lABF[ind],(S3$lbf_variable[1,ind]+S3$lbf_variable[2,ind]))



DD=D4
SS=S4

z=(DD$beta/sqrt(DD$varbeta))
abf=approx.bf.estimates(z=z,V=DD$varbeta,sdY=DD$sdY)
head(abf)
cor(abf$lABF,t(SS$lbf_variable))

plot(abf$lABF,SS$lbf_variable[1,])

ind=which(abf$lABF>=10)
cor(abf$lABF[ind],t(S3$lbf_variable[,ind]))
cor(abf$lABF[ind],(S3$lbf_variable[1,ind]+S3$lbf_variable[2,ind]))



######
pp_to_logbf_scale <- function(pp, pi=1e-4,abf) {
  
  n <- ncol(pp)
  
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
  i=1
  for (i in 1:nrow(pp)){
    v=logbf0[i,]
    #ind=which.min(pp[i,])
    #delta=v[ind]-abf[ind]
    ind=order(v)[1:10]
    delta=mean(v[ind])-mean(abf[ind])
    v=v-delta
    logbf0[i,]=v
  }
  
  return(logbf0)
}

#####
DD=D3
SS=runsusie(DD)

z=(DD$beta/sqrt(DD$varbeta))
abf=approx.bf.estimates(z=z,V=DD$varbeta,sdY=DD$sdY)
ppp=logbf_to_pp(SS$lbf_variable,pi=0.002)
bf_v=pp_to_logbf_scale(ppp,pi=0.002,abf=abf$lABF)
lm(bf_v[1,]~SS$lbf_variable[1,])

test=function(d1,d2){
  s1=runsusie(d1)
  s2=runsusie(d2)
  gold=coloc.susie(s1,s2)
  
  DD=d1
  SS=runsusie(DD)
  z=(DD$beta/sqrt(DD$varbeta))
  abf=approx.bf.estimates(z=z,V=DD$varbeta,sdY=DD$sdY)
  ppp=logbf_to_pp(SS$lbf_variable,pi=0.002)
  bf_v=pp_to_logbf_scale(ppp,pi=0.002,abf=abf$lABF)
  
  lst=c("lbf_variable","sets")
  s1_test=s1[lst]
  class(s1_test)="susie"
  s1_test$lbf_variable=bf_v
 
  DD=d2
  SS=runsusie(DD)
  z=(DD$beta/sqrt(DD$varbeta))
  abf=approx.bf.estimates(z=z,V=DD$varbeta,sdY=DD$sdY)
  ppp=logbf_to_pp(SS$lbf_variable,pi=0.002)
  bf_v=pp_to_logbf_scale(ppp,pi=0.002,abf=abf$lABF) 
  s2_test=s2[lst]
  class(s2_test)="susie"
  s2_test$lbf_variable=bf_v
  
  res_test=coloc.susie(s1_test,s2_test)
  print("gold")
  print(gold$summary)
  print("test")
  print(res_test$summary)
}
test(D1,D2)
test(D1,D3)
test(D1,D4)
test(D2,D3)
test(D2,D4)
test(D3,D4)




