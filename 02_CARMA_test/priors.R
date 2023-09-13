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




