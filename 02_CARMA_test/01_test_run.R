setwd("~/Projects/Sanger_OT_susie_coloc/02_CARMA_test/")
load("APOE_CAMRA_data.RData")

ind=sample(1:481,200)
ind=1:nrow(sumstat)

z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstat$Z[ind]/3
ld.list[[1]]<-as.matrix(ld)[ind,ind]
lambda.list[[1]]<-1

library("CARMA")
library(Matrix)
source("CARMA_fixed_sigma.R")
CARMA.results<-CARMA_fixed_sigma(z.list,ld.list,lambda.list=lambda.list,
                                 outlier.switch=T)

CARMA.results=CARMA.results[[1]]

CSs=CARMA.results$`Credible set`[[2]]

prior0=p_gamma_l_1(m=481,maxk=10)/p_gamma_l_0(m=481,maxk=10)
pip=CARMA.results$all.C.list[[3]]
logBFs=log((pip/(1-pip))/prior0)







