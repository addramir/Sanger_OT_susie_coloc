setwd("~/Projects/Sanger_OT_susie_coloc/02_CARMA_test/")
load("APOE_CAMRA_data.RData")

ind=sample(1:481,481)
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
source("priors.R")
CARMA.results<-CARMA_fixed_sigma(z.list,ld.list,lambda.list=lambda.list,
                                 outlier.switch=T)

CARMA.results=CARMA.results[[1]]

CSs=CARMA.results$`Credible set`[[2]]

prior0=p_gamma_l_1(m=481,maxk=10)/p_gamma_l_0(m=481,maxk=10)
pip=CARMA.results$all.C.list[[3]]
logBFs=log((pip/(1-pip))/prior0)

#####
library(coloc)
data(coloc_test_data)
attach(coloc_test_data)

z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-D3$beta/sqrt(D3$varbeta)
ld.list[[1]]<-D3$LD
lambda.list[[1]]<-1

CARMA.results<-CARMA_fixed_sigma(z.list,ld.list,lambda.list=lambda.list,
                                 outlier.switch=T)
CARMA.results=CARMA.results[[1]]

#S3=runsusie(D3)
library(susieR)
SSr=susie_rss(z = z.list[[1]],R=ld.list[[1]],n=1000)


BFs=CARMA.results$all.C.list[[1]][[1]]
MS=CARMA.results$all.C.list[[1]][[2]]

df=CARMA.results
NUM_CS=length(df$`Credible set`[[2]])
i=1
for (i in NUM_CS){
  cs_snps=df$`Credible set`[[2]][[i]]
  cred=apply(MS[,cs_snps],MAR=1,sum)
  ind=which(cred>=1)
  log(sum(exp(BFs[ind]))/(length(ind)*500))
}


sum(exp(aa[ which(model.space[,i]==1)]))/prob.sum



