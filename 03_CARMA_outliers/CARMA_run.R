

source("C:/Projects/Sanger_OT_susie_coloc/02_CARMA_test/CARMA_orig.R")
source("C:/Projects/Sanger_OT_susie_coloc/02_CARMA_test/CARMA_orig_Rcpp.R")

library(coloc)
data(coloc_test_data)
attach(coloc_test_data)

z.list<-list()
ld.list<-list()
lambda.list<-list()
z=D3$beta/sqrt(D3$varbeta)
z.list[[1]]<-z
ld.list[[1]]<-D3$LD
lambda.list[[1]]<-1

ind=sample(1:500,20)
z.list[[1]]<-z[ind]
ld.list[[1]]<-D3$LD[ind,ind]

library(Matrix)


w.list=NULL
output.labels='.'
label.list=NULL
effect.size.prior='Spike-slab'
rho.index=0.99
BF.index=10
EM.dist='Logistic'
Max.Model.Dim=2e+5
all.iter=3
all.inner.iter=10
input.alpha=0
epsilon.threshold=1e-5
printing.log=F
num.causal=10
y.var=1
tau=0.04
outlier.switch=T
outlier.BF.index=1/3.2
prior.prob.computation='Logistic'



CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,
                                 outlier.switch=T)
CARMA.results=CARMA.results[[1]]