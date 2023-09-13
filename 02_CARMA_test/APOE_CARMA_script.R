library(CARMA)
library(data.table)
library(magrittr)
library(dplyr)
library(devtools)
library(R.utils)
library(Matrix)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(glmnet)

setwd('~/Desktop')
sumstat <- fread('/Users/hn9/Documents/GitHub/carmapy/tests/APOE_locus_sumstats.txt.gz',
                 sep = "\t", header = T, check.names = F, data.table = F,
                 stringsAsFactors = F)

ld <- fread('/Users/hn9/Documents/GitHub/carmapy/tests/APOE_locus_ld.txt.gz',
            sep = "\t", header = T, check.names = F, data.table = F,
            stringsAsFactors = F)


##### set up the input files for CARMA
sumstat$Z <- as.numeric(sumstat$Z)
#ld <- as.matrix(ld)
z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstat$Z
ld.list[[1]]<-as.matrix(ld)
lambda.list[[1]]<-1

CARMA.results<-CARMA_fixed_sigma(z.list,ld.list,lambda.list=lambda.list,
                                      outlier.switch=T)

sumstat.result = sumstat %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
  for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
    sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
  }
}
###### write the GWAS summary statistics with PIP and CS
fwrite(x = sumstat.result,
       file = "apoe_carma.txt.gz",
       sep = "\t", quote = F, na = "NA", row.names = F, col.names = T,
       compress = "gzip")