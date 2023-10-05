library(CARMA)

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

CARMA.results<-CARMA_fixed_sigma(z.list,ld.list,lambda.list=lambda.list,
                                 outlier.switch=T)
CARMA.results=CARMA.results[[1]]

c0=CARMA.results
c0$`Credible set`
#S3=runsusie(D3)
library(susieR)
SSr=susie_rss(z = z.list[[1]],R=ld.list[[1]],n=1000)
s0=SSr



#####
z=D3$beta/sqrt(D3$varbeta)

z[89]=-z[89]

z.list[[1]]<-z

CARMA.results<-CARMA_fixed_sigma(z.list,ld.list,lambda.list=lambda.list,
                                 outlier.switch=T)
CARMA.results=CARMA.results[[1]]
CARMA.results$`Credible set`

SSr=susie_rss(z = z.list[[1]],R=ld.list[[1]],n=1000)
apply(SSr$alpha,MAR=1,FUN=function(x){y=1:length(x); ind=order(x,decreasing = TRUE); x[ind[1:2]]})

#####
z=D3$beta/sqrt(D3$varbeta)

lll=D3$LD
la=lll[105,]
la=la[-105]
za=z[105]

z=z[-105]
lll=lll[-105,-105]

z=c(z,za,za+10)
z.list[[1]]<-z

new_l=array(1,c(501,501))
new_l[1:499,1:499]=lll
new_l[501,1:499]=new_l[500,1:499]=la
new_l[1:499,501]=new_l[1:499,500]=la

ld.list[[1]]<-new_l


CARMA.results<-CARMA_fixed_sigma(z.list,ld.list,lambda.list=lambda.list,
                                 outlier.switch=T)
CARMA.results=CARMA.results[[1]]
CARMA.results$`Credible set`
CARMA.results$Outliers

SSr=susie_rss(z = z.list[[1]],R=ld.list[[1]],n=1000)
apply(SSr$alpha,MAR=1,FUN=function(x){y=1:length(x); ind=order(x,decreasing = TRUE); x[ind[1:5]]})

ind_o=sort(CARMA.results$Outliers[,1])
SSr=susie_rss(z = z.list[[1]][-ind_o],R=ld.list[[1]][-ind_o,-ind_o],n=1000)
apply(SSr$alpha,MAR=1,FUN=function(x){y=1:length(x); ind=order(x,decreasing = TRUE); x[ind[1:5]]})




#####
z=D3$beta/sqrt(D3$varbeta)

ind=sample(1:500,10)
z[ind]=z[ind]+rnorm(length(ind),sd=2)

z.list[[1]]<-z

ld.list[[1]]<-D3$LD
CARMA.results<-CARMA_fixed_sigma(z.list,ld.list,lambda.list=lambda.list,
                                 outlier.switch=T)
CARMA.results=CARMA.results[[1]]
CARMA.results$`Credible set`
CARMA.results$Outliers

SSr=susie_rss(z = z.list[[1]],R=ld.list[[1]],n=1000)
apply(SSr$alpha,MAR=1,FUN=function(x){y=1:length(x); ind=order(x,decreasing = TRUE); x[ind[1:2]]})


ind_o=sort(CARMA.results$Outliers[,1])
SSr=susie_rss(z = z.list[[1]][-ind_o],R=ld.list[[1]][-ind_o,-ind_o],n=1000)
apply(SSr$alpha,MAR=1,FUN=function(x){y=1:length(x); ind=order(x,decreasing = TRUE); x[ind[1:5]]})


#####


