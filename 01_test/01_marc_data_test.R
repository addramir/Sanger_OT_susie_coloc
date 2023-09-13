library(coloc)
library(data.table)

setwd("~/Projects/susie_coloc_data/01_marc_test/")
t1=fread("CARMA_finemapping_cd4_test.txt",data.table=FALSE)
t2=fread("finemapping_cd4_test_v2.txt",data.table=FALSE)

head(t1)
head(t2)

table(t1$feature_id)
table(t2$feature_id)

table(t1$QTL%in%t2$QTL)

ind=match(t1$QTL,t2$QTL)
table(t2$QTL[ind]==t1$QTL)
t2_excl=t2[-ind,]
t2=t2[ind,]

table(t2$beta==t1$beta)

cor(t2$SusieRss_pip,t1$CARMA_PIP)
plot(t2$SusieRss_pip,t1$CARMA_PIP)

table(t1$CARMA_PIP>0.95)
table(t2$SusieRss_pip>0.95)
table(t2_excl$SusieRss_pip>0.95)

ind=which(t2$SusieRss_pip>0.1 | t1$CARMA_PIP>0.1)
cor(t2$SusieRss_pip[ind],t1$CARMA_PIP[ind])
plot(t2$SusieRss_pip[ind],t1$CARMA_PIP[ind])

colnames(t1)
colnames(t2)

