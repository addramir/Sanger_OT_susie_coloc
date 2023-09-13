setwd("~/Projects/susie_coloc_data/fg_example/")
t1=fread("finemapping_full_finngen_R8_AB1_AMOEBIASIS.FINEMAP.snp",data.table=F)
t2=fread("finemapping_full_finngen_R8_AB1_AMOEBIASIS.SUSIE.snp",data.table=F)

table(t1$rsid==t2$rsid)

plot(t1$prob,t2$cs_specific_prob)


plot(t1$log10bf,t2$sd1)


prior0=p_gamma_l_1(m=nrow(t1),maxk=10)/p_gamma_l_0(m=nrow(t1),maxk=10)
pip=t1$prob
logBFs=log10((pip/(1-pip))/prior0)

plot(t1$log10bf,logBFs)
summary(t1$log10bf-logBFs)

logBFs=log10(exp(finemap_logBF(pip,maxk=10)))


t1=fread("finemapping_full_finngen_R8_AB1_EBV.FINEMAP.snp",data.table=F)
t2=fread("finemapping_full_finngen_R8_AB1_EBV.SUSIE.snp",data.table=F)

table(t1$rsid==t2$rsid)

plot(t1$prob,t2$cs_specific_prob)


plot(t1$log10bf,t2$sd1)

prior0=p_gamma_l_1(m=nrow(t1),maxk=10)/p_gamma_l_0(m=nrow(t1),maxk=10)
pip=t1$prob
logBFs=log10((pip/(1-pip))/prior0)

plot(t1$log10bf,logBFs)
summary(t1$log10bf-logBFs)
