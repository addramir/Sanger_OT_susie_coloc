setwd("Projects/susie_coloc_data/finn_saige_cs/")

library(data.table)

df=fread("FINN_cs.csv",data.table=F)

df=df[df$is99_credset=="TRUE",]

dim(df)
length(df$lead_variant_id)
length(unique(df$lead_variant_id))
lst=unique(df$lead_variant_id)

out=array(NA,c(length(lst),2))
i=1
for (i in 1:length(lst)){
  lv=lst[i]
  dff=df[df$lead_variant_id==lv,]
  out[i,1]=nrow(dff)
  out[i,2]=sum(dff$postprob)
}

summary(out[,2])
ind=which(out[,2]<=1 & out[,2]>0.95)
out_finn=out[ind,]
summary(out_finn[,2])


######

df=fread("SAIGE_cs.csv",data.table=F)

df=df[df$is99_credset=="TRUE",]

dim(df)
length(df$lead_variant_id)
length(unique(df$lead_variant_id))
lst=unique(df$lead_variant_id)

out=array(NA,c(length(lst),2))
i=1
for (i in 1:length(lst)){
  lv=lst[i]
  dff=df[df$lead_variant_id==lv,]
  out[i,1]=nrow(dff)
  out[i,2]=sum(dff$postprob)
}

summary(out[,2])
ind=which(out[,2]<=1 & out[,2]>0.95)
out_saige=out[ind,]
summary(out_saige[,2])

length(out_saige[,1])
summary(out_saige[,1])
length(out_finn[,1])
summary(out_finn[,1])
t.test(out_saige[,1],out_finn[,1])









