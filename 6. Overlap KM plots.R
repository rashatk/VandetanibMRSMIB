###-------SIGNATURE SCORES
#load metabric log intensity expression data 
metab=read.table("data_mrna_agilent_microarray.txt", header=T, sep="\t")
#matrix preparation and mean duplicates
df=subset(metab, duplicated(metab$Hugo_Symbol))
means=do.call(rbind,lapply(lapply(split(df,df$Hugo_Symbol),`[`,2:ncol(df)),colMeans))
meanrows=rownames(means)
library(dplyr)
singles2=filter(metab, !(Hugo_Symbol %in% meanrows))
rownames(singles2)=singles2$Hugo_Symbol
singles2=singles2[-c(1)]
means=data.frame(means)
names(means)
merged_df=rbind(singles2,means)
merged_df=merged_df[-c(1)]
#load clinical data and subset ER status
clin_pt=read.table("data_clinical_patient.txt", header=T, sep="\t")
ER_mb=clin_pt[c(1,7)]
rownames(ER_mb)=ER_mb$PATIENT_ID
ER_mb=ER_mb[-c(1)]
names(merged_df)=gsub("[.]", "-", names(merged_df))
ER_mb=clin_pt[c(1,7)]
rownames(ER_mb)=ER_mb$PATIENT_ID
ER_mb=ER_mb[-c(1)]
ERpos=subset(ER_mb, ER_mb$ER_IHC=="Positve")
pts=rownames(ERpos)
ERdata=subset(merged_df, names(merged_df)%in% pts)

###Median centroid score calculation
#subset by signature among ER+ patients

#----ALL_down
mtx=na.omit(ERdata[ALL_down,]) ##specify signature
subsetted=data.frame(t(mtx))
subsetted$score=apply(mtx,2,median)
centroid=subsetted[c("score")]

quantile(centroid$score)
centroid$median=ifelse(centroid$score>6.895341,"High","Low")
quantile(centroid$score, probs=c(0.3333,0.6666))
centroid$tertile=ifelse(centroid$score>8.009411,"High", ifelse(centroid$score<7.535284, "Low", "Mid"))
rownames(centroid) <- gsub("[.]", "-", rownames(centroid))
centpts=rownames(centroid)

##---divide clinical data for survival curves
names(clin_pt)
surv=clin_pt[c(1,14,15,23,24)]
rownames(surv)=surv$PATIENT_ID
surv=surv[-c(1)]
survptz=subset(surv, rownames(surv) %in% centpts)
pts2=rownames(survptz)
centroid=centroid %>% arrange(factor(rownames(centroid), levels = pts2))

##merge columns and omit those without a status
merged=cbind(centroid,survptz)

unique(merged$OS_STATUS)
unique(merged$RFS_STATUS)

#----selecting overall survival data
library(survminer)
library(survival)

survcolumns <- merged
survivv <- as.matrix(survcolumns)
unique(survcolumns$OS_STATUS)
survv <-sub("0:LIVING", "0", survivv)
survv2 <-sub("1:DECEASED", "1", survv)
survivaldf <- as.data.frame(survv2)
OS_MONTHS <- survivaldf$OS_MONTHS
OS_STATUS <- survivaldf$OS_STATUS
monthz <- as.numeric(OS_MONTHS)
status <- as.numeric(OS_STATUS)

#plot by median cutoff
survivaldf$median2=factor(survivaldf$median, levels=c('High','Low'))
fit <- survfit(Surv(monthz, status) ~ median2, data = survivaldf)
res.cox <- coxph(Surv(monthz, status) ~ median2, data = survivaldf)
summary(res.cox, conf.int=FALSE)
nrow(survivaldf)
ggsurvplot(fit, data = survivaldf, pval=TRUE, pval.method=TRUE, linetype = 1, size=1.5, title="Overlapping downregulated genes \n in METABRIC patients", subtitle="HR=0.74", break.x.by = 25, xlab = "Time (months)", ylab="Overall survival probability", conf.int = FALSE,  xlim=c(0, 200), legend.labs =c("High", "Low"), palette =c("#ca4781","#E7B800"), risk.table = FALSE, tables.theme = theme_cleantable(), tables.y.text = FALSE, font.subtitle = c(10, "black"),font.x=10, font.y=13, font.main=c(12,"bold"),font.tickslab=c(10),censor.size=0, font.legend=c(10,"bold"), pval.size=4, pval.method.coord=c(7, 0.24), pval.method.size=4, legend = c(0.14, 0.45))


#----selecting RFS data
survcolumns <- merged
survivv <- as.matrix(survcolumns)
unique(survcolumns$RFS_STATUS)
survv <-sub("0:Not Recurred", "0", survivv)
survv2 <-sub("1:Recurred", "1", survv)
survivaldf <- as.data.frame(survv2)
RFS_MONTHS <- survivaldf$RFS_MONTHS
RFS_STATUS <- survivaldf$RFS_STATUS
monthz <- as.numeric(RFS_MONTHS)
status <- as.numeric(RFS_STATUS)

#plot by median cutoff
survivaldf$median2=factor(survivaldf$median, levels=c('High','Low'))
fit <- survfit(Surv(monthz, status) ~ median2, data = survivaldf)
res.cox <- coxph(Surv(monthz, status) ~ median2, data = survivaldf)
summary(res.cox, conf.int=FALSE)
nrow(survivaldf)
ggsurvplot(fit, data = survivaldf, pval=TRUE, pval.method=TRUE, linetype = 1, size=1.5, title="Overlapping downregulated genes \n in METABRIC patients", subtitle="HR=0.60", break.x.by = 25, ylab="Recurrence free survival probability",xlab = "Time (months)", conf.int = FALSE,  xlim=c(0, 200), legend.labs =c("High", "Low"), palette =c("#ca4781", "#E7B800"), risk.table = FALSE, tables.theme = theme_cleantable(), tables.y.text = FALSE, font.subtitle = c(10, "black"),font.x=10, font.y=13, font.main=c(12,"bold"),font.tickslab=c(10),censor.size=0, font.legend=c(10,"bold"), pval.size=4, pval.method.coord=c(7, 0.24), pval.method.size=4, legend = c(0.14, 0.45))

