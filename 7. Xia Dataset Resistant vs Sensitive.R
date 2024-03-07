##Endocrine therapy treated tumors database analysis
library(NOISeq)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)

##Data preparation
groups=read.csv("scotland_groupdata.csv",header=TRUE,row.names=1)#Samples list prepared on excel from sheet provided in supplementary data by xia et al, with duplicates averaged and annotated by _m 
groups$double=ifelse(groups$subtype=="iR","Intrinsic",ifelse(groups$subtype=="aR", "Adaptive", ifelse(groups$subtype=="Sensitive Pre-adaptive","PMP",ifelse(groups$subtype=="Sensitive","Sensitive", ""))))
groups$subtype=ifelse(groups$subtype=="iR","Intrinsic",ifelse(groups$subtype=="aR", "Adaptive", ifelse(groups$subtype=="Sensitive Pre-adaptive","Sensitive",ifelse(groups$subtype=="Sensitive","Sensitive", ""))))

#Load gene matrix with mean biological replicates
matrix=read.table("meanbio_scotland.txt", header=T, sep="\t")#Gene matrix was prepared according to duplicate samples by calculating means of duplicate samples based on matrix from xia et al
names(matrix) <- gsub("[.]", "-", names(matrix))
names(matrix) <- gsub("X", "", names(matrix))

#Data filtering, normalization and log2 transformation
mat=matrix[rowSums(matrix) > 0,]
data=uqua(mat, long = 1000, lc = 0, k = 0)
my_data=data.frame(data)
dataz=log2(my_data)
names(dataz) <- gsub("[.]", "-", names(dataz))
names(dataz) <- gsub("X", "", names(dataz))


##--------Analysis
#Gene set scoring
overlapdown=read.csv("Overlaps_mRNAseq_DOWN_05082023.csv",header=TRUE,row.names = 1)
overlapup=read.csv("Overlaps_mRNAseq_UP_05082023.csv",header=TRUE,row.names = 1)
down_genes=overlapdown$ALL_down
up_genes=overlapup$ALL_up

###Median score calculation
#----down_genes
mtx=subset(dataz, rownames(dataz) %in% down_genes)
subsetted=mtx[c(1),]
subsetted=data.frame(t(subsetted))
subsetted$median=apply(mtx,2,median)
subsetted=subsetted[c(-1)]
centroid=subsetted
centroid$catt=groups$type
centroid$group=groups$subtype
centroid$double=groups$double

library(ggplot2)
library(ggpubr)

unique(centroid$double)
unique(centroid$catt)

x=subset(centroid, centroid$catt=="Resistant")
median(x$median)
x=subset(centroid, centroid$catt=="Sensitive")
median(x$median)

ggplot(centroid, aes(x = catt, y = median, fill=catt))+geom_violin(trim=FALSE)+ scale_fill_manual(name= "Sample Response",values=c("grey40", "grey"))+xlab("") + theme_classic()+  ylab("VANDETANIB_OVERLAP_DOWN Score") + geom_signif(comparisons = list(c("Sensitive", "Resistant")), map_signif_level=TRUE, size=0.5,textsize=5,step_increase = 0.05, tip_length = 0.01, test="wilcox.test" , y_position = 8)+ theme(text=element_text(size=15),axis.text.x = element_text(size=12, face="bold"),axis.text.y = element_text(size=12),axis.title = element_text(size=12,face="bold"),legend.title = element_text(face="bold", size=14),legend.text = element_text(size=14),legend.position = "none")+scale_x_discrete(labels=c('Resistant', 'Sensitive'))+stat_summary(fun="mean", geom="point", size=2, color="black")+ geom_boxplot(width=0.1)
bact <- compare_means(median~catt, data = centroid,method = 'wilcox.test')
bact