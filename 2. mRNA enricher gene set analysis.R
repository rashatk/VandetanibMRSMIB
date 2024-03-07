#Enricher analysis
setwd("C:/Users/rasha/Downloads/FinalMRNASEQcodes-20231028T190637Z-001/FinalMRNASEQcodes")

#Load reference datasets
library(clusterProfiler)
GO_file = "MERGED_GSEA.breast_PEROU_SPANHEIMER.RET.VAND.gmt"
data=read.gmt(GO_file)


#Select significantly upregulated and downregulated genes
df <- read.csv("MCF7_DES_03312023.csv", header=TRUE, row.names=1)
df.top=df[((df$padj < 0.05)),]
upz <- rownames(subset(df.top, df.top$log2FoldChange>0.5))
downz <- rownames(subset(df.top, df.top$log2FoldChange< -0.5))

#Enrichment analysis
genes <- upz
egmt <- enricher(genes, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod = "BH", TERM2GENE=data)
file=data.frame(egmt)
rownames(file)=NULL
file2 <- file[-c(2)]
file2$GR <- sub("\\/.*","\\/", file2$BgRatio, perl=TRUE)
file2$GR = stringr::str_replace(file2$GR, "/" , "")
file2$GR = as.numeric(file2$Count)/as.numeric(file2$GR)
up = file2

genes <- downz
egmt <- enricher(genes, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod = "BH", TERM2GENE=data)
file=data.frame(egmt)
rownames(file)=NULL
file2 <- file[-c(2)]
file2$GR <- sub("\\/.*","\\/", file2$BgRatio, perl=TRUE)
file2$GR = stringr::str_replace(file2$GR, "/" , "")
file2$GR = as.numeric(file2$Count)/as.numeric(file2$GR)
down = file2

##select gene sets
setz=down #run the below analysis for "down" then replace with "up" and re-run
setz=subset(setz,setz$p.adjust<0.05)
setz1=rbind(setz[grepl("E2|ESTROGEN|ESTRADIOL|Estrogen|estrogen|estradiol|ESR1|ESR2|_ESR_|ESRRA|ANDROGEN", setz$ID),]) 
setz1=setz1[order(setz1$GR),]
setz1=unique(setz1)
setz2=rbind(setz[grepl("TAMOXIFEN|RESISTANCE|SERM|FULVESTRANT|ENDOCRINE|Tamoxifen|tamoxifen|resist", setz$ID),])
setz2=setz2[order(setz2$GR),]
setz2=unique(setz2)
setz3=rbind(setz[grepl("MAPK|ERK1|ERK2|MAP2K|MEK", setz$ID),])#MAPK pathways
setz3=setz3[order(setz3$GR),]
setz3=unique(setz3)
setz4=rbind(setz[grepl("PI3K|AKT|MTOR|akt|pi3k", setz$ID),])#PI3K AKT pathways
setz4=setz4[order(setz4$GR),]
setz4=unique(setz4)
setz5=rbind(setz[grepl("KRAS|ras|Ras", setz$ID),])#KRAS pathways
setz5=setz5[order(setz5$GR),]
setz5=unique(setz5)
setz6=rbind(setz[grepl("EGFR|egfr", setz$ID),])#EGFR pathways
setz6=setz6[order(setz6$GR),]
setz6=unique(setz6)
setz7=rbind(setz[grepl("VEGF", setz$ID),])#VEGF pathways
setz7=setz7[order(setz7$GR),]
setz7=unique(setz7)
setz8=rbind(setz[grepl("RET", setz$ID),])#RET pathways
setz8=setz8[order(setz8$GR),]
setz8=unique(setz8)
setz9=rbind(setz[grepl("ERBB|Erbb", setz$ID),]) #ERRB2 pathways
setz9=setz9[order(setz9$GR),]
setz9=unique(setz9)
setz10=rbind(setz[grepl("MYC|myc", setz$ID),])#MYC pathways
setz10=setz10[order(setz10$GR),]
setz10=unique(setz10)
setz11=rbind(setz[grepl("P53|APOPTOSIS|CYCLE|G2M|PHASE|DEATH|CYCLIN|CCND|RB|cyclin|Cyclin|p53|proliferat|Proliferat|HDAC", setz$ID),]) # p53 AND CYCLE PATHWAYS
setz11=setz11[order(setz11$GR),]
setz11=unique(setz11)
setz12=setz[grepl("FGF", setz$ID),]
setz12=unique(setz12)
setz13=setz[grepl("LUMINAL|BASAL|NORMAL|SUBTYPE|INVASIVE|luminal|Luminal|EMT|Claudin", setz$ID),] ##phenotype
setz13=setz13[order(setz13$GR),]
setz13=unique(setz13)
setz14=setz[grepl("JAK|STAT|JNK|JUN|FOS", setz$ID),]
setz14=unique(setz14)
setz15=setz[grepl("TNF|NFKB|IKBKB", setz$ID),]
setz15=unique(setz15)
setz16=setz[grepl("WNT|CATENIN", setz$ID),]
setz16=unique(setz16)
setz1$category="Estrogen & ESR1"
setz2$category="AE"
setz3$category="MAPK"
setz4$category="AKT"
setz5$category="KRAS"
setz6$category="EGFR"
setz7$category="VEGF"
setz8$category="RET"
setz9$category="ERBB2"
setz10$category="MYC"
setz11$category="Proliferation"
setz12$category="FGFR"
setz13$category="Phenotype"
setz14$category="JAK/STAT"
setz15$category="NFKB"
setz16$category="WNT"
dffz_down=rbind(setz1,setz2,setz3,setz4,setz5,setz6,setz7,setz9,setz10,setz11,setz12,setz13,setz14,setz15,setz16) #run after "down" as setz
dffz_down$group="down"

setz=up
setz=subset(setz,setz$p.adjust<0.05)
setz1=rbind(setz[grepl("E2|ESTROGEN|ESTRADIOL|Estrogen|estrogen|estradiol|ESR1|ESR2|_ESR_|ESRRA|ANDROGEN", setz$ID),]) 
setz1=setz1[order(setz1$GR),]
setz1=unique(setz1)
setz2=rbind(setz[grepl("TAMOXIFEN|RESISTANCE|SERM|FULVESTRANT|ENDOCRINE|Tamoxifen|tamoxifen|resist", setz$ID),])
setz2=setz2[order(setz2$GR),]
setz2=unique(setz2)
setz3=rbind(setz[grepl("MAPK|ERK1|ERK2|MAP2K|MEK", setz$ID),])#MAPK pathways
setz3=setz3[order(setz3$GR),]
setz3=unique(setz3)
setz4=rbind(setz[grepl("PI3K|AKT|MTOR|akt|pi3k", setz$ID),])#PI3K AKT pathways
setz4=setz4[order(setz4$GR),]
setz4=unique(setz4)
setz5=rbind(setz[grepl("KRAS|ras|Ras", setz$ID),])#KRAS pathways
setz5=setz5[order(setz5$GR),]
setz5=unique(setz5)
setz6=rbind(setz[grepl("EGFR|egfr", setz$ID),])#EGFR pathways
setz6=setz6[order(setz6$GR),]
setz6=unique(setz6)
setz7=rbind(setz[grepl("VEGF", setz$ID),])#VEGF pathways
setz7=setz7[order(setz7$GR),]
setz7=unique(setz7)
setz8=rbind(setz[grepl("RET", setz$ID),])#RET pathways
setz8=setz8[order(setz8$GR),]
setz8=unique(setz8)
setz9=rbind(setz[grepl("ERBB|Erbb", setz$ID),]) #ERRB2 pathways
setz9=setz9[order(setz9$GR),]
setz9=unique(setz9)
setz10=rbind(setz[grepl("MYC|myc", setz$ID),])#MYC pathways
setz10=setz10[order(setz10$GR),]
setz10=unique(setz10)
setz11=rbind(setz[grepl("P53|APOPTOSIS|CYCLE|G2M|PHASE|DEATH|CYCLIN|CCND|cyclin
                        |Cyclin|p53|proliferat|Proliferat|HDAC|RB", setz$ID),]) # p53 AND CYCLE PATHWAYS
setz11=setz11[order(setz11$GR),]
setz11=unique(setz11)
setz12=setz[grepl("FGF", setz$ID),]
setz12=unique(setz12)
setz13=setz[grepl("LUMINAL|BASAL|NORMAL|SUBTYPE|INVASIVE|
                  luminal|Luminal|EMT|Claudin", setz$ID),] ##phenotype
setz13=setz13[order(setz13$GR),]
setz13=unique(setz13)
setz14=setz[grepl("JAK|STAT|JNK|JUN|FOS", setz$ID),]
setz14=unique(setz14)
setz15=setz[grepl("TNF|NFKB|IKBKB", setz$ID),]
setz15=unique(setz15)
setz16=setz[grepl("WNT|CATENIN", setz$ID),]
setz16=unique(setz16)
setz1$category="Estrogen & ESR1"
setz2$category="AE"
setz3$category="MAPK"
setz4$category="AKT"
setz5$category="KRAS"
setz6$category="EGFR"
setz7$category="VEGF"
setz8$category="RET"
setz9$category="ERBB2"
setz10$category="MYC"
setz11$category="Proliferation"
setz12$category="FGFR"
setz13$category="Phenotype"
setz14$category="JAK/STAT"
setz15$category="NFKB"
setz16$category="WNT"
dffz_up=rbind(setz1,setz2,setz3,setz4,setz5,setz6,setz7,setz9,setz10,setz11,setz12,setz13,setz14,setz15,setz16) #run after "up" as setz
dffz_up$group="up"

#modify directionality of plotting by allocating negative GR for downregulated sets and positive GR for upregulated sets
newdown=dffz_down
newdown$GRmod=-dffz_down$GR
newup=dffz_up
newup$GRmod=dffz_up$GR
dffz=rbind(newdown,newup)
write.csv(dffz,"EnricherOutput_MCF7.csv")

#assign factors for grouped plotting
dffz$Enr_f = factor(dffz$category, levels=c('Estrogen & ESR1','AE','MAPK','AKT','EGFR','VEGF','ERBB2','KRAS','MYC','Proliferation','Phenotype','JNK', 'NFKB','WNT', 'FGFR', 'JAK/STAT'))
dffz$group=ifelse(dffz$group=="up","UPREGULATED","DOWNREGULATED")
dffz$Enr_group=factor(dffz$group, levels=c('DOWNREGULATED','UPREGULATED'))

#plot_everything
library(ggplot2)
ggplot(data = dffz, aes(x = GR, y = ID,fill=p.adjust)) + geom_bar(stat="identity")+facet_grid(Enr_f~., scales = "free", space = "free", switch="y") + theme(axis.title.y=element_blank(),legend.position = "bottom",text=element_text(size=13))+xlab("Gene Ratio by Set")

#filter out any sets that pop-up unrelated to breast cancer
dffz=dffz[!grepl("SYNAPTIC|ARSENIC|GOBP_ENDOCRINE|GOBP_REGULATION_OF_ENDOCRINE|NEUROENDOCRINE|NFE2L2|PROSTATE|LUNG|KIDNEY|E2F3|MEISSNER|HEPATOBLASTOMA|ESOPHAGEAL|NEUROTRANSMITTER|INFANT|GLIOBLASTOMA|PLATEAU|BRUINS|GTPASE|DOXORUBICIN|MARTINEZ|INTESTINE|LTE2|PANCREAS|FETAL|DEPRIVATION|AGING|KRAS.DF|DOCETACEL|E2F1|HP_|Lung|ANATOMICAL|NCX|METHYLATED|YOSHIMURA|MTOR_UP.|CARBOXYLIC|LINKAGE|AKT_UP|VEGF|P53_DN|HAMAI|KRAS.BREAST|PEREZ|JNK_DN|JAK2_DN|UP.N4.|YANG|VANTVEER",dffz$ID),]

ggplot(data = dffz, aes(x = GRmod, y = ID,fill=p.adjust)) + geom_bar(stat="identity")+facet_grid(Enr_f~., scales = "free", space = "free", switch="y") + theme(axis.title.y=element_blank(),legend.position = "bottom",text=element_text(size=13))+xlab("Gene Ratio by Set")

#plot subsets to discern gene sets of interest
ggplot(data = dffz[grepl("ESR", dffz$category),], aes(x = GRmod, y = ID,fill=p.adjust)) + geom_bar(stat="identity")+facet_grid(Enr_f~., scales = "free", space = "free", switch="y") + theme(axis.title.y=element_blank(),legend.position = "bottom",text=element_text(size=13))+xlab("Gene Ratio by Set")


###-------without grouping, order by gene ratio, omitted certain gene sets and designated directionality of sets to delinate duplicate sets in one direction per category
dffz=dffz[grepl("BOWIE_RESPONSE_TO_TAMOXIFEN|STEIN_ESTROGEN_RESPONSE_NOT_VIA_ESRRA|SMID_BREAST_CANCER_ERBB2_UP|FRASOR_TAMOXIFEN_RESPONSE_UP|MASSARWEH_RESPONSE_TO_ESTRADIOL|BHAT_ESR1_TARGETS_NOT_VIA_AKT1_DN|MCF7.E2.repressed.genes|CREIGHTON_ENDOCRINE_THERAPY_RESISTANCE_5|Taube_EMT_down|WP_MAPK_SIGNALING_PATHWAY|KEGG_MAPK_SIGNALING_PATHWAY|EGFR_UP.V1_UP|WP_EGFR_SIGNALING_PATHWAY|Duke_Module16_pi3k|HALLMARK_MYC_TARGETS_UP|DANG_MYC_TARGETS_UP|Duke_Module18_ras|LAMB_CCND1_TARGETS|Wirapati_Proliferation|FRASOR_RESPONSE_TO_SERM_OR_FULVESTRANT_DN|Proliferation_Cluster",dffz$ID),]
up=subset(dffz, dffz$group=="UPREGULATED")
up$ID <- paste0(up$ID, "__up")
down=subset(dffz, dffz$group=="DOWNREGULATED")
down$ID <- paste0(down$ID, "__down")
dffz=rbind(up,down)
library(ggpubr)
ggplot(data = dffz[order(dffz$GRmod),], aes(x = GRmod, y = reorder(ID,GRmod),fill=p.adjust))+geom_bar(stat="identity")+xlab("dGR")+theme_pubclean()+theme(axis.title.y=element_blank(),legend.position = "right",text=element_text(size=13,face="bold"))+scale_fill_continuous(name = "FDR")
#modify names/design using adobe illustrator

###----new plot layout
dffz$ID = gsub("__up", "", dffz$ID)
dffz$ID = gsub("__down", "", dffz$ID)

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(extrafont)
library(patchwork)
dffz$ID
dffz$ID = gsub("Duke_Module16_pi3k_Mike_PMID.20335537", "Duke_Module16_pi3k_Mike_PMID", dffz$ID)
dffz$ID = gsub("Duke_Module18_ras_Mike_PMID.20335537", "Duke_Module18_ras_Mike_PMID", dffz$ID)
dffz$ID = gsub("Proliferation_Cluster_BMC.Med.Genomics.2011_PMID.21214954", "Proliferation_Cluster", dffz$ID)
dffz$ID = gsub("Wirapati_Proliferation_Breast.Cancer.Res.2008_PMID.18662380", "Wirapati_Proliferation", dffz$ID)
dffz$ID = gsub("Taube_EMT_down_PNAS.2010_PMID.20713713", "Taube_EMT_down", dffz$ID)

a=ggplot(data = dffz[grepl("UP",dffz$group),], aes(x = GRmod, y = reorder(ID,GRmod),fill=p.adjust))+ggtitle("MCF7")+geom_bar(stat="identity")+xlab("dGR")+theme_pubclean()+theme(axis.title.y=element_blank(),legend.position = "right",text=element_text(size=10),plot.title = element_text(face="bold"),axis.text.y = element_text(size=13),legend.title=element_text(size=10),legend.text = element_text(size=8),legend.key.height = unit(0.5,'cm'),legend.key.width=unit(0.5,'cm'))+scale_fill_gradient(low="#ca7272",high="#cf3a3a",name="upFDR")

b=ggplot(data = dffz[grepl("DOWN",dffz$group),], aes(x = GRmod, y = reorder(ID,GRmod),fill=p.adjust))+geom_bar(stat="identity")+xlab("dGR")+theme_pubclean()+theme(axis.title.y=element_blank(),legend.position = "right",text=element_text(size=12),axis.text.y = element_text(size=13),legend.title=element_text(size=10),legend.text = element_text(size=8),legend.key.height = unit(0.5,'cm'),legend.key.width=unit(0.5,'cm'))+scale_fill_gradient(low="#7c97fb",high="#2a57fc",name="downFDR")
layout <- "
A
B
B
"
a + b + plot_layout(design=layout)
plot=a + b + plot_layout(design=layout)
ggsave("MCF7_newenricher.pdf",width=8.5,height=8,dpi=300)


##---MCF7TAM
##Select significantly upregulated and downregulated genes
df <- read.csv("MCF7TAM_DES_03312023.csv", header=TRUE, row.names=1)
df.top=df[((df$padj < 0.05)),]
upz <- rownames(subset(df.top, df.top$log2FoldChange>0.5))
downz <- rownames(subset(df.top, df.top$log2FoldChange< -0.5))

#Enrichment analysis
genes <- upz
egmt <- enricher(genes, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod = "BH", TERM2GENE=data)
file=data.frame(egmt)
rownames(file)=NULL
file2 <- file[-c(2)]
file2$GR <- sub("\\/.*","\\/", file2$BgRatio, perl=TRUE)
file2$GR = stringr::str_replace(file2$GR, "/" , "")
file2$GR = as.numeric(file2$Count)/as.numeric(file2$GR)
up = file2

genes <- downz
egmt <- enricher(genes, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod = "BH", TERM2GENE=data)
file=data.frame(egmt)
rownames(file)=NULL
file2 <- file[-c(2)]
file2$GR <- sub("\\/.*","\\/", file2$BgRatio, perl=TRUE)
file2$GR = stringr::str_replace(file2$GR, "/" , "")
file2$GR = as.numeric(file2$Count)/as.numeric(file2$GR)
down = file2

##select gene sets
setz=down #run the below analysis for "down" then replace with "up" and re-run
setz=subset(setz,setz$p.adjust<0.05)
setz1=rbind(setz[grepl("E2|ESTROGEN|ESTRADIOL|Estrogen|estrogen|estradiol|ESR1|ESR2|_ESR_|ESRRA|ANDROGEN", setz$ID),]) 
setz1=setz1[order(setz1$GR),]
setz1=unique(setz1)
setz2=rbind(setz[grepl("TAMOXIFEN|RESISTANCE|SERM|FULVESTRANT|ENDOCRINE|Tamoxifen|tamoxifen|resist", setz$ID),])
setz2=setz2[order(setz2$GR),]
setz2=unique(setz2)
setz3=rbind(setz[grepl("MAPK|ERK1|ERK2|MAP2K|MEK", setz$ID),])#MAPK pathways
setz3=setz3[order(setz3$GR),]
setz3=unique(setz3)
setz4=rbind(setz[grepl("PI3K|AKT|MTOR|akt|pi3k", setz$ID),])#PI3K AKT pathways
setz4=setz4[order(setz4$GR),]
setz4=unique(setz4)
setz5=rbind(setz[grepl("KRAS|ras|Ras", setz$ID),])#KRAS pathways
setz5=setz5[order(setz5$GR),]
setz5=unique(setz5)
setz6=rbind(setz[grepl("EGFR|egfr", setz$ID),])#EGFR pathways
setz6=setz6[order(setz6$GR),]
setz6=unique(setz6)
setz7=rbind(setz[grepl("VEGF", setz$ID),])#VEGF pathways
setz7=setz7[order(setz7$GR),]
setz7=unique(setz7)
setz8=rbind(setz[grepl("RET", setz$ID),])#RET pathways
setz8=setz8[order(setz8$GR),]
setz8=unique(setz8)
setz9=rbind(setz[grepl("ERBB|Erbb", setz$ID),]) #ERRB2 pathways
setz9=setz9[order(setz9$GR),]
setz9=unique(setz9)
setz10=rbind(setz[grepl("MYC|myc", setz$ID),])#MYC pathways
setz10=setz10[order(setz10$GR),]
setz10=unique(setz10)
setz11=rbind(setz[grepl("P53|APOPTOSIS|CYCLE|G2M|PHASE|DEATH|CYCLIN|CCND|RB|cyclin|Cyclin|p53|proliferat|Proliferat|HDAC", setz$ID),]) # p53 AND CYCLE PATHWAYS
setz11=setz11[order(setz11$GR),]
setz11=unique(setz11)
setz12=setz[grepl("FGF", setz$ID),]
setz12=unique(setz12)
setz13=setz[grepl("LUMINAL|BASAL|NORMAL|SUBTYPE|INVASIVE|luminal|Luminal|EMT|Claudin", setz$ID),] ##phenotype
setz13=setz13[order(setz13$GR),]
setz13=unique(setz13)
setz14=setz[grepl("JAK|STAT|JNK|JUN|FOS", setz$ID),]
setz14=unique(setz14)
setz15=setz[grepl("TNF|NFKB|IKBKB", setz$ID),]
setz15=unique(setz15)
setz16=setz[grepl("WNT|CATENIN", setz$ID),]
setz16=unique(setz16)
setz1$category="Estrogen & ESR1"
setz2$category="AE"
setz3$category="MAPK"
setz4$category="AKT"
setz5$category="KRAS"
setz6$category="EGFR"
setz7$category="VEGF"
setz8$category="RET"
setz9$category="ERBB2"
setz10$category="MYC"
setz11$category="Proliferation"
setz12$category="FGFR"
setz13$category="Phenotype"
setz14$category="JAK/STAT"
setz15$category="NFKB"
setz16$category="WNT"
dffz_down=rbind(setz1,setz2,setz3,setz4,setz5,setz6,setz7,setz9,setz10,setz11,setz12,setz13,setz14,setz15,setz16) #run after "down" as setz
dffz_down$group="down"

setz=up
setz=subset(setz,setz$p.adjust<0.05)
setz1=rbind(setz[grepl("E2|ESTROGEN|ESTRADIOL|Estrogen|estrogen|estradiol|ESR1|ESR2|_ESR_|ESRRA|ANDROGEN", setz$ID),]) 
setz1=setz1[order(setz1$GR),]
setz1=unique(setz1)
setz2=rbind(setz[grepl("TAMOXIFEN|RESISTANCE|SERM|FULVESTRANT|ENDOCRINE|Tamoxifen|tamoxifen|resist", setz$ID),])
setz2=setz2[order(setz2$GR),]
setz2=unique(setz2)
setz3=rbind(setz[grepl("MAPK|ERK1|ERK2|MAP2K|MEK", setz$ID),])#MAPK pathways
setz3=setz3[order(setz3$GR),]
setz3=unique(setz3)
setz4=rbind(setz[grepl("PI3K|AKT|MTOR|akt|pi3k", setz$ID),])#PI3K AKT pathways
setz4=setz4[order(setz4$GR),]
setz4=unique(setz4)
setz5=rbind(setz[grepl("KRAS|ras|Ras", setz$ID),])#KRAS pathways
setz5=setz5[order(setz5$GR),]
setz5=unique(setz5)
setz6=rbind(setz[grepl("EGFR|egfr", setz$ID),])#EGFR pathways
setz6=setz6[order(setz6$GR),]
setz6=unique(setz6)
setz7=rbind(setz[grepl("VEGF", setz$ID),])#VEGF pathways
setz7=setz7[order(setz7$GR),]
setz7=unique(setz7)
setz8=rbind(setz[grepl("RET", setz$ID),])#RET pathways
setz8=setz8[order(setz8$GR),]
setz8=unique(setz8)
setz9=rbind(setz[grepl("ERBB|Erbb", setz$ID),]) #ERRB2 pathways
setz9=setz9[order(setz9$GR),]
setz9=unique(setz9)
setz10=rbind(setz[grepl("MYC|myc", setz$ID),])#MYC pathways
setz10=setz10[order(setz10$GR),]
setz10=unique(setz10)
setz11=rbind(setz[grepl("P53|APOPTOSIS|CYCLE|G2M|PHASE|DEATH|CYCLIN|CCND|cyclin
                        |Cyclin|p53|proliferat|Proliferat|HDAC|RB", setz$ID),]) # p53 AND CYCLE PATHWAYS
setz11=setz11[order(setz11$GR),]
setz11=unique(setz11)
setz12=setz[grepl("FGF", setz$ID),]
setz12=unique(setz12)
setz13=setz[grepl("LUMINAL|BASAL|NORMAL|SUBTYPE|INVASIVE|
                  luminal|Luminal|EMT|Claudin", setz$ID),] ##phenotype
setz13=setz13[order(setz13$GR),]
setz13=unique(setz13)
setz14=setz[grepl("JAK|STAT|JNK|JUN|FOS", setz$ID),]
setz14=unique(setz14)
setz15=setz[grepl("TNF|NFKB|IKBKB", setz$ID),]
setz15=unique(setz15)
setz16=setz[grepl("WNT|CATENIN", setz$ID),]
setz16=unique(setz16)
setz1$category="Estrogen & ESR1"
setz2$category="AE"
setz3$category="MAPK"
setz4$category="AKT"
setz5$category="KRAS"
setz6$category="EGFR"
setz7$category="VEGF"
setz8$category="RET"
setz9$category="ERBB2"
setz10$category="MYC"
setz11$category="Proliferation"
setz12$category="FGFR"
setz13$category="Phenotype"
setz14$category="JAK/STAT"
setz15$category="NFKB"
setz16$category="WNT"
dffz_up=rbind(setz1,setz2,setz3,setz4,setz5,setz6,setz7,setz9,setz10,setz11,setz12,setz13,setz14,setz15,setz16) #run after "up" as setz
dffz_up$group="up"

#modify directionality of plotting by allocating negative GR for downregulated sets and positive GR for upregulated sets
newdown=dffz_down
newdown$GRmod=-dffz_down$GR
newup=dffz_up
newup$GRmod=dffz_up$GR
dffz=rbind(newdown,newup)
write.csv(dffz,"EnricherOutput_MCF7TAM.csv")

#assign factors for grouped plotting
dffz$Enr_f = factor(dffz$category, levels=c('Estrogen & ESR1','AE','MAPK','AKT','EGFR','VEGF','ERBB2','KRAS','MYC','Proliferation','Phenotype','JNK', 'NFKB','WNT', 'FGFR', 'JAK/STAT'))
dffz$group=ifelse(dffz$group=="up","UPREGULATED","DOWNREGULATED")
dffz$Enr_group=factor(dffz$group, levels=c('DOWNREGULATED','UPREGULATED'))

#plot_everything
library(ggplot2)
ggplot(data = dffz, aes(x = GR, y = ID,fill=p.adjust)) + geom_bar(stat="identity")+facet_grid(Enr_f~., scales = "free", space = "free", switch="y") + theme(axis.title.y=element_blank(),legend.position = "bottom",text=element_text(size=13))+xlab("Gene Ratio by Set")

#filter out any sets that pop-up unrelated to breast cancer
dffz=dffz[!grepl("SYNAPTIC|ARSENIC|GOBP_ENDOCRINE|GOBP_REGULATION_OF_ENDOCRINE|NEUROENDOCRINE|NFE2L2|PROSTATE|LUNG|KIDNEY|E2F3|MEISSNER|HEPATOBLASTOMA|ESOPHAGEAL|NEUROTRANSMITTER|INFANT|GLIOBLASTOMA|PLATEAU|BRUINS|GTPASE|DOXORUBICIN|MARTINEZ|INTESTINE|LTE2|PANCREAS|FETAL|DEPRIVATION|AGING|KRAS.DF|DOCETACEL|E2F1|HP_|Lung|ANATOMICAL|NCX|METHYLATED|YOSHIMURA|MTOR_UP.|CARBOXYLIC|LINKAGE|AKT_UP|VEGF|P53_DN|HAMAI|KRAS.BREAST|PEREZ|JNK_DN|JAK2_DN|UP.N4.|YANG|VANTVEER",dffz$ID),]

ggplot(data = dffz, aes(x = GRmod, y = ID,fill=p.adjust)) + geom_bar(stat="identity")+facet_grid(Enr_f~., scales = "free", space = "free", switch="y") + theme(axis.title.y=element_blank(),legend.position = "bottom",text=element_text(size=13))+xlab("Gene Ratio by Set")

#plot subsets to discern gene sets of interest
ggplot(data = dffz[grepl("ESR", dffz$category),], aes(x = GRmod, y = ID,fill=p.adjust)) + geom_bar(stat="identity")+facet_grid(Enr_f~., scales = "free", space = "free", switch="y") + theme(axis.title.y=element_blank(),legend.position = "bottom",text=element_text(size=13))+xlab("Gene Ratio by Set")


###-------without grouping, order by gene ratio, omitted certain gene sets and designated directionality of sets to delinate duplicate sets in one direction per category
dffz=dffz[grepl("FRASOR_TAMOXIFEN_RESPONSE_UP|FRASOR_RESPONSE_TO_ESTRADIOL_UP|REACTOME_NUCLEAR_SIGNALING_BY_ERBB4|MASSARWEH_rESPONSE_TO_ESTRADIOL|BHAT_ESR1_TARGETS_VIA_AKT1_UP|STOSSI_RESPONSE_TO_ESTRADIOL|MM_Erbb2.like|CREIGHTON_ENDOCRINE_THERAPY_RESISTANCE_5|HALLMARK_KRAS_SIGNALING_UP|Duke_Module16_pi3k|CYCLIN_D1_KE_.V1_UP|PID_MYC_ACTIV_PATHWAY|Duke_Module05_egfr|PID_MYC_REPRESS_PATHWAY|POOLA_INVASIVE_BREAST_CANCER_UP|SWEET_KRAS_TARGETS_DN|Wirapati_Proliferation|STEIN_ESRRA_TARGETS_RESPONSIVE_TO_ESTROGEN_DN|FRASOR_RESPONSE_TO_SERM_OR_FULVESTRANT_DN|Proliferation_Cluster|Prosigna_Proliferation|MAPK_activity_Score",dffz$ID),]


up=subset(dffz, dffz$group=="UPREGULATED")
up$ID <- paste0(up$ID, "__up")
down=subset(dffz, dffz$group=="DOWNREGULATED")
down$ID <- paste0(down$ID, "__down")
dffz=rbind(up,down)
library(ggpubr)
ggplot(data = dffz[order(dffz$GRmod),], aes(x = GRmod, y = reorder(ID,GRmod),fill=p.adjust))+geom_bar(stat="identity")+xlab("dGR")+theme_pubclean()+theme(axis.title.y=element_blank(),legend.position = "right",text=element_text(size=13,face="bold"))+scale_fill_continuous(name = "FDR")
#modify names/design using adobe illustrator

###----new plot layout
dffz$ID = gsub("__up", "", dffz$ID)
dffz$ID = gsub("__down", "", dffz$ID)

library(ggplot2)
library(ggpubr)
library(patchwork)
dffz$ID
dffz$ID = gsub("Proliferation_Cluster_BMC.Med.Genomics.2011_PMID.21214954", "Proliferation_Cluster", dffz$ID)
dffz$ID = gsub("Wirapati_Proliferation_Breast.Cancer.Res.2008_PMID.18662380", "Wirapati_Proliferation", dffz$ID)
dffz$ID = gsub("Prosigna_Proliferation_18_BMC.Med.Genomics.2015_PMID.26297356", "Prosigna_Proliferation", dffz$ID)
dffz$ID = gsub("MAPK_activity_Score_MPAS_PrecisionOncology.2018_PMID.TBD", "MAPK_activity_Score", dffz$ID)
dffz$ID = gsub("MM_Erbb2.like_Genome.Biol.2013_PMID.24220145", "MM_Erbb2.like", dffz$ID)
dffz$ID = gsub("Duke_Module05_egfr_Mike_PMID.20335537", "Duke_Module05_egfr_Mike", dffz$ID)
dffz$ID = gsub("Duke_Module16_pi3k_Mike_PMID.20335537", "Duke_Module16_pi3k_Mike", dffz$ID)

a=ggplot(data = dffz[grepl("UP",dffz$group),], aes(x = GRmod, y = reorder(ID,GRmod),fill=p.adjust))+ggtitle("MCF7TAM")+geom_bar(stat="identity")+xlab("dGR")+theme_pubclean()+theme(axis.title.y=element_blank(),legend.position = "right",text=element_text(size=10),plot.title = element_text(face="bold"),axis.text.y = element_text(size=13),legend.title=element_text(size=10),legend.text = element_text(size=8),legend.key.height = unit(0.5,'cm'),legend.key.width=unit(0.5,'cm'))+scale_fill_gradient(low="#ca7272",high="#cf3a3a",name="upFDR")

b=ggplot(data = dffz[grepl("DOWN",dffz$group),], aes(x = GRmod, y = reorder(ID,GRmod),fill=p.adjust))+geom_bar(stat="identity")+xlab("dGR")+theme_pubclean()+theme(axis.title.y=element_blank(),legend.position = "right",text=element_text(size=12),axis.text.y = element_text(size=13),legend.title=element_text(size=10),legend.text = element_text(size=8),legend.key.height = unit(0.5,'cm'),legend.key.width=unit(0.5,'cm'))+scale_fill_gradient(low="#7c97fb",high="#2a57fc",name="downFDR")
layout <- "
A
B
B
B
"
a + b + plot_layout(design=layout)
plot=a + b + plot_layout(design=layout)
ggsave("MCF7TAM_newenricher.pdf",width=9,height=8,dpi=300)

MCF7TAMup=dffz[grepl("UP",dffz$group),]
MCF7TAMup$line="MCF7TAM"
MCF7TAMdown=dffz[grepl("DOWN",dffz$group),]
MCF7TAMdown$line="MCF7TAM"
MCF7TAMsets=rbind(MCF7TAMup,MCF7TAMdown)

