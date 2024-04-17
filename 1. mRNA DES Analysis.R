library(ComplexHeatmap)
library(fabricatr)
library(dplyr)
library(DESeq2)
library(viridis)
library(edgeR)
library(EnsDb.Hsapiens.v79)
library(NOISeq)
library(EnhancedVolcano)

setwd("C:/Users/rasha/Downloads/FinalMRNASEQcodes-20231028T190637Z-001/FinalMRNASEQcodes")
input=read.csv("rsem.merged.gene_counts.csv",header=TRUE,row.names=1)
input=input[-c(1)]
names(input)

#prepare design matrix
samples=read.csv("sampledesign.csv",header=TRUE)
mod=model.matrix(~0+treatment+cell.line, data=samples)

#divide tam and vand
new=input
names(new)
control_cols <- new[c(19:21,10:12,4:6,34:36,25:27)]
tam_cols <- new[c(13:15,28:30)]
vand_cols <- new[c(16:18,7:9,1:3,31:33,22:24)]

dfCV <- cbind(control_cols, vand_cols)
dfCT <- cbind(control_cols, tam_cols)
dfTV <- cbind(tam_cols, vand_cols)

####ANALYSIS of either combination set
##vandetanib divisions
CV_MCF7 <- dfCV[c(1:3,16:18)]
CV_MCF7EXE <-dfCV[c(7:9,22:24)]
CV_T47DTAM <-dfCV[c(13:15,28:30)]

#after correlation analysis, two samples were of low quality/outliers and were omitted 
CV_MCF7TAM <-dfCV[c(4:6,20:21)]
CV_T47D <- dfCV[c(11:12,25:27)]

#DES analysis
###for each line, specify set and line details

####MCF7
set <- CV_MCF7
groups <- data.frame(samples=colnames(CV_MCF7),category = c(rep("control", 3), rep("treated",3)))
coldata=groups

#data filtering and normalization
df=set[rowSums(set) > 0,]
library(edgeR)
keep <- filterByExpr(df, design = mod)
df <- df[keep,]
library(NOISeq)
data=uqua(df, long = 1000, lc = 0, k = 0)
data=as.data.frame(data)

#run DES
library(DESeq2)
roundedcounts=round(data)
groups$category=factor(groups$category)
dds <- DESeqDataSetFromMatrix(countData = roundedcounts,
                              colData = groups,
                              design = ~ category)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$category <- relevel(dds$category, ref = "control")
dds <- DESeq(dds)
res <- results(dds)
DESeq_output <- data.frame(res)
DES_MCF7=na.omit(DESeq_output)
df=DES_MCF7
df.top <- df[(df$padj < 0.05) & (abs(df$log2FoldChange) > 0.5),]

#gene symbols for top genes and rearrange
library("EnsDb.Hsapiens.v79")
ensembl.genes <- rownames(df.top)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
df.top=subset(df.top, rownames(df.top) != setdiff(rownames(df.top), geneIDs1$GENEID))
df.top$SYMBOL=geneIDs1$SYMBOL
df.top$ENSEMBL=rownames(df.top)

#matrix for heatmap plotting
rlog_out <- vst(dds, blind=FALSE) #get normalized count data from dds object
mat<-assay(rlog_out)[rownames(df.top), coldata$samples] #sig genes x samples
colnames(mat) <- coldata$samples
row.names(mat)=df.top$SYMBOL
row.names(df.top)=df.top$SYMBOL
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mat)
num_keep <- 20
Genez <- rbind(tail(df.top[order(df.top$log2FoldChange),], num_keep), (head(df.top[order(df.top$log2FoldChange),], num_keep)))
rows_keep <- rownames(Genez)
matrixx <- mat.scaled[rows_keep,]

#select top significantly regulated genes
MCF7_sigup=rownames(df.top[(df.top$padj < 0.05) & (df.top$log2FoldChange > 0.5),])
MCF7_sigdown=rownames(df.top[(df.top$padj < 0.05) & (df.top$log2FoldChange < -0.5),])
write.csv(df.top, "MCF7_DES_03312023.csv")
write.csv(df,"MCF7_alloutputDES.csv")

#allsignificants
MCF7_allsigup=rownames(df[(df$padj < 0.05) & (df$log2FoldChange > 0),])
MCF7_allsigdown=rownames(df[(df$padj < 0.05) & (df$log2FoldChange < 0),])
MCF7_increased=rownames(df[(df$padj < 0.05) & (df$log2FoldChange > 0.5),])
MCF7_decreased=rownames(df[(df$padj < 0.05) & (df$log2FoldChange < -0.5),])

#prepare tree file
tree=df[(df$padj < 0.05),]
tree=tree[c(2)]
tree$direction=ifelse(tree$log2FoldChange>0, "+1","-1")
tree$absolute=abs(tree$log2FoldChange)
tree=tree[c(2,3)]
write.csv(tree,"MCF7_sigupdown_03312023.csv")


##grouping samples and plotting heatmap
#define sample groups
grps <- groups[c(2)]
grps
grps<- ifelse(grps == "control" , "Vehicle", "Vandetanib") 
grps
colnames(grps)=NULL
rownames(grps)=NULL
grps
GRP1 <- as.matrix(grps)
ann <- data.frame(GRP1)
colnames(ann) <- c('Condition')
colours <- list('Condition' = c("Vehicle" = "#c6dbef", "Vandetanib" =  "#6baed6"))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_height = 0.5,
                            gap = unit(0.5, 'mm'),
                            show_annotation_name = FALSE)


library(ComplexHeatmap)
library(viridis)
Heatmap(matrixx, 
        cluster_rows = F,
        name="Expression",
        cluster_columns = F,
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 8),
        col=viridis(n=6,option="C"),
        column_title_gp = gpar(fontsize=7),
        column_title_side = "bottom",
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_max_width = unit(10,"in"),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        column_order = order(as.numeric(gsub("column", "", names(data)))),
        row_order = order(as.numeric(gsub("row", "", rownames(Genez)))),
        top_annotation = colAnn,
        heatmap_legend_param = list(legend_direction = "horizontal"),
        height = unit(120, "mm"),
        width = unit(20, "mm"))

###volcano plot
##volcano plot
library(EnhancedVolcano)
sc <-df

#select significant genes and genes of interest
subset_sigreg=sc[(abs(sc$log2FoldChange)>0.5) & (sc$padj < 0.05),]
ensembl.genes <- rownames(subset_sigreg)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
listo=geneIDs1$GENEID
sub=subset(subset_sigreg, rownames(subset_sigreg) %in% listo)
head(geneIDs1)
head(sub)
individual=subset(sc, !rownames(sc) %in% listo)
individual$SYMBOL="NotMapped"
sub$SYMBOL=geneIDs1$SYMBOL
merged=rbind(sub,individual)
merged$ENSEMBL=rownames(merged)

genezadd=c("RET","EGFR")
genez=subset(sub,abs(sub$log2FoldChange)>2)
genez=genez$SYMBOL
genez=unlist(strsplit(genez,","))
genez=c(genezadd,genez)
sub2=subset(merged, merged$SYMBOL %in%genez)

#assign grouping
keyvals <- ifelse(merged$log2FoldChange < -0.5 & merged$padj < 0.05, 'blue', ifelse(merged$log2FoldChange > 0.5 & merged$padj < 0.05, 'red','grey'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'grey'] <- 'Intermediate'
names(keyvals)[keyvals == 'blue'] <- 'Downregulated'
names(keyvals)[keyvals == 'red'] <- 'Upregulated'

EnhancedVolcano(merged,
                lab = merged$SYMBOL,
                x = 'log2FoldChange',
                y = 'padj',
                title = NULL,
                subtitle="",
                subtitleLabSize = 9,
                axisLabSize = 10,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3,
                labSize = 3,
                colCustom=keyvals,
                colAlpha = 0.2,
                boxedLabels=FALSE,
                selectLab = sub2$SYMBOL,
                ylab = bquote(~-Log[10] ~ FDR),
                drawConnectors = FALSE,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendPosition = 'right',
                captionLabSize = 0,
                legendLabSize = 11,
                legendIconSize = 3.0,
                max.overlaps = 20,
                cutoffLineType = 'dashed',
                cutoffLineCol = 'grey70',
                cutoffLineWidth = 0.4,
                border='full',
                borderWidth = 0.6,
                borderColour = 'black')
#9x11


####MCF7TAM
set <- CV_MCF7TAM
groups <- data.frame(samples=colnames(CV_MCF7TAM),category = c(rep("control", 3), rep("treated",2)))
coldata=groups

#data filtering and normalization
df=set[rowSums(set) > 0,]
library(edgeR)
keep <- filterByExpr(df, design = mod)
df <- df[keep,]
library(NOISeq)
data=uqua(df, long = 1000, lc = 0, k = 0)
data=as.data.frame(data)

#run DES
library(DESeq2)
roundedcounts=round(data)
groups$category=factor(groups$category)
dds <- DESeqDataSetFromMatrix(countData = roundedcounts,
                              colData = groups,
                              design = ~ category)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$category <- relevel(dds$category, ref = "control")
dds <- DESeq(dds)
res <- results(dds)
DESeq_output <- data.frame(res)
df <- na.omit(DESeq_output)
DES_MCF7TAM=na.omit(DESeq_output)
df=DES_MCF7TAM
df.top <- df[(df$padj < 0.05) & (abs(df$log2FoldChange) > 0.5),]


#gene symbols for top genes and rearrange
library("EnsDb.Hsapiens.v79")
ensembl.genes <- rownames(df.top)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
df.top=subset(df.top, rownames(df.top) != setdiff(rownames(df.top), geneIDs1$GENEID))
df.top$SYMBOL=geneIDs1$SYMBOL
df.top$ENSEMBL=rownames(df.top)


#matrix for heatmap plotting
rlog_out <- vst(dds, blind=FALSE) #get normalized count data from dds object
mat<-assay(rlog_out)[rownames(df.top), coldata$samples] #sig genes x samples
colnames(mat) <- coldata$samples
row.names(mat)=df.top$SYMBOL
row.names(df.top)=df.top$SYMBOL
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mat)
num_keep <- 20
Genez <- rbind(tail(df.top[order(df.top$log2FoldChange),], num_keep), (head(df.top[order(df.top$log2FoldChange),], num_keep)))
Genez <- Genez[order(Genez$log2FoldChange, decreasing = TRUE),]
rows_keep <- rownames(Genez)
matrixx <- mat.scaled[rows_keep,]

subset(df.top,rownames(df.top) == "RET")

#select top significantly regulated genes
MCF7TAM_sigup=rownames(df.top[(df.top$padj < 0.05) & (df.top$log2FoldChange > 0.5),])
MCF7TAM_sigdown=rownames(df.top[(df.top$padj < 0.05) & (df.top$log2FoldChange < -0.5),])
write.csv(df.top, "MCF7TAM_DES_03312023.csv")
write.csv(df,"MCF7TAM_alloutputDES.csv")

#allsignificants
MCF7TAM_allsigup=rownames(df[(df$padj < 0.05) & (df$log2FoldChange > 0),])
MCF7TAM_allsigdown=rownames(df[(df$padj < 0.05) & (df$log2FoldChange < 0),])

MCF7TAM_increased=rownames(df[(df$padj < 0.05) & (df$log2FoldChange > 0.5),])
MCF7TAM_decreased=rownames(df[(df$padj < 0.05) & (df$log2FoldChange < -0.5),])

#prepare tree file
tree=df[(df$padj < 0.05),]
tree=tree[c(2)]
tree$direction=ifelse(tree$log2FoldChange>0, "+1","-1")
tree$absolute=abs(tree$log2FoldChange)
tree=tree[c(2,3)]
write.csv(tree,"MCF7TAM_sigupdown_03312023.csv")


##grouping samples and plotting heatmap
#define sample groups
grps <- groups[c(2)]
grps
grps<- ifelse(grps == "control" , "Vehicle", "Vandetanib") 
grps
colnames(grps)=NULL
rownames(grps)=NULL
grps
GRP1 <- as.matrix(grps)
ann <- data.frame(GRP1)
colnames(ann) <- c('Condition')
colours <- list('Condition' = c("Vehicle" = "#c6dbef", "Vandetanib" =  "#6baed6"))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_height = 0.5,
                            gap = unit(0.5, 'mm'),
                            show_annotation_name = FALSE)


library(ComplexHeatmap)
library(viridis)
Heatmap(matrixx, 
        cluster_rows = F,
        name="Expression",
        cluster_columns = F,
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 8),
        column_title_gp = gpar(fontsize=7),
        column_title_side = "bottom",
        col=viridis(n=5,option="C"),
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_max_width = unit(10,"in"),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        column_order = order(as.numeric(gsub("column", "", names(data)))),
        row_order = order(as.numeric(gsub("row", "", rownames(Genez)))),
        top_annotation = colAnn,
        heatmap_legend_param = list(legend_direction = "horizontal"),
        height = unit(120, "mm"),
        width = unit(15, "mm"))

###volcano plot
##volcano plot
library(EnhancedVolcano)
sc <-df

#select significant genes and genes of interest
subset_sigreg=sc[(abs(sc$log2FoldChange)>0.5) & (sc$padj < 0.05),]
ensembl.genes <- rownames(subset_sigreg)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
listo=geneIDs1$GENEID
sub=subset(subset_sigreg, rownames(subset_sigreg) %in% listo)
head(geneIDs1)
head(sub)
individual=subset(sc, !rownames(sc) %in% listo)
individual$SYMBOL="NotMapped"
sub$SYMBOL=geneIDs1$SYMBOL
merged=rbind(sub,individual)
merged$ENSEMBL=rownames(merged)

genezadd=c("RET","EGFR")
genez=subset(sub,abs(sub$log2FoldChange)>1.5)
genez=genez$SYMBOL
genez=unlist(strsplit(genez,","))
genez=c(genezadd,genez)
sub2=subset(merged, merged$SYMBOL %in%genez)

#assign grouping
keyvals <- ifelse(merged$log2FoldChange < -0.5 & merged$padj < 0.05, 'blue', ifelse(merged$log2FoldChange > 0.5 & merged$padj < 0.05, 'red','grey'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'grey'] <- 'Intermediate'
names(keyvals)[keyvals == 'blue'] <- 'Downregulated'
names(keyvals)[keyvals == 'red'] <- 'Upregulated'

EnhancedVolcano(merged,
                lab = merged$SYMBOL,
                x = 'log2FoldChange',
                y = 'padj',
                title = NULL,
                subtitle="",
                subtitleLabSize = 9,
                axisLabSize = 10,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3,
                labSize = 3,
                colCustom=keyvals,
                colAlpha = 0.2,
                boxedLabels=FALSE,
                selectLab = sub2$SYMBOL,
                ylab = bquote(~-Log[10] ~ FDR),
                drawConnectors = FALSE,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendPosition = 'right',
                captionLabSize = 0,
                legendLabSize = 11,
                legendIconSize = 3.0,
                max.overlaps = 20,
                cutoffLineType = 'dashed',
                cutoffLineCol = 'grey70',
                cutoffLineWidth = 0.4,
                border='full',
                borderWidth = 0.6,
                borderColour = 'black')
#9x11


####MCF7EXE
set <- CV_MCF7EXE[c(1:4,6)]
groups <- data.frame(samples=colnames(set),category = c(rep("control", 3), rep("treated",2)))
coldata=groups

#data filtering and normalization
df=set[rowSums(set) > 0,]
library(edgeR)
keep <- filterByExpr(df, design = mod)
df <- df[keep,]
library(NOISeq)
data=uqua(df, long = 1000, lc = 0, k = 0)
data=as.data.frame(data)

#run DES
library(DESeq2)
roundedcounts=round(data)
groups$category=factor(groups$category)
dds <- DESeqDataSetFromMatrix(countData = roundedcounts,
                              colData = groups,
                              design = ~ category)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$category <- relevel(dds$category, ref = "control")
dds <- DESeq(dds)
res <- results(dds)
DESeq_output <- data.frame(res)
df <- na.omit(DESeq_output)
DES_MCF7EXE=na.omit(DESeq_output)
df=DES_MCF7EXE
df.top <- df[(df$padj < 0.05) & (abs(df$log2FoldChange) > 0.5),]


#gene symbols for top genes and rearrange
library("EnsDb.Hsapiens.v79")
ensembl.genes <- rownames(df.top)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
df.top=subset(df.top, rownames(df.top) != setdiff(rownames(df.top), geneIDs1$GENEID))
df.top$SYMBOL=geneIDs1$SYMBOL
df.top$ENSEMBL=rownames(df.top)

#matrix for heatmap plotting
rlog_out <- vst(dds, blind=FALSE) #get normalized count data from dds object
mat<-assay(rlog_out)[rownames(df.top), coldata$samples] #sig genes x samples
colnames(mat) <- coldata$samples
row.names(mat)=df.top$SYMBOL
row.names(df.top)=df.top$SYMBOL
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mat)
num_keep <- 20
Genez <- rbind(tail(df.top[order(df.top$log2FoldChange),], num_keep), (head(df.top[order(df.top$log2FoldChange),], num_keep)))
Genez <- Genez[order(Genez$log2FoldChange, decreasing = TRUE),]
rows_keep <- rownames(Genez)
matrixx <- mat.scaled[rows_keep,]


#select top significantly regulated genes
MCF7EXE_sigup=rownames(df.top[(df.top$padj < 0.05) & (df.top$log2FoldChange > 0.5),])
MCF7EXE_sigdown=rownames(df.top[(df.top$padj < 0.05) & (df.top$log2FoldChange < -0.5),])
write.csv(df.top, "MCF7EXE_DES_03312023.csv")
write.csv(df,"MCF7EXE_alloutputDES.csv")

#allsignificants
MCF7EXE_allsigup=rownames(df[(df$padj < 0.05) & (df$log2FoldChange > 0),])
MCF7EXE_allsigdown=rownames(df[(df$padj < 0.05) & (df$log2FoldChange < 0),])

MCF7EXE_increased=rownames(df[(df$padj < 0.05) & (df$log2FoldChange > 0.5),])
MCF7EXE_decreased=rownames(df[(df$padj < 0.05) & (df$log2FoldChange < -0.5),])

#prepare tree file
tree=df[(df$padj < 0.05),]
tree=tree[c(2)]
tree$direction=ifelse(tree$log2FoldChange>0, "+1","-1")
tree$absolute=abs(tree$log2FoldChange)
tree=tree[c(2,3)]
write.csv(tree,"MCF7EXE_sigupdown_03312023.csv")

##grouping samples and plotting heatmap
#define sample groups
grps <- groups[c(2)]
grps
grps<- ifelse(grps == "control" , "Vehicle", "Vandetanib") 
grps
colnames(grps)=NULL
rownames(grps)=NULL
grps
GRP1 <- as.matrix(grps)
ann <- data.frame(GRP1)
colnames(ann) <- c('Condition')
colours <- list('Condition' = c("Vehicle" = "#c6dbef", "Vandetanib" =  "#6baed6"))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_height = 0.5,
                            gap = unit(0.5, 'mm'),
                            show_annotation_name = FALSE)


library(ComplexHeatmap)
library(viridis)
Heatmap(matrixx, 
        cluster_rows = F,
        name="Expression",
        cluster_columns = F,
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 8),
        col=viridis(n=6,option="C"),
        column_title_gp = gpar(fontsize=7),
        column_title_side = "bottom",
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_max_width = unit(10,"in"),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        column_order = order(as.numeric(gsub("column", "", names(data)))),
        row_order = order(as.numeric(gsub("row", "", rownames(Genez)))),
        top_annotation = colAnn,
        heatmap_legend_param = list(legend_direction = "horizontal"),
        height = unit(120, "mm"),
        width = unit(20, "mm"))


###volcano plot
##volcano plot
library(EnhancedVolcano)
sc <-df

#select significant genes and genes of interest
subset_sigreg=sc[(abs(sc$log2FoldChange)>0.5) & (sc$padj < 0.05),]
ensembl.genes <- rownames(subset_sigreg)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
listo=geneIDs1$GENEID
sub=subset(subset_sigreg, rownames(subset_sigreg) %in% listo)
head(geneIDs1)
head(sub)
individual=subset(sc, !rownames(sc) %in% listo)
individual$SYMBOL="NotMapped"
sub$SYMBOL=geneIDs1$SYMBOL
merged=rbind(sub,individual)
merged$ENSEMBL=rownames(merged)

genezadd=c("RET","EGFR")
genez=subset(sub,abs(sub$log2FoldChange)>1.5)
genez=genez$SYMBOL
genez=unlist(strsplit(genez,","))
genez=c(genezadd,genez)
sub2=subset(merged, merged$SYMBOL %in%genez)

#assign grouping
keyvals <- ifelse(merged$log2FoldChange < -0.5 & merged$padj < 0.05, 'blue', ifelse(merged$log2FoldChange > 0.5 & merged$padj < 0.05, 'red','grey'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'grey'] <- 'Intermediate'
names(keyvals)[keyvals == 'blue'] <- 'Downregulated'
names(keyvals)[keyvals == 'red'] <- 'Upregulated'

EnhancedVolcano(merged,
                lab = merged$SYMBOL,
                x = 'log2FoldChange',
                y = 'padj',
                title = NULL,
                subtitle="",
                subtitleLabSize = 9,
                axisLabSize = 10,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3,
                labSize = 3,
                colCustom=keyvals,
                colAlpha = 0.2,
                boxedLabels=FALSE,
                selectLab = sub2$SYMBOL,
                ylab = bquote(~-Log[10] ~ FDR),
                drawConnectors = FALSE,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendPosition = 'right',
                captionLabSize = 0,
                legendLabSize = 11,
                legendIconSize = 3.0,
                max.overlaps = 20,
                cutoffLineType = 'dashed',
                cutoffLineCol = 'grey70',
                cutoffLineWidth = 0.4,
                border='full',
                borderWidth = 0.6,
                borderColour = 'black')

#9x11

