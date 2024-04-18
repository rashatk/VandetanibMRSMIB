#Plot all significant genes
dftop_MCF7=read.csv("MCF7_alloutputDES.csv",header=TRUE,row.names=1)
dftop_MCF7TAM=read.csv("MCF7TAM_alloutputDES.csv",header=TRUE,row.names=1)
dftop_MCF7EXE=read.csv("MCF7EXE_alloutputDES.csv",header=TRUE,row.names=1)


#plot all overlapping genes
down_mcf7=subset(dftop_MCF7, rownames(dftop_MCF7) %in% xALL_down)
down_mcf7$line="MCF7"
ensembl.genes <- rownames(down_mcf7)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
down_mcf7$SYMBOL=geneIDs1$SYMBOL
rownames(down_mcf7)=NULL
up_mcf7=subset(dftop_MCF7, rownames(dftop_MCF7) %in% xALL_up)
up_mcf7$line="MCF7"
ensembl.genes <- rownames(up_mcf7)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
up_mcf7$SYMBOL=geneIDs1$SYMBOL
rownames(up_mcf7)=NULL

down_mcf7tam=subset(dftop_MCF7TAM, rownames(dftop_MCF7TAM) %in% xALL_down)
down_mcf7tam$line="MCF7TAM"
ensembl.genes <- rownames(down_mcf7tam)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
down_mcf7tam$SYMBOL=geneIDs1$SYMBOL
rownames(down_mcf7tam)=NULL
up_mcf7tam=subset(dftop_MCF7TAM, rownames(dftop_MCF7TAM) %in% xALL_up)
up_mcf7tam$line="MCF7TAM"
ensembl.genes <- rownames(up_mcf7tam)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
up_mcf7tam$SYMBOL=geneIDs1$SYMBOL
rownames(up_mcf7tam)=NULL

down_mcf7exe=subset(dftop_MCF7EXE, rownames(dftop_MCF7EXE) %in% xALL_down)
down_mcf7exe$line="MCF7EXE"
ensembl.genes <- rownames(down_mcf7exe)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
down_mcf7exe$SYMBOL=geneIDs1$SYMBOL
rownames(down_mcf7exe)=NULL
up_mcf7exe=subset(dftop_MCF7EXE, rownames(dftop_MCF7EXE) %in% xALL_up)
up_mcf7exe$line="MCF7EXE"
ensembl.genes <- rownames(up_mcf7exe)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
up_mcf7exe$SYMBOL=geneIDs1$SYMBOL
rownames(up_mcf7exe)=NULL


merged=rbind(down_mcf7, down_mcf7tam,down_mcf7exe,up_mcf7,up_mcf7tam,up_mcf7exe)
merged$direction=ifelse(merged$log2FoldChange>0,"UP","DOWN")
DOWN=subset(merged, merged$direction=="DOWN")
UP=subset(merged,merged$direction=="UP")
merged=rbind(DOWN,UP)
merged$SYMBOL <- factor(merged$SYMBOL, levels=unique(as.character(merged$SYMBOL)) )


ggplot(merged, aes(x=line,y=SYMBOL,fill=log2FoldChange))+geom_tile()+theme(axis.text.x = element_text(angle = 90))+ylab("")+xlab("")+ theme(axis.line = element_blank(),axis.ticks = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white"),text = element_text(size=13), axis.text.x = element_text(vjust = 1, hjust = 1),legend.position="right",legend.text = element_text(size=10), legend.title =element_text(size=10))+scale_fill_gradient2(low="blue",high="red",name="Log2FC")



#Plot only notable genes after searching
#ALL *******
library("EnsDb.Hsapiens.v79")
ensembl.genes <- c("SELENBP1","AURKA","TPX2","SUN2","ESRP1","CCND1","ECT2","CDC20","TOP2A","CCNB1","LY6E","PLK1","SHMT2","IQGAP3","FTL","BAG1","WBP2")
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
genelist=geneIDs1$GENEID

mcf7=subset(dftop_MCF7, rownames(dftop_MCF7) %in% genelist)
ge=ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(mcf7), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
mcf7$genes=ge$SYMBOL
mcf7$line="MCF7"
rownames(mcf7)=NULL

mcf7tam=subset(dftop_MCF7TAM, rownames(dftop_MCF7TAM) %in% genelist)
ge=ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(mcf7tam), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
mcf7tam$genes=ge$SYMBOL
mcf7tam$line="MCF7TAM"
rownames(mcf7tam)=NULL

mcf7exe=subset(dftop_MCF7EXE, rownames(dftop_MCF7EXE) %in% genelist)
ge=ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(mcf7exe), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
mcf7exe$genes=ge$SYMBOL
mcf7exe$line="MCF7EXE"
rownames(mcf7exe)=NULL

merged=rbind(mcf7exe,mcf7tam,mcf7)
merged$direction=ifelse(merged$log2FoldChange>0,"UP","DOWN")
DOWN=subset(merged, merged$direction=="DOWN")
UP=subset(merged,merged$direction=="UP")
merged=rbind(DOWN,UP)
merged=as.data.frame(merged)
names(merged)
merged=merged[c(2,7,8,9)]
merged$genes <- factor(merged$genes, levels=unique(as.character(merged$genes)) )


ggplot(merged, aes(x = line,y=genes,fill=log2FoldChange))+geom_tile()+theme(axis.text.x = element_text(angle = 90,hjust=1))+ylab("")+xlab("")+ theme(axis.line = element_blank(),axis.text.x = element_text(size=12,face = "bold"),axis.ticks = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white"),text = element_text(size=15,), legend.text = element_text(size=10),legend.title =element_text(size=10,face="bold"))+scale_fill_gradient2(low="blue",high="red",name="Log2FC")
