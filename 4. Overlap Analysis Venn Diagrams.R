##--repeat for each cell line and run below lines for each
DOWN_MCF7EXE=MCF7EXE_sigdown
UP_MCF7EXE=MCF7EXE_sigup
  
DOWN_MCF7TAM=MCF7TAM_sigdown
UP_MCF7TAM=MCF7TAM_sigup
  
DOWN_MCF7=MCF7_sigdown
UP_MCF7=MCF7_sigup

#venn diagram
library(VennDiagram)
library(RColorBrewer)
library(eulerr)

#down
listo <- list(MCF7EXE = DOWN_MCF7EXE, MCF7TAM = DOWN_MCF7TAM, MCF7=DOWN_MCF7)
listo
col <- c('#a3e4e6', '#fddd8c', '#efa5c1')
fit1 <- euler(listo)
fit1
plot(fit1,
     quantities = list(cex = 1.4, alpha = 0.6),
     labels = list(cex = 1.5, alpha = 0.9),
     fills = list(fill = col, alpha = 0.6),
     edges = NULL)

venn.diagram(x = list(DOWN_MCF7TAM,DOWN_MCF7EXE,DOWN_MCF7), main="Downregulated with Vandetanib",main.fontfamily = "sans",main.cex = 0.6,main.fontface = "bold",category.names = c("MCF7TAM" , "MCF7EXE" , "MCF7"),filename = 'DownVenn05082023.png',
  output=TRUE,imagetype="png",height = 600 , width = 600 , resolution = 300, compression = "lzw",lwd = 2,lty = 'blank', fill = col,cex = .6,fontface = "bold",  fontfamily = "sans",cat.cex = 0.6,cat.default.pos = "outer",cat.pos = c(-27, 27, 135),cat.dist = c(0.055, 0.055, 0.085),cat.fontfamily = "sans",rotation = 1)


library(VennDiagram)
overlap=calculate.overlap(x=listo)
overlap
df <- data.frame(lapply(overlap, function(x) {
  x <- unlist(x)
  length(x) <- max(lengths(overlap))
  return(x)
}))

na.omit(length(unique(df$a5)))
na.omit(length(unique(df$a2)))
na.omit(length(unique(df$a4)))
na.omit(length(unique(df$a6)))
na.omit(length(unique(df$a1)))
na.omit(length(unique(df$a3)))
na.omit(length(unique(df$a7)))

ALL_down=na.omit(unique(df$a5))
EXE_TAM_down=na.omit(unique(df$a2))
M_EXE_down=na.omit(unique(df$a4))
M_TAM_down=na.omit(unique(df$a6))
Unique_EXE_down=na.omit(unique(df$a1))
Unique_TAM_down=na.omit(unique(df$a3))
Unique_M_down=na.omit(unique(df$a7))

names(df)=c("ALL_down","EXE_TAM_down","M_EXE_down","M_TAM_down","Unique_EXE_down","Unique_TAM_down","Unique_M_down")
library("EnsDb.Hsapiens.v79")
ensembl.genes <- ALL_down
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
xALL_down=geneIDs1$GENEID
ensembl.genes <- EXE_TAM_down
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
xEXE_TAM_down=geneIDs1$GENEID
ensembl.genes <- M_EXE_down
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
xM_EXE_down=geneIDs1$GENEID
ensembl.genes <- M_TAM_down
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
xM_TAM_down=geneIDs1$GENEID
ensembl.genes <- Unique_EXE_down
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
xUnique_EXE_down=geneIDs1$GENEID
ensembl.genes <- Unique_TAM_down
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
xUnique_TAM_down=geneIDs1$GENEID
ensembl.genes <- Unique_M_down
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
xUnique_M_down=geneIDs1$GENEID
list2=list(ALL_down,EXE_TAM_down,M_TAM_down,M_EXE_down,Unique_M_down,Unique_TAM_down,Unique_EXE_down)
df2 <- data.frame(lapply(list2, function(x) {
  x <- unlist(x)
  length(x) <- max(lengths(list2))
  return(x)
}))
names(df2)=c("ALL_down","EXE_TAM_down","M_TAM_down","M_EXE_down","Unique_M_down","Unique_TAM_down","Unique_EXE_down")
write.csv(df2,"Overlaps_mRNAseq_DOWN_05082023.csv")


#up
listo <- list(MCF7EXE = UP_MCF7EXE, MCF7TAM = UP_MCF7TAM, MCF7=UP_MCF7)
listo
col <- c('#a3e4e6', '#fddd8c', '#efa5c1')
fit1 <- euler(listo)
fit1
plot(fit1,
     quantities = list(cex = 1.4, alpha = 0.6),
     labels = list(cex = 1.5, alpha = 0.9),
     fills = list(fill = col, alpha = 0.6),
     edges = NULL)

venn.diagram(x = list(UP_MCF7TAM,UP_MCF7EXE,UP_MCF7), main="Upregulated with Vandetanib",main.fontfamily = "sans",main.cex = 0.6,main.fontface = "bold",category.names = c("MCF7TAM" , "MCF7EXE" , "MCF7"),filename = 'UpVenn05082023.png',output=TRUE,imagetype="png",height = 600 , width = 600 , resolution = 300, compression = "lzw",lwd = 2,lty = 'blank', fill = col,cex = .6,fontface = "bold",  fontfamily = "sans",cat.cex = .6,cat.default.pos = "outer",cat.pos = c(-27, 27, 135),cat.dist = c(0.055, 0.055, 0.085),cat.fontfamily = "sans",rotation = 1)

library(VennDiagram)
overlap=calculate.overlap(x=listo)
overlap
df <- data.frame(lapply(overlap, function(x) {
  x <- unlist(x)
  length(x) <- max(lengths(overlap))
  return(x)
}))

na.omit(length(unique(df$a5)))
na.omit(length(unique(df$a2)))
na.omit(length(unique(df$a4)))
na.omit(length(unique(df$a6)))
na.omit(length(unique(df$a1)))
na.omit(length(unique(df$a3)))
na.omit(length(unique(df$a7)))

ALL_up=na.omit(unique(df$a5))
EXE_TAM_up=na.omit(unique(df$a2))
M_EXE_up=na.omit(unique(df$a4))
M_TAM_up=na.omit(unique(df$a6))
Unique_EXE_up=na.omit(unique(df$a1))
Unique_TAM_up=na.omit(unique(df$a3))
Unique_M_up=na.omit(unique(df$a7))

names(df)=c("ALL_up","EXE_TAM_up","M_EXE_up","M_TAM_up","Unique_EXE_up","Unique_TAM_up","Unique_M_up")
library("EnsDb.Hsapiens.v79")
ensembl.genes <- ALL_up
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
xALL_up=geneIDs1$GENEID
ensembl.genes <- EXE_TAM_up
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
xEXE_TAM_up=geneIDs1$GENEID
ensembl.genes <- M_EXE_up
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
xM_EXE_up=geneIDs1$GENEID
ensembl.genes <- M_TAM_up
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
xM_TAM_up=geneIDs1$GENEID
ensembl.genes <- Unique_EXE_up
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
xUnique_EXE_up=geneIDs1$GENEID
ensembl.genes <- Unique_TAM_up
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
xUnique_TAM_up=geneIDs1$GENEID
ensembl.genes <- Unique_M_up
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
xUnique_M_up=geneIDs1$GENEID
list2=list(ALL_up,EXE_TAM_up,M_TAM_up,M_EXE_up,Unique_M_up,Unique_TAM_up,Unique_EXE_up)
df2 <- data.frame(lapply(list2, function(x) {
  x <- unlist(x)
  length(x) <- max(lengths(list2))
  return(x)
}))
names(df2)=c("ALL_up","EXE_TAM_up","M_TAM_up","M_EXE_up","Unique_M_up","Unique_TAM_up","Unique_EXE_up")
write.csv(df2,"Overlaps_mRNAseq_UP_05082023.csv")
