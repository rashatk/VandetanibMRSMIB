install.packages("reshape2")
library(reshape2)
library(ComplexHeatmap)

#merge line gene sets
nrow(MCF7sets)
nrow(MCF7TAMsets)
nrow(MCF7EXEsets)
mergedo=rbind(MCF7sets, MCF7TAMsets, MCF7EXEsets)
write.csv(mergedo,"mergedo")
nrow(mergedo)

#load
mergedo=read.csv("mergedo",header=TRUE,row.names = 1)

##attempt1
dff=mergedo[c(1,10,12,15)]
riko=dcast(dff, ID +category~ line, value.var = "GRmod", fun.aggregate = mean)
riko=riko[order (riko$category), ]
riko=riko[ ! duplicated(riko$ID), ]
rownames(riko)=riko$ID

design=riko
riko=riko[-c(1,2)]
riko$MCF7[is.nan(riko$MCF7)]<-NA
riko$MCF7TAM[is.nan(riko$MCF7TAM)]<-NA
riko$MCF7EXE[is.nan(riko$MCF7EXE)]<-NA
riko2=as.matrix(riko)

grps <- design$category
grps
colnames(grps)=NULL
rownames(grps)=NULL
grps
GRP1 <- as.matrix(grps)
ann <- data.frame(GRP1)
colnames(ann) <- c('Gene Set Category')
colours <- list('Gene Set Category' = c("AE" = "#c6dbef", "AKT" =  "#6baed6", "EGFR" =  "navy","ERBB2" =  "pink","Estrogen & ESR1" =  "orange","KRAS" =  "red","MAPK" =  "skyblue","MYC" =  "grey","NFKB" =  "purple","Phenotype" =  "green","Proliferation" =  "tomato"))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'row',
                            col = colours,
                            annotation_height = 0.5,
                            gap = unit(0.5, 'mm'),
                            show_annotation_name = FALSE)

#heatmap
Heatmap(riko2, 
        na_col = "white", 
        cluster_rows = F,
        name= "dGR",
        row_names_gp = gpar(fontsize = 10),
        row_names_max_width = unit(10,"in"),
        column_names_gp = gpar(fontsize=9),
        height = unit(200, "mm"),
        width = unit(15, "mm"),
        right_annotation = colAnn)



#favorable
favorablesets=rownames(riko2)

#subsetmassive
MCF7subsetmassive=read.csv("EnricherOutput_MCF7.csv",row.names = 1,header = TRUE)
MCF7TAMsubsetmassive=read.csv("EnricherOutput_MCF7TAM.csv",row.names = 1,header = TRUE)
MCF7EXEsubsetmassive=read.csv("EnricherOutput_MCF7EXE.csv",row.names = 1,header = TRUE)

MCF7zz=subset(MCF7subsetmassive, MCF7subsetmassive$ID %in% favorablesets)
MCF7zz$line="MCF7"
MCF7TAMzz=subset(MCF7TAMsubsetmassive, MCF7TAMsubsetmassive$ID %in% favorablesets)
MCF7TAMzz$line="MCF7TAM"
MCF7EXEzz=subset(MCF7EXEsubsetmassive, MCF7EXEsubsetmassive$ID %in% favorablesets)
MCF7EXEzz$line="MCF7EXE"
mergedo=rbind(MCF7zz, MCF7TAMzz, MCF7EXEzz)

##attempt1
colnames(mergedo)
dff=mergedo[c(1,10,12,13)]
riko=dcast(dff, ID +category~ line, value.var = "GRmod", fun.aggregate = mean)
riko=riko[order (riko$category), ]
riko=riko[ ! duplicated(riko$ID), ]
rownames(riko)=riko$ID

design=riko
riko=riko[-c(1,2)]
riko$MCF7[is.nan(riko$MCF7)]<-NA
riko$MCF7TAM[is.nan(riko$MCF7TAM)]<-NA
riko$MCF7EXE[is.nan(riko$MCF7EXE)]<-NA
riko2=as.matrix(riko)


grps <- design$category
grps
colnames(grps)=NULL
rownames(grps)=NULL
grps <- gsub('AE', 'Anti-Estrogen', grps)
grps <- gsub('NFKB', 'zz', grps)
grps <- gsub('Phenotype', 'zz', grps)
grps <- gsub('Proliferation', 'zz', grps)
grps <- gsub('zz', 'Phenotype & Proliferation', grps)
grps <- gsub('ERBB2', 'ERBB', grps)
grps
GRP1 <- as.matrix(grps)
ann <- data.frame(GRP1)
colnames(ann) <- c('Gene Set Category')
colours <- list('Gene Set Category' = c("Anti-Estrogen" = "skyblue", "AKT" =  "#6baed6", "EGFR" =  "dodgerblue","ERBB" =  "lavenderblush2","Estrogen & ESR1" =  "pink","KRAS" =  "deeppink3","MAPK" =  "slateblue1","MYC" =  "purple3", "Phenotype & Proliferation" =  "darkseagreen"))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'row',
                            col = colours,
                            annotation_height = 0.5,
                            gap = unit(0.5, 'mm'),
                            show_annotation_name = FALSE)

#heatmap
Heatmap(riko2, 
        na_col = "white", 
        cluster_rows = F,
        name= "dGR",
        row_names_gp = gpar(fontsize = 10),
        row_names_max_width = unit(10,"in"),
        column_names_gp = gpar(fontsize=11),
        height = unit(130, "mm"),
        width = unit(15, "mm"),
        right_annotation = colAnn)
