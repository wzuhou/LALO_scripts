##############
### LALO heatmap
##############
#Correlation plot
{ library(DESeq2)
  library(edgeR)
  library(dplyr)
  library(limma)
  library("RColorBrewer")
  library("ggplot2")}
##########################

exprSet=assay(dds)
group_list=colData(dds)
#rm(list = ls())

# ##############
# #Top genes
# ##############
# TopComp<-results(dds, contrast=c("group", '3', '2'))
# TopComp
# head(TopComp$pvalue)
# TopComp2<-TopComp[order(TopComp$padj),]
# test <- head(TopComp2,50)
# choose_gene=head(rownames(TopComp2),50) ## 50 maybe better
# head(TopComp2,50)

###############

TS="Gonad"
# TS="Heart"
TS="Hypothalamus"
 TS="Liver"
TISSUE=c("Gonad","Heart","Hypothalamus","Liver")#


  ###############
  group_list=colData(dds)
  exprSet=assay(dds)
  dim(exprSet)
  exprSet <- as.matrix(exprSet)
  
  GeneOfInterest<-c('<MyGeneOfInterest>') # e.g., c("FKBP5",'ME1','TKFC','PDHA2','CS')
  
  GeneOfInterest_scaffold<-c('<MyGeneOfInterest>') # e.g., c("FKBP5",'ME1','TKFC','PDHA2','CS')
  
  ann_colors = list(
    Tissue = c("Gonad"="#FF0000" , "Heart"="#00A08A","Hypothalamus"="#F98400", "Liver"="#5BBCD6"),
    LH_Sub_Stage = c(Arrival = "#0072B2",  Incubation = "#F0E442" ),
    Weather=c(Extreme_Spring="#F8AFA8",No_Storm="#74A089",Storm="#FDDDA0"),
    Conditions=c('1'="#81A88D",'2'="#FDD262",'3'="#D3DDDC",'4'="#C7B19C"))
  
  Order <- row.names(sampleTable[order(sampleTable$Conditions),])
  sortkey <- factor(row.names(sampleTable[order(sampleTable$Tissue),]), levels = Order)
  exprSet <- as.matrix(exprSet)
  exprSet <- exprSet[,sortkey]
  choose_matrix=exprSet[row.names(exprSet)%in%GeneOfInterest_scaffold,]
  row.names(choose_matrix) <- str_split(row.names(choose_matrix), "[.]",simplify = T)[,1]
  choose_matrix <- choose_matrix[,Order]
  choose_matrix=t(scale(t(choose_matrix)))
  choose_matrix
  colnames(choose_matrix)
  
  #Correlation plot
  group_list
  tmp=data.frame(group_list[,c(3,11,14,20)])
  
  names(tmp) <- c("Tissue","LH_Sub_Stage","Weather","Conditions")
  tmp <- tmp[order(tmp$Conditions),]
  all(rownames(tmp) == colnames(choose_matrix))
  
  pheatmap::pheatmap(choose_matrix,
                     #treeheight_row=0,
                     color = wesanderson::wes_palette("Zissou1", 100,type="continuous"),
                     #cutree_cols = 4,
                     annotation_colors = ann_colors,
                     treeheight_col=0,
                     main=paste(TS,"RNA",sep="_"),
                     cluster_rows=T, 
                     cluster_cols=F,
                     fontsize_col= 9,
                     width=9,
                     show_colnames=F,
                     filename = paste0("./RNA_Seq_LALO/GENE_expression_plot/Heatmap/",TS,"_RNA_GeneOfInterest.png"),
                     annotation_col= tmp
                     #labels_col=tmp$Conditions#labels_col=ID$breed
  )
###END###
  
  
  #dev.off()
}

