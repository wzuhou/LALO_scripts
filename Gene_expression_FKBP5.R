{library(dplyr)
  library(ggplot2)
  library("wesanderson")
  library(DESeq2)
  library(edgeR)
  library(limma)
  library("RColorBrewer")
}

#Visualize 
#########################
#log2 expression level
#########################
{GENE="FKBP5"
gene_name="FKBP5"
TS="Hypothalamus"
GeneOfInterest_scaffold #if exist
}

{  try(expr=GENE_all <- plotCounts(dds, gene=GENE, intgroup=c("ID","Tissue","Conditions"),cex=1.2,normalized= F,returnData=T),silent = F)
#GENE_all <- cbind(GENE_all,str_split(row.names(GENE_all),"_",simplify = T)[,1:3])
head(GENE_all)
GENE_TS<- GENE_all[GENE_all$Tissue==paste0(TS),]
}

# ##ALL Tissues
  GENE_all$Tissue<-str_replace(GENE_all$Tissue,'Gonad' , "Testis")
  GENE_all$Tissue<-str_replace(GENE_all$Tissue,'ADRE' , "Adrenal")
  GENE_all$Tissue<-str_replace(GENE_all$Tissue,'FAT' , "Fat")
  GENE_all$Tissue<-str_replace(GENE_all$Tissue,'PIT' , "Pituitary")
  Myorder = c("Testis" ,
              "Heart" ,
              "Hypothalamus",
              "Liver" ,
              "Adrenal",
              "Pituitary",
              "Fat") 
  
GENE_all$count <- log2(GENE_all$count)#Log2 transform #Dont run multiple times, make sure the value are not double calculated.
  
{w1 <- ggplot(GENE_all,aes(x=Conditions,y=count,fill=factor(Tissue,levels = Myorder )))+
    geom_boxplot(alpha=0.7,color="black")+#,position = "jitter"
    #geom_point(size=2,alpha=0.7,color="black",position = position_jitter(width =0.1,height = 0))+#
    #scale_fill_manual(values =c("#00AFBB", "#E7B800", "#FC4E07",'darkblue'))+
    scale_fill_manual(values =c(wes_palette("Darjeeling1",n=5,type ="discrete"),wes_palette("Darjeeling2",n=3,type ="discrete")))+
    labs(title = paste0(gene_name))+ # ,subtitle =paste0("All conditions")
    #geom_smooth(count~ Conditions, span = 0.8)+
    facet_wrap(~ factor(Tissue,levels = Myorder ),drop=F,scales="free",ncol = 4)+ #Tissue
    theme_bw()+
    xlab(NULL) +
    ylab("log2 expression" ) + #"log2 expression count" "TPM"
    theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
          plot.subtitle = element_text(hjust = 0.5,size = 15,face = "bold"),
          axis.text=element_text(size=15),
          axis.title=element_text(size=15),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position = c(0.88, 0.25)) +
    guides(fill = guide_legend(title ="Tissues",label.hjust=0.01,label.vjust=0.5,ncol = 1))
  w1}
 ggsave(paste0("M:/ROSLIN/RNA_2LALO/GENE_expression_plot/Gene","_",gene_name,".jpg"),w1,width = 8,height = 8,dpi=500)
 ggsave(paste0("M:/ROSLIN/RNA_2LALO/GENE_expression_plot/Gene","_",gene_name,".PDF"),w1,width = 9,height =9)
  ###END###
