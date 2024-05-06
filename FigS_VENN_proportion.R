setwd('./RNA_Seq_LALO/VENN/')
# Load library
library(VennDiagram)
library(ggvenn)

### read in
compare="4_VS_3"
Mytitle='Snowstorm'
{
set1 <- read.table(paste0("./RNA_Seq_LALO/Gonad/GENE_Gonad_Conditions_",compare,"_Volcano_ALL.txt"),header=F) #Gonad
set2 <- read.table(paste0("./RNA_Seq_LALO/Hypothalamus/GENE_Hypothalamus_Conditions_",compare,"_Volcano_ALL.txt"),header=F) #Hyp
set3 <- read.table(paste0("./RNA_Seq_LALO/Liver/GENE_Liver_Conditions_",compare,"_Volcano_ALL.txt"),header=F) #Liver
set4 <- read.table(paste0("./RNA_Seq_LALO/Heart/GENE_Heart_Conditions_",compare,"_Volcano_ALL.txt"),header=F) #Heart
}

##########################
###Regular venn plot
###ggven
##########################
library(ggvenn)
x1 <- list(
  Gonad = set1$V1, 
  Hypothalamus = set2$V1, 
  Liver = set3$V1,
  Heart = set4$V1
)
#ID<-names(x1)
IDs<-lengths(x1)
ID<-list()
for (i in 1:length(IDs)){
  #print(i)
  N<-print(paste0(names(IDs)[i],"\n",IDs[i]))
  ID<-append(ID,N)
}
names(x1)<-ID

df <-
  ggvenn(
    x1,
    show_percentage = F ,
    show_elements = F, #gene list
    fill_color = c("#00AFBB","#FC4E07",'darkblue',"#E7B800"),
    fill_alpha = 0.3,
    stroke_size = 0.5,
    text_size = 3, #number inside
    set_name_size = 2 #name of group
  )+
  #scale_x_continuous(expand = expansion(mult = .2))+
  scale_y_continuous(expand = expansion(mult = .2))+
  ggtitle( Mytitle) +
  theme(plot.title = element_text(size=10,hjust = 0.5))
  
df 
#the label of gonad will be cut, unclear why, solution: https://gaospecial.github.io/ggVennDiagram/articles/using-ggVennDiagram.html#:~:text=(ggplot2)-,long%20category%20names,-If%20you%20use
ggsave(paste0(compare,'.pdf'),df,width = 3,height = 3) #,dpi = 600

tail(list_to_data_frame(x1))
x2<-list_to_data_frame(x1)

dev.off()

###########
### Fat and liver
###########
compare="43_liver_VS_fat"
{
  set1 <- read.table(paste0("./RNA_Seq_LALO/Liver/GENE_Liver_Conditions_4_VS_3_Volcano_UP.txt"),header=F) #Liver up
  set2 <- read.table(paste0("./RNA_Seq_LALO/Liver/GENE_Liver_Conditions_4_VS_3_Volcano_DOWN.txt"),header=F) #Liver down
  ####Different dir
  set3 <- read.table(paste0("./RNA_2LALO/FAT/GENE_FAT_Conditions_4_VS_3_Volcano_UP.txt"),header=F) #Fat up
  set4 <- read.table(paste0("./RNA_2LALO/FAT/GENE_FAT_Conditions_4_VS_3_Volcano_DOWN.txt"),header=F) #Fat down
}

library(ggvenn)
x3 <- list(
  Liver_down = set2$V1,
  Liver_up = set1$V1, 
  
  Fat_up = set3$V1,
  Fat_down = set4$V1
)

IDs<-lengths(x3)
ID<-list()
for (i in 1:length(IDs)){
  #print(i)
  N<-print(paste0(names(IDs)[i],"\n",IDs[i]))
  ID<-append(ID,N)
}
names(x3) <-ID

df <-
  ggvenn(
    x3,
    show_percentage = F ,
    show_elements = F, #gene list
    fill_color = c("azure3","azure",'darkgoldenrod1',"darkgoldenrod4"),
    fill_alpha = 0.5,
    stroke_size = 0.5,
    text_size = 3, #number inside
    set_name_size = 2 #name of group
  )+
  #scale_x_continuous(expand = expansion(mult = .2))+
  scale_y_continuous(expand = expansion(mult = .2))
df 
#the label of gonad will be cut, unclear why, solution: https://gaospecial.github.io/ggVennDiagram/articles/using-ggVennDiagram.html#:~:text=(ggplot2)-,long%20category%20names,-If%20you%20use
ggsave(paste0(compare,'.png'),df,width = 3,height = 3,dpi = 600)

dev.off()

################
### 21 vs 43###
################
library(eulerr)
#genes <- paste("gene",1:1000,sep="")
compare="Gonad"
{
  set1 <- read.table(paste0("./RNA_Seq_LALO/",compare,"/GENE_",compare,"_Conditions_1_VS_2","_Volcano_ALL.txt"),header=F) #Liver 12
  set2 <- read.table(paste0("./RNA_Seq_LALO/",compare,"/GENE_",compare,"_Conditions_4_VS_3","_Volcano_ALL.txt"),header=F) #Liver 43
  # set3 <- read.table(paste0("./RNA_2LALO/FAT/GENE_FAT_Conditions_4_VS_3_Volcano_UP.txt"),header=F) #Fat up
  # set4 <- read.table(paste0("./RNA_2LALO/FAT/GENE_FAT_Conditions_4_VS_3_Volcano_DOWN.txt"),header=F) #Fat down
}

x4 <- list(
  Extreme_spring = set1, 
  Snowstorm = set2
)

###nVennR-proportional-atypical
library(nVennR)
obj<-plotVenn(x4, setColors=c("skyblue", "pink1"), borderWidth = 0) #, "mediumorchid",'gray50'
df2<-listVennRegions(obj)
df3<-lapply(df2,length)

df4<-stringr::str_split(names(df3),'[[(]]',simplify = T) [,2]
df4<-stringr::str_replace(df4,'[[)]]' , "" )
df4<-stringr::str_replace_all(df4,',' , "&" )
df4<-stringr::str_replace_all(df4,' ' , "" )

df3<-as.data.frame(df3)
names(df3)<-df4
df3


for (i in 1:length(df3)){
  GRP=colnames(df3)[i]
  VALUE=df3[i]
  #print(GRP)
  cat(paste0('"',GRP,'"',"=",VALUE,','))
}

#####manually copy-paste#####
fit1<-euler(c(
  "Snowstorm"=99,"Extreme_spring"=10
  ))
p1<-plot(fit1,
         #quantities = TRUE,
         quantities=list(cex=1.5),
         # fill = "transparent",
         fills = list(fill = c("#FC4E07",'darkblue',"#00AFBB","#E7B800"), alpha = 0.5),
         lty = 1, #1:4,
         legend = T,
         #labels = list(font =2,cex=0.8), #labels = list(font =2,cex=0.8,pos = 4,y =c(3,-3 ),just = c( "bottom")),
         main = list(label=compare,y=-0.2,cex=1.2))
p1
dev.off()

pdf(file=paste0("21_vs_43_",compare,'.pdf'),
    width=5, height=4) #, res=600
p1
dev.off()

####end#####
