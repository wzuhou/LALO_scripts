library(dplyr)
library(tidyr)
library(stringr)
library(writexl)

TISSUE=c("Gonad","Heart","Hypothalamus","Liver")#
#Making the dataframe used for interation between conditions

dat <- crossing(x=1:4,y=1:4) #Case=y
dat2 <- dat %>% 
  filter(!row_number()==2)%>%
  group_by(grp = paste0(pmin(x, y), pmax(x, y))) %>%
  filter(!x==y) %>%
  filter(row_number() == min(row_number())) %>%
  ungroup() %>%
  filter(!grp%in%c("13","24",'14'))%>%
  select(-grp)%>%
  mutate(Compare=c('spring','LHS','storm'))

dat2
outsheets=NULL
#----------------------------------------------
#----------------------------------------------
  if (exists('TISSUE')) {
    ITERA=4
  } else{
    ITERA=1  }
  
  for (j in seq(ITERA)){  #seq(ITERA)
    #j=1
    try(expr=eval(parse(text=paste0("TS <- TISSUE[",j,"]"))),silent = T)
    try(expr=print(TS),silent = T)
    wd=paste0("./RNA_Seq_LALO/",TS)
    setwd(wd)

    for (i in 1:length(dat2$x)){
    #  i=3
      Base=dat2$x[i] #3 #dat2$x[i]
      Case=dat2$y[i] #4 #dat2$y[i]
    n=paste(TS,"Conditions",Case,"VS",Base,sep="_")
    n
    UP <- read.table(paste0(n,"_Volcano_UP.txt"),header=T,sep="\t")
    DOWN <- read.table(paste0(n,"_Volcano_DOWN.txt"),header=T,sep="\t")
    ALL <- rbind(UP,DOWN)
    ALL <- ALL%>%arrange(padj)
    ALL$Gene_name <- str_split(ALL$Gene, "\\.",simplify = T)[,1]
    ALL <- ALL[,c('Gene',"Gene_name","log2FoldChange",'padj')]
    
    SHEET=paste(toupper(substr(TS,1,3)),dat2$Compare[i],sep='_')
    sheets<-eval(parse(text = paste0('list(', SHEET,' = ALL)')))
    outsheets <- c( outsheets,sheets)
    }
  }


#----------------------------------------------
# RNA 2
#----------------------------------------------

TISSUE=c("PIT","ADRE","FAT")#
ITERA=3

for (j in seq(ITERA)){  #seq(ITERA)
  
  try(expr=eval(parse(text=paste0("TS <- TISSUE[",j,"]"))),silent = T)
  try(expr=print(TS),silent = T)
  wd=paste0("./RNA_2LALO/",TS)
  setwd(wd)
  
  i=3
  Base=dat2$x[i] #3 #dat2$x[i]
  Case=dat2$y[i] #4 #dat2$y[i]
  n=paste(TS,Case,"VS",Base,sep="_") #different input name
  n
  UP <- read.table(paste0(n,"_Volcano_UP.txt"),header=T,sep="\t")
  DOWN <- read.table(paste0(n,"_Volcano_DOWN.txt"),header=T,sep="\t")
  ALL <- rbind(UP,DOWN)
  ALL <- ALL%>%arrange(padj)
  ALL$Gene_name <- str_split(ALL$Gene, "\\.",simplify = T)[,1]
  ALL <- ALL[,c('Gene',"Gene_name","log2FoldChange",'padj')]

  SHEET=paste(toupper(substr(TS,1,3)),dat2$Compare[i],sep='_')
  sheets<-eval(parse(text = paste0('list(', SHEET,' = ALL)')))
  outsheets <- c( outsheets,sheets)
}


write_xlsx(outsheets, "./Paper1_LALO/Table_DEG.xlsx")



