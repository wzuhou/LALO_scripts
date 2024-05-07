#R-volcano plot for RNAseq results
{library(stringr)
  library(ggplot2)
  library(ggrepel)}

DEseq_DEG=DEseq_DEG
n=paste0(prefix)
#Input
#exprSet,DEseq_DEG,n='DEseq2'
## two columns of DEG needed, log2FoldChange and pvalue

#######################
#####Cut-off threshold
#logFC_cutoff <- with(DEseq_DEG,mean(abs( log2FoldChange)) + 2*sd(abs(log2FoldChange)) ) #logFC_cutoff=0.58
DEseq_DEG$change = as.factor(ifelse(DEseq_DEG$padj < 0.05 & abs(DEseq_DEG$log2FoldChange) > logFC_cutoff,
                                   ifelse(DEseq_DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT'))
table(DEseq_DEG$change)

#####Change gene name suffix
DEseq_DEG$Gene <- str_split(row.names(DEseq_DEG), "[.]",simplify = T)[,1] 
#Or
#DEseq_DEG$Gene <- row.names(DEseq_DEG)

write.table(row.names(subset(DEseq_DEG,change=="UP")),
            paste0("GENE_",n,"_Volcano_UP.txt"),row.names=F,quote=F,sep="\t",col.names=F)
write.table(row.names(subset(DEseq_DEG,change=="DOWN")),
            paste0("GENE_",n,"_Volcano_DOWN.txt"),row.names=F,quote=F,sep="\t",col.names=F)
write.table(row.names(subset(DEseq_DEG,change=="UP" |change=="DOWN")),
            paste0("GENE_",n,"_Volcano_ALL.txt"),row.names=F,quote=F,sep="\t",col.names=F)

this_tile <- paste0(n ,
                     # '\nCutoff for logFC is ',round(logFC_cutoff,3),
                     '\nUP gene: ',nrow(DEseq_DEG[DEseq_DEG$change =='UP',]) ,
                    '\nDOWN gene: ',nrow(DEseq_DEG[DEseq_DEG$change =='DOWN',]))

###Optional: highlight & label
# Add highlight and annotation information
DEseq_DEG<-DEseq_DEG %>%
  mutate( is_annotate=ifelse(Gene %in% GeneOfInterest &row_number()<=10, "yes", "no")) #if in the list and top 10 rows
#&str_detect(Gene, "^(?!MCM[0-9]*)")
table(DEseq_DEG$is_annotate)
#DEseq_DEG%>%  filter(is_annotate =='yes')

############################### 
###Ploting
paste0('g',i,j,k)
assign(paste0('g',i,j,k),
       ggplot(data=DEseq_DEG, aes(x=log2FoldChange, y=-log10(padj), color=change)) +
         geom_point(alpha=0.8, size=1.8) +
         ####################OPTIONAL
         # Add highlighted points
         #geom_point(data=subset(DEseq_DEG, is_highlight=="yes"), color="firebrick4", size=2,alpha=0.8) +
         # Add label using ggrepel to avoid overlapping
         geom_label_repel(data=subset(DEseq_DEG, is_annotate=="yes"), aes(label=Gene), size=2.5,max.overlaps =15 ) +
         ####################OPTIONAL
         theme_set(theme_set(theme_classic(base_size=15)))+
         #xlab("log2 fold change") + ylab("-log10 FDR P-value") +
         ggtitle( this_tile ) +
         scale_colour_manual(values = c('#3B9AB2','black','#EBCC2A'))+
         theme_classic()+
         expand_limits(x = 0, y = 0)+
         theme(
           #text = element_text(size = 22),
           plot.title = element_text(size=10,hjust = 0.5), #hjust adjust the position with 1 aligning to right
           axis.title.x=element_blank(),axis.title.y=element_blank(),
           legend.position = "none")
)
# ##when you need to made changes to the figure
 for (j in seq(1,4)){
  for (i in seq(1,3)){
     print(paste0('g',i,j))
 assign(paste0('g',i,j),
       eval(parse(text=paste0('g',i,j)))+  
  theme(axis.title.x=element_blank(),axis.title.y=element_blank()))
}}

#eval(parse(text=paste0('g',i,j,k)))
#save grid arrange
ggsave(eval(parse(text=paste0('g',i,j,k))) ,filename = paste0("DEseq2_",n,'_volcano.png'),width = 4,height = 4)

###END###


