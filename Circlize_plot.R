#######draw_v_U_D(exprSet,DEseq_DEG,paste0(prefix))
{library(stringr)
  library(ggplot2)
  library(ggrepel)
}

#need_DEG=DEseq_DEG
#n=paste0(prefix)

draw_v_U_D <- function(exprSet,need_DEG,n='DEseq2'){
## we only need two columns of DEG, which are log2FoldChange and pvalue
# ## heatmap
# library(pheatmap)
# choose_gene=head(rownames(need_DEG),50) ## 50 maybe better
# choose_matrix=exprSet[choose_gene,]
# choose_matrix=t(scale(t(choose_matrix)))
# pheatmap(choose_matrix,filename = paste0(n,'_need_DEG_top50_heatmap.png'))
#logFC_cutoff <- with(need_DEG,mean(abs( log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
logFC_cutoff=0.58

need_DEG$change = as.factor(ifelse(need_DEG$padj < 0.05 & abs(need_DEG$log2FoldChange) > logFC_cutoff,
                                   ifelse(need_DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
)
#table(need_DEG$change)


need_DEG$Gene <- str_split(row.names(need_DEG), "[.]",simplify = T)[,1]#####Change gene name
#need_DEG$Gene <- row.names(need_DEG)

write.table(subset(need_DEG,change=="UP"),
            paste0(n,"_Volcano_UP.txt"),row.names=F,quote=F,sep="\t")
write.table(subset(need_DEG,change=="DOWN"),
            paste0(n,"_Volcano_DOWN.txt"),row.names=F,quote=F,sep="\t")
write.table(row.names(subset(need_DEG,change=="UP")),
            paste0("GENE_",n,"_Volcano_UP.txt"),row.names=F,quote=F,sep="\t",col.names=F)
write.table(row.names(subset(need_DEG,change=="DOWN")),
            paste0("GENE_",n,"_Volcano_DOWN.txt"),row.names=F,quote=F,sep="\t",col.names=F)
write.table(row.names(subset(need_DEG,change=="UP" |change=="DOWN")),
            paste0("GENE_",n,"_Volcano_ALL.txt"),row.names=F,quote=F,sep="\t",col.names=F)

this_tile <- paste0(n ,
                     # '\nCutoff for logFC is ',round(logFC_cutoff,3),
                     '\nUP gene: ',nrow(need_DEG[need_DEG$change =='UP',]) ,
                    '\nDOWN gene: ',nrow(need_DEG[need_DEG$change =='DOWN',])
)

#########Optional: highlight & label############
# Add highlight and annotation information
need_DEG<-need_DEG %>%
  #mutate( is_highlight=ifelse( -log10(padj)>max(-log10(padj))*0.95 ,"yes", "no")) %>%
  #mutate( is_annotate=ifelse(Gene %in% GeneOfInterest &row_number()<=5, "yes", "no")) #if in the list and top 10 rows
  #mutate( is_annotate=ifelse(Gene %in% GeneOfInterest&-log10(padj)>max(-log10(padj))*0.95 , "yes", "no")) #if in the list and top 5% Pval
  mutate( is_annotate=ifelse(Gene %in% GeneOfInterest & change != 'NOT' & str_detect(Gene, "^(?!LOC[0-9]*)") & row_number()<=20, "yes", "no")) 
#liver 43  & rownames(need_DEG)!='NPAS2.1' 
#&str_detect(Gene, "^(?!MCM[0-9]*)")
# mutate( is_annotate=ifelse(row_number()<=5 & change != 'NOT' & str_detect(Gene, "^(?!LOC[0-9]*)") &str_detect(Gene, "^(?!Scaffold.*)") | Gene=='FKBP5' , "yes", "no")) #if top 5
#mutate( is_annotate=ifelse(-log10(padj)>max(-log10(padj))*0.95 , "yes", "no")) #if pval top 0.05

#table(need_DEG$is_annotate)
#need_DEG%>%  filter(is_annotate =='yes')

#########Making plots#####################
paste0('g',i,j,k)
assign(paste0('g',i,j,k),
       ggplot(data=need_DEG, aes(x=log2FoldChange, y=-log10(padj), color=change)) +
         geom_point(alpha=0.8, size=1.8) +
         ####################OPTIONAL
         # Add highlighted points
         #geom_point(data=subset(need_DEG, is_highlight=="yes"), color="firebrick4", size=2,alpha=0.8) +
         # Add label using ggrepel to avoid overlapping
         geom_label_repel(data=subset(need_DEG, is_annotate=="yes"), aes(label=Gene), size=2.5,max.overlaps =15 ) +
         ####################OPTIONAL
         theme_set(theme_set(theme_classic(base_size=15)))+
         #xlab("log2 fold change") + ylab("-log10 FDR P-value") +
         ggtitle( this_tile ) +
         #scale_colour_manual(values = c('dodgerblue4','black','firebrick4')) ## corresponding to the levels(res$change)+
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
# for (j in seq(1,4)){
#   for (i in seq(1,3)){
#     print(paste0('g',i,j))
# assign(paste0('g',i,j),
#        eval(parse(text=paste0('g',i,j)))+  
#   theme(axis.title.x=element_blank(),axis.title.y=element_blank()))
# }}

#eval(parse(text=paste0('g',i,j,k)))
ggsave(eval(parse(text=paste0('g',i,j,k))) ,filename = paste0("DEseq2_",n,'_volcano.png'),width = 4,height = 4)
########
#save grid arrange
# my_list<-vector()
# for (k in c(1,3,2)){
#     for (l in seq(1:4)){
#     my_list<-c(my_list,paste0('g',k,l))
#   }
# }
# my_list=c(my_list)
# glist <-paste0(my_list,collapse=',')
# library('gridExtra')
# grid.arrange(g11,g12,g13,g14,g31,g32,g33,g34,g21,g22,g23,g24, ncol = 4) ## display plot
#ggsave(filename = "C:/Users/zwu33/Downloads/ROSLIN/WCS/WCS_RNA/RNA_Seq_LALO/DEGs_all_label.png", arrangeGrob(g11,g12,g13,g14,g31,g32,g33,g34,g21,g22,g23,g24,ncol = 4),dpi = 600,width = 8,height = 8)  ## save plot

#save each
#for (g in my_list){
#ggsave(eval(parse(text=paste0(g))),filename = paste0("C:/Users/zwu33/Downloads/ROSLIN/RNA_Seq_LALO/DEGs_all_label_fig/",g,".png"),width = 4,height = 4)
#}

#plot_grid(g11,g12,g13,g14,g31,g32,g33,g34,g21,g22,g23,g24, align = "v", ncol = 4)

# library(patchwork)
# (plots = wrap_plots(g11,g12,g13,g14,g31,g32,g33,g34,g21,g22,g23,g24, ncol = 4)) +
#    # plot_annotation(title = "log2 fold change",
#     #                subtitle = "-log10 FDR P-value") &
#     theme(text = element_text(size = 5),
#           axis.title.x=element_blank())
          #plot.title = element_text(vjust = -110, hjust = 0.50),
          #plot.subtitle = element_text(vjust = -55, hjust = -0.01, size = 12))
ggsave(eval(parse(text=paste0('g',i,j,k))) ,filename = paste0("DEseq2_",n,'_volcano.png'),width = 4,height = 4)
#ggsave(g,filename = paste0("Selected_DEseq2_Vol_",n,'_volcano_Label_Highlight.png'),width = 4,height = 4)
#ggsave(g,filename = 'M:/ROSLIN/RNA_Seq_LALO/Legend.png',width = 4.5,height = 4.5)
}

