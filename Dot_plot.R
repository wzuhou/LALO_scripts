#!/usr/bin/env Rscript
## Make Dot Plot with Percent Divergence on color scale
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
library('wesanderson')

{setwd('./LALO_genome/Genome_quality/DOT_plot/')
  h_lines=TRUE# TRUE#
  keep_ref=30
  similarity=TRUE #"show % identity (-s): "
  on_target=TRUE #"show % identity for on-target alignments only (-t): "
  interactive}

#a=paste0(as.character(1:30),collapse=',')
a='1,1A,2,3,4,4A,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30' #ZF1
a='Z' #ZF2
a= "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30" #Gal

{  setwd('./LALO_genome/Genome_quality/DOT_plot/')
    input_filename='fixed_ZF_i30_m.filter.coords' 
    output_filename='fixed_ZF_i30_m.filter.coords_1_22'
    min_query_aln=350000 #Gal: 30000# sex: 100000 #ZF: 300000 
    min_align= 14500 #Gal: 1000#sex: 1000 #ZF: 20000 #ZF(22):15600
    refIDs=a #paste0(as.character(1:30),collapse=',') #NULL #cat(a)
    }
}

subset(alignments,refID=='30')
#Before plot
{ 
# read in alignments
alignments = read.table(input_filename, stringsAsFactors = F, skip = 5)
alignments = alignments[,-c(3,6,9,11,14)]

# set column names
colnames(alignments) = c("refStart","refEnd","queryStart","queryEnd","lenAlnRef","lenAlnQuery","percentID","percentAlnCovRef","percentAlnCovQuery","refID","queryID")

cat(paste0("\nNumber of alignments: ", nrow(alignments),"\n"))
cat(paste0("Number of query sequences: ", length(unique(alignments$queryID)),"\n"))

#Modified the scaffold name
alignments$queryID2 <- gsub('Scaffold_Scaffold_', '', alignments$queryID)
alignments$queryID2 <- gsub('Scaffold_', '', alignments$queryID2)
alignments$queryID2 <- as.factor(as.character(alignments$queryID2))
alignments$queryID <-alignments$queryID2
sort(as.numeric(unique(alignments$refID)))

# sort by ref chromosome sizes, keep top X chromosomes OR keep specified IDs
if(is.null(refIDs)){
  chromMax = tapply(alignments$refEnd, alignments$refID, max)
  if(is.null(keep_ref)){
    keep_ref = length(chromMax)
  }
  refIDsToKeepOrdered = names(sort(chromMax, decreasing = T)[1:keep_ref])
  alignments = alignments[which(alignments$refID %in% refIDsToKeepOrdered),]
  
} else {
  refIDsToKeepOrdered = unlist(strsplit(refIDs, ","))
  alignments = alignments[which(alignments$refID %in% refIDsToKeepOrdered),]
}

# filter queries by alignment length, for now include overlapping intervals
queryLenAgg = tapply(alignments$lenAlnQuery, alignments$queryID, sum)
alignments = alignments[which(alignments$queryID %in% names(queryLenAgg)[which(queryLenAgg > min_query_aln)]),]

# filter alignment by length
alignments = alignments[which(alignments$lenAlnQuery > min_align),]

# re-filter queries by alignment length, for now include overlapping intervals
queryLenAgg = tapply(alignments$lenAlnQuery, alignments$queryID, sum)
alignments = alignments[which(alignments$queryID %in% names(queryLenAgg)[which(queryLenAgg > min_query_aln)]),]

cat(paste0("\nAfter filtering... Number of alignments: ", nrow(alignments),"\n"))
cat(paste0("After filtering... Number of query sequences: ", length(unique(alignments$queryID)),"\n\n"))

# sort df on ref
alignments$refID = factor(alignments$refID, levels = refIDsToKeepOrdered) # set order of refID
alignments = alignments[with(alignments,order(refID,refStart)),]
chromMax = tapply(alignments$refEnd, alignments$refID, max)

# # sort df on query
# alignments$queryID = factor(alignments$queryID, levels = refIDsToKeepOrdered) # set order of refID
# alignments = alignments[with(alignments,order(queryID,refStart)),]
# chromMax = tapply(alignments$refEnd, alignments$queryID, max)

# make new ref alignments for dot plot
if(length(levels(alignments$refID)) > 1){
  alignments$refStart2 = alignments$refStart + sapply(as.character(alignments$refID), function(x) ifelse(x == names((chromMax))[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax)) - 1]) )
  alignments$refEnd2 = alignments$refEnd +     sapply(as.character(alignments$refID), function(x) ifelse(x == names((chromMax))[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax)) - 1]) )
} else {
  alignments$refStart2 = alignments$refStart
  alignments$refEnd2 = alignments$refEnd
  
}

## queryID sorting step 1/2
# sort levels of factor 'queryID' based on longest alignment
alignments$queryID = factor(alignments$queryID, levels=unique(as.character(alignments$queryID))) 
queryMaxAlnIndex = tapply(alignments$lenAlnQuery,
                          alignments$queryID,
                          which.max,
                          simplify = F)
alignments$queryID = factor(alignments$queryID, levels = unique(as.character(alignments$queryID))[order(mapply(
  function(x, i)
    alignments$refStart2[which(i == alignments$queryID)][x],
  queryMaxAlnIndex,
  names(queryMaxAlnIndex)
))])

## queryID sorting step 2/2
## sort levels of factor 'queryID' based on longest aggregrate alignmentst to refID's
# per query ID, get aggregrate alignment length to each refID 
queryLenAggPerRef = sapply((levels(alignments$queryID)), function(x) tapply(alignments$lenAlnQuery[which(alignments$queryID == x)], alignments$refID[which(alignments$queryID == x)], sum) )
if(length(levels(alignments$refID)) > 1){
  queryID_Ref = apply(queryLenAggPerRef, 2, function(x) rownames(queryLenAggPerRef)[which.max(x)])
} else {queryID_Ref = sapply(queryLenAggPerRef, function(x) names(queryLenAggPerRef)[which.max(x)])}
# set order for queryID
alignments$queryID = factor(alignments$queryID, levels = (levels(alignments$queryID))[order(match(queryID_Ref, levels(alignments$refID)))])

#  flip query starts stops to forward if most align are in reverse complement
queryRevComp = tapply(alignments$queryEnd - alignments$queryStart, alignments$queryID, function(x) sum(x)) < 0
queryRevComp = names(queryRevComp)[which(queryRevComp)]
queryMax = tapply(c(alignments$queryEnd, alignments$queryStart), c(alignments$queryID,alignments$queryID), max)
names(queryMax) = levels(alignments$queryID)
alignments$queryStart[which(alignments$queryID %in% queryRevComp)] = queryMax[match(as.character(alignments$queryID[which(alignments$queryID %in% queryRevComp)]), names(queryMax))] - alignments$queryStart[which(alignments$queryID %in% queryRevComp)] + 1
alignments$queryEnd[which(alignments$queryID %in% queryRevComp)] = queryMax[match(as.character(alignments$queryID[which(alignments$queryID %in% queryRevComp)]), names(queryMax))] - alignments$queryEnd[which(alignments$queryID %in% queryRevComp)] + 1

## make new query alignments for dot plot
# subtract queryStart and Ends by the minimum alignment coordinate + 1
queryMin = tapply(c(alignments$queryEnd, alignments$queryStart), c(alignments$queryID,alignments$queryID), min)
names(queryMin) = levels(alignments$queryID)
alignments$queryStart = as.numeric(alignments$queryStart - queryMin[match(as.character(alignments$queryID),names(queryMin))] + 1)
alignments$queryEnd = as.numeric(alignments$queryEnd - queryMin[match(as.character(alignments$queryID),names(queryMin))] + 1)

queryMax = tapply(c(alignments$queryEnd, alignments$queryStart), c(alignments$queryID,alignments$queryID), max)
names(queryMax) = levels(alignments$queryID)
alignments$queryStart2 = alignments$queryStart + sapply(as.character(alignments$queryID), function(x) ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1]) )
alignments$queryEnd2 = alignments$queryEnd + sapply(as.character(alignments$queryID), function(x) ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1]) )

# get mean percent ID per contig
#   calc percent ID based on on-target alignments only
if(on_target & length(levels(alignments$refID)) > 1){
  alignments$queryTarget = queryID_Ref[match(as.character(alignments$queryID), names(queryID_Ref))]
  alignmentsOnTarget = alignments[which(as.character(alignments$refID) == alignments$queryTarget),]
  scaffoldIDmean = tapply(alignmentsOnTarget$percentID, alignmentsOnTarget$queryID, mean)
  alignments$percentIDmean = as.numeric(scaffoldIDmean[match(as.character(alignments$queryID), names(scaffoldIDmean))])
  alignments$percentIDmean[which(as.character(alignments$refID) != alignments$queryTarget)] = NA
} else{
  scaffoldIDmean = tapply(alignments$percentID, alignments$queryID, mean)
  alignments$percentIDmean = as.numeric(scaffoldIDmean[match(as.character(alignments$queryID), names(scaffoldIDmean))])
}
}

#My: change the scaffold id to numbers
# alignments$queryID2 <- gsub('Scaffold_Scaffold_', '', alignments$queryID)
# alignments$queryID2 <- gsub('Scaffold_', '', alignments$queryID2)
# alignments$queryID2 <- as.factor(as.character(alignments$queryID2))

###############
#####plot
yTickMarks = tapply(alignments$queryEnd2, alignments$queryID, max)
options(warn = -1) # turn off warnings
#Pal=wes_palette("Zissou1",n=100, type = "seq")
if (similarity) {
  gp = ggplot(alignments) +
    geom_point(
      mapping = aes(x = refStart2, y = queryStart2, color = percentIDmean),
      size = 0.009
    ) +
    geom_point(
      mapping = aes(x = refEnd2, y = queryEnd2, color = percentIDmean),
      size = 0.009
    ) +
    geom_segment(
      aes(
        x = refStart2,
        xend = refEnd2,
        y = queryStart2,
        yend = queryEnd2,
        color = percentIDmean,
        #size=5, lineend = "butt",
        text = sprintf(
          'Query ID: %s<br>Query Start Pos: %s<br>Query End Pos: %s<br>Target ID: %s<br>Target Start Pos: %s<br>Target End Pos: %s<br>Length: %s kb',
          queryID,
          queryStart,
          queryEnd,
          refID,
          refStart,
          refEnd,
          round(lenAlnQuery / 1000, 1)
        )
      )
    ) +
    scale_x_continuous(breaks = cumsum(as.numeric(chromMax)), expand = c(0, 0),guide = guide_axis(n.dodge=3),
                       labels = levels(alignments$refID)) +
    theme_bw() +
    theme(text = element_text(size = 8)) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title = element_text(size = 20),
      axis.text.x = element_text(size = 5),
      axis.text.y = element_text(size = 10),#, angle = 15
      legend.key.size = unit(0.3, 'cm'), #change legend key size
      legend.key.height = unit(1, 'cm'), #change legend key height
      legend.key.width = unit(0.3, 'cm'), #change legend key width
      legend.title = element_text(size=5), #change legend title font size
      legend.text = element_text(size=5) #change legend text font size
    ) +
    scale_y_continuous(breaks = yTickMarks, guide = guide_axis(n.dodge=3),expand = c(0, 0), labels = substr(levels(alignments$queryID), start = 1, stop = 20)) +
    { if(h_lines){ geom_hline(yintercept = yTickMarks,
                                  color = "grey60",
                                  size = .1) }} +
    
    #scale_color_distiller(palette = "Spectral") + #,limits = c(85, 99)
    #scale_color_distiller(palette = 'GnBu',direction = 1)+
    scale_colour_viridis_c(direction = -1)+
    labs(color = "Mean Identity%")+ #, 
         # title = paste0(   paste0("Post-filtering number of alignments: ", nrow(alignments),"\t\t\t\t"),
         #                   paste0("minimum alignment length (-m): ", min_align,"\n"),
         #                   paste0("Post-filtering number of queries: ", length(unique(alignments$queryID)),"\t\t\t\t\t\t\t\t"),
         #                   paste0("minimum query aggregate alignment length (-q): ", min_query_aln)
         # )) +
    xlab("T. guttata") +#"T. guttata" #"G. gallus"
    ylab("C. lapponicus")
} 

#}
gp
ggsave(filename = paste0(output_filename, "_final.png"), width = 7, height = 5, units = "in", dpi = 550, limitsize = F)
unique(alignments$queryID)

cat(paste0(   paste0("Post-filtering number of alignments: ", nrow(alignments),"\t"),
              paste0("minimum alignment length (-m): ", min_align,"\n"),
              paste0("Post-filtering number of queries: ", length(unique(alignments$queryID)),"\t"),
              paste0("minimum query aggregate alignment length (-q): ", min_query_aln)),
    paste0('LALO scaffold number: ', sort(unique(alignments$queryID))),sep='\n'
)
sort(as.numeric(as.character(unique(alignments$queryID))))

###END###
