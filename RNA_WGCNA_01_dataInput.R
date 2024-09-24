#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

# Display the current working directory
getwd();
setwd('./WGCNA')

# Load the WGCNA package
library(WGCNA);
library(data.table)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
femData = fread("LALO_TPM_genecounts.txt",header=T,stringsAsFactors = T,sep='\t');
femData = fread("LALO_vsdata_genecounts.txt",header=T,stringsAsFactors = T,sep='\t');

# Take a quick look at what is in the data set:
dim(femData);
names(femData);

femData[1:4,1:4]
#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================
nsample=48
if (ncol(femData)== nsample+6){
  #featurecouns input
geneinfo=femData[, c(1:6)]
datExpr0 = as.data.frame(t(femData[, -c(1:6)]));
names(datExpr0) = femData$Geneid;
rownames(datExpr0) = names(femData)[-c(1:6)];
write.csv(geneinfo,'LALO_tpm_Gene_Annotation.csv',quote = F,row.names = F)
}else if (ncol(femData)==nsample+1){
  #vst input
  geneinfo=femData[,49]
datExpr0 = as.data.frame(t(femData[,-49]));
names(datExpr0) = femData$Geneid;
rownames(datExpr0) = names(femData)[-c(49)];
write.csv(geneinfo,'LALO_vst_Gene_Annotation.csv',quote = F)
}

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
# plot
pdf(file = "Plots/sampleClustering_revised.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# Plot a line to show the cut
abline(h = 40, col = "red");
dev.off()
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 0, minSize = 10)
table(clust)
# clust 1 or 0 contains the samples we want to keep.
keepSamples = (clust==0)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr) #so one sample was removed

cat(nGenes,nSamples)
#=====================================================================================
#
#  Code chunk 7
#Phenotype data
#=====================================================================================


traitData = read.table("./WGCNA/LALO_trait.txt", header = TRUE,stringsAsFactors = T)
dim(traitData)
names(traitData)

#Filter trait
#traitData<-traitData[c( 'ID',"Tissue","LH_Sub_Stage","Conditions","Year","cort0","Fat_score","Weight" )]

# remove columns that hold information we do not need.
# allTraits = traitData[, -c(31, 16)];
# allTraits = allTraits[, c(2, 11:36) ];
allTraits = traitData
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.

femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$ID);
datTraits = allTraits[traitRows, -1];#remove the ID 
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage();


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
library(dplyr)
datTraits0=datTraits
datTraits0 <- datTraits0 %>% 
  mutate(Tissue_num=ifelse(Tissue=='Gonad',1,
                           ifelse(Tissue=='Heart',2,
                                  ifelse(Tissue=='Hypothalamus',3,
                                         ifelse(Tissue=='Liver',4,0)))))
###mutate colums for tissues
names(datTraits0)
datTraits=datTraits0[,c(6,9,10,13,14,22)]
datTraits=mutate_all(datTraits, function(x) as.numeric(as.character(x)))#convert the whole table to numeric
traitColors = numbers2colors(datTraits, signed = FALSE,colors = RColorBrewer::brewer.pal(name='BuPu',n=9));
#traitColors = numbers2colors(datTraits, signed = FALSE,colors = lisa::lisa_palette("KatsushikaHokusai", 1000, "continuous"));
traitColors = numbers2colors(datTraits, signed = FALSE,colors =rev(vangogh::vangogh_palette("Landscape",100,"continuous")))
# Plot the sample dendrogram and the colors underneath.
pdf(file = "Plots/vst_sampleClustering_pheno_revised.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


save(datExpr, datTraits, file = "MyLALO-revised-01-dataInput.RData")

                     
                     #END#
