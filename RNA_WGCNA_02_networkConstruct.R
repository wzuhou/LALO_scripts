#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
 
setwd('./WGCNA')
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "MyLALO-revised-01-dataInput.RData" );#"MyLALO-01-dataInput.RData"
#The variable lnames contains the names of loaded variables.
lnames
set.seed(15)

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:13), seq(from = 14, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType='signed')
sft$powerEstimate

# Plot the results:
pdf('Plots/soft_thresholding_power_vst_signed_revised.pdf',9,5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# We use the following power for the power adjacency function.
#beta1=12
beta1=11
Connectivity=softConnectivity(datExpr,power=beta1)-1
par(mfrow=c(1,1))
scaleFreePlot(Connectivity, main=paste("soft threshold, power=",beta1), truncated=F);
#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


net = blockwiseModules(datExpr, power = 11, #sft$powerEstimate, or 5 or 13
                       networkType='signed',
                       #TOMType = "signed", 
                       minModuleSize = 30,
                       #maxBlockSize = 6000,
                       #reassignThreshold = 0, 
                       mergeCutHeight = 0.15,
                       deepSplit = 3, #how sensitive module detection should be to module splitting
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       #saveTOMFileBase = "LALO_vst_TOM", 
                       verbose = 3)

#To see how many modules were identified and what the module sizes are
table(net$colors) 

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

pdf('Plots/vst_dendrogram_module_color_revised.pdf',12,9)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
unique(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], 
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    addGuide = TRUE,
                    dendroLabels = FALSE,
                    hang = 0.03,
                    guideHang = 0.05
                    )
dev.off()
## assign all of the gene to their corresponding module 
## hclust for the genes.

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "LALO_vst_revised-02-networkConstruction-auto.RData")


#END#
