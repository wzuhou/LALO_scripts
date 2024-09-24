#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

# Display the current working directory
getwd();
setwd('./WGCNA')
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
# workingDir = ".";
# setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "MyLALO_revised-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "LALO_vst_revised-02-networkConstruction-auto.RData"); 
lnames
DATA='vst'
set.seed(10)

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0) 
moduleTraitCor = cor(MEs, datTraits, use = "p");##Modulecolor x pheno
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#moduleTraitPvalue[row.names(moduleTraitPvalue)=='MElightyellow','Tissue_num'] #check a specific correlation level
#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================



# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

sizeGrWindow(10,6)
tiff(paste0('Plots/',DATA,'_Module-trait relationships_revised','.tiff'),15,10,units = "in",res=350)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed (50),
               #colors =colorRampPalette(brewer.pal(11,'RdBu')),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

################################################


PHENO='cort0'
# names(datTraits)==PHENO
phenotype_measurement  = as.data.frame(datTraits[,names(datTraits)==PHENO]);
names(phenotype_measurement) = paste0(PHENO)
y=phenotype_measurement 
GS1=as.numeric(cor(y,datExpr, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance,
                          moduleColors, mean, na.rm=T)

# Plot it
sizeGrWindow(8,7)
pdf(paste0('Plots/',DATA,'_',PHENO,'_Gene_sig_across_module_revised.pdf'),8,7)
par(mfrow = c(1,1),cex=0.6)
plotModuleSignificance(GeneSignificance,moduleColors) #y is average value of gene significant
dev.off()


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Define variable weight containing the weight column of datTrait
phenotype_measurement = as.data.frame(datTraits$cort0);
names(phenotype_measurement) = "cort0"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
datExpr[1:4,1:4]
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, phenotype_measurement, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(phenotype_measurement), sep="");
names(GSPvalue) = paste("p.GS.", names(phenotype_measurement), sep="");

GSPvalue[row.names(GSPvalue)=='FMN2.2',]
GSPvalue[row.names(GSPvalue)=='FKBP5',]
geneModuleMembership[row.names(geneModuleMembership)=='FKBP5',]##in general all very low

geneTraitSignificance_F = as.data.frame(cor(datExpr[names(datExpr)=='FKBP5'], phenotype_measurement, use = "p"));
GSPvalue_F = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
library(ggplot2)
library(dplyr)
# Add the regression line
gene_Trait_plot <- data.frame(FKBP5=datExpr[names(datExpr)=='FKBP5'],
  Trait=datTraits$cort0,
  Tissue_num=datTraits$Tissue_num
  )

gene_Trait_plot <- gene_Trait_plot %>% 
  mutate(Tissue=ifelse(Tissue_num==1,'Gonad',
                       ifelse(Tissue_num==2,'Heart',
                              ifelse(Tissue_num==3,'Hypothalamus',
                                     ifelse(Tissue_num==4,'Liver',0)))))
gene_Trait_plot$Tissue=as.factor(gene_Trait_plot$Tissue)
p1 <- ggplot(gene_Trait_plot, aes(x=Trait, y=FKBP5,color=Tissue)) + 
  geom_point(alpha=0.4)+
  geom_smooth(method=glm, se=F,fullrange=TRUE)+# Remove the confidence interval , se=FALSE
  #scale_fill_manual(values = alpha(c("#00AFBB", "#E7B800", "#FC4E07",'darkblue'),0.1))+
  scale_color_manual(values = alpha(c("#00AFBB", "#E7B800", "#FC4E07",'darkblue'),0.7))+
  theme_classic()
p1
ggsave('Plots/FKBP5_cor_cort.tiff',width = 5,height = 1.8,dpi=350,p1)

#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


names(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================

##list genes in specific modules
COLOR='black'
mg <- stringr::str_split(names(datExpr)[moduleColors==COLOR],'[[.]]',simplify = T)[,1]

write.table(mg,paste0('./WGCNA/Plots/pheno_cor/list_gene_',COLOR,'.txt'),row.names=F,quote = F,col.names = F)
table(moduleColors)
#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================

DATA='vst'###?
annot = read.csv(file = paste0("LALO_",DATA,"_Gene_Annotation.csv"));
dim(annot);
dim(datExpr)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$Geneid)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================

# Create the starting data frame
geneInfo0 = data.frame(Geneid = probes,
                      geneSymbol = annot$Geneid[probes2annot],
                      #LocusLinkID = annot$LocusLinkID[probes2annot],
                      moduleColor = moduleColors,
                      geneTraitSignificance,
                      GSPvalue)
row.names(geneInfo0)
row.names(geneInfo0)<-probes
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, phenotype_measurement, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.cort0));
geneInfo = geneInfo0[geneOrder, ]

FKBP5 <- geneInfo[geneInfo$Geneid =='FKBP5',]
geneInfo[geneInfo$Geneid =='FKBP5',1:5]
FKBP5[,stringr::str_detect(names(FKBP5),paste0("^\\Q", "p.", "\\E")) & FKBP5[1,]<= 0.05]
FKBP5[,stringr::str_detect(names(FKBP5),paste0("^\\Q", "p.", "\\E"))]
#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================

DATA
write.csv(geneInfo, file = paste0("LALO_",DATA,"_geneInfo_pheno_revised.csv"))

geneOrder2 = order(geneInfo0$p.GS.cort0,geneInfo0$moduleColor,decreasing = F);
geneInfo2 = geneInfo0[geneOrder2, ]
geneInfo3 <- geneInfo2[geneInfo2$p.GS.cort0<0.05,1:5]
write.csv(geneInfo3, file = paste0("LALO_",DATA,"_geneInfo_pheno_cort0.05.csv"))
#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

#Select modules
module = "black"
module='steelblue'
column = match(module, modNames);
moduleGenes = moduleColors==module;
PHENO
pdf(paste0('Plots/',DATA,'_',PHENO,'_',module,'_cor_Module_gene.pdf'),5,3)
#sizeGrWindow(6,6);
par(mfrow = c(1,1),cex=0.6);
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), #correlation coef with ME
                   abs(geneTraitSignificance[moduleGenes, 1]),#correlation coef
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for pheno",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1, cex.lab = 1, cex.axis = 1, col = module)

dev.off()

moduleTraitCor[row.names(moduleTraitCor)==paste0('ME',module),PHENO]
moduleTraitPvalue[row.names(moduleTraitPvalue)==paste0('ME',module),PHENO]

geneInfo[geneInfo$Geneid =='PENK',1:5]
geneInfo[geneInfo$Geneid =='FAM105B',1:5]
module='blue'
PHENO='cort0'

