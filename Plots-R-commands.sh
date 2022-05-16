### Created by E.Zarza in April - May 2022 to analyse ITS data for manuscript: Fungal diversity in shade-coffee plantations in Soconusco, Mexico
### Eugenia Zarza, Alejandra López-Pastrana, Anne Damon, Griselda Karina Guillén-Navarro, Luz Verónica García-Fajardo.
### El Colegio de la Frontera Sur - Tapachula - Mexico
### Analyses were executed on Ubuntu 18. It requires R v4,and other libraries indicated bellow.
### This is not intended to run as script. Commands should be executed one by one.

#Create a conda environment and install R v.4 and other libraries

conda create --name R4 python=3.8
conda activate R4
conda install -c conda-forge r-rcpp
conda install -c bioconda bioconductor-dada2
conda install -c bioconda bioconductor-phyloseq
conda install -c conda-forge r-phangorn
conda install -c conda-forge r-tidyverse

#open R
R
#R version 4.0.5 (2021-03-31) -- "Shake and Throw"

#and  installed ampvis with:
install.packages("remotes")
remotes::install_github("MadsAlbertsen/ampvis2")

#load package
library(ampvis2)

#set working directory
setwd("/home/eugenia/R-fungi-peerJ")

#import otu table for Fungi, and metadata
 fungi <- amp_load(
  otutable = "Only_Fungi_OTU_table.tsv",
  metadata = "Fungi-assigned-metadata.csv"
)
fungi
#ampvis2 object with 3 elements.
#Summary of OTU table:
#     Samples         OTUs  Total#Reads    Min#Reads    Max#Reads Median#Reads
#          48         3783      2260671         2402       126471      39191.5
#   Avg#Reads
#    47097.31
#
#Assigned taxonomy:
#     Kingdom       Phylum        Class        Order       Family        Genus
#  3783(100%) 2977(78.69%) 2552(67.46%) 2352(62.17%) 2016(53.29%) 1879(49.67%)
#     Species
#1632(43.14%)
#
#Metadata variables: 5
# SampleID, Type, Site, Group1, Group2

#Keep only samples with more than 10,000 reads
fungi10K <- amp_subset_samples(
  fungi,
  minreads = 10000
)
#3 samples and 78 OTUs have been filtered
#Before: 48 samples and 3783 OTUs
#After: 45 samples and 3705 OTUs

#estimate alpha diversity
alphadiversityresult <- amp_alphadiv(fungi10K,
 measure = c("shannon", "simpson"),
 rarefy = 20000
)
#The following sample(s) have not been rarefied (less than 20000 reads):
#BOP22, BOP43, BOP109, BOP121, HOP32, HOP85, HOP100, HOP106

#export result with rarefaction at 20K Reads
 write.table(alphadiversityresult, "fungi10K_shannon_simpson.txt")

#extract shannon diversity column
alphadiversityresult[ , c("SampleID", "Shannon")]
shannon_result<-alphadiversityresult[ , c("SampleID", "Shannon")]
simpson_result<-alphadiversityresult[ , c("SampleID", "Simpson")]
write.table(shannon_result, "fungi10K_shannon_result.txt")
write.table(simpson_result, "fungi10K_simpson_result.txt")

#estimate alpha diversity with rarefaction at 12.408K, the lowest number of reads in samples with more than 10K reads
alphadiverare12K <- amp_alphadiv(fungi10K,
measure = c("shannon", "simpson"),
rarefy = 12408
)
#export result 12K
 write.table(alphadiverare12K, "fungi_rare12K_shannon_simpson.txt")

 #extract shannon diversity column 12K
 alphadiverare12K[ , c("SampleID", "Shannon")]
 shannon_result12K<-alphadiverare12K[ , c("SampleID", "Shannon")]
 simpson_result12K<-alphadiverare12K[ , c("SampleID", "Simpson")]
 write.table(shannon_result12K, "fungi_rare12K_shannon_result.txt")
 write.table(simpson_result12K, "fungi_rare12K_simpson_result.txt")

 shannon12K_result<-read.table("/home/eugenia/R-fungi-peerJ/fungi_rare12K_shannon_result.txt", header=TRUE, row.names=1)
#changed treatmen names because did not accept '-' in names
 metadata<-read.table("/home/eugenia/R-fungi-peerJ/FungiAssignedMetadata.csv", header=TRUE, row.names=1)

 ### Install amplicon with R console within the ampvis conda environment to create alpha diversity boxplots
install.packages("remotes")
remotes::install_github("microbiota/amplicon")
library(amplicon)
setwd("/home/eugenia/R-fungi-peerJ/")
shannon_result<-read.table("/home/eugenia/R-fungi-peerJ/fungi10K_shannon_result.txt", header=TRUE, row.names=1)
metadata<-read.table("/home/eugenia/R-fungi-peerJ/FungiAssignedMetadata.csv", header=TRUE, row.names=1)
#boxplots for alpha diversity results, rarefied at 20K, were produced with the following commands
#by Site
png(
"shannon_site2.png",
width     = 3.25,
height    = 3.25,
units     = "in",
res       = 1200,
pointsize = 4
)
par(
mar      = c(5, 5, 2, 2),
cex.axis = 2,
cex.lab  = 2
)
alpha_boxplot(shannon_result, metadata, "Shannon", "Site")
dev.off()
#by plant type
png("shannon_type.png", width = 3.25, height = 3.25, units= "in", res = 1200, pointsize = 4)
par(mar = c(5, 5, 2, 2),cex.axis = 2, cex.lab = 2)
alpha_boxplot(shannon_result, metadata, "Shannon", "Type")
 dev.off()

#by groups (Trunk high, trunk low, branch, twig)
png("shannon_group1.png", width = 3.25, height = 3.25, units= "in", res = 1200, pointsize = 4)
par(mar = c(5, 5, 2, 2),cex.axis = 2, cex.lab = 2)
alpha_boxplot(shannon_result, metadata, "Shannon", "Group1")
dev.off()
#by groups (plant+position)
png("shannon_group2.png", width = 3.25, height = 3.25, units= "in", res = 1200, pointsize = 4)
par(mar = c(5, 5, 2, 2),cex.axis = 2, cex.lab = 2)
alpha_boxplot(shannon_result, metadata, "Shannon", "Group2")
dev.off()

#now the same procedure but with simpson index
simpson_result<-read.table("/home/eugenia/R-fungi-peerJ/fungi10K_simpson_result.txt", header=TRUE, row.names=1)
#by Site
png("simpson_site.png", width = 3.25, height = 3.25, units = "in", res = 1200,pointsize = 4)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab  = 2)
alpha_boxplot(simpson_result, metadata, "Simpson", "Site")
dev.off()
#by plant type
png("simpson_type.png", width = 3.25, height = 3.25, units= "in", res = 1200, pointsize = 4)
par(mar = c(5, 5, 2, 2),cex.axis = 2, cex.lab = 2)
alpha_boxplot(simpson_result, metadata, "Simpson", "Type")
 dev.off()
#by groups (Trunk high, trunk low, branch, twig)
png("simpson_group1.png", width = 3.25, height = 3.25, units= "in", res = 1200, pointsize = 4)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
alpha_boxplot(simpson_result, metadata, "Simpson", "Group1")
dev.off()
#by groups (plant+position)
png("simpson_group2.png", width = 3.25, height = 3.25, units= "in", res = 1200, pointsize = 4)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
alpha_boxplot(simpson_result, metadata, "Simpson", "Group2")
dev.off()

#now import results from analysis with Shannon index rarefied at 12K
shannon12K_result<-read.table("/home/eugenia/R-fungi-peerJ/fungi_rare12K_shannon_result.txt", header=TRUE, row.names=1)

#boxplots for alpha diversity results, rarefied at 12K, were produced with the following commands
#by Site
png("shannon12K_site.png", width = 3.25, height = 3.25, units = "in", res = 1200, pointsize = 3)
par(mar = c(5, 5, 2, 2),cex.axis = 2, cex.lab  = 2)
alpha_boxplot(shannon12K_result, metadata, "Shannon", "Site")
dev.off()
#by plant type
png("shannon12K_type.png", width = 3.25, height = 3.25, units= "in", res = 1200, pointsize = 4)
par(mar = c(5, 5, 2, 2),cex.axis = 2, cex.lab = 2)
alpha_boxplot(shannon12K_result, metadata, "Shannon", "Type")
 dev.off()
#by groups (Trunk high, trunk low, branch, twig)
png("shannon12K_group1.png", width = 3.25, height = 3.25, units= "in", res = 1200, pointsize = 4)
par(mar = c(5, 5, 2, 2),cex.axis = 2, cex.lab = 2)
alpha_boxplot(shannon12K_result, metadata, "Shannon", "Group1")
dev.off()
#by groups (plant+position)
png("shannon12K_group2.png", width = 3.25, height = 3.25, units= "in", res = 1200, pointsize = 4)
par(mar = c(5, 5, 2, 2),cex.axis = 2, cex.lab = 2)
alpha_boxplot(shannon12K_result, metadata, "Shannon", "Group2")
dev.off()

#now import results from analysis with Simpson index rarefied at 12K
simpson12K_result<-read.table("/home/eugenia/R-fungi-peerJ/fungi_rare12K_simpson_result.txt", header=TRUE, row.names=1)

#boxplots for alpha diversity results, rarefied at 12K, were produced with the following commands
#by Site
png("simpson12K_site.png", width = 3.25, height = 3.25, units = "in", res = 1200, pointsize = 4)
par(mar = c(5, 5, 2, 2),cex.axis = 2, cex.lab  = 2)
alpha_boxplot(simpson12K_result, metadata, "Simpson", "Site")
dev.off()
#by plant type
png("simpson12K_type.png", width = 3.25, height = 3.25, units= "in", res = 1200, pointsize = 4)
par(mar = c(5, 5, 2, 2),cex.axis = 2, cex.lab = 2)
alpha_boxplot(simpson12K_result, metadata, "Simpson", "Type")
 dev.off()
#by groups (Trunk high, trunk low, branch, twig)
png("simpson12K_group1.png", width = 3.25, height = 3.25, units= "in", res = 1200, pointsize = 4)
par(mar = c(5, 5, 2, 2),cex.axis = 2, cex.lab = 2)
alpha_boxplot(simpson12K_result, metadata, "Simpson", "Group1")
dev.off()
#by groups (plant+position)
png("simpson12K_group2.png", width = 3.25, height = 3.25, units= "in", res = 1200, pointsize = 4)
par(mar = c(5, 5, 2, 2),cex.axis = 2, cex.lab = 2)
alpha_boxplot(simpson12K_result, metadata, "Simpson", "Group2")
dev.off()

#Ordination with library amplicon
setwd("/home/eugenia/R-fungi-peerJ/")
fungi_asvs<-read.table("/home/eugenia/R-fungi-peerJ/Fungi-assigned-OTU-table.csv", header=TRUE, fill = TRUE, row.names=1)
BetaDiv(otu=fungi_asvs, map=metadata, group="Site", dist="bray", method="NMDS", Micromet="adonis", pvalue.cutoff=0.05)
#output R in ampvis environment
#ID           stat        p
#1 B_VS_H adonis:R 0.036 p: 0.001

#[[5]]
#[1] "adonis:R 0.036 p: 0.001"

png("NMDS_bray_Sitev3.png",width = 4, height = 4, units= "in", res = 1200, pointsize = 1)
BetaDiv(otu=fungi_asvs, map=metadata, group="Site", dist="bray", method="NMDS", Micromet="adonis", pvalue.cutoff=0.05)
dev.off()
#ID           stat        p
#1 B_VS_H adonis:R 0.036 p: 0.001
#[1] "adonis:R 0.036 p: 0.001"

png("NMDS_bray_Type.png",width = 4, height = 4, units= "in", res = 1200, pointsize = 1)
BetaDiv(otu=fungi_asvs, map=metadata, group="Type", dist="bray", method="NMDS", Micromet="adonis", pvalue.cutoff=0.05)
dev.off()
#ID           stat        p
# Coffee_VS_Tree adonis:R 0.026 p: 0.002

png("NMDS_bray_Group1.png",width = 4, height = 4, units= "in", res = 1200, pointsize = 1)
BetaDiv(otu=fungi_asvs, map=metadata, group="Group1", dist="bray", method="NMDS", Micromet="adonis", pvalue.cutoff=0.05)
dev.off()
#ID           stat        p
#1      Twig_VS_TrunkLow adonis:R 0.055 p: 0.001
#2     Twig_VS_TrunkHigh adonis:R 0.052 p: 0.005
#3        Twig_VS_Branch adonis:R 0.048  p: 0.08
#4 TrunkLow_VS_TrunkHigh adonis:R 0.041 p: 0.707
#5    TrunkLow_VS_Branch adonis:R 0.046 p: 0.073
#6   TrunkHigh_VS_Branch adonis:R 0.042 p: 0.654

#[1] "adonis:R 0.069 p: 0.02"

png("NMDS_bray_Group2.png",width = 4, height = 4, units= "in", res = 1200, pointsize = 1)
BetaDiv(otu=fungi_asvs, map=metadata, group="Group2", dist="bray", method="NMDS", Micromet="adonis", pvalue.cutoff=0.05)
dev.off()
#ID           stat        p
#1       Coffee_Twig_VS_Coffee_TrunkLow  adonis:R 0.08 p: 0.001
#2             Coffee_Twig_VS_Tree_Twig adonis:R 0.103 p: 0.105
#3      Coffee_Twig_VS_Coffee_TrunkHigh  adonis:R 0.08 p: 0.012
#4           Coffee_Twig_VS_Tree_Branch adonis:R 0.103 p: 0.117
#5        Coffee_Twig_VS_Tree_TrunkHigh adonis:R 0.101 p: 0.135
#6         Coffee_Twig_VS_Tree_TrunkLow adonis:R 0.102 p: 0.094
#7         Coffee_Twig_VS_Coffee_Branch adonis:R 0.075 p: 0.103
#8         Coffee_TrunkLow_VS_Tree_Twig adonis:R 0.106 p: 0.005
#9  Coffee_TrunkLow_VS_Coffee_TrunkHigh adonis:R 0.066 p: 0.533
#10      Coffee_TrunkLow_VS_Tree_Branch adonis:R 0.101 p: 0.013
#11   Coffee_TrunkLow_VS_Tree_TrunkHigh adonis:R 0.082 p: 0.619
#12    Coffee_TrunkLow_VS_Tree_TrunkLow  adonis:R 0.09 p: 0.546
#13    Coffee_TrunkLow_VS_Coffee_Branch adonis:R 0.068 p: 0.213
#14       Tree_Twig_VS_Coffee_TrunkHigh adonis:R 0.105 p: 0.044
#15            Tree_Twig_VS_Tree_Branch adonis:R 0.143 p: 0.407
#16         Tree_Twig_VS_Tree_TrunkHigh adonis:R 0.154 p: 0.117
#17          Tree_Twig_VS_Tree_TrunkLow  adonis:R 0.16 p: 0.027
#18          Tree_Twig_VS_Coffee_Branch adonis:R 0.105 p: 0.003
#19     Coffee_TrunkHigh_VS_Tree_Branch adonis:R 0.103  p: 0.05
#20  Coffee_TrunkHigh_VS_Tree_TrunkHigh adonis:R 0.098 p: 0.134
#21   Coffee_TrunkHigh_VS_Tree_TrunkLow   adonis:R 0.1 p: 0.035
#22   Coffee_TrunkHigh_VS_Coffee_Branch adonis:R 0.067 p: 0.442
#23       Tree_Branch_VS_Tree_TrunkHigh  adonis:R 0.14 p: 0.598
#24        Tree_Branch_VS_Tree_TrunkLow adonis:R 0.149 p: 0.189
#25        Tree_Branch_VS_Coffee_Branch adonis:R 0.101 p: 0.019
#26     Tree_TrunkHigh_VS_Tree_TrunkLow adonis:R 0.138     p: 1
#27     Tree_TrunkHigh_VS_Coffee_Branch adonis:R 0.094 p: 0.131
#28      Tree_TrunkLow_VS_Coffee_Branch adonis:R 0.094 p: 0.141

#[1] "adonis:R 0.16 p: 0.005"

#boxplot for FUNGUIld Results https://r-coder.com/boxplot-en-r/
#The file Fungi_taxonomic_guild_assignmen sheet in S2_ITS_Sequence_count_Taxon_Guild_assignment.xlsx was edited:
#Summarized information of ASVs relative read counts per microsite, and split into tab delimeted files per each guild
#Each guild was imported in R to create Boxplots

#Create vector for colours
myrainbow <- c("#FF0000","#FFBF00","#80FF00","green4","#00FFFF","#0040FF","#8000FF","#FF00BF")

#Create graphs for guild ' Pathotroph-Saprotroph'
#Import table
PatSap <-read.table("pat-sap.csv", header=TRUE, row.names=1)
#Stack counts per column, each of them representing a microsite
PatSap_stackNA <- stack(PatSap)
#As we will transform counts to log scale, 0 were converted to NAs, to avoid inf values
PatSap_stackNA$values[(PatSap_stackNA$values == 0)] <- NA
    head(PatSap_stackNA)
#save as png
png("PatSap.png", width = 8, height = 4, units= "in", res = 600, pointsize = 4)
#set graph parameters
par(mar = c(8, 11, 5, 2), cex.axis = 3, cex.lab  = 2, cex.main = 6)
#plot
        boxplot(PatSap_stackNA$values ~ PatSap_stackNA$ind,
        log = "y", col = myrainbow, las = 2,
        main = "Pathotroph-Saprotroph",
        xlab = "",  # Etiqueta eje X Pathotroph-Saprotroph
        ylab = "") #Relative read count (log)

dev.off()

#Create graphs for guild 'Pathotroph-Saprotroph-Symbiotroph'
#Import table
PatSapSym <-read.table("pat-sap-sym.csv", header=TRUE, row.names=1)
#Stack counts per column, each of them representing a microsite
PatSapSym_stackNA <- stack(PatSapSym)
#As we will transform counts to log scale, 0 were converted to NAs, to avoid inf values
PatSapSym_stackNA$values[(PatSapSym_stackNA$values == 0)] <- NA
  head(PatSapSym_stackNA)
#save as png
png("PatSapSym.png", width = 8, height = 4, units= "in", res = 600, pointsize = 4)
#set graph parameters
par(mar = c(8, 11, 5, 2), cex.axis = 3, cex.lab  = 2, cex.main = 6)
#plot
        boxplot(PatSapSym_stackNA$values ~ PatSapSym_stackNA$ind,
        log = "y", col = myrainbow, las = 2,
        main = "Pathotroph-Saprotroph-Symbiotroph",
        xlab = "",  # Etiqueta eje X Pathotroph-Saprotroph-Symbiotroph
        ylab = "") #Relative read count (log)

dev.off()

#Create graphs for guild 'Pathotroph'
#Import table
Pato <-read.table("pato.csv", header=TRUE, row.names=1)
#Stack counts per column, each of them representing a microsite
Pato_stackNA <- stack(Pato)
#As we will transform counts to log scale, 0 were converted to NAs, to avoid inf values
Pato_stackNA$values[(Pato_stackNA$values == 0)] <- NA
    head(Pato_stackNA)
#save as png
png("Pato.png", width = 8, height = 4, units= "in", res = 600, pointsize = 4)
#set graph parameters
par(mar = c(8, 11, 5, 2), cex.axis = 3, cex.lab  = 2, cex.main = 6)
#plot
        boxplot(Pato_stackNA$values ~ Pato_stackNA$ind,
        log = "y", col = myrainbow, las = 2,
        main = "Pathotroph",
        xlab = "",  # Etiqueta eje X Pathotroph
        ylab = "") #Relative read count (log)

dev.off()

#Create graphs for guild 'Pathogen-Saprotroph-Symbiotroph'
#Import table
PathogenSapSym <-read.table("pathogen-sap-sym.csv", header=TRUE, row.names=1)
#Stack counts per column, each of them representing a microsite
PathogenSapSym_stackNA <- stack(PathogenSapSym)
#As we will transform counts to log scale, 0 were converted to NAs, to avoid inf values
PathogenSapSym_stackNA$values[(PathogenSapSym_stackNA$values == 0)] <- NA
    head(PathogenSapSym_stackNA)
#save as png
png("PathogenSapSym.png", width = 8, height = 4, units= "in", res = 600, pointsize = 4)
#set graph parameters
par(mar = c(8, 11, 5, 2), cex.axis = 3, cex.lab  = 2, cex.main = 6)
#plot
        boxplot(PathogenSapSym_stackNA$values ~ PathogenSapSym_stackNA$ind,
        log = "y", col = myrainbow, las = 2,
        main = "Pathogen-Saprotroph-Symbiotroph",
        xlab = "",  # Etiqueta eje X Pathogen-Saprotroph-Symbiotroph
        ylab = "") #Relative read count (log)

dev.off()


#Create graphs for guild 'Pathotroph-Symbiotroph'
#Import table
PatSym <-read.table("pat-sym.csv", header=TRUE, row.names=1)
#Stack counts per column, each of them representing a microsite
PatSym_stackNA <- stack(PatSym)
#As we will transform counts to log scale, 0 were converted to NAs, to avoid inf values
PatSym_stackNA$values[(PatSym_stackNA$values == 0)] <- NA
    head(PatSym_stackNA)
#save as png
png("PatSym.png", width = 8, height = 4, units= "in", res = 600, pointsize = 4)
#set graph parameters
par(mar = c(8, 11, 5, 2), cex.axis = 3, cex.lab  = 2, cex.main = 6)
#plot
        boxplot(PatSym_stackNA$values ~ PatSym_stackNA$ind,
        log = "y", col = myrainbow, las = 2,
        main = "Pathotroph-Symbiotroph",
        xlab = "",  # Etiqueta eje X Pathotroph-Symbiotroph
        ylab = "") #Relative read count (log)

dev.off()

#Create graphs for guild 'Saprotroph'
#Import table
Sap <-read.table("sap.csv", header=TRUE, row.names=1)
#Stack counts per column, each of them representing a microsite
Sap_stackNA <- stack(Sap)
#As we will transform counts to log scale, 0 were converted to NAs, to avoid inf values
Sap_stackNA$values[(Sap_stackNA$values == 0)] <- NA
    head(Sap_stackNA)
#save as png
png("Sap.png", width = 8, height = 4, units= "in", res = 600, pointsize = 4)
#set graph parameters
par(mar = c(8, 11, 5, 2), cex.axis = 3, cex.lab  = 2, cex.main = 6)
#plot
        boxplot(Sap_stackNA$values ~ Sap_stackNA$ind,
        log = "y", col = myrainbow, las = 2,
        main = "Saprotroph",
        xlab = "",  # Etiqueta eje X Saprotroph
        ylab = "") #Relative read count (log)

dev.off()

#Create graphs for guild 'Saprotroph-Pathotroph-Symbiotroph'
#Import table
SapPatSym <-read.table("sap-pat-sym.csv", header=TRUE, row.names=1)
#Stack counts per column, each of them representing a microsite
SapPatSym_stackNA <- stack(SapPatSym)
#As we will transform counts to log scale, 0 were converted to NAs, to avoid inf values
SapPatSym_stackNA$values[(SapPatSym_stackNA$values == 0)] <- NA
    head(SapPatSym_stackNA)
#save as png
png("SapPatSym.png", width = 8, height = 4, units= "in", res = 600, pointsize = 4)
#set graph parameters
par(mar = c(8, 11, 5, 2), cex.axis = 3, cex.lab  = 2, cex.main = 6)
#plot
        boxplot(SapPatSym_stackNA$values ~ SapPatSym_stackNA$ind,
        log = "y", col = myrainbow, las = 2,
        main = "Saprotroph-Pathotroph-Symbiotroph",
        xlab = "",  # Etiqueta eje X Saprotroph-Pathotroph-Symbiotroph
        ylab = "") #Relative read count (log)

dev.off()

#Create graphs for guild 'Saprotroph-Symbiotroph'
#Import table
SapSym <-read.table("sap-sym.csv", header=TRUE, row.names=1)
#Stack counts per column, each of them representing a microsite
SapSym_stackNA <- stack(SapSym)
#As we will transform counts to log scale, 0 were converted to NAs, to avoid inf values
SapSym_stackNA$values[(SapSym_stackNA$values == 0)] <- NA
    head(SapSym_stackNA)
#save as png
png("SapSym.png", width = 8, height = 4, units= "in", res = 600, pointsize = 4)
#set graph parameters
par(mar = c(8, 11, 5, 2), cex.axis = 3, cex.lab  = 2, cex.main = 6)
#plot
        boxplot(SapSym_stackNA$values ~ SapSym_stackNA$ind,
        log = "y", col = myrainbow, las = 2,
        main = "Saprotroph-Symbiotroph",
        xlab = "",  # Etiqueta eje X Saprotroph-Symbiotroph
        ylab = "") #Relative read count (log)

dev.off()

#Create graphs for guild 'Symbiotroph'
#Import table
Sym <-read.table("sym.csv", header=TRUE, row.names=1)
#Stack counts per column, each of them representing a microsite
Sym_stackNA <- stack(Sym)
#As we will transform counts to log scale, 0 were converted to NAs, to avoid inf values
Sym_stackNA$values[(Sym_stackNA$values == 0)] <- NA
    head(Sym_stackNA)
#save as png
png("Sym.png", width = 8, height = 4, units= "in", res = 600, pointsize = 4)
#set graph parameters
par(mar = c(8, 11, 5, 2), cex.axis = 3, cex.lab  = 2, cex.main = 6)
#plot
        boxplot(Sym_stackNA$values ~ Sym_stackNA$ind,
        log = "y", col = myrainbow, las = 2,
        main = "Symbiotroph",
        xlab = "",  # Etiqueta eje X Symbiotroph
        ylab = "") #Relative read count (log)

dev.off()

###################Create heatmap at class level ######
#Steps performed in QIIME2 and bash
### export table - creates a biom artifact
qiime tools export \
  --input-path filtered-Fungi-table-level3.qza \
  --output-path exported-filtered-tables
### change name
mv exported-filtered-tables/feature-table.biom exported-filtered-tables/filtered-Fungi-table-level3.biom
### export biom to tsv
biom convert -i exported-filtered-tables/filtered-Fungi-table-level3.biom -o filtered-Fungi-table-level3.tsv --to-tsv

### get first line with header
grep 'OTU ID' filtered-Fungi-table-level3.tsv | sed 's/\#OTU ID/Class/g'> filtered-Fungi-table-class.tsv

###get only class name and frequencies per sample
awk 'NR>3' filtered-Fungi-table-level3.tsv | awk -F ";" '{print $3}' | grep -v 'c__'$'\t' | grep -v 'c__unidentified' | grep -v '^__' | sed -e 's/c__//g' -e 's/\[//g' -e 's/\]/_b/g' -e 's/_cls_Incertae_sedis//g'>> filtered-Fungi-table-class.tsv

### this table needs to be transposed (https://www.thelinuxrain.com/articles/transposing-rows-and-columns-3-methods) to follow R heatmap tutorial. The program datamash was installed and executed.
cat filtered-Fungi-table-class.tsv | datamash transpose > filtered-Fungi-table-class.transposed.tsv

###work in R in a new conda environment
#follow tutorial https://www.molecularecologist.com/2013/08/making-heatmaps-with-r-for-microbiome-analysis/
R
#R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
#Install Heatplus
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Heatplus")
library(ggplot2)
library(Heatplus)
library(vegan)
library(RColorBrewer)

setwd("/home/eugenia/R-fungi-peerJ")

## create vector to color samples according to microsite
## "Coffee twig", 	red	"#FF0000",
## “Coffee trunk low”, 	orange	"#FFBF00",
## “Coffee trunk high”, 	Chartreuse	"#80FF00",
## “Coffee branch”, 	green4	"green4",
## “Tree twig”, 	cyan	"#00FFFF",
## “Tree trunk low”, 	navyblue	"#0040FF",
## “Tree trunk high”, 	purple	"#8000FF",
## “Tree branch”	magenta	"#FF00BF"


#colors defined with the rainbow palette, except for green4
micrositeRB <-c("#FF0000", "#FFBF00", "#00FFFF", "#80FF00", "#FF00BF", "#8000FF", "#0040FF", "#FFBF00", "#FF0000", "green4", "#80FF00", "#FFBF00", "#00FFFF", "#FF00BF", "#8000FF", "#0040FF", "green4", "#FF0000", "green4", "#80FF00", "#FFBF00", "#FF0000", "green4", "#80FF00", "#FF0000", "#FFBF00", "#00FFFF", "#80FF00", "#FF00BF", "#8000FF", "#0040FF", "#FFBF00", "#FF0000", "green4", "#80FF00", "#FFBF00", "#00FFFF", "#FF00BF", "#8000FF", "#0040FF", "green4", "#FF0000", "green4", "#80FF00", "#FFBF00", "#FF0000", "green4",  "#80FF00")

### colorRampPalette is in the RColorBrewer package.
### This creates a colour palette that shades from khaki1 to red in RGB space with 50 unique colours
scalekhakired <- colorRampPalette(c("khaki1", "red"), space = "rgb", bias=0.2)(50)

### import data
class.data <- read.csv("filtered-Fungi-table-class.transposed.tsv", sep="\t", header=TRUE)
### strip off the sample ids and convert them to row names so that the data matrix contains only sequence count data.
### the header for sample names was "Class" due to previous table transposing
row.names(class.data) <- class.data$Class
class.data <- class.data[, -1]
### transform the raw counts of reads to proportions within a sample:
class.prop <- class.data/rowSums(class.data)
class.prop[1:3, 1:3]
# calculate the Bray-Curtis dissimilarity matrix on the full dataset:
class.dist <- vegdist(class.prop, method = "bray")
### Do average linkage hierarchical clustering. Other options are 'complete' or 'single'.
### You'll need to choose the one that best fits the needs of your situation and your data.
row.clus <- hclust(class.dist, "aver")
###  add a column dendrogram to cluster the classes that occur more often together.
### you have to transpose the dataset to get the genera as rows
class.dist.t <- vegdist(t(class.prop), method = "bray")
col.clus <- hclust(class.dist.t, "aver")
### make the heatmap with Rowv = as.dendrogram(row.clus)

#save as a png file
png("class_hm3.png", width = 15, height = 10, units= "in", res = 600, pointsize = 18)
#set graph parameters, margins parameters are the reference for legend location
par(mar = c(8, 6.5, 0.5, 4), cex.axis = 4, cex.lab = 4, cex.main = 6)
#plot heatmap
heatmap(as.matrix(class.prop), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scalekhakired, RowSideColors = micrositeRB, margins = c(10, 4))

#followed https://r-coder.com/leyenda-r/ #x = "topleft"
legend(x = "topleft", cex=0.65, #position and size of legend
      legend = c("Coffee twig", "Coffee trunk low", "Coffee trunk high", "Coffee branch", "Tree twig", "Tree trunk low", "Tree trunk high", "Tree branch"),
      fill = c("#FF0000", "#FFBF00", "#80FF00", "green4", "#00FFFF", "#0040FF", "#8000FF", "#FF00BF"), #follows rainbow palette
      border = "black")
dev.off()

#Analysis with indicspecies, created new environment with bioconda
#Installed devtools with conda 4.10.3
R
#R version 4.1.3 (2022-03-10) -- "One Push-Up"
#install indicspecies https://emf-creaf.github.io/indicspecies/
devtools::install_github("emf-creaf/indicspecies")
install.packages("permute")
#followed tutorial https://cran.r-project.org/web/packages/indicspecies/vignettes/IndicatorSpeciesAnalysis.html#bibliography
setwd("/home/eugenia/R-fungi-peerJ/")
library(indicspecies)

cat Fungi-assigned-OTU-table.csv | datamash transpose > Fungi-assigned-OTU.transposed.tsv

Fungi_OTUs<-read.table("/home/eugenia/R-fungi-peerJ/Fungi-assigned-OTU.transposed.tsv", header=TRUE, row.names=1)
groups <- c("CoffeeTwig", "CoffeeTrunkLow", "TreesTwig", "CoffeeTrunkHigh", "TreesBranch", "TreesTrunkHigh", "TreesTrunkLow", "CoffeeTrunkLow", "CoffeeTwig", "CoffeeBranch", "CoffeeTrunkHigh", "CoffeeTrunkLow", "TreesTwig", "TreesBranch", "TreesTrunkHigh", "TreesTrunkLow", "CoffeeBranch", "CoffeeTwig", "CoffeeBranch", "CoffeeTrunkHigh", "CoffeeTrunkLow", "CoffeeTwig", "CoffeeBranch", "CoffeeTrunkHigh", "CoffeeTwig", "CoffeeTrunkLow", "TreesTwig", "CoffeeTrunkHigh", "TreesBranch", "TreesTrunkHigh", "TreesTrunkLow", "CoffeeTrunkLow", "CoffeeTwig", "CoffeeBranch", "CoffeeTrunkHigh", "CoffeeTrunkLow", "TreesTwig", "TreesBranch", "TreesTrunkHigh", "TreesTrunkLow", "CoffeeBranch", "CoffeeTwig", "CoffeeBranch", "CoffeeTrunkHigh", "CoffeeTrunkLow", "CoffeeTwig", "CoffeeBranch", "CoffeeTrunkHigh")

indval = multipatt(Fungi_OTUs, groups,
                   control = how(nperm=999))

                   #extract widespread taxa
                   widespred <- indval$sign
                   write.table(widespred, file = "widespread.csv")
