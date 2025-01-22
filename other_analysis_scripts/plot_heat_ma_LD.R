setwd("/Users/rdk_mac_work/Dropbox/RishiMac3Work/DC174_revision//")

library(LDheatmap)
require(snpStats)
library(ggplot2)
library(viridis)
viridis.palette <- rev(viridis(n = 100, alpha = 1))

#WHOLE VCF
#load data
z <- read.csv(file = "Square_LD_wholeVCF.ld.gz", sep = "\t", header = FALSE)
snp_file <- read.csv(file = "Square_LD_wholeVCF.bim", sep = "\t", header = FALSE)
SNPIDs <- snp_file$V4

#make into matrix
LD <-as.matrix(z)

tiff(paste("./LD_plot_wholeVCF", ".tiff", sep = ""), height=30, width=30, units="in", res=300)
LDheatmap(LD,genetic.distances=SNPIDs, distances = "physical",
          color=viridis.palette,
          LDmeasure="r", 
          title="Pairwise LD", 
          add.map=T, add.key=T,
          geneMapLabelX=NULL, geneMapLabelY=NULL, geneMapLocation = 0.1)
dev.off()

#now do reduced plot to see if it helps us interpret:
z2 <- z

z2[z2 < .4] <- as.numeric("NaN")

#make into matrix
LD <-as.matrix(z2)

tiff(paste("./LD_plot_R2_0.4_plus_wholeVCF", ".tiff", sep = ""), height=30, width=30, units="in", res=300)
LDheatmap(LD,genetic.distances=SNPIDs, distances = "physical",
          color=viridis.palette,
          LDmeasure="r", 
          title="Pairwise LD", 
          add.map=T, add.key=T,
          geneMapLabelX=NULL, geneMapLabelY=NULL, geneMapLocation = 0.1)
dev.off()

######################################
#CK
#load data
z <- read.csv(file = "Square_LD_CK_VCF.ld.gz", sep = "\t", header = FALSE)
snp_file <- read.csv(file = "Square_LD_CK_VCF.bim", sep = "\t", header = FALSE)
SNPIDs <- snp_file$V4

#make into matrix
LD <-as.matrix(z)

tiff(paste("./LD_plot_CK_VCF", ".tiff", sep = ""), height=30, width=30, units="in", res=300)
LDheatmap(LD,genetic.distances=SNPIDs, distances = "physical",
          color=viridis.palette,
          LDmeasure="r", 
          title="Pairwise LD", 
          add.map=T, add.key=T,
          geneMapLabelX=NULL, geneMapLabelY=NULL, geneMapLocation = 0.1)
dev.off()

#now do reduced plot to see if it helps us interpret:
z2 <- z

z2[z2 < .4] <- as.numeric("NaN")

#make into matrix
LD <-as.matrix(z2)

tiff(paste("./LD_plot_R2_0.4_plus_CK_VCF", ".tiff", sep = ""), height=30, width=30, units="in", res=300)
LDheatmap(LD,genetic.distances=SNPIDs, distances = "physical",
          color=viridis.palette,
          LDmeasure="r", 
          title="Pairwise LD", 
          add.map=T, add.key=T,
          geneMapLabelX=NULL, geneMapLabelY=NULL, geneMapLocation = 0.1)
dev.off()

######################################
#CO
#load data
z <- read.csv(file = "Square_LD_CO_VCF.ld.gz", sep = "\t", header = FALSE)
snp_file <- read.csv(file = "Square_LD_CO_VCF.bim", sep = "\t", header = FALSE)
SNPIDs <- snp_file$V4

#make into matrix
LD <-as.matrix(z)

tiff(paste("./LD_plot_CO_VCF", ".tiff", sep = ""), height=30, width=30, units="in", res=300)
LDheatmap(LD,genetic.distances=SNPIDs, distances = "physical",
          color=viridis.palette,
          LDmeasure="r", 
          title="Pairwise LD", 
          add.map=T, add.key=T,
          geneMapLabelX=NULL, geneMapLabelY=NULL, geneMapLocation = 0.1)
dev.off()

#now do reduced plot to see if it helps us interpret:
z2 <- z

z2[z2 < .4] <- as.numeric("NaN")

#make into matrix
LD <-as.matrix(z2)

tiff(paste("./LD_plot_R2_0.4_plus_CO_VCF", ".tiff", sep = ""), height=30, width=30, units="in", res=300)
LDheatmap(LD,genetic.distances=SNPIDs, distances = "physical",
          color=viridis.palette,
          LDmeasure="r", 
          title="Pairwise LD", 
          add.map=T, add.key=T,
          geneMapLabelX=NULL, geneMapLabelY=NULL, geneMapLocation = 0.1)
dev.off()

######################################
#KO
#load data
z <- read.csv(file = "Square_LD_KO_VCF.ld.gz", sep = "\t", header = FALSE)
snp_file <- read.csv(file = "Square_LD_KO_VCF.bim", sep = "\t", header = FALSE)
SNPIDs <- snp_file$V4

#make into matrix
LD <-as.matrix(z)

tiff(paste("./LD_plot_KO_VCF", ".tiff", sep = ""), height=30, width=30, units="in", res=300)
LDheatmap(LD,genetic.distances=SNPIDs, distances = "physical",
          color=viridis.palette,
          LDmeasure="r", 
          title="Pairwise LD", 
          add.map=T, add.key=T,
          geneMapLabelX=NULL, geneMapLabelY=NULL, geneMapLocation = 0.1)
dev.off()

#now do reduced plot to see if it helps us interpret:
z2 <- z

z2[z2 < .4] <- as.numeric("NaN")

#make into matrix
LD <-as.matrix(z2)

tiff(paste("./LD_plot_R2_0.4_plus_KO_VCF", ".tiff", sep = ""), height=30, width=30, units="in", res=300)
LDheatmap(LD,genetic.distances=SNPIDs, distances = "physical",
          color=viridis.palette,
          LDmeasure="r", 
          title="Pairwise LD", 
          add.map=T, add.key=T,
          geneMapLabelX=NULL, geneMapLabelY=NULL, geneMapLocation = 0.1)
dev.off()

