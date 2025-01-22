#De-kayne et al. 

library(ggplot2)
library("ggrepel")
library(tidyverse)
library(ggplot2)

#files in /data/martin/genomics/analyses/Danaus_popgen/DC174/analysis/admixture/chr15_onlySVs
#set directory and load background file
setwd("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/PCAs_final/")
background <- read.csv("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/background/DC174_background.csv", header = T, sep = ",")
bacakground_uniq_col <- unique(background$colour)
bacakground_uniq_country <- unique(background$country)


morph_cols <- c(alcippus="#5db463", orientis="#2143d1", klugii="#ffac07", chrysippus="#bc4754", intermediate="black", unknown="gray")


#DC174 <- read.csv("./SVs_only/DC174_data.csv", as.is=T, header = FALSE)

background$decider <- "gray"

for(i in 1:nrow(background)){
  GC <- as.character(background$backgroundColour)[i]
  FT <- as.character(background$forewingBand)[i]
  if(GC == "light" & FT == "present"){
    background$decider[i] <- "#bc4754"
  }
  if(GC == "light" & FT == "absent"){
    background$decider[i] <- "#ffac07"
  }
  if(GC == "dark" & FT == "present"){
    background$decider[i] <- "#2143d1"
  }
  if(GC == "intermediate"){
    background$decider[i] <- "black"
  }
  if(FT == "partial"){
    background$decider[i] <- "black"
  }
}

eigenvec <- "./SVs_only/admix_in_chr15_onlySVs.eigenvec"
eigenval <- "./SVs_only/admix_in_chr15_onlySVs.eigenval"

pca <- read_table(eigenvec, col_names = FALSE)
eigenval <- scan(eigenval)

pca <- pca[,-1]

# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

colour_plot <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  col_name <- subset(background, as.character(background$indiv) == as.character(pca$ind[i]))
  colour_plot[i] <- as.character(col_name$decider)
}

pca <- as.data.frame(data.frame(pca, colour_plot))

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)
tiff("./SVs_only/SV_only_PCA_PC1PC2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75)
text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
#legend('topleft', 
#  legend = c(bacakground_uniq_country), 
#  col = c(bacakground_uniq_col),
#  pch = 16,
#  pt.cex = 1,
#  cex = 0.5) 
dev.off()

tiff("./SVs_only/SV_only_PCA_PC1PC2_nolab.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75)
dev.off()

tiff("./SVs_only/SV_only_PCA_PC1PC3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75)
text(pca$PC1, pca$PC3, pca$ind, cex = 0.5, pos = 3)
#legend('topright', 
# legend = c(bacakground_uniq_country), 
# col = c(bacakground_uniq_col),
# pch = 16,
# pt.cex = 1,
# cex = 0.5) 
dev.off()

tiff("./SVs_only/SV_only_PCA_PC1PC3_nolab.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75)
dev.off()

tiff("./SVs_only/SV_only_PCA_PC2PC3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC2, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75)
text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)
#legend('topright', 
# legend = c(bacakground_uniq_country), 
# col = c(bacakground_uniq_col),
# pch = 16,
# pt.cex = 1,
# cex = 0.5) 
dev.off()

tiff("./SVs_only/SV_only_PCA_PC2PC3_nolab.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC2, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75)
dev.off()

##SUB_SAHARAN
eigenvec <- "./SVs_only/admix_in_chr15_onlySVs_sahara.eigenvec"
eigenval <- "./SVs_only/admix_in_chr15_onlySVs_sahara.eigenval"

pca <- read_table(eigenvec, col_names = FALSE)
eigenval <- scan(eigenval)

pca <- pca[,-1]

# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

colour_plot <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  col_name <- subset(background, as.character(background$indiv) == as.character(pca$ind[i]))
  colour_plot[i] <- as.character(col_name$decider)
}

pca <- as.data.frame(data.frame(pca, colour_plot))

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)
tiff("./SVs_only/SV_only_PCA_SUB_PC1PC2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75)
text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
#legend('topleft', 
#  legend = c(bacakground_uniq_country), 
#  col = c(bacakground_uniq_col),
#  pch = 16,
#  pt.cex = 1,
#  cex = 0.5) 
dev.off()

tiff("./SVs_only/SV_only_PCA_SUB_PC1PC2_nolab.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75)
dev.off()

tiff("./SVs_only/SV_only_PCA_SUB_PC1PC3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75)
text(pca$PC1, pca$PC3, pca$ind, cex = 0.5, pos = 3)
#legend('topright', 
# legend = c(bacakground_uniq_country), 
# col = c(bacakground_uniq_col),
# pch = 16,
# pt.cex = 1,
# cex = 0.5) 
dev.off()

tiff("./SVs_only/SV_only_PCA_SUB_PC1PC3_nolab.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75)
dev.off()

tiff("./SVs_only/SV_only_PCA_SUB_PC2PC3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC2, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75)
text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)
#legend('topright', 
# legend = c(bacakground_uniq_country), 
# col = c(bacakground_uniq_col),
# pch = 16,
# pt.cex = 1,
# cex = 0.5) 
dev.off()

tiff("./SVs_only/SV_only_PCA_SUB_PC2PC3_nolab.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC2, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75)
dev.off()

#####################SUB_SAHARAN FULL GENOME (NO CHR15)
eigenvec <- "./SVs_only/11_nochr15_sub_unfilt_out.eigenvec"
eigenval <- "./SVs_only/11_nochr15_sub_unfilt_out.eigenval"

pca <- read_table(eigenvec, col_names = FALSE)
eigenval <- scan(eigenval)

pca <- pca[,-1]

# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

colour_plot <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  col_name <- subset(background, as.character(background$indiv) == as.character(pca$ind[i]))
  colour_plot[i] <- as.character(col_name$decider)
}

pca <- as.data.frame(data.frame(pca, colour_plot))

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)
tiff("./SVs_only/FULL_GENO_PCA_SUB_PC1PC2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75)
text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
#legend('topleft', 
#  legend = c(bacakground_uniq_country), 
#  col = c(bacakground_uniq_col),
#  pch = 16,
#  pt.cex = 1,
#  cex = 0.5) 
dev.off()

tiff("./SVs_only/FULL_GENO_PCA_SUB_PC1PC2_nolab.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75)
dev.off()

tiff("./SVs_only/FULL_GENO_PCA_SUB_PC1PC3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75)
text(pca$PC1, pca$PC3, pca$ind, cex = 0.5, pos = 3)
#legend('topright', 
# legend = c(bacakground_uniq_country), 
# col = c(bacakground_uniq_col),
# pch = 16,
# pt.cex = 1,
# cex = 0.5) 
dev.off()

tiff("./SVs_only/FULL_GENO_PCA_SUB_PC1PC3_nolab.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75)
dev.off()

tiff("./SVs_only/FULL_GENO_PCA_SUB_PC2PC3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC2, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75)
text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)
#legend('topright', 
# legend = c(bacakground_uniq_country), 
# col = c(bacakground_uniq_col),
# pch = 16,
# pt.cex = 1,
# cex = 0.5) 
dev.off()

tiff("./SVs_only/FULL_GENO_PCA_SUB_PC2PC3_nolab.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC2, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75)
dev.off()
