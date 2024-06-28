#set directory and load background file
setwd("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/PCAs_final/")
background <- read.csv("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/background/DC174_background.csv", header = T, sep = ",")
bacakground_uniq_col <- unique(background$colour)
bacakground_uniq_country <- unique(background$country)

#########
# all coloured by country 
#########

################# 1. 1 all indivs ############################### 
# plot the full 174 dataset (without outgroups)
#load libraries
library(ggplot2)
library("ggrepel")
library(tidyverse)
library(ggplot2)

pca_list <- c(1:9)
titles <- c("all_excl_chr15", "chr15_full", "chr15_nosvs", "chr15_region1", "chr15_region2", "chr15_region4", "chr15_region3", "chr15_region1.1", "chr15_region1.2")

for(i in 1:9){
  pca_numb <- pca_list[i]
  title <- titles[i]
  if(i == 1){
    eigenvec <- paste(pca_numb, "_all_filt_out.eigenvec", sep = "")
    eigenval <- paste(pca_numb, "_all_filt_out.eigenval", sep = "")
  }
  if(i > 1){
    eigenvec <- paste(pca_numb, "_chr15_filt_out.eigenvec", sep = "")
    eigenval <- paste(pca_numb, "_chr15_filt_out.eigenval", sep = "")
  }
  
  pca <- read_table(eigenvec, col_names = FALSE)
  eigenval <- scan(eigenval)
  
  pca <- pca[,-1]
  
  # set names
  names(pca)[1] <- "ind"
  names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
  
  colour_plot <- rep(NA, length(pca$ind))
  for (i in 1:nrow(pca)){
    col_name <- subset(background, as.character(background$indiv) == as.character(pca$ind[i]))
    colour_plot[i] <- as.character(col_name$region)
  }
  
  pca <- as.data.frame(data.frame(pca, colour_plot))
  
  pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
  
  # make plot
  a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
  a + ylab("Percentage variance explained") + theme_light()
  
  # calculate the cumulative sum of the percentage variance explained
  cumsum(pve$pve)
  tiff(paste("plots/", title, "_PC1PC2_", ".tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC1, 
       pca$PC2, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
       ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
  #legend('topleft', 
       #  legend = c(bacakground_uniq_country), 
       #  col = c(bacakground_uniq_col),
       #  pch = 16,
       #  pt.cex = 1,
       #  cex = 0.5) 
  dev.off()
  
  tiff(paste("plots/", title, "_PC1PC3_", ".tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC1, 
       pca$PC3, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
       ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC1, pca$PC3, pca$ind, cex = 0.5, pos = 3)
  #legend('topright', 
        # legend = c(bacakground_uniq_country), 
        # col = c(bacakground_uniq_col),
        # pch = 16,
        # pt.cex = 1,
        # cex = 0.5) 
  dev.off()
  
  tiff(paste("plots/", title, "_PC2PC3_", ".tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC2, 
       pca$PC3, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
       ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)
  #legend('topright', 
        # legend = c(bacakground_uniq_country), 
        # col = c(bacakground_uniq_col),
        # pch = 16,
        # pt.cex = 1,
        # cex = 0.5) 
  dev.off()
}

################# 1. 2 all indivs - not filtered ############################### 

pca_list <- c(1:9)
titles <- c("all_excl_chr15", "chr15_full", "chr15_nosvs", "chr15_region1", "chr15_region2", "chr15_region4", "chr15_region3", "chr15_region1.1", "chr15_region1.2")

for(i in 1:9){
  pca_numb <- pca_list[i]
  title <- titles[i]
  if(i == 1){
    eigenvec <- paste(pca_numb, "_all_unfilt_out.eigenvec", sep = "")
    eigenval <- paste(pca_numb, "_all_unfilt_out.eigenval", sep = "")
  }
  if(i > 1){
    eigenvec <- paste(pca_numb, "_chr15_unfilt_out.eigenvec", sep = "")
    eigenval <- paste(pca_numb, "_chr15_unfilt_out.eigenval", sep = "")
  }
  
  pca <- read_table(eigenvec, col_names = FALSE)
  eigenval <- scan(eigenval)
  
  pca <- pca[,-1]
  
  # set names
  names(pca)[1] <- "ind"
  names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
  
  colour_plot <- rep(NA, length(pca$ind))
  for (i in 1:nrow(pca)){
    col_name <- subset(background, as.character(background$indiv) == as.character(pca$ind[i]))
    colour_plot[i] <- as.character(col_name$region)
  }
  
  pca <- as.data.frame(data.frame(pca, colour_plot))
  
  pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
  
  # make plot
  a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
  a + ylab("Percentage variance explained") + theme_light()
  
  # calculate the cumulative sum of the percentage variance explained
  cumsum(pve$pve)
  tiff(paste("plots/", title, "_PC1PC2_unfilt", ".tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC1, 
       pca$PC2, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
       ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
  #legend('topleft', 
  #  legend = c(bacakground_uniq_country), 
  #  col = c(bacakground_uniq_col),
  #  pch = 16,
  #  pt.cex = 1,
  #  cex = 0.5) 
  dev.off()
  
  tiff(paste("plots/", title, "_PC1PC3_unfilt", ".tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC1, 
       pca$PC3, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
       ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC1, pca$PC3, pca$ind, cex = 0.5, pos = 3)
  #legend('topright', 
  # legend = c(bacakground_uniq_country), 
  # col = c(bacakground_uniq_col),
  # pch = 16,
  # pt.cex = 1,
  # cex = 0.5) 
  dev.off()
  
  tiff(paste("plots/", title, "_PC2PC3_unfilt", ".tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC2, 
       pca$PC3, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
       ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)
  #legend('topright', 
  # legend = c(bacakground_uniq_country), 
  # col = c(bacakground_uniq_col),
  # pch = 16,
  # pt.cex = 1,
  # cex = 0.5) 
  dev.off()
}

#########
#now all coloured by phenotype - background colour
#########

#do hindwing colour
background$backgroundColour_colour <- c()
for (i in 1:nrow(background)){
  if (background$backgroundColour[i]=="light"){
    background$backgroundColour_colour[i] <- "burlywood1"}
  if (background$backgroundColour[i]=="intermediate"){
    background$backgroundColour_colour[i] <- "darkorange"}
  if (background$backgroundColour[i]=="dark"){
    background$backgroundColour_colour[i] <- "darkorange4"}
  if (background$backgroundColour[i]=="missing"){
    background$backgroundColour_colour[i] <- "hotpink"}
}

################# 1. 3 all indivs ############################### 
# plot the full 90 dataset (without outgroups)
#load libraries

pca_list <- c(1:9)
titles <- c("all_excl_chr15", "chr15_full", "chr15_nosvs", "chr15_region1", "chr15_region2", "chr15_region4", "chr15_region3", "chr15_region1.1", "chr15_region1.2")

for(i in 1:9){
  pca_numb <- pca_list[i]
  title <- titles[i]
  if(i == 1){
    eigenvec <- paste(pca_numb, "_all_filt_out.eigenvec", sep = "")
    eigenval <- paste(pca_numb, "_all_filt_out.eigenval", sep = "")
  }
  if(i > 1){
    eigenvec <- paste(pca_numb, "_chr15_filt_out.eigenvec", sep = "")
    eigenval <- paste(pca_numb, "_chr15_filt_out.eigenval", sep = "")
  }
  
  pca <- read_table(eigenvec, col_names = FALSE)
  eigenval <- scan(eigenval)
  
  pca <- pca[,-1]
  
  # set names
  names(pca)[1] <- "ind"
  names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
  
  colour_plot <- rep(NA, length(pca$ind))
  for (i in 1:nrow(pca)){
    col_name <- subset(background, as.character(background$indiv) == as.character(pca$ind[i]))
    colour_plot[i] <- as.character(col_name$backgroundColour_colour)
  }
  
  pca <- as.data.frame(data.frame(pca, colour_plot))
  
  pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
  
  # make plot
  a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
  a + ylab("Percentage variance explained") + theme_light()
  
  # calculate the cumulative sum of the percentage variance explained
  cumsum(pve$pve)
  tiff(paste("plots/phenotype_coloured/background/", title, "_PC1PC2_", "background.tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC1, 
       pca$PC2, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
       ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
  #legend('topleft', 
  #  legend = c(bacakground_uniq_country), 
  #  col = c(bacakground_uniq_col),
  #  pch = 16,
  #  pt.cex = 1,
  #  cex = 0.5) 
  dev.off()
  
  tiff(paste("plots/phenotype_coloured/background/", title, "_PC1PC3_", "background.tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC1, 
       pca$PC3, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
       ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC1, pca$PC3, pca$ind, cex = 0.5, pos = 3)
  #legend('topright', 
  # legend = c(bacakground_uniq_country), 
  # col = c(bacakground_uniq_col),
  # pch = 16,
  # pt.cex = 1,
  # cex = 0.5) 
  dev.off()
  
  tiff(paste("plots/phenotype_coloured/background/", title, "_PC2PC3_", "background.tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC2, 
       pca$PC3, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
       ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)
  #legend('topright', 
  # legend = c(bacakground_uniq_country), 
  # col = c(bacakground_uniq_col),
  # pch = 16,
  # pt.cex = 1,
  # cex = 0.5) 
  dev.off()
}


################# 1. 4 all indivs - not filtered ############################### 

pca_list <- c(1:9)
titles <- c("all_excl_chr15", "chr15_full", "chr15_nosvs", "chr15_region1", "chr15_region2", "chr15_region4", "chr15_region3", "chr15_region1.1", "chr15_region1.2")

for(i in 1:9){
  pca_numb <- pca_list[i]
  title <- titles[i]
  if(i == 1){
    eigenvec <- paste(pca_numb, "_all_unfilt_out.eigenvec", sep = "")
    eigenval <- paste(pca_numb, "_all_unfilt_out.eigenval", sep = "")
  }
  if(i > 1){
    eigenvec <- paste(pca_numb, "_chr15_unfilt_out.eigenvec", sep = "")
    eigenval <- paste(pca_numb, "_chr15_unfilt_out.eigenval", sep = "")
  }
  
  pca <- read_table(eigenvec, col_names = FALSE)
  eigenval <- scan(eigenval)
  
  pca <- pca[,-1]
  
  # set names
  names(pca)[1] <- "ind"
  names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
  
  colour_plot <- rep(NA, length(pca$ind))
  for (i in 1:nrow(pca)){
    col_name <- subset(background, as.character(background$indiv) == as.character(pca$ind[i]))
    colour_plot[i] <- as.character(col_name$backgroundColour_colour)
  }
  
  pca <- as.data.frame(data.frame(pca, colour_plot))
  
  pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
  
  # make plot
  a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
  a + ylab("Percentage variance explained") + theme_light()
  
  # calculate the cumulative sum of the percentage variance explained
  cumsum(pve$pve)
  tiff(paste("plots/phenotype_coloured/background/", title, "_PC1PC2_unfilt_", "background.tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC1, 
       pca$PC2, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
       ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
  #legend('topleft', 
  #  legend = c(bacakground_uniq_country), 
  #  col = c(bacakground_uniq_col),
  #  pch = 16,
  #  pt.cex = 1,
  #  cex = 0.5) 
  dev.off()
  
  tiff(paste("plots/phenotype_coloured/background/", title, "_PC1PC3_unfilt_", "background.tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC1, 
       pca$PC3, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
       ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC1, pca$PC3, pca$ind, cex = 0.5, pos = 3)
  #legend('topright', 
  # legend = c(bacakground_uniq_country), 
  # col = c(bacakground_uniq_col),
  # pch = 16,
  # pt.cex = 1,
  # cex = 0.5) 
  dev.off()
  
  tiff(paste("plots/phenotype_coloured/background/", title, "_PC2PC3_unfilt_", "background.tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC2, 
       pca$PC3, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
       ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)
  #legend('topright', 
  # legend = c(bacakground_uniq_country), 
  # col = c(bacakground_uniq_col),
  # pch = 16,
  # pt.cex = 1,
  # cex = 0.5) 
  dev.off()
}



#########
#now all coloured by phenotype - forewing band
#########

background$forewingBand_colour <- c()
for (i in 1:nrow(background)){
  if (background$forewingBand[i]=="absent"){
    background$forewingBand_colour[i] <- "darkorange"}
  if (background$forewingBand[i]=="partial"){
    background$forewingBand_colour[i] <- "grey"}
  if (background$forewingBand[i]=="present"){
    background$forewingBand_colour[i] <- "black"}
  if (background$forewingBand[i]=="missing"){
    background$forewingBand_colour[i] <- "hotpink"}
}

################# 1. 9 all indivs ############################### 
# plot the full 90 dataset (without outgroups)
#load libraries

pca_list <- c(1:9)
titles <- c("all_excl_chr15", "chr15_full", "chr15_nosvs", "chr15_region1", "chr15_region2", "chr15_region4", "chr15_region3", "chr15_region1.1", "chr15_region1.2")

for(i in 1:9){
  pca_numb <- pca_list[i]
  title <- titles[i]
  if(i == 1){
    eigenvec <- paste(pca_numb, "_all_filt_out.eigenvec", sep = "")
    eigenval <- paste(pca_numb, "_all_filt_out.eigenval", sep = "")
  }
  if(i > 1){
    eigenvec <- paste(pca_numb, "_chr15_filt_out.eigenvec", sep = "")
    eigenval <- paste(pca_numb, "_chr15_filt_out.eigenval", sep = "")
  }
  
  pca <- read_table(eigenvec, col_names = FALSE)
  eigenval <- scan(eigenval)
  
  pca <- pca[,-1]
  
  # set names
  names(pca)[1] <- "ind"
  names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
  
  colour_plot <- rep(NA, length(pca$ind))
  for (i in 1:nrow(pca)){
    col_name <- subset(background, as.character(background$indiv) == as.character(pca$ind[i]))
    colour_plot[i] <- as.character(col_name$forewingBand_colour)
  }
  
  pca <- as.data.frame(data.frame(pca, colour_plot))
  
  pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
  
  # make plot
  a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
  a + ylab("Percentage variance explained") + theme_light()
  
  # calculate the cumulative sum of the percentage variance explained
  cumsum(pve$pve)
  tiff(paste("plots/phenotype_coloured/forewing/", title, "_PC1PC2_", "forewing_band.tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC1, 
       pca$PC2, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
       ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
  #legend('topleft', 
  #  legend = c(bacakground_uniq_country), 
  #  col = c(bacakground_uniq_col),
  #  pch = 16,
  #  pt.cex = 1,
  #  cex = 0.5) 
  dev.off()
  
  tiff(paste("plots/phenotype_coloured/forewing/", title, "_PC1PC3_", "forewing_band.tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC1, 
       pca$PC3, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
       ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC1, pca$PC3, pca$ind, cex = 0.5, pos = 3)
  #legend('topright', 
  # legend = c(bacakground_uniq_country), 
  # col = c(bacakground_uniq_col),
  # pch = 16,
  # pt.cex = 1,
  # cex = 0.5) 
  dev.off()
  
  tiff(paste("plots/phenotype_coloured/forewing/", title, "_PC2PC3_", "forewing_band.tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC2, 
       pca$PC3, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
       ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)
  #legend('topright', 
  # legend = c(bacakground_uniq_country), 
  # col = c(bacakground_uniq_col),
  # pch = 16,
  # pt.cex = 1,
  # cex = 0.5) 
  dev.off()
}

################# 1. 11 all indivs - not filtered ############################### 

pca_list <- c(1:9)
titles <- c("all_excl_chr15", "chr15_full", "chr15_nosvs", "chr15_region1", "chr15_region2", "chr15_region4", "chr15_region3", "chr15_region1.1", "chr15_region1.2")

for(i in 1:9){
  pca_numb <- pca_list[i]
  title <- titles[i]
  if(i == 1){
    eigenvec <- paste(pca_numb, "_all_unfilt_out.eigenvec", sep = "")
    eigenval <- paste(pca_numb, "_all_unfilt_out.eigenval", sep = "")
  }
  if(i > 1){
    eigenvec <- paste(pca_numb, "_chr15_unfilt_out.eigenvec", sep = "")
    eigenval <- paste(pca_numb, "_chr15_unfilt_out.eigenval", sep = "")
  }
  
  pca <- read_table(eigenvec, col_names = FALSE)
  eigenval <- scan(eigenval)
  
  pca <- pca[,-1]
  
  # set names
  names(pca)[1] <- "ind"
  names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
  
  colour_plot <- rep(NA, length(pca$ind))
  for (i in 1:nrow(pca)){
    col_name <- subset(background, as.character(background$indiv) == as.character(pca$ind[i]))
    colour_plot[i] <- as.character(col_name$forewingBand_colour)
  }
  
  pca <- as.data.frame(data.frame(pca, colour_plot))
  
  pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
  
  # make plot
  a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
  a + ylab("Percentage variance explained") + theme_light()
  
  # calculate the cumulative sum of the percentage variance explained
  cumsum(pve$pve)
  tiff(paste("plots/phenotype_coloured/forewing/", title, "_PC1PC2_unfilt_", "forewing_band.tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC1, 
       pca$PC2, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
       ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
  #legend('topleft', 
  #  legend = c(bacakground_uniq_country), 
  #  col = c(bacakground_uniq_col),
  #  pch = 16,
  #  pt.cex = 1,
  #  cex = 0.5) 
  dev.off()
  
  tiff(paste("plots/phenotype_coloured/forewing/", title, "_PC1PC3_unfilt_", "forewing_band.tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC1, 
       pca$PC3, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
       ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC1, pca$PC3, pca$ind, cex = 0.5, pos = 3)
  #legend('topright', 
  # legend = c(bacakground_uniq_country), 
  # col = c(bacakground_uniq_col),
  # pch = 16,
  # pt.cex = 1,
  # cex = 0.5) 
  dev.off()
  
  tiff(paste("plots/phenotype_coloured/forewing/", title, "_PC2PC3_unfilt_", "forewing_band.tiff", sep = ""), height=8, width=8, units="in", res=300, compression="lzw")
  plot(pca$PC2, 
       pca$PC3, 
       col = as.character(pca$colour_plot), 
       pch=16,
       lwd = 2,
       xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
       ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
       main = title)
  text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)
  #legend('topright', 
  # legend = c(bacakground_uniq_country), 
  # col = c(bacakground_uniq_col),
  # pch = 16,
  # pt.cex = 1,
  # cex = 0.5) 
  dev.off()
}


###################extra#######
#now want PC2/3 vs colour

pca <- read_table("9_chr15_unfilt_out.eigenvec", col_names = FALSE)
eigenval <- scan("9_chr15_unfilt_out.eigenval")
  
pca <- pca[,-1]
  
  # set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
  
colour_plot <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  col_name <- subset(background, as.character(background$indiv) == as.character(pca$ind[i]))
  colour_plot[i] <- as.character(col_name$backgroundColour_colour)
}
  
pca <- as.data.frame(data.frame(pca, colour_plot))
  
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
  
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
  

light_df <- subset(pca, pca$colour_plot == "burlywood1")
int_df <- subset(pca, pca$colour_plot == "darkorange")
dark_df <- subset(pca, pca$colour_plot == "darkorange4")

all_df <- as.data.frame(rbind(light_df, int_df, dark_df))

par(mfrow=c(3,1), mai = c(0.4, 0.6, 0.2, 0.1))
plot(as.integer(as.factor(all_df$colour_plot)), all_df$PC1, 
     type = 'p', main = "region 1.2-PC1",
     xaxt = "n", ylab = "PC1")
axis(1, at = c(1, 2, 3), labels = c("light", "int", "dark"))

abline(lm(all_df$PC1~as.integer(as.factor(all_df$colour_plot))))
a_rsq <- summary(lm(all_df$PC1~as.integer(as.factor(all_df$colour_plot))))$r.squared
a <- summary(lm(all_df$PC1~as.integer(as.factor(all_df$colour_plot))))$coefficients[,4]  
a[[2]]
legend(2.5,0.05, legend = c(paste("r-squared=", round(a_rsq, 3), sep = ""), paste("p-value=", round(a[[2]], 6), sep = "")))

plot(as.integer(as.factor(all_df$colour_plot)), all_df$PC2, 
     type = 'p', main = "region 1.2-PC2",
     xaxt = "n", ylab = "PC2")
axis(1, at = c(1, 2, 3), labels = c("light", "int", "dark"))
abline(lm(all_df$PC2~as.integer(as.factor(all_df$colour_plot))))
b_rsq <- summary(lm(all_df$PC2~as.integer(as.factor(all_df$colour_plot))))$r.squared
b <- summary(lm(all_df$PC2~as.integer(as.factor(all_df$colour_plot))))$coefficients[,4]  
b[[2]]
legend(2.5,-0.2, legend = c(paste("r-squared=", round(b_rsq, 5), sep = ""), paste("p-value=", round(b[[2]], 6), sep = "")))


plot(as.integer(as.factor(all_df$colour_plot)), all_df$PC3, 
     type = 'p', main = "region 1.2-PC3",
     xaxt = "n", ylab = "PC3")
axis(1, at = c(1, 2, 3), labels = c("light", "int", "dark"))
abline(lm(all_df$PC3~as.integer(as.factor(all_df$colour_plot))))
c_rsq <- summary(lm(all_df$PC3~as.integer(as.factor(all_df$colour_plot))))$r.squared
c <- summary(lm(all_df$PC3~as.integer(as.factor(all_df$colour_plot))))$coefficients[,4]  
c[[2]]
legend(2.5,-0.10, legend = c(paste("r-squared=", round(c_rsq, 5), sep = ""), paste("p-value=", round(c[[2]], 6), sep = "")))



par(mfrow=c(3,1), mai = c(0.5, 0.3, 0.2, 0.1))
boxplot(light_df$PC1, int_df$PC1, dark_df$PC1, names= (c("light", "int", "dark")), main = "region 1.2 - PC1")
boxplot(light_df$PC2, int_df$PC2, dark_df$PC2, names= (c("light", "int", "dark")), main = "region 1.2 - PC2")
boxplot(light_df$PC3, int_df$PC3, dark_df$PC3, names= (c("light", "int", "dark")), main = "region 1.2 - PC3")

library(vioplot)
vioplot(light_df$PC1, int_df$PC1, dark_df$PC1, main = "region 1.2 - PC1")
vioplot(light_df$PC2, int_df$PC2, dark_df$PC2, main = "region 1.2 - PC2")
vioplot(light_df$PC3, int_df$PC3, dark_df$PC3, main = "region 1.2 - PC3")


#prepare sub-saharan input
a <- as.data.frame(subset(background, background$country == "Tunisia"))
b <- as.data.frame(subset(background, background$country == "Italy"))
c <- as.data.frame(subset(background, background$country == "missing"))
d <- as.data.frame(subset(background, background$country == "Philippines"))
e <- as.data.frame(subset(background, background$country == "Spain"))
f <- as.data.frame(subset(background, background$country == "StHelena"))
samples <- as.data.frame(rbind(a,b,c,d,e,f))
samples <- as.data.frame(samples$indiv)



pca <- read_table("10_chr15_sub_unfilt_out.eigenvec", col_names = FALSE)
eigenval <- scan("10_chr15_sub_unfilt_out.eigenval")

pca <- pca[,-1]

# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

#want to colour samples
#RSA <- "#2143d1"
#Nigeria and Ghana <- "#bc4754"
#hybrid zone <- Rwanda + SM16K (NRB) and SM19M (MPL) samples
colour_plot <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  col_name <- subset(background, as.character(background$indiv) == as.character(pca$ind[i]))
  col_type <- as.character(col_name$region)
  if(col_type == "royalblue2"){
    colour_plot[i] <- "#2143d1"
  }
  if(col_type == "springgreen4"){
    colour_plot[i] <- "#bc4754"
  }
  if(col_type == "tomato"){
    colour_plot[i] <- "#ffac07"
  }
  if(col_type == "darkgrey"){
    colour_plot[i] <- "black"
  }
  if(col_type == "tomato4"){
    colour_plot[i] <- "darkgrey"
  }
}

pca <- as.data.frame(data.frame(pca, colour_plot))

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)
tiff("plots/sub_saharan_chr15/sub_saharan_PC1PC2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
     main = "PC1 vs PC2")
text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
#legend('topleft', 
#  legend = c(bacakground_uniq_country), 
#  col = c(bacakground_uniq_col),
#  pch = 16,
#  pt.cex = 1,
#  cex = 0.5) 
dev.off()

tiff("plots/sub_saharan_chr15/sub_saharan_PC1PC3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "PC1 vs PC3")
text(pca$PC1, pca$PC3, pca$ind, cex = 0.5, pos = 3)
#legend('topright', 
# legend = c(bacakground_uniq_country), 
# col = c(bacakground_uniq_col),
# pch = 16,
# pt.cex = 1,
# cex = 0.5) 
dev.off()

tiff("plots/sub_saharan_chr15/sub_saharan_PC2PC3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC2, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "PC2 vs PC3")
text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)
#legend('topright', 
# legend = c(bacakground_uniq_country), 
# col = c(bacakground_uniq_col),
# pch = 16,
# pt.cex = 1,
# cex = 0.5) 
dev.off()

tiff("plots/sub_saharan_chr15/sub_saharan_PC1PC2_nolabs.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
     main = "PC1 vs PC2")
dev.off()

tiff("plots/sub_saharan_chr15/sub_saharan_PC1PC3_nolabs.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "PC1 vs PC3")
dev.off()

tiff("plots/sub_saharan_chr15/sub_saharan_PC2PC3_nolabs.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC2, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "PC2 vs PC3")
dev.off()

#prepare sub-saharan input
pca <- read_table("11_nochr15_sub_unfilt_out.eigenvec", col_names = FALSE)
eigenval <- scan("11_nochr15_sub_unfilt_out.eigenval")

pca <- pca[,-1]

# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

#want to colour samples
#RSA <- "#2143d1"
#Nigeria and Ghana <- "#bc4754"
#hybrid zone <- Rwanda + SM16K (NRB) and SM19M (MPL) samples
colour_plot <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  col_name <- subset(background, as.character(background$indiv) == as.character(pca$ind[i]))
  col_type <- as.character(col_name$region)
  if(col_type == "royalblue2"){
    colour_plot[i] <- "#2143d1"
  }
  if(col_type == "springgreen4"){
    colour_plot[i] <- "#bc4754"
  }
  if(col_type == "tomato"){
    colour_plot[i] <- "#ffac07"
  }
  if(col_type == "darkgrey"){
    colour_plot[i] <- "black"
  }
  if(col_type == "tomato4"){
    colour_plot[i] <- "darkgrey"
  }
}

pca <- as.data.frame(data.frame(pca, colour_plot))

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)
tiff("plots/sub_saharan_chr15/sub_saharan_PC1PC2_fullgeno.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
     main = "PC1 vs PC2")
text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
#legend('topleft', 
#  legend = c(bacakground_uniq_country), 
#  col = c(bacakground_uniq_col),
#  pch = 16,
#  pt.cex = 1,
#  cex = 0.5) 
dev.off()

tiff("plots/sub_saharan_chr15/sub_saharan_PC1PC3_fullgeno.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "PC1 vs PC3")
text(pca$PC1, pca$PC3, pca$ind, cex = 0.5, pos = 3)
#legend('topright', 
# legend = c(bacakground_uniq_country), 
# col = c(bacakground_uniq_col),
# pch = 16,
# pt.cex = 1,
# cex = 0.5) 
dev.off()

tiff("plots/sub_saharan_chr15/sub_saharan_PC2PC3_fullgeno.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC2, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "PC2 vs PC3")
text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)
#legend('topright', 
# legend = c(bacakground_uniq_country), 
# col = c(bacakground_uniq_col),
# pch = 16,
# pt.cex = 1,
# cex = 0.5) 
dev.off()

tiff("plots/sub_saharan_chr15/sub_saharan_PC1PC2_nolabs_fullgeno.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
     main = "PC1 vs PC2")
dev.off()

tiff("plots/sub_saharan_chr15/sub_saharan_PC1PC3_nolabs_fullgeno.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "PC1 vs PC3")
dev.off()

tiff("plots/sub_saharan_chr15/sub_saharan_PC2PC3_nolabs_fullgeno.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC2, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=16,
     lwd = 2,
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "PC2 vs PC3")
dev.off()
