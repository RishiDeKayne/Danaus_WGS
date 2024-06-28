#plot phasing 
setwd("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/phasing//")
background <- read.csv("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/background/DC174_background.csv", header = T, sep = ",")
bacakground_uniq_col <- unique(background$colour)
bacakground_uniq_country <- unique(background$country)

phase_stats <- read.csv("full_whatshap_stats.txt", sep = " ", head = FALSE)
phase_stats <- phase_stats[,-5]
phase_stats <- phase_stats[,1:5]
colnames(phase_stats) <- c("indiv", "contig", "blocks", "longest_block", "proportion")

#get contig names
chromosomes <- read.csv("../gemma/Dchry2.2.fa.fai", head = FALSE, sep = "\t")
rename <- read.csv("../gemma/rename_file.txt", head = FALSE, sep = " ")
rename$name <- paste(rename$V7, rename$V8)
rename$name <- gsub(" ", "", rename$name)

#add colour column to full data for plotting
phase_stats$colour <- c()

for (i in 1:nrow(phase_stats)){
  individual <- as.character(phase_stats$indiv)[i]
  col_df <- subset(background, background$indiv == individual)
  phase_stats$colour[i] <- col_df$colour
}

stat_df <- as.data.frame(rename$V1)
colnames(stat_df) <- c("contig")
stat_df$mean_block_number <- c()
stat_df$max_block_proportion <- c()
stat_df$mean_block_proportion <- c()

for (i in 1:length(rename$name)){
  chr <- as.character(rename$V1)[i]
  contig_df <- subset(phase_stats, phase_stats$contig == chr)
  plot(contig_df$blocks, pch = 16, col = contig_df$colour, main = paste("# blocks - ", chr, sep = ""), xlab = "individual", ylab = "number of phased blocks")
  plot(contig_df$proportion, pch = 16, col = contig_df$colour, main = paste("longest block % - ", chr, sep = ""), xlab = "individual", ylab = "% of contig covered by longest phased block")
  stat_df$mean_block_number[i] <- mean(contig_df$blocks)
  stat_df$max_block_proportion[i] <- max(contig_df$proportion)
  stat_df$mean_block_proportion[i] <- mean(contig_df$proportion)
}
