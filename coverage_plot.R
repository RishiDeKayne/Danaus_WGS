#coverage of DC174  
data <- read.csv(file = "Dropbox/RishiMac2/Danaus/DC174/bam_raw_coverage_full_plot.txt", sep = ",", header = F)
hist(data$V2, breaks = 10)

mean(data$V2)

colnames(data) <- c("indiv_bam", "coverage")
