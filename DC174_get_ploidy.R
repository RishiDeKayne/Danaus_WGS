#script to add sex info for haploid genotying of contig 1.1
setwd("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/")
background <- read.csv("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/background/DC174_background.csv", header = T, sep = ",")

indivs <- read.csv("./VCF_indiv_list.txt", header = F)

indivs$ploidy <- c()
for (i in 1:length(indivs$V1)){
  indiv_name <- as.character(indivs$V1[i])
  background_sub <- subset(background, as.character(background$indiv) == indiv_name)
  sex <- background_sub$sex
  if(sex == "M"){
    indivs$ploidy[i] <- 2
  }
  if(sex == "F"){
    indivs$ploidy[i] <- 1
  }
}

write.table(indivs, file = "sample.ploidy.txt", col.names = F, row.names = F, quote = F, sep = "\t")
