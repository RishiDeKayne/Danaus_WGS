#Code to make supplementary table of background info
setwd("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/background/")
background <- read.csv("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/background/DC174_background.csv", header = T, sep = ",")
extra_info <- read.csv("SM_Samples.xlsx - SM Samples.csv", na.strings=c("NA","GONE"))
extra_info$ID <- gsub("-", "", extra_info$ID)

newbackground <- background
newbackground$accession <- "TBD"
newbackground$white <- "unknown"
newbackground$background <- "unknown"
newbackground$forewing <- "unknown"

accession_df <- subset(extra_info, nchar(as.character(extra_info$Accession)) > 0)

for(i in 1:nrow(background)){
  individual_name <- as.character(background$indiv[i])
  if(individual_name %in% accession_df$ID){
    indiv_subdf <- subset(extra_info, extra_info$ID == individual_name)
    newbackground$accession[i] <- indiv_subdf$Accession
  }
}

white_df <- subset(extra_info, nchar(as.character(extra_info$hindwingWhite)) > 0)

for(i in 1:nrow(background)){
  individual_name <- as.character(background$indiv[i])
  if(individual_name %in% white_df$ID){
    indiv_subdf <- subset(extra_info, extra_info$ID == individual_name)
    newbackground$white[i] <- indiv_subdf$hindwingWhite
  }
}

ground_df <- subset(extra_info, nchar(as.character(extra_info$groundColour)) > 0)

for(i in 1:nrow(background)){
  individual_name <- as.character(background$indiv[i])
  if(individual_name %in% ground_df$ID){
    indiv_subdf <- subset(extra_info, extra_info$ID == individual_name)
    newbackground$background[i] <- indiv_subdf$groundColour
  }
}

forewing_df <- subset(extra_info, nchar(as.character(extra_info$forewingTip)) > 0)

for(i in 1:nrow(background)){
  individual_name <- as.character(background$indiv[i])
  if(individual_name %in% forewing_df$ID){
    indiv_subdf <- subset(extra_info, extra_info$ID == individual_name)
    newbackground$forewing[i] <- indiv_subdf$forewingTip
  }
}

#supp file
supp_file <- as.data.frame(cbind(newbackground$indiv, newbackground$country, newbackground$region_narrow,
                   newbackground$lat, newbackground$long, newbackground$sex, newbackground$hindwingWhite,
                   newbackground$forewingBand, newbackground$backgroundColour, newbackground$accession))

colnames(supp_file) <- c("ID", "Country", "Region", "Latitute (decimal)", "Longitude (decimal)",
                         "Sex", "White Hindwing", "Forewing Band", "Background Colour", "NCBI Accession")


write.table(supp_file, file = "/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/background/supplementary_table_file.csv", quote = F, sep = ",", row.names =F)
