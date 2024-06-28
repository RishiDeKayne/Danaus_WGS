#gemma input
setwd("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/PCAs/")
background <- read.csv("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174//Background/DC174_background.csv", header = T)
bacakground_uniq_col <- unique(background$colour)
bacakground_uniq_country <- unique(background$country)

#####BACKGROUND####
#background phenotype
background$backgroundColour_valuer <- c()
for (i in 1:nrow(background)){
  if (background$backgroundColour[i]=="light"){
    background$backgroundColour_valuer[i] <- 1}
  if (background$backgroundColour[i]=="intermediate"){
    background$backgroundColour_valuer[i] <- 2}
  if (background$backgroundColour[i]=="dark"){
    background$backgroundColour_valuer[i] <- 3}
  if (background$backgroundColour[i]=="missing"){
    background$backgroundColour_valuer[i] <- "missing"}
}

#write file
background_colour_values <- as.data.frame(cbind(background$indiv, background$backgroundColour_valuer))
background_colour_values <- subset(background_colour_values, background_colour_values$V2 != "missing")
write.table(background_colour_values, "../gemma/background_colour_values.csv", col.names = FALSE, row.names=FALSE,sep="\t", quote = FALSE)

#####FOREWING####
#forewing band phenotype
background$forewing_band_valuer <- c()
for (i in 1:nrow(background)){
  if (background$forewingBand[i]=="absent"){
    background$forewing_band_valuer[i] <- 1}
  if (background$forewingBand[i]=="partial"){
    background$forewing_band_valuer[i] <- 2}
  if (background$forewingBand[i]=="present"){
    background$forewing_band_valuer[i] <- 3}
  if (background$forewingBand[i]=="missing"){
    background$forewing_band_valuer[i] <- "missing"}
}

#write file
forewing_band_values <- as.data.frame(cbind(background$indiv, background$forewing_band_valuer))
forewing_band_values <- subset(forewing_band_values, forewing_band_values$V2 != "missing")
write.table(forewing_band_values, "../gemma/forewing_band_values.csv", col.names = FALSE, row.names=FALSE,sep="\t", quote = FALSE)

#####WHITE####
#white  phenotype
background$white_valuer <- c()
for (i in 1:nrow(background)){
  if (background$hindwingWhite[i]=="absent"){
    background$white_valuer[i] <- 1}
  if (background$hindwingWhite[i]=="partial"){
    background$white_valuer[i] <- 2}
  if (background$hindwingWhite[i]=="present"){
    background$white_valuer[i] <- 3}
  if (background$hindwingWhite[i]=="missing"){
    background$white_valuer[i] <- "missing"}
}

#write file
white_values <- as.data.frame(cbind(background$indiv, background$white_valuer))
white_values <- subset(white_values, white_values$V2 != "missing")
write.table(white_values, "../gemma/white_values.csv", col.names = FALSE, row.names=FALSE,sep="\t", quote = FALSE)
