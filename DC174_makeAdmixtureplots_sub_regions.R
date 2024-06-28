#set directory and load background file
setwd("~/Dropbox/RishiMac2/Danaus/DC174/admixture/regions/")
background <- read.csv("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174//Background/DC174_background.csv", header = T)

###########################part 3: region1.1####
fam <- read.csv("admix_in_chr15_1.1.fam", sep = " ", header = FALSE)

background$numb <- 1:174

background$country_numb <- c()
for(i in 1:nrow(background)){
  if(background$country[i] == "Fuerteventura"){
    background$country_numb[i] <- 1
  }
  if(background$country[i] == "Italy"){
    background$country_numb[i] <- 2
  }
  if(background$country[i] == "Tunisia"){
    background$country_numb[i] <- 3
  }
  if(background$country[i] == "Nigeria"){
    background$country_numb[i] <- 4
  }
  if(background$country[i] == "Ghana"){
    background$country_numb[i] <- 5
  }
  if(background$country[i] == "Kenya"){
    background$country_numb[i] <- 6
  }
  if(background$country[i] == "Rwanda"){
    background$country_numb[i] <- 7
  }
  if(background$country[i] == "RSA"){
    background$country_numb[i] <- 8
  }
  if(background$country[i] == "StHelena"){
    background$country_numb[i] <- 9
  }
  if(background$country[i] == "Philippines"){
    background$country_numb[i] <- 10
  }
}

background_ordered <- background[order(background$country_numb),]
background_ordered$newnumb <- 1:174
background_newordered <- background_ordered[order(background_ordered$numb),]

Fuerta_indivs <- subset(background_newordered, background_newordered$country == "Spain")
Italy_indivs <- subset(background_newordered, background_newordered$country == "Italy")
Tunisia_indivs <- subset(background_newordered, background_newordered$country == "Tunisia")
Nigeria_indivs <- subset(background_newordered, background_newordered$country == "Nigeria")
Ghana_indivs <- subset(background_newordered, background_newordered$country == "Ghana")
Kenya_indivs <- subset(background_newordered, background_newordered$country == "Kenya")
Rwanda_indivs <- subset(background_newordered, background_newordered$country == "Rwanda")
RSA_indivs <- subset(background_newordered, background_newordered$country == "RSA")
StHelena_indivs <- subset(background_newordered, background_newordered$country == "StHelena")
Philippines_indivs <- subset(background_newordered, background_newordered$country == "Philippines")
missing <- subset(background_newordered, background_newordered$country == "missing")

popsum <- sum(nrow(Fuerta_indivs), nrow(Italy_indivs), nrow(Tunisia_indivs), nrow(Nigeria_indivs),
    nrow(Ghana_indivs), nrow(Kenya_indivs), nrow(Rwanda_indivs), nrow(RSA_indivs), 
    nrow(StHelena_indivs), nrow(Philippines_indivs), nrow(missing))

ad2 <- read.table("admix_in_chr15_1.1.2.Q")
rownames(ad2) <- fam$V1
ad2$order <- background_newordered$newnumb
ad2_or <- ad2[order(ad2$V1),]
ad2_o <- ad2[order(ad2$order),]
ad2_of <- ad2_o[,1:2]
barplot(t(as.matrix(ad2_of)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2", cex.names = 0.3)

tiff("all_chr15_1.1_K2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad2_of_Fuerta <- ad2_of[which(rownames(ad2_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Fuerta)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Fuerta", cex.names = 0.3)
ad2_of_Italy <- ad2_of[which(rownames(ad2_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Italy)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Italy", cex.names = 0.3)
ad2_of_Tunisia <- ad2_of[which(rownames(ad2_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Tunisia)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Tunis", cex.names = 0.3)

ad2_of_Nigeria <- ad2_of[which(rownames(ad2_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Nigeria)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Nigera", cex.names = 0.3)
ad2_of_Ghana <- ad2_of[which(rownames(ad2_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Ghana)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Ghana", cex.names = 0.3)
ad2_of_Kenya <- ad2_of[which(rownames(ad2_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Kenya)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Kenya", cex.names = 0.3)
ad2_of_Rwanda <- ad2_of[which(rownames(ad2_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Rwanda)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Rwanda", cex.names = 0.3)
ad2_of_RSA <- ad2_of[which(rownames(ad2_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad2_of_RSA)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - RSA", cex.names = 0.3)

ad2_of_StHelena <- ad2_of[which(rownames(ad2_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad2_of_StHelena)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - StHelena", cex.names = 0.3)
ad2_of_Philippines <- ad2_of[which(rownames(ad2_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Philippines)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Phili", cex.names = 0.3)
dev.off()

ad3 <- read.table("admix_in_chr15_1.1.3.Q")
rownames(ad3) <- fam$V1
ad3$order <- background_newordered$newnumb
ad3_or <- ad3[order(ad3$V2),]
ad3_or <- ad3_or[order(ad3_or$V3),]
ad3_or <- ad3_or[order(ad3_or$V1),]
ad3_o <- ad3_or[order(ad3_or$order),]
ad3_of <- ad3_o[,1:3]
barplot(t(as.matrix(ad3_of)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3", cex.names = 0.3)

tiff("all_chr15_1.1_K3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad3_of_Fuerta <- ad3_of[which(rownames(ad3_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Fuerta", cex.names = 0.3)
ad3_of_Italy <- ad3_of[which(rownames(ad3_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Italy", cex.names = 0.3)
ad3_of_Tunisia <- ad3_of[which(rownames(ad3_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Tunis", cex.names = 0.3)

ad3_of_Nigeria <- ad3_of[which(rownames(ad3_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Nigera", cex.names = 0.3)
ad3_of_Ghana <- ad3_of[which(rownames(ad3_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Ghana", cex.names = 0.3)
ad3_of_Kenya <- ad3_of[which(rownames(ad3_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Kenya", cex.names = 0.3)
ad3_of_Rwanda <- ad3_of[which(rownames(ad3_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Rwanda", cex.names = 0.3)
ad3_of_RSA <- ad3_of[which(rownames(ad3_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad3_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - RSA", cex.names = 0.3)

ad3_of_StHelena <- ad3_of[which(rownames(ad3_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad3_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - StHelena", cex.names = 0.3)
ad3_of_Philippines <- ad3_of[which(rownames(ad3_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Phili", cex.names = 0.3)
dev.off()

ad4 <- read.table("admix_in_chr15_1.1.4.Q")
rownames(ad4) <- fam$V1
ad4$order <- background_newordered$newnumb
ad4_or <- ad4[order(ad4$V1),]
ad4_or <- ad4_or[order(ad4_or$V2),]
ad4_or <- ad4_or[order(ad4_or$V3),]
ad4_or <- ad4_or[order(ad4_or$V4),]

ad4_o <- ad4[order(ad4$order),]
ad4_of <- ad4_o[,1:4]
barplot(t(as.matrix(ad4_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4", cex.names = 0.3)

tiff("all_chr15_1.1_K4.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad4_of_Fuerta <- ad4_of[which(rownames(ad4_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Fuerta", cex.names = 0.3)
ad4_of_Italy <- ad4_of[which(rownames(ad4_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Italy", cex.names = 0.3)
ad4_of_Tunisia <- ad4_of[which(rownames(ad4_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Tunis", cex.names = 0.3)

ad4_of_Nigeria <- ad4_of[which(rownames(ad4_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Nigera", cex.names = 0.3)
ad4_of_Ghana <- ad4_of[which(rownames(ad4_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Ghana", cex.names = 0.3)
ad4_of_Kenya <- ad4_of[which(rownames(ad4_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Kenya", cex.names = 0.3)
ad4_of_Rwanda <- ad4_of[which(rownames(ad4_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Rwanda", cex.names = 0.3)
ad4_of_RSA <- ad4_of[which(rownames(ad4_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad4_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - RSA", cex.names = 0.3)

ad4_of_StHelena <- ad4_of[which(rownames(ad4_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad4_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - StHelena", cex.names = 0.3)
ad4_of_Philippines <- ad4_of[which(rownames(ad4_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Phili", cex.names = 0.3)
dev.off()

ad5 <- read.table("admix_in_chr15_1.1.5.Q")
rownames(ad5) <- fam$V1
ad5$order <- background_newordered$newnumb
ad5_or <- ad5[order(ad5$V1),]
ad5_o <- ad5[order(ad5$order),]
ad5_of <- ad5_o[,1:5]
barplot(t(as.matrix(ad5_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5", cex.names = 0.3)

tiff("all_chr15_1.1_K5.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad5_of_Fuerta <- ad5_of[which(rownames(ad5_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Fuerta", cex.names = 0.3)
ad5_of_Italy <- ad5_of[which(rownames(ad5_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Italy", cex.names = 0.3)
ad5_of_Tunisia <- ad5_of[which(rownames(ad5_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Tunis", cex.names = 0.3)

ad5_of_Nigeria <- ad5_of[which(rownames(ad5_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Nigera", cex.names = 0.3)
ad5_of_Ghana <- ad5_of[which(rownames(ad5_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Ghana", cex.names = 0.3)
ad5_of_Kenya <- ad5_of[which(rownames(ad5_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Kenya", cex.names = 0.3)
ad5_of_Rwanda <- ad5_of[which(rownames(ad5_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Rwanda", cex.names = 0.3)
ad5_of_RSA <- ad5_of[which(rownames(ad5_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad5_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - RSA", cex.names = 0.3)

ad5_of_StHelena <- ad5_of[which(rownames(ad5_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad5_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - StHelena", cex.names = 0.3)
ad5_of_Philippines <- ad5_of[which(rownames(ad5_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Phili", cex.names = 0.3)
dev.off()

ad6 <- read.table("admix_in_chr15_1.1.6.Q")
rownames(ad6) <- fam$V1
ad6$order <- background_newordered$newnumb
ad6_or <- ad6[order(ad6$V1),]
ad6_o <- ad6[order(ad6$order),]
ad6_of <- ad6_o[,1:6]
barplot(t(as.matrix(ad6_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6", cex.names = 0.3)


tiff("all_chr15_1.1_K6.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad6_of_Fuerta <- ad6_of[which(rownames(ad6_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Fuerta", cex.names = 0.3)
ad6_of_Italy <- ad6_of[which(rownames(ad6_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Italy", cex.names = 0.3)
ad6_of_Tunisia <- ad6_of[which(rownames(ad6_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Tunis", cex.names = 0.3)

ad6_of_Nigeria <- ad6_of[which(rownames(ad6_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Nigera", cex.names = 0.3)
ad6_of_Ghana <- ad6_of[which(rownames(ad6_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Ghana", cex.names = 0.3)
ad6_of_Kenya <- ad6_of[which(rownames(ad6_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Kenya", cex.names = 0.3)
ad6_of_Rwanda <- ad6_of[which(rownames(ad6_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Rwanda", cex.names = 0.3)
ad6_of_RSA <- ad6_of[which(rownames(ad6_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad6_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - RSA", cex.names = 0.3)

ad6_of_StHelena <- ad6_of[which(rownames(ad6_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad6_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - StHelena", cex.names = 0.3)
ad6_of_Philippines <- ad6_of[which(rownames(ad6_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Phili", cex.names = 0.3)
dev.off()

ad7 <- read.table("admix_in_chr15_1.1.7.Q")
rownames(ad7) <- fam$V1
ad7$order <- background_newordered$newnumb
ad7_or <- ad7[order(ad7$V3),]
ad7_or <- ad7_or[order(ad7_or$V2),]
ad7_or <- ad7_or[order(ad7_or$V1),]
ad7_o <- ad7[order(ad7$order),]
ad7_of <- ad7_o[,1:7]
barplot(t(as.matrix(ad7_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7", cex.names = 0.3)

tiff("all_chr15_1.1_K7.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad7_of_Fuerta <- ad7_of[which(rownames(ad7_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Fuerta", cex.names = 0.3)
ad7_of_Italy <- ad7_of[which(rownames(ad7_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Italy", cex.names = 0.3)
ad7_of_Tunisia <- ad7_of[which(rownames(ad7_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Tunis", cex.names = 0.3)

ad7_of_Nigeria <- ad7_of[which(rownames(ad7_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Nigera", cex.names = 0.3)
ad7_of_Ghana <- ad7_of[which(rownames(ad7_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Ghana", cex.names = 0.3)
ad7_of_Kenya <- ad7_of[which(rownames(ad7_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Kenya", cex.names = 0.3)
ad7_of_Rwanda <- ad7_of[which(rownames(ad7_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Rwanda", cex.names = 0.3)
ad7_of_RSA <- ad7_of[which(rownames(ad7_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad7_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - RSA", cex.names = 0.3)

ad7_of_StHelena <- ad7_of[which(rownames(ad7_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad7_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - StHelena", cex.names = 0.3)
ad7_of_Philippines <- ad7_of[which(rownames(ad7_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Phili", cex.names = 0.3)
dev.off()

ad8 <- read.table("admix_in_chr15_1.1.8.Q")
rownames(ad8) <- fam$V1
ad8$order <- background_newordered$newnumb
ad8_or <- ad8[order(ad8$V1),]
ad8_o <- ad8[order(ad8$order),]
ad8_of <- ad8_o[,1:8]

barplot(t(as.matrix(ad8_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8", cex.names = 0.3)

tiff("all_chr15_1.1_K8.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad8_of_Fuerta <- ad8_of[which(rownames(ad8_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Fuerta", cex.names = 0.3)
ad8_of_Italy <- ad8_of[which(rownames(ad8_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Italy", cex.names = 0.3)
ad8_of_Tunisia <- ad8_of[which(rownames(ad8_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Tunis", cex.names = 0.3)

ad8_of_Nigeria <- ad8_of[which(rownames(ad8_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Nigera", cex.names = 0.3)
ad8_of_Ghana <- ad8_of[which(rownames(ad8_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Ghana", cex.names = 0.3)
ad8_of_Kenya <- ad8_of[which(rownames(ad8_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Kenya", cex.names = 0.3)
ad8_of_Rwanda <- ad8_of[which(rownames(ad8_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Rwanda", cex.names = 0.3)
ad8_of_RSA <- ad8_of[which(rownames(ad8_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad8_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - RSA", cex.names = 0.3)

ad8_of_StHelena <- ad8_of[which(rownames(ad8_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad8_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - StHelena", cex.names = 0.3)
ad8_of_Philippines <- ad8_of[which(rownames(ad8_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Phili", cex.names = 0.3)
dev.off()

ad9 <- read.table("admix_chr15_in.9.Q")
rownames(ad9) <- fam$V1
ad9$order <- background_newordered$newnumb
ad9_or <- ad9[order(ad9$V1),]
ad9_o <- ad9[order(ad9$order),]
ad9_of <- ad9_o[,1:9]

barplot(t(as.matrix(ad9_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9", cex.names = 0.3)

tiff("all_chr15_1.1_K9.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad9_of_Fuerta <- ad9_of[which(rownames(ad9_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Fuerta", cex.names = 0.3)
ad9_of_Italy <- ad9_of[which(rownames(ad9_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Italy", cex.names = 0.3)
ad9_of_Tunisia <- ad9_of[which(rownames(ad9_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Tunis", cex.names = 0.3)

ad9_of_Nigeria <- ad9_of[which(rownames(ad9_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Nigera", cex.names = 0.3)
ad9_of_Ghana <- ad9_of[which(rownames(ad9_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Ghana", cex.names = 0.3)
ad9_of_Kenya <- ad9_of[which(rownames(ad9_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Kenya", cex.names = 0.3)
ad9_of_Rwanda <- ad9_of[which(rownames(ad9_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Rwanda", cex.names = 0.3)
ad9_of_RSA <- ad9_of[which(rownames(ad9_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad9_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - RSA", cex.names = 0.3)

ad9_of_StHelena <- ad9_of[which(rownames(ad9_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad9_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - StHelena", cex.names = 0.3)
ad9_of_Philippines <- ad9_of[which(rownames(ad9_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Phili", cex.names = 0.3)
dev.off()

ad10 <- read.table("admix_in_chr15_1.1.10.Q")
rownames(ad10) <- fam$V1
ad10$order <- background_newordered$newnumb
ad10_or <- ad10[order(ad10$V3),]
ad10_or <- ad10_or[order(ad10_or$V2),]
ad10_or <- ad10_or[order(ad10_or$V1),]
ad10_o <- ad10[order(ad10$order),]
ad10_of <- ad10_o[,1:10]

tiff("all_chr15_1.1_K10.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad10_of_Fuerta <- ad10_of[which(rownames(ad10_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Fuerta", cex.names = 0.3)
ad10_of_Italy <- ad10_of[which(rownames(ad10_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Italy", cex.names = 0.3)
ad10_of_Tunisia <- ad10_of[which(rownames(ad10_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Tunis", cex.names = 0.3)

ad10_of_Nigeria <- ad10_of[which(rownames(ad10_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Nigera", cex.names = 0.3)
ad10_of_Ghana <- ad10_of[which(rownames(ad10_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Ghana", cex.names = 0.3)
ad10_of_Kenya <- ad10_of[which(rownames(ad10_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Kenya", cex.names = 0.3)
ad10_of_Rwanda <- ad10_of[which(rownames(ad10_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Rwanda", cex.names = 0.3)
ad10_of_RSA <- ad10_of[which(rownames(ad10_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad10_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - RSA", cex.names = 0.3)

ad10_of_StHelena <- ad10_of[which(rownames(ad10_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad10_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - StHelena", cex.names = 0.3)
ad10_of_Philippines <- ad10_of[which(rownames(ad10_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Phili", cex.names = 0.3)
dev.off()

# cv error plot:
k_vals <- c(2:10)
cv_vals <- c(0.21159, 0.20095, 0.19305, 0.18889, 0.17838, 0.17958, 0.19436, 0.18897, 0.19431)

plot(k_vals, cv_vals, type = 'l', lty = 2, ylab = "CV error", xlab = "K", main = "reg 1.1 - best k=6")
points(k_vals, cv_vals, pch = 16, col= "black")
points(6, 0.17838, pch = 16, col= "red")


###########################part 4: region1.2######

ad2 <- read.table("admix_in_chr15_1.2.2.Q")
rownames(ad2) <- fam$V1
ad2$order <- background_newordered$newnumb
ad2_or <- ad2[order(ad2$V1),]
ad2_o <- ad2[order(ad2$order),]
ad2_of <- ad2_o[,1:2]
barplot(t(as.matrix(ad2_of)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2", cex.names = 0.3)

tiff("all_chr15_1.2_K2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad2_of_Fuerta <- ad2_of[which(rownames(ad2_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Fuerta)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Fuerta", cex.names = 0.3)
ad2_of_Italy <- ad2_of[which(rownames(ad2_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Italy)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Italy", cex.names = 0.3)
ad2_of_Tunisia <- ad2_of[which(rownames(ad2_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Tunisia)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Tunis", cex.names = 0.3)

ad2_of_Nigeria <- ad2_of[which(rownames(ad2_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Nigeria)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Nigera", cex.names = 0.3)
ad2_of_Ghana <- ad2_of[which(rownames(ad2_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Ghana)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Ghana", cex.names = 0.3)
ad2_of_Kenya <- ad2_of[which(rownames(ad2_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Kenya)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Kenya", cex.names = 0.3)
ad2_of_Rwanda <- ad2_of[which(rownames(ad2_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Rwanda)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Rwanda", cex.names = 0.3)
ad2_of_RSA <- ad2_of[which(rownames(ad2_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad2_of_RSA)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - RSA", cex.names = 0.3)

ad2_of_StHelena <- ad2_of[which(rownames(ad2_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad2_of_StHelena)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - StHelena", cex.names = 0.3)
ad2_of_Philippines <- ad2_of[which(rownames(ad2_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Philippines)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Phili", cex.names = 0.3)
dev.off()

ad3 <- read.table("admix_in_chr15_1.2.3.Q")
rownames(ad3) <- fam$V1
ad3$order <- background_newordered$newnumb
ad3_or <- ad3[order(ad3$V2),]
ad3_or <- ad3_or[order(ad3_or$V3),]
ad3_or <- ad3_or[order(ad3_or$V1),]
ad3_o <- ad3_or[order(ad3_or$order),]
ad3_of <- ad3_o[,1:3]
barplot(t(as.matrix(ad3_of)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3", cex.names = 0.3)

tiff("all_chr15_1.2_K3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad3_of_Fuerta <- ad3_of[which(rownames(ad3_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Fuerta", cex.names = 0.3)
ad3_of_Italy <- ad3_of[which(rownames(ad3_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Italy", cex.names = 0.3)
ad3_of_Tunisia <- ad3_of[which(rownames(ad3_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Tunis", cex.names = 0.3)

ad3_of_Nigeria <- ad3_of[which(rownames(ad3_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Nigera", cex.names = 0.3)
ad3_of_Ghana <- ad3_of[which(rownames(ad3_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Ghana", cex.names = 0.3)
ad3_of_Kenya <- ad3_of[which(rownames(ad3_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Kenya", cex.names = 0.3)
ad3_of_Rwanda <- ad3_of[which(rownames(ad3_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Rwanda", cex.names = 0.3)
ad3_of_RSA <- ad3_of[which(rownames(ad3_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad3_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - RSA", cex.names = 0.3)

ad3_of_StHelena <- ad3_of[which(rownames(ad3_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad3_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - StHelena", cex.names = 0.3)
ad3_of_Philippines <- ad3_of[which(rownames(ad3_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Phili", cex.names = 0.3)
dev.off()

ad4 <- read.table("admix_in_chr15_1.2.4.Q")
rownames(ad4) <- fam$V1
ad4$order <- background_newordered$newnumb
ad4_or <- ad4[order(ad4$V1),]
ad4_or <- ad4_or[order(ad4_or$V2),]
ad4_or <- ad4_or[order(ad4_or$V3),]
ad4_or <- ad4_or[order(ad4_or$V4),]

ad4_o <- ad4[order(ad4$order),]
ad4_of <- ad4_o[,1:4]
barplot(t(as.matrix(ad4_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4", cex.names = 0.3)

tiff("all_chr15_1.2_K4.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad4_of_Fuerta <- ad4_of[which(rownames(ad4_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Fuerta", cex.names = 0.3)
ad4_of_Italy <- ad4_of[which(rownames(ad4_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Italy", cex.names = 0.3)
ad4_of_Tunisia <- ad4_of[which(rownames(ad4_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Tunis", cex.names = 0.3)

ad4_of_Nigeria <- ad4_of[which(rownames(ad4_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Nigera", cex.names = 0.3)
ad4_of_Ghana <- ad4_of[which(rownames(ad4_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Ghana", cex.names = 0.3)
ad4_of_Kenya <- ad4_of[which(rownames(ad4_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Kenya", cex.names = 0.3)
ad4_of_Rwanda <- ad4_of[which(rownames(ad4_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Rwanda", cex.names = 0.3)
ad4_of_RSA <- ad4_of[which(rownames(ad4_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad4_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - RSA", cex.names = 0.3)

ad4_of_StHelena <- ad4_of[which(rownames(ad4_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad4_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - StHelena", cex.names = 0.3)
ad4_of_Philippines <- ad4_of[which(rownames(ad4_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Phili", cex.names = 0.3)
dev.off()

ad5 <- read.table("admix_in_chr15_1.2.5.Q")
rownames(ad5) <- fam$V1
ad5$order <- background_newordered$newnumb
ad5_or <- ad5[order(ad5$V1),]
ad5_o <- ad5[order(ad5$order),]
ad5_of <- ad5_o[,1:5]
barplot(t(as.matrix(ad5_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5", cex.names = 0.3)

tiff("all_chr15_1.2_K5.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad5_of_Fuerta <- ad5_of[which(rownames(ad5_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Fuerta", cex.names = 0.3)
ad5_of_Italy <- ad5_of[which(rownames(ad5_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Italy", cex.names = 0.3)
ad5_of_Tunisia <- ad5_of[which(rownames(ad5_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Tunis", cex.names = 0.3)

ad5_of_Nigeria <- ad5_of[which(rownames(ad5_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Nigera", cex.names = 0.3)
ad5_of_Ghana <- ad5_of[which(rownames(ad5_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Ghana", cex.names = 0.3)
ad5_of_Kenya <- ad5_of[which(rownames(ad5_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Kenya", cex.names = 0.3)
ad5_of_Rwanda <- ad5_of[which(rownames(ad5_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Rwanda", cex.names = 0.3)
ad5_of_RSA <- ad5_of[which(rownames(ad5_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad5_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - RSA", cex.names = 0.3)

ad5_of_StHelena <- ad5_of[which(rownames(ad5_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad5_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - StHelena", cex.names = 0.3)
ad5_of_Philippines <- ad5_of[which(rownames(ad5_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Phili", cex.names = 0.3)
dev.off()

ad6 <- read.table("admix_in_chr15_1.2.6.Q")
rownames(ad6) <- fam$V1
ad6$order <- background_newordered$newnumb
ad6_or <- ad6[order(ad6$V1),]
ad6_o <- ad6[order(ad6$order),]
ad6_of <- ad6_o[,1:6]
barplot(t(as.matrix(ad6_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6", cex.names = 0.3)


tiff("all_chr15_1.2_K6.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad6_of_Fuerta <- ad6_of[which(rownames(ad6_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Fuerta", cex.names = 0.3)
ad6_of_Italy <- ad6_of[which(rownames(ad6_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Italy", cex.names = 0.3)
ad6_of_Tunisia <- ad6_of[which(rownames(ad6_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Tunis", cex.names = 0.3)

ad6_of_Nigeria <- ad6_of[which(rownames(ad6_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Nigera", cex.names = 0.3)
ad6_of_Ghana <- ad6_of[which(rownames(ad6_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Ghana", cex.names = 0.3)
ad6_of_Kenya <- ad6_of[which(rownames(ad6_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Kenya", cex.names = 0.3)
ad6_of_Rwanda <- ad6_of[which(rownames(ad6_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Rwanda", cex.names = 0.3)
ad6_of_RSA <- ad6_of[which(rownames(ad6_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad6_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - RSA", cex.names = 0.3)

ad6_of_StHelena <- ad6_of[which(rownames(ad6_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad6_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - StHelena", cex.names = 0.3)
ad6_of_Philippines <- ad6_of[which(rownames(ad6_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Phili", cex.names = 0.3)
dev.off()

ad7 <- read.table("admix_in_chr15_1.2.7.Q")
rownames(ad7) <- fam$V1
ad7$order <- background_newordered$newnumb
ad7_or <- ad7[order(ad7$V3),]
ad7_or <- ad7_or[order(ad7_or$V2),]
ad7_or <- ad7_or[order(ad7_or$V1),]
ad7_o <- ad7[order(ad7$order),]
ad7_of <- ad7_o[,1:7]
barplot(t(as.matrix(ad7_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7", cex.names = 0.3)

tiff("all_chr15_1.2_K7.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad7_of_Fuerta <- ad7_of[which(rownames(ad7_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Fuerta", cex.names = 0.3)
ad7_of_Italy <- ad7_of[which(rownames(ad7_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Italy", cex.names = 0.3)
ad7_of_Tunisia <- ad7_of[which(rownames(ad7_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Tunis", cex.names = 0.3)

ad7_of_Nigeria <- ad7_of[which(rownames(ad7_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Nigera", cex.names = 0.3)
ad7_of_Ghana <- ad7_of[which(rownames(ad7_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Ghana", cex.names = 0.3)
ad7_of_Kenya <- ad7_of[which(rownames(ad7_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Kenya", cex.names = 0.3)
ad7_of_Rwanda <- ad7_of[which(rownames(ad7_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Rwanda", cex.names = 0.3)
ad7_of_RSA <- ad7_of[which(rownames(ad7_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad7_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - RSA", cex.names = 0.3)

ad7_of_StHelena <- ad7_of[which(rownames(ad7_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad7_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - StHelena", cex.names = 0.3)
ad7_of_Philippines <- ad7_of[which(rownames(ad7_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Phili", cex.names = 0.3)
dev.off()

ad8 <- read.table("admix_in_chr15_1.2.8.Q")
rownames(ad8) <- fam$V1
ad8$order <- background_newordered$newnumb
ad8_or <- ad8[order(ad8$V1),]
ad8_o <- ad8[order(ad8$order),]
ad8_of <- ad8_o[,1:8]

barplot(t(as.matrix(ad8_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8", cex.names = 0.3)

tiff("all_chr15_1.2_K8.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad8_of_Fuerta <- ad8_of[which(rownames(ad8_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Fuerta", cex.names = 0.3)
ad8_of_Italy <- ad8_of[which(rownames(ad8_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Italy", cex.names = 0.3)
ad8_of_Tunisia <- ad8_of[which(rownames(ad8_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Tunis", cex.names = 0.3)

ad8_of_Nigeria <- ad8_of[which(rownames(ad8_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Nigera", cex.names = 0.3)
ad8_of_Ghana <- ad8_of[which(rownames(ad8_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Ghana", cex.names = 0.3)
ad8_of_Kenya <- ad8_of[which(rownames(ad8_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Kenya", cex.names = 0.3)
ad8_of_Rwanda <- ad8_of[which(rownames(ad8_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Rwanda", cex.names = 0.3)
ad8_of_RSA <- ad8_of[which(rownames(ad8_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad8_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - RSA", cex.names = 0.3)

ad8_of_StHelena <- ad8_of[which(rownames(ad8_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad8_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - StHelena", cex.names = 0.3)
ad8_of_Philippines <- ad8_of[which(rownames(ad8_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Phili", cex.names = 0.3)
dev.off()

ad9 <- read.table("admix_in_chr15_1.2.9.Q")
rownames(ad9) <- fam$V1
ad9$order <- background_newordered$newnumb
ad9_or <- ad9[order(ad9$V1),]
ad9_o <- ad9[order(ad9$order),]
ad9_of <- ad9_o[,1:9]

barplot(t(as.matrix(ad9_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9", cex.names = 0.3)

tiff("all_chr15_1.2_K9.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad9_of_Fuerta <- ad9_of[which(rownames(ad9_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Fuerta", cex.names = 0.3)
ad9_of_Italy <- ad9_of[which(rownames(ad9_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Italy", cex.names = 0.3)
ad9_of_Tunisia <- ad9_of[which(rownames(ad9_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Tunis", cex.names = 0.3)

ad9_of_Nigeria <- ad9_of[which(rownames(ad9_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Nigera", cex.names = 0.3)
ad9_of_Ghana <- ad9_of[which(rownames(ad9_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Ghana", cex.names = 0.3)
ad9_of_Kenya <- ad9_of[which(rownames(ad9_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Kenya", cex.names = 0.3)
ad9_of_Rwanda <- ad9_of[which(rownames(ad9_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Rwanda", cex.names = 0.3)
ad9_of_RSA <- ad9_of[which(rownames(ad9_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad9_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - RSA", cex.names = 0.3)

ad9_of_StHelena <- ad9_of[which(rownames(ad9_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad9_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - StHelena", cex.names = 0.3)
ad9_of_Philippines <- ad9_of[which(rownames(ad9_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Phili", cex.names = 0.3)
dev.off()

ad10 <- read.table("admix_in_chr15_1.2.10.Q")
rownames(ad10) <- fam$V1
ad10$order <- background_newordered$newnumb
ad10_or <- ad10[order(ad10$V3),]
ad10_or <- ad10_or[order(ad10_or$V2),]
ad10_or <- ad10_or[order(ad10_or$V1),]
ad10_o <- ad10[order(ad10$order),]
ad10_of <- ad10_o[,1:10]

tiff("all_chr15_1.2_K10.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad10_of_Fuerta <- ad10_of[which(rownames(ad10_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Fuerta", cex.names = 0.3)
ad10_of_Italy <- ad10_of[which(rownames(ad10_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Italy", cex.names = 0.3)
ad10_of_Tunisia <- ad10_of[which(rownames(ad10_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Tunis", cex.names = 0.3)

ad10_of_Nigeria <- ad10_of[which(rownames(ad10_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Nigera", cex.names = 0.3)
ad10_of_Ghana <- ad10_of[which(rownames(ad10_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Ghana", cex.names = 0.3)
ad10_of_Kenya <- ad10_of[which(rownames(ad10_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Kenya", cex.names = 0.3)
ad10_of_Rwanda <- ad10_of[which(rownames(ad10_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Rwanda", cex.names = 0.3)
ad10_of_RSA <- ad10_of[which(rownames(ad10_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad10_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - RSA", cex.names = 0.3)

ad10_of_StHelena <- ad10_of[which(rownames(ad10_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad10_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - StHelena", cex.names = 0.3)
ad10_of_Philippines <- ad10_of[which(rownames(ad10_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Phili", cex.names = 0.3)
dev.off()


# cv error plot:
k_vals <- c(2:10)
cv_vals <- c(0.22629, 0.15140, 0.12404, 0.11721, 0.10469, 0.10818, 0.10502, 0.10656, 0.10820)
plot(k_vals, cv_vals, type = 'l', lty = 2, ylab = "CV error", xlab = "K", main = "reg 1.2 - best k=6")
points(k_vals, cv_vals, pch = 16, col= "black")
points(6, 0.10469, pch = 16, col= "red")


###########################part 5: region2)######
ad2 <- read.table("admix_in_chr15_2.2.Q")
rownames(ad2) <- fam$V1
ad2$order <- background_newordered$newnumb
ad2_or <- ad2[order(ad2$V1),]
ad2_o <- ad2[order(ad2$order),]
ad2_of <- ad2_o[,1:2]
barplot(t(as.matrix(ad2_of)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2", cex.names = 0.3)

tiff("all_chr15_2_K2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad2_of_Fuerta <- ad2_of[which(rownames(ad2_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Fuerta)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Fuerta", cex.names = 0.3)
ad2_of_Italy <- ad2_of[which(rownames(ad2_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Italy)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Italy", cex.names = 0.3)
ad2_of_Tunisia <- ad2_of[which(rownames(ad2_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Tunisia)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Tunis", cex.names = 0.3)

ad2_of_Nigeria <- ad2_of[which(rownames(ad2_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Nigeria)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Nigera", cex.names = 0.3)
ad2_of_Ghana <- ad2_of[which(rownames(ad2_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Ghana)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Ghana", cex.names = 0.3)
ad2_of_Kenya <- ad2_of[which(rownames(ad2_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Kenya)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Kenya", cex.names = 0.3)
ad2_of_Rwanda <- ad2_of[which(rownames(ad2_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Rwanda)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Rwanda", cex.names = 0.3)
ad2_of_RSA <- ad2_of[which(rownames(ad2_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad2_of_RSA)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - RSA", cex.names = 0.3)

ad2_of_StHelena <- ad2_of[which(rownames(ad2_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad2_of_StHelena)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - StHelena", cex.names = 0.3)
ad2_of_Philippines <- ad2_of[which(rownames(ad2_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Philippines)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Phili", cex.names = 0.3)
dev.off()

ad3 <- read.table("admix_in_chr15_2.3.Q")
rownames(ad3) <- fam$V1
ad3$order <- background_newordered$newnumb
ad3_or <- ad3[order(ad3$V2),]
ad3_or <- ad3_or[order(ad3_or$V3),]
ad3_or <- ad3_or[order(ad3_or$V1),]
ad3_o <- ad3_or[order(ad3_or$order),]
ad3_of <- ad3_o[,1:3]
barplot(t(as.matrix(ad3_of)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3", cex.names = 0.3)

tiff("all_chr15_2_K3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad3_of_Fuerta <- ad3_of[which(rownames(ad3_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Fuerta", cex.names = 0.3)
ad3_of_Italy <- ad3_of[which(rownames(ad3_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Italy", cex.names = 0.3)
ad3_of_Tunisia <- ad3_of[which(rownames(ad3_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Tunis", cex.names = 0.3)

ad3_of_Nigeria <- ad3_of[which(rownames(ad3_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Nigera", cex.names = 0.3)
ad3_of_Ghana <- ad3_of[which(rownames(ad3_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Ghana", cex.names = 0.3)
ad3_of_Kenya <- ad3_of[which(rownames(ad3_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Kenya", cex.names = 0.3)
ad3_of_Rwanda <- ad3_of[which(rownames(ad3_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Rwanda", cex.names = 0.3)
ad3_of_RSA <- ad3_of[which(rownames(ad3_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad3_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - RSA", cex.names = 0.3)

ad3_of_StHelena <- ad3_of[which(rownames(ad3_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad3_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - StHelena", cex.names = 0.3)
ad3_of_Philippines <- ad3_of[which(rownames(ad3_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Phili", cex.names = 0.3)
dev.off()

ad4 <- read.table("admix_in_chr15_2.4.Q")
rownames(ad4) <- fam$V1
ad4$order <- background_newordered$newnumb
ad4_or <- ad4[order(ad4$V1),]
ad4_or <- ad4_or[order(ad4_or$V2),]
ad4_or <- ad4_or[order(ad4_or$V3),]
ad4_or <- ad4_or[order(ad4_or$V4),]

ad4_o <- ad4[order(ad4$order),]
ad4_of <- ad4_o[,1:4]
barplot(t(as.matrix(ad4_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4", cex.names = 0.3)

tiff("all_chr15_2_K4.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad4_of_Fuerta <- ad4_of[which(rownames(ad4_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Fuerta", cex.names = 0.3)
ad4_of_Italy <- ad4_of[which(rownames(ad4_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Italy", cex.names = 0.3)
ad4_of_Tunisia <- ad4_of[which(rownames(ad4_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Tunis", cex.names = 0.3)

ad4_of_Nigeria <- ad4_of[which(rownames(ad4_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Nigera", cex.names = 0.3)
ad4_of_Ghana <- ad4_of[which(rownames(ad4_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Ghana", cex.names = 0.3)
ad4_of_Kenya <- ad4_of[which(rownames(ad4_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Kenya", cex.names = 0.3)
ad4_of_Rwanda <- ad4_of[which(rownames(ad4_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Rwanda", cex.names = 0.3)
ad4_of_RSA <- ad4_of[which(rownames(ad4_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad4_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - RSA", cex.names = 0.3)

ad4_of_StHelena <- ad4_of[which(rownames(ad4_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad4_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - StHelena", cex.names = 0.3)
ad4_of_Philippines <- ad4_of[which(rownames(ad4_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Phili", cex.names = 0.3)
dev.off()

ad5 <- read.table("admix_in_chr15_2.5.Q")
rownames(ad5) <- fam$V1
ad5$order <- background_newordered$newnumb
ad5_or <- ad5[order(ad5$V1),]
ad5_o <- ad5[order(ad5$order),]
ad5_of <- ad5_o[,1:5]
barplot(t(as.matrix(ad5_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5", cex.names = 0.3)

tiff("all_chr15_2_K5.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad5_of_Fuerta <- ad5_of[which(rownames(ad5_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Fuerta", cex.names = 0.3)
ad5_of_Italy <- ad5_of[which(rownames(ad5_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Italy", cex.names = 0.3)
ad5_of_Tunisia <- ad5_of[which(rownames(ad5_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Tunis", cex.names = 0.3)

ad5_of_Nigeria <- ad5_of[which(rownames(ad5_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Nigera", cex.names = 0.3)
ad5_of_Ghana <- ad5_of[which(rownames(ad5_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Ghana", cex.names = 0.3)
ad5_of_Kenya <- ad5_of[which(rownames(ad5_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Kenya", cex.names = 0.3)
ad5_of_Rwanda <- ad5_of[which(rownames(ad5_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Rwanda", cex.names = 0.3)
ad5_of_RSA <- ad5_of[which(rownames(ad5_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad5_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - RSA", cex.names = 0.3)

ad5_of_StHelena <- ad5_of[which(rownames(ad5_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad5_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - StHelena", cex.names = 0.3)
ad5_of_Philippines <- ad5_of[which(rownames(ad5_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Phili", cex.names = 0.3)
dev.off()

ad6 <- read.table("admix_in_chr15_2.6.Q")
rownames(ad6) <- fam$V1
ad6$order <- background_newordered$newnumb
ad6_or <- ad6[order(ad6$V1),]
ad6_o <- ad6[order(ad6$order),]
ad6_of <- ad6_o[,1:6]
barplot(t(as.matrix(ad6_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6", cex.names = 0.3)


tiff("all_chr15_2_K6.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad6_of_Fuerta <- ad6_of[which(rownames(ad6_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Fuerta", cex.names = 0.3)
ad6_of_Italy <- ad6_of[which(rownames(ad6_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Italy", cex.names = 0.3)
ad6_of_Tunisia <- ad6_of[which(rownames(ad6_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Tunis", cex.names = 0.3)

ad6_of_Nigeria <- ad6_of[which(rownames(ad6_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Nigera", cex.names = 0.3)
ad6_of_Ghana <- ad6_of[which(rownames(ad6_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Ghana", cex.names = 0.3)
ad6_of_Kenya <- ad6_of[which(rownames(ad6_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Kenya", cex.names = 0.3)
ad6_of_Rwanda <- ad6_of[which(rownames(ad6_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Rwanda", cex.names = 0.3)
ad6_of_RSA <- ad6_of[which(rownames(ad6_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad6_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - RSA", cex.names = 0.3)

ad6_of_StHelena <- ad6_of[which(rownames(ad6_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad6_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - StHelena", cex.names = 0.3)
ad6_of_Philippines <- ad6_of[which(rownames(ad6_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Phili", cex.names = 0.3)
dev.off()

ad7 <- read.table("admix_in_chr15_2.7.Q")
rownames(ad7) <- fam$V1
ad7$order <- background_newordered$newnumb
ad7_or <- ad7[order(ad7$V3),]
ad7_or <- ad7_or[order(ad7_or$V2),]
ad7_or <- ad7_or[order(ad7_or$V1),]
ad7_o <- ad7[order(ad7$order),]
ad7_of <- ad7_o[,1:7]
barplot(t(as.matrix(ad7_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7", cex.names = 0.3)

tiff("all_chr15_2_K7.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad7_of_Fuerta <- ad7_of[which(rownames(ad7_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Fuerta", cex.names = 0.3)
ad7_of_Italy <- ad7_of[which(rownames(ad7_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Italy", cex.names = 0.3)
ad7_of_Tunisia <- ad7_of[which(rownames(ad7_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Tunis", cex.names = 0.3)

ad7_of_Nigeria <- ad7_of[which(rownames(ad7_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Nigera", cex.names = 0.3)
ad7_of_Ghana <- ad7_of[which(rownames(ad7_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Ghana", cex.names = 0.3)
ad7_of_Kenya <- ad7_of[which(rownames(ad7_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Kenya", cex.names = 0.3)
ad7_of_Rwanda <- ad7_of[which(rownames(ad7_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Rwanda", cex.names = 0.3)
ad7_of_RSA <- ad7_of[which(rownames(ad7_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad7_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - RSA", cex.names = 0.3)

ad7_of_StHelena <- ad7_of[which(rownames(ad7_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad7_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - StHelena", cex.names = 0.3)
ad7_of_Philippines <- ad7_of[which(rownames(ad7_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Phili", cex.names = 0.3)
dev.off()

ad8 <- read.table("admix_in_chr15_2.8.Q")
rownames(ad8) <- fam$V1
ad8$order <- background_newordered$newnumb
ad8_or <- ad8[order(ad8$V1),]
ad8_o <- ad8[order(ad8$order),]
ad8_of <- ad8_o[,1:8]

barplot(t(as.matrix(ad8_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8", cex.names = 0.3)

tiff("all_chr15_2_K8.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad8_of_Fuerta <- ad8_of[which(rownames(ad8_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Fuerta", cex.names = 0.3)
ad8_of_Italy <- ad8_of[which(rownames(ad8_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Italy", cex.names = 0.3)
ad8_of_Tunisia <- ad8_of[which(rownames(ad8_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Tunis", cex.names = 0.3)

ad8_of_Nigeria <- ad8_of[which(rownames(ad8_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Nigera", cex.names = 0.3)
ad8_of_Ghana <- ad8_of[which(rownames(ad8_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Ghana", cex.names = 0.3)
ad8_of_Kenya <- ad8_of[which(rownames(ad8_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Kenya", cex.names = 0.3)
ad8_of_Rwanda <- ad8_of[which(rownames(ad8_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Rwanda", cex.names = 0.3)
ad8_of_RSA <- ad8_of[which(rownames(ad8_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad8_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - RSA", cex.names = 0.3)

ad8_of_StHelena <- ad8_of[which(rownames(ad8_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad8_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - StHelena", cex.names = 0.3)
ad8_of_Philippines <- ad8_of[which(rownames(ad8_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Phili", cex.names = 0.3)
dev.off()

ad9 <- read.table("admix_in_chr15_2.9.Q")
rownames(ad9) <- fam$V1
ad9$order <- background_newordered$newnumb
ad9_or <- ad9[order(ad9$V1),]
ad9_o <- ad9[order(ad9$order),]
ad9_of <- ad9_o[,1:9]

barplot(t(as.matrix(ad9_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9", cex.names = 0.3)

tiff("all_chr15_2_K9.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad9_of_Fuerta <- ad9_of[which(rownames(ad9_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Fuerta", cex.names = 0.3)
ad9_of_Italy <- ad9_of[which(rownames(ad9_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Italy", cex.names = 0.3)
ad9_of_Tunisia <- ad9_of[which(rownames(ad9_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Tunis", cex.names = 0.3)

ad9_of_Nigeria <- ad9_of[which(rownames(ad9_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Nigera", cex.names = 0.3)
ad9_of_Ghana <- ad9_of[which(rownames(ad9_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Ghana", cex.names = 0.3)
ad9_of_Kenya <- ad9_of[which(rownames(ad9_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Kenya", cex.names = 0.3)
ad9_of_Rwanda <- ad9_of[which(rownames(ad9_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Rwanda", cex.names = 0.3)
ad9_of_RSA <- ad9_of[which(rownames(ad9_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad9_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - RSA", cex.names = 0.3)

ad9_of_StHelena <- ad9_of[which(rownames(ad9_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad9_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - StHelena", cex.names = 0.3)
ad9_of_Philippines <- ad9_of[which(rownames(ad9_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Phili", cex.names = 0.3)
dev.off()

ad10 <- read.table("admix_in_chr15_2.10.Q")
rownames(ad10) <- fam$V1
ad10$order <- background_newordered$newnumb
ad10_or <- ad10[order(ad10$V3),]
ad10_or <- ad10_or[order(ad10_or$V2),]
ad10_or <- ad10_or[order(ad10_or$V1),]
ad10_o <- ad10[order(ad10$order),]
ad10_of <- ad10_o[,1:10]

tiff("all_chr15_2_K10.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad10_of_Fuerta <- ad10_of[which(rownames(ad10_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Fuerta", cex.names = 0.3)
ad10_of_Italy <- ad10_of[which(rownames(ad10_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Italy", cex.names = 0.3)
ad10_of_Tunisia <- ad10_of[which(rownames(ad10_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Tunis", cex.names = 0.3)

ad10_of_Nigeria <- ad10_of[which(rownames(ad10_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Nigera", cex.names = 0.3)
ad10_of_Ghana <- ad10_of[which(rownames(ad10_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Ghana", cex.names = 0.3)
ad10_of_Kenya <- ad10_of[which(rownames(ad10_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Kenya", cex.names = 0.3)
ad10_of_Rwanda <- ad10_of[which(rownames(ad10_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Rwanda", cex.names = 0.3)
ad10_of_RSA <- ad10_of[which(rownames(ad10_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad10_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - RSA", cex.names = 0.3)

ad10_of_StHelena <- ad10_of[which(rownames(ad10_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad10_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - StHelena", cex.names = 0.3)
ad10_of_Philippines <- ad10_of[which(rownames(ad10_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Phili", cex.names = 0.3)
dev.off()

# cv error plot:
k_vals <- c(2:10)
cv_vals <- c(0.24583, 0.16874, 0.14242, 0.13033, 0.12113, 0.12320, 0.12075, 0.12064, 0.13057)
plot(k_vals, cv_vals, type = 'l', lty = 2, ylab = "CV error", xlab = "K", main = "reg 2 - best k=9")
points(k_vals, cv_vals, pch = 16, col= "black")
points(9, 0.12064, pch = 16, col= "red")

###########################part 6: region4)######
ad2 <- read.table("admix_in_chr15_4.2.Q")
rownames(ad2) <- fam$V1
ad2$order <- background_newordered$newnumb
ad2_or <- ad2[order(ad2$V1),]
ad2_o <- ad2[order(ad2$order),]
ad2_of <- ad2_o[,1:2]
barplot(t(as.matrix(ad2_of)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2", cex.names = 0.3)

tiff("all_chr15_4_K2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad2_of_Fuerta <- ad2_of[which(rownames(ad2_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Fuerta)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Fuerta", cex.names = 0.3)
ad2_of_Italy <- ad2_of[which(rownames(ad2_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Italy)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Italy", cex.names = 0.3)
ad2_of_Tunisia <- ad2_of[which(rownames(ad2_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Tunisia)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Tunis", cex.names = 0.3)

ad2_of_Nigeria <- ad2_of[which(rownames(ad2_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Nigeria)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Nigera", cex.names = 0.3)
ad2_of_Ghana <- ad2_of[which(rownames(ad2_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Ghana)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Ghana", cex.names = 0.3)
ad2_of_Kenya <- ad2_of[which(rownames(ad2_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Kenya)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Kenya", cex.names = 0.3)
ad2_of_Rwanda <- ad2_of[which(rownames(ad2_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Rwanda)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Rwanda", cex.names = 0.3)
ad2_of_RSA <- ad2_of[which(rownames(ad2_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad2_of_RSA)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - RSA", cex.names = 0.3)

ad2_of_StHelena <- ad2_of[which(rownames(ad2_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad2_of_StHelena)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - StHelena", cex.names = 0.3)
ad2_of_Philippines <- ad2_of[which(rownames(ad2_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Philippines)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Phili", cex.names = 0.3)
dev.off()

ad3 <- read.table("admix_in_chr15_4.3.Q")
rownames(ad3) <- fam$V1
ad3$order <- background_newordered$newnumb
ad3_or <- ad3[order(ad3$V2),]
ad3_or <- ad3_or[order(ad3_or$V3),]
ad3_or <- ad3_or[order(ad3_or$V1),]
ad3_o <- ad3_or[order(ad3_or$order),]
ad3_of <- ad3_o[,1:3]
barplot(t(as.matrix(ad3_of)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3", cex.names = 0.3)

tiff("all_chr15_4_K3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad3_of_Fuerta <- ad3_of[which(rownames(ad3_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Fuerta", cex.names = 0.3)
ad3_of_Italy <- ad3_of[which(rownames(ad3_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Italy", cex.names = 0.3)
ad3_of_Tunisia <- ad3_of[which(rownames(ad3_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Tunis", cex.names = 0.3)

ad3_of_Nigeria <- ad3_of[which(rownames(ad3_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Nigera", cex.names = 0.3)
ad3_of_Ghana <- ad3_of[which(rownames(ad3_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Ghana", cex.names = 0.3)
ad3_of_Kenya <- ad3_of[which(rownames(ad3_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Kenya", cex.names = 0.3)
ad3_of_Rwanda <- ad3_of[which(rownames(ad3_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Rwanda", cex.names = 0.3)
ad3_of_RSA <- ad3_of[which(rownames(ad3_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad3_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - RSA", cex.names = 0.3)

ad3_of_StHelena <- ad3_of[which(rownames(ad3_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad3_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - StHelena", cex.names = 0.3)
ad3_of_Philippines <- ad3_of[which(rownames(ad3_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Phili", cex.names = 0.3)
dev.off()

ad4 <- read.table("admix_in_chr15_4.4.Q")
rownames(ad4) <- fam$V1
ad4$order <- background_newordered$newnumb
ad4_or <- ad4[order(ad4$V1),]
ad4_or <- ad4_or[order(ad4_or$V2),]
ad4_or <- ad4_or[order(ad4_or$V3),]
ad4_or <- ad4_or[order(ad4_or$V4),]

ad4_o <- ad4[order(ad4$order),]
ad4_of <- ad4_o[,1:4]
barplot(t(as.matrix(ad4_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4", cex.names = 0.3)

tiff("all_chr15_4_K4.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad4_of_Fuerta <- ad4_of[which(rownames(ad4_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Fuerta", cex.names = 0.3)
ad4_of_Italy <- ad4_of[which(rownames(ad4_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Italy", cex.names = 0.3)
ad4_of_Tunisia <- ad4_of[which(rownames(ad4_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Tunis", cex.names = 0.3)

ad4_of_Nigeria <- ad4_of[which(rownames(ad4_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Nigera", cex.names = 0.3)
ad4_of_Ghana <- ad4_of[which(rownames(ad4_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Ghana", cex.names = 0.3)
ad4_of_Kenya <- ad4_of[which(rownames(ad4_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Kenya", cex.names = 0.3)
ad4_of_Rwanda <- ad4_of[which(rownames(ad4_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Rwanda", cex.names = 0.3)
ad4_of_RSA <- ad4_of[which(rownames(ad4_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad4_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - RSA", cex.names = 0.3)

ad4_of_StHelena <- ad4_of[which(rownames(ad4_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad4_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - StHelena", cex.names = 0.3)
ad4_of_Philippines <- ad4_of[which(rownames(ad4_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Phili", cex.names = 0.3)
dev.off()

ad5 <- read.table("admix_in_chr15_4.5.Q")
rownames(ad5) <- fam$V1
ad5$order <- background_newordered$newnumb
ad5_or <- ad5[order(ad5$V1),]
ad5_o <- ad5[order(ad5$order),]
ad5_of <- ad5_o[,1:5]
barplot(t(as.matrix(ad5_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5", cex.names = 0.3)

tiff("all_chr15_4_K5.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad5_of_Fuerta <- ad5_of[which(rownames(ad5_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Fuerta", cex.names = 0.3)
ad5_of_Italy <- ad5_of[which(rownames(ad5_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Italy", cex.names = 0.3)
ad5_of_Tunisia <- ad5_of[which(rownames(ad5_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Tunis", cex.names = 0.3)

ad5_of_Nigeria <- ad5_of[which(rownames(ad5_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Nigera", cex.names = 0.3)
ad5_of_Ghana <- ad5_of[which(rownames(ad5_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Ghana", cex.names = 0.3)
ad5_of_Kenya <- ad5_of[which(rownames(ad5_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Kenya", cex.names = 0.3)
ad5_of_Rwanda <- ad5_of[which(rownames(ad5_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Rwanda", cex.names = 0.3)
ad5_of_RSA <- ad5_of[which(rownames(ad5_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad5_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - RSA", cex.names = 0.3)

ad5_of_StHelena <- ad5_of[which(rownames(ad5_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad5_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - StHelena", cex.names = 0.3)
ad5_of_Philippines <- ad5_of[which(rownames(ad5_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Phili", cex.names = 0.3)
dev.off()

ad6 <- read.table("admix_in_chr15_4.6.Q")
rownames(ad6) <- fam$V1
ad6$order <- background_newordered$newnumb
ad6_or <- ad6[order(ad6$V1),]
ad6_o <- ad6[order(ad6$order),]
ad6_of <- ad6_o[,1:6]
barplot(t(as.matrix(ad6_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6", cex.names = 0.3)


tiff("all_chr15_4_K6.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad6_of_Fuerta <- ad6_of[which(rownames(ad6_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Fuerta", cex.names = 0.3)
ad6_of_Italy <- ad6_of[which(rownames(ad6_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Italy", cex.names = 0.3)
ad6_of_Tunisia <- ad6_of[which(rownames(ad6_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Tunis", cex.names = 0.3)

ad6_of_Nigeria <- ad6_of[which(rownames(ad6_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Nigera", cex.names = 0.3)
ad6_of_Ghana <- ad6_of[which(rownames(ad6_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Ghana", cex.names = 0.3)
ad6_of_Kenya <- ad6_of[which(rownames(ad6_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Kenya", cex.names = 0.3)
ad6_of_Rwanda <- ad6_of[which(rownames(ad6_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Rwanda", cex.names = 0.3)
ad6_of_RSA <- ad6_of[which(rownames(ad6_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad6_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - RSA", cex.names = 0.3)

ad6_of_StHelena <- ad6_of[which(rownames(ad6_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad6_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - StHelena", cex.names = 0.3)
ad6_of_Philippines <- ad6_of[which(rownames(ad6_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Phili", cex.names = 0.3)
dev.off()

ad7 <- read.table("admix_in_chr15_4.7.Q")
rownames(ad7) <- fam$V1
ad7$order <- background_newordered$newnumb
ad7_or <- ad7[order(ad7$V3),]
ad7_or <- ad7_or[order(ad7_or$V2),]
ad7_or <- ad7_or[order(ad7_or$V1),]
ad7_o <- ad7[order(ad7$order),]
ad7_of <- ad7_o[,1:7]
barplot(t(as.matrix(ad7_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7", cex.names = 0.3)

tiff("all_chr15_4_K7.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad7_of_Fuerta <- ad7_of[which(rownames(ad7_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Fuerta", cex.names = 0.3)
ad7_of_Italy <- ad7_of[which(rownames(ad7_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Italy", cex.names = 0.3)
ad7_of_Tunisia <- ad7_of[which(rownames(ad7_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Tunis", cex.names = 0.3)

ad7_of_Nigeria <- ad7_of[which(rownames(ad7_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Nigera", cex.names = 0.3)
ad7_of_Ghana <- ad7_of[which(rownames(ad7_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Ghana", cex.names = 0.3)
ad7_of_Kenya <- ad7_of[which(rownames(ad7_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Kenya", cex.names = 0.3)
ad7_of_Rwanda <- ad7_of[which(rownames(ad7_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Rwanda", cex.names = 0.3)
ad7_of_RSA <- ad7_of[which(rownames(ad7_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad7_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - RSA", cex.names = 0.3)

ad7_of_StHelena <- ad7_of[which(rownames(ad7_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad7_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - StHelena", cex.names = 0.3)
ad7_of_Philippines <- ad7_of[which(rownames(ad7_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Phili", cex.names = 0.3)
dev.off()

ad8 <- read.table("admix_in_chr15_4.8.Q")
rownames(ad8) <- fam$V1
ad8$order <- background_newordered$newnumb
ad8_or <- ad8[order(ad8$V1),]
ad8_o <- ad8[order(ad8$order),]
ad8_of <- ad8_o[,1:8]

barplot(t(as.matrix(ad8_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8", cex.names = 0.3)

tiff("all_chr15_4_K8.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad8_of_Fuerta <- ad8_of[which(rownames(ad8_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Fuerta", cex.names = 0.3)
ad8_of_Italy <- ad8_of[which(rownames(ad8_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Italy", cex.names = 0.3)
ad8_of_Tunisia <- ad8_of[which(rownames(ad8_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Tunis", cex.names = 0.3)

ad8_of_Nigeria <- ad8_of[which(rownames(ad8_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Nigera", cex.names = 0.3)
ad8_of_Ghana <- ad8_of[which(rownames(ad8_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Ghana", cex.names = 0.3)
ad8_of_Kenya <- ad8_of[which(rownames(ad8_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Kenya", cex.names = 0.3)
ad8_of_Rwanda <- ad8_of[which(rownames(ad8_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Rwanda", cex.names = 0.3)
ad8_of_RSA <- ad8_of[which(rownames(ad8_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad8_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - RSA", cex.names = 0.3)

ad8_of_StHelena <- ad8_of[which(rownames(ad8_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad8_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - StHelena", cex.names = 0.3)
ad8_of_Philippines <- ad8_of[which(rownames(ad8_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Phili", cex.names = 0.3)
dev.off()

ad9 <- read.table("admix_in_chr15_4.9.Q")
rownames(ad9) <- fam$V1
ad9$order <- background_newordered$newnumb
ad9_or <- ad9[order(ad9$V1),]
ad9_o <- ad9[order(ad9$order),]
ad9_of <- ad9_o[,1:9]

barplot(t(as.matrix(ad9_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9", cex.names = 0.3)

tiff("all_chr15_4_K9.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad9_of_Fuerta <- ad9_of[which(rownames(ad9_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Fuerta", cex.names = 0.3)
ad9_of_Italy <- ad9_of[which(rownames(ad9_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Italy", cex.names = 0.3)
ad9_of_Tunisia <- ad9_of[which(rownames(ad9_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Tunis", cex.names = 0.3)

ad9_of_Nigeria <- ad9_of[which(rownames(ad9_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Nigera", cex.names = 0.3)
ad9_of_Ghana <- ad9_of[which(rownames(ad9_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Ghana", cex.names = 0.3)
ad9_of_Kenya <- ad9_of[which(rownames(ad9_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Kenya", cex.names = 0.3)
ad9_of_Rwanda <- ad9_of[which(rownames(ad9_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Rwanda", cex.names = 0.3)
ad9_of_RSA <- ad9_of[which(rownames(ad9_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad9_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - RSA", cex.names = 0.3)

ad9_of_StHelena <- ad9_of[which(rownames(ad9_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad9_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - StHelena", cex.names = 0.3)
ad9_of_Philippines <- ad9_of[which(rownames(ad9_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Phili", cex.names = 0.3)
dev.off()

ad10 <- read.table("admix_in_chr15_4.10.Q")
rownames(ad10) <- fam$V1
ad10$order <- background_newordered$newnumb
ad10_or <- ad10[order(ad10$V3),]
ad10_or <- ad10_or[order(ad10_or$V2),]
ad10_or <- ad10_or[order(ad10_or$V1),]
ad10_o <- ad10[order(ad10$order),]
ad10_of <- ad10_o[,1:10]

tiff("all_chr15_4_K10.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad10_of_Fuerta <- ad10_of[which(rownames(ad10_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Fuerta", cex.names = 0.3)
ad10_of_Italy <- ad10_of[which(rownames(ad10_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Italy", cex.names = 0.3)
ad10_of_Tunisia <- ad10_of[which(rownames(ad10_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Tunis", cex.names = 0.3)

ad10_of_Nigeria <- ad10_of[which(rownames(ad10_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Nigera", cex.names = 0.3)
ad10_of_Ghana <- ad10_of[which(rownames(ad10_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Ghana", cex.names = 0.3)
ad10_of_Kenya <- ad10_of[which(rownames(ad10_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Kenya", cex.names = 0.3)
ad10_of_Rwanda <- ad10_of[which(rownames(ad10_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Rwanda", cex.names = 0.3)
ad10_of_RSA <- ad10_of[which(rownames(ad10_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad10_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - RSA", cex.names = 0.3)

ad10_of_StHelena <- ad10_of[which(rownames(ad10_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad10_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - StHelena", cex.names = 0.3)
ad10_of_Philippines <- ad10_of[which(rownames(ad10_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Phili", cex.names = 0.3)
dev.off()

# cv error plot:
k_vals <- c(2:10)
cv_vals <- c(0.25524, 0.23387, 0.21533, 0.20582, 0.20382, 0.20439, 0.20421, 0.20115, 0.21060)
plot(k_vals, cv_vals, type = 'l', lty = 2, ylab = "CV error", xlab = "K", main = "reg 4 - best k=9")
points(k_vals, cv_vals, pch = 16, col= "black")
points(9, 0.20115, pch = 16, col= "red")

######cv-error plot####
tiff("cv-error_regions_indiv.tiff", height=6, width=6, units="in", res=400, compression="lzw")
par(mfrow=c(2,2))
# cv error plot reg 1.1:
k_vals <- c(2:10)
cv_vals <- c(0.21159, 0.20095, 0.19305, 0.18889, 0.17838, 0.17958, 0.19436, 0.18897, 0.19431)

plot(k_vals, cv_vals, type = 'l', lty = 2, ylab = "CV error", xlab = "K", main = "reg 1.1 - best k=6")
points(k_vals, cv_vals, pch = 16, col= "black")
points(6, 0.17838, pch = 16, col= "red")

# cv error plot reg 1.2:
k_vals <- c(2:10)
cv_vals <- c(0.22629, 0.15140, 0.12404, 0.11721, 0.10469, 0.10818, 0.10502, 0.10656, 0.10820)
plot(k_vals, cv_vals, type = 'l', lty = 2, ylab = "CV error", xlab = "K", main = "reg 1.2 - best k=6")
points(k_vals, cv_vals, pch = 16, col= "black")
points(6, 0.10469, pch = 16, col= "red")

# cv error plot reg 2:
k_vals <- c(2:10)
cv_vals <- c(0.24583, 0.16874, 0.14242, 0.13033, 0.12113, 0.12320, 0.12075, 0.12064, 0.13057)
plot(k_vals, cv_vals, type = 'l', lty = 2, ylab = "CV error", xlab = "K", main = "reg 2 - best k=9")
points(k_vals, cv_vals, pch = 16, col= "black")
points(9, 0.12064, pch = 16, col= "red")

# cv error plot reg 4:
k_vals <- c(2:10)
cv_vals <- c(0.25524, 0.23387, 0.21533, 0.20582, 0.20382, 0.20439, 0.20421, 0.20115, 0.21060)
plot(k_vals, cv_vals, type = 'l', lty = 2, ylab = "CV error", xlab = "K", main = "reg 4 - best k=9")
points(k_vals, cv_vals, pch = 16, col= "black")
points(9, 0.20115, pch = 16, col= "red")

dev.off()

#################plot reference indivs
klugii_list <- c("SM15W61", "SM15W66", "SM15W69", "SM15W72", "SM15W74", "SM18W01")
klugii <- background_newordered[which(background_newordered$indiv %in% klugii_list),]

chrys_list <- c("SM16N01", "SM16N04", "SM16N05", "SM16N06", "SM16N20", "SM16N37")
chrys <- background_newordered[which(background_newordered$indiv %in% chrys_list),]

ori_list <- c("SM21TS03", "SM21TS05", "SM21TS16", "SM21TS33", "SM21TS34", "SM21TS35", "SM21TS36")
ori <- background_newordered[which(background_newordered$indiv %in% ori_list),]

weird_list <- c("SM21SP06", "SM21SP09")
weird <- background_newordered[which(background_newordered$indiv %in% weird_list),]

#region 1.1 k=6
ad6 <- read.table("admix_in_chr15_1.1.6.Q")
rownames(ad6) <- fam$V1
ad6$order <- background_newordered$newnumb
ad6_or <- ad6[order(ad6$V1),]
ad6_o <- ad6[order(ad6$order),]
ad6_of <- ad6_o[,1:6]
barplot(t(as.matrix(ad6_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6", cex.names = 0.3)


tiff("all_chr15_1.1_K6_chr_klu_ori_weird.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,2))
ad6_of_chr <- ad6_of[which(rownames(ad6_of) %in% chrys$indiv),]
barplot(t(as.matrix(ad6_of_chr)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - chrys", cex.names = 0.3)
ad6_of_klu <- ad6_of[which(rownames(ad6_of) %in% klugii$indiv),]
barplot(t(as.matrix(ad6_of_klu)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - klugii", cex.names = 0.3)
ad6_of_ori <- ad6_of[which(rownames(ad6_of) %in% ori$indiv),]
barplot(t(as.matrix(ad6_of_ori)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - ori", cex.names = 0.3)
ad6_of_weird <- ad6_of[which(rownames(ad6_of) %in% weird$indiv),]
barplot(t(as.matrix(ad6_of_weird)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - weird", cex.names = 0.3)
dev.off()

#region 1.2 k=6
ad6 <- read.table("admix_in_chr15_1.2.6.Q")
rownames(ad6) <- fam$V1
ad6$order <- background_newordered$newnumb
ad6_or <- ad6[order(ad6$V1),]
ad6_o <- ad6[order(ad6$order),]
ad6_of <- ad6_o[,1:6]
barplot(t(as.matrix(ad6_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6", cex.names = 0.3)


tiff("all_chr15_1.2_K6_chr_klu_ori_weird.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,2))
ad6_of_chr <- ad6_of[which(rownames(ad6_of) %in% chrys$indiv),]
barplot(t(as.matrix(ad6_of_chr)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - chrys", cex.names = 0.3)
ad6_of_klu <- ad6_of[which(rownames(ad6_of) %in% klugii$indiv),]
barplot(t(as.matrix(ad6_of_klu)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - klugii", cex.names = 0.3)
ad6_of_ori <- ad6_of[which(rownames(ad6_of) %in% ori$indiv),]
barplot(t(as.matrix(ad6_of_ori)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - ori", cex.names = 0.3)
ad6_of_weird <- ad6_of[which(rownames(ad6_of) %in% weird$indiv),]
barplot(t(as.matrix(ad6_of_weird)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - weird", cex.names = 0.3)
dev.off()

#region 2 k=9
ad9 <- read.table("admix_in_chr15_2.9.Q")
rownames(ad9) <- fam$V1
ad9$order <- background_newordered$newnumb
ad9_or <- ad9[order(ad9$V1),]
ad9_o <- ad9[order(ad9$order),]
ad9_of <- ad9_o[,1:9]
barplot(t(as.matrix(ad9_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9", cex.names = 0.3)


tiff("all_chr15_2_K9_chr_klu_ori_weird.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,2))
ad9_of_chr <- ad9_of[which(rownames(ad9_of) %in% chrys$indiv),]
barplot(t(as.matrix(ad9_of_chr)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - chrys", cex.names = 0.3)
ad9_of_klu <- ad9_of[which(rownames(ad9_of) %in% klugii$indiv),]
barplot(t(as.matrix(ad9_of_klu)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - klugii", cex.names = 0.3)
ad9_of_ori <- ad9_of[which(rownames(ad9_of) %in% ori$indiv),]
barplot(t(as.matrix(ad9_of_ori)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - ori", cex.names = 0.3)
ad9_of_weird <- ad9_of[which(rownames(ad9_of) %in% weird$indiv),]
barplot(t(as.matrix(ad9_of_weird)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - weird", cex.names = 0.3)
dev.off()

#region 4 k=9
ad9 <- read.table("admix_in_chr15_4.9.Q")
rownames(ad9) <- fam$V1
ad9$order <- background_newordered$newnumb
ad9_or <- ad9[order(ad9$V1),]
ad9_o <- ad9[order(ad9$order),]
ad9_of <- ad9_o[,1:9]
barplot(t(as.matrix(ad9_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9", cex.names = 0.3)


tiff("all_chr15_4_K9_chr_klu_ori_weird.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,2))
ad9_of_chr <- ad9_of[which(rownames(ad9_of) %in% chrys$indiv),]
barplot(t(as.matrix(ad9_of_chr)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - chrys", cex.names = 0.3)
ad9_of_klu <- ad9_of[which(rownames(ad9_of) %in% klugii$indiv),]
barplot(t(as.matrix(ad9_of_klu)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - klugii", cex.names = 0.3)
ad9_of_ori <- ad9_of[which(rownames(ad9_of) %in% ori$indiv),]
barplot(t(as.matrix(ad9_of_ori)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - ori", cex.names = 0.3)
ad9_of_weird <- ad9_of[which(rownames(ad9_of) %in% weird$indiv),]
barplot(t(as.matrix(ad9_of_weird)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - weird", cex.names = 0.3)
dev.off()

##########################
#CHR15 SVs only

#plot cv error:
vals <- c(0.24417, 0.19517, 0.18068, 0.16924, 0.16244, 0.16756, 0.16557, 0.16193, 0.16216, 0.16540, 0.17455)
CV_error <- as.data.frame(2:(length(vals)+1))
colnames(CV_error) <- "k"
CV_error$CV <- vals

tiff("./CHR_15_SV_ONLYCV_plot.tiff", height=4, width=4, units="in", res=300, compression="lzw")
plot(CV_error$k, CV_error$CV, pch = 16, ylab = "CV error from admixture run", xlab = "K")
lines(CV_error$k, CV_error$CV, col = "grey50", lty = 3)
dev.off()

fam <- read.csv("chr15_onlySVs/admix_in_chr15_onlySVs.fam", sep = " ", header = FALSE)

background$numb <- 1:174

background$country_numb <- c()
for(i in 1:nrow(background)){
  if(background$country[i] == "Fuerteventura"){
    background$country_numb[i] <- 1
  }
  if(background$country[i] == "Italy"){
    background$country_numb[i] <- 2
  }
  if(background$country[i] == "Tunisia"){
    background$country_numb[i] <- 3
  }
  if(background$country[i] == "Nigeria"){
    background$country_numb[i] <- 4
  }
  if(background$country[i] == "Ghana"){
    background$country_numb[i] <- 5
  }
  if(background$country[i] == "Kenya"){
    background$country_numb[i] <- 6
  }
  if(background$country[i] == "Rwanda"){
    background$country_numb[i] <- 7
  }
  if(background$country[i] == "RSA"){
    background$country_numb[i] <- 8
  }
  if(background$country[i] == "StHelena"){
    background$country_numb[i] <- 9
  }
  if(background$country[i] == "Philippines"){
    background$country_numb[i] <- 10
  }
}

background_ordered <- background[order(background$country_numb),]
background_ordered$newnumb <- 1:174
background_newordered <- background_ordered[order(background_ordered$numb),]

Fuerta_indivs <- subset(background_newordered, background_newordered$country == "Spain")
Italy_indivs <- subset(background_newordered, background_newordered$country == "Italy")
Tunisia_indivs <- subset(background_newordered, background_newordered$country == "Tunisia")
Nigeria_indivs <- subset(background_newordered, background_newordered$country == "Nigeria")
Ghana_indivs <- subset(background_newordered, background_newordered$country == "Ghana")
Kenya_indivs <- subset(background_newordered, background_newordered$country == "Kenya")
Rwanda_indivs <- subset(background_newordered, background_newordered$country == "Rwanda")
RSA_indivs <- subset(background_newordered, background_newordered$country == "RSA")
StHelena_indivs <- subset(background_newordered, background_newordered$country == "StHelena")
Philippines_indivs <- subset(background_newordered, background_newordered$country == "Philippines")
missing <- subset(background_newordered, background_newordered$country == "missing")

popsum <- sum(nrow(Fuerta_indivs), nrow(Italy_indivs), nrow(Tunisia_indivs), nrow(Nigeria_indivs),
              nrow(Ghana_indivs), nrow(Kenya_indivs), nrow(Rwanda_indivs), nrow(RSA_indivs), 
              nrow(StHelena_indivs), nrow(Philippines_indivs), nrow(missing))

ad2 <- read.table("chr15_onlySVs/admix_in_chr15_onlySVs.2.Q")
rownames(ad2) <- fam$V1
ad2$order <- background_newordered$newnumb
ad2_or <- ad2[order(ad2$V1),]
ad2_o <- ad2[order(ad2$order),]
ad2_of <- ad2_o[,1:2]
barplot(t(as.matrix(ad2_of)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2", cex.names = 0.3)

tiff("all_chr15_SV_ONLY_K2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad2_of_Fuerta <- ad2_of[which(rownames(ad2_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Fuerta)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Fuerta", cex.names = 0.3)
ad2_of_Italy <- ad2_of[which(rownames(ad2_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Italy)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Italy", cex.names = 0.3)
ad2_of_Tunisia <- ad2_of[which(rownames(ad2_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Tunisia)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Tunis", cex.names = 0.3)

ad2_of_Nigeria <- ad2_of[which(rownames(ad2_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Nigeria)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Nigera", cex.names = 0.3)
ad2_of_Ghana <- ad2_of[which(rownames(ad2_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Ghana)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Ghana", cex.names = 0.3)
ad2_of_Kenya <- ad2_of[which(rownames(ad2_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Kenya)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Kenya", cex.names = 0.3)
ad2_of_Rwanda <- ad2_of[which(rownames(ad2_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Rwanda)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Rwanda", cex.names = 0.3)
ad2_of_RSA <- ad2_of[which(rownames(ad2_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad2_of_RSA)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - RSA", cex.names = 0.3)

ad2_of_StHelena <- ad2_of[which(rownames(ad2_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad2_of_StHelena)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - StHelena", cex.names = 0.3)
ad2_of_Philippines <- ad2_of[which(rownames(ad2_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad2_of_Philippines)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2 - Phili", cex.names = 0.3)
dev.off()

ad3 <- read.table("chr15_onlySVs/admix_in_chr15_onlySVs.3.Q")
rownames(ad3) <- fam$V1
ad3$order <- background_newordered$newnumb
ad3_or <- ad3[order(ad3$V2),]
ad3_or <- ad3_or[order(ad3_or$V3),]
ad3_or <- ad3_or[order(ad3_or$V1),]
ad3_o <- ad3_or[order(ad3_or$order),]
ad3_of <- ad3_o[,1:3]
barplot(t(as.matrix(ad3_of)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3", cex.names = 0.3)

tiff("all_chr15_SV_ONLY_K3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad3_of_Fuerta <- ad3_of[which(rownames(ad3_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Fuerta", cex.names = 0.3)
ad3_of_Italy <- ad3_of[which(rownames(ad3_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Italy", cex.names = 0.3)
ad3_of_Tunisia <- ad3_of[which(rownames(ad3_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Tunis", cex.names = 0.3)

ad3_of_Nigeria <- ad3_of[which(rownames(ad3_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Nigera", cex.names = 0.3)
ad3_of_Ghana <- ad3_of[which(rownames(ad3_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Ghana", cex.names = 0.3)
ad3_of_Kenya <- ad3_of[which(rownames(ad3_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Kenya", cex.names = 0.3)
ad3_of_Rwanda <- ad3_of[which(rownames(ad3_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Rwanda", cex.names = 0.3)
ad3_of_RSA <- ad3_of[which(rownames(ad3_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad3_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - RSA", cex.names = 0.3)

ad3_of_StHelena <- ad3_of[which(rownames(ad3_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad3_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - StHelena", cex.names = 0.3)
ad3_of_Philippines <- ad3_of[which(rownames(ad3_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad3_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3 - Phili", cex.names = 0.3)
dev.off()

ad4 <- read.table("chr15_onlySVs/admix_in_chr15_onlySVs.4.Q")
rownames(ad4) <- fam$V1
ad4$order <- background_newordered$newnumb
ad4_or <- ad4[order(ad4$V1),]
ad4_or <- ad4_or[order(ad4_or$V2),]
ad4_or <- ad4_or[order(ad4_or$V3),]
ad4_or <- ad4_or[order(ad4_or$V4),]

ad4_o <- ad4[order(ad4$order),]
ad4_of <- ad4_o[,1:4]
barplot(t(as.matrix(ad4_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4", cex.names = 0.3)

tiff("all_chr15_SV_ONLY_K4.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad4_of_Fuerta <- ad4_of[which(rownames(ad4_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Fuerta", cex.names = 0.3)
ad4_of_Italy <- ad4_of[which(rownames(ad4_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Italy", cex.names = 0.3)
ad4_of_Tunisia <- ad4_of[which(rownames(ad4_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Tunis", cex.names = 0.3)

ad4_of_Nigeria <- ad4_of[which(rownames(ad4_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Nigera", cex.names = 0.3)
ad4_of_Ghana <- ad4_of[which(rownames(ad4_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Ghana", cex.names = 0.3)
ad4_of_Kenya <- ad4_of[which(rownames(ad4_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Kenya", cex.names = 0.3)
ad4_of_Rwanda <- ad4_of[which(rownames(ad4_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Rwanda", cex.names = 0.3)
ad4_of_RSA <- ad4_of[which(rownames(ad4_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad4_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - RSA", cex.names = 0.3)

ad4_of_StHelena <- ad4_of[which(rownames(ad4_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad4_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - StHelena", cex.names = 0.3)
ad4_of_Philippines <- ad4_of[which(rownames(ad4_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad4_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4 - Phili", cex.names = 0.3)
dev.off()

ad5 <- read.table("chr15_onlySVs/admix_in_chr15_onlySVs.5.Q")
rownames(ad5) <- fam$V1
ad5$order <- background_newordered$newnumb
ad5_or <- ad5[order(ad5$V1),]
ad5_o <- ad5[order(ad5$order),]
ad5_of <- ad5_o[,1:5]
barplot(t(as.matrix(ad5_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5", cex.names = 0.3)

tiff("all_chr15_SV_ONLY_K5.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad5_of_Fuerta <- ad5_of[which(rownames(ad5_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Fuerta", cex.names = 0.3)
ad5_of_Italy <- ad5_of[which(rownames(ad5_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Italy", cex.names = 0.3)
ad5_of_Tunisia <- ad5_of[which(rownames(ad5_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Tunis", cex.names = 0.3)

ad5_of_Nigeria <- ad5_of[which(rownames(ad5_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Nigera", cex.names = 0.3)
ad5_of_Ghana <- ad5_of[which(rownames(ad5_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Ghana", cex.names = 0.3)
ad5_of_Kenya <- ad5_of[which(rownames(ad5_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Kenya", cex.names = 0.3)
ad5_of_Rwanda <- ad5_of[which(rownames(ad5_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Rwanda", cex.names = 0.3)
ad5_of_RSA <- ad5_of[which(rownames(ad5_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad5_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - RSA", cex.names = 0.3)

ad5_of_StHelena <- ad5_of[which(rownames(ad5_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad5_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - StHelena", cex.names = 0.3)
ad5_of_Philippines <- ad5_of[which(rownames(ad5_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad5_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5 - Phili", cex.names = 0.3)
dev.off()

ad6 <- read.table("chr15_onlySVs/admix_in_chr15_onlySVs.6.Q")
rownames(ad6) <- fam$V1
ad6$order <- background_newordered$newnumb
ad6_or <- ad6[order(ad6$V1),]
ad6_o <- ad6[order(ad6$order),]
ad6_of <- ad6_o[,1:6]
barplot(t(as.matrix(ad6_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6", cex.names = 0.3)


tiff("all_chr15_SV_ONLY_K6.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad6_of_Fuerta <- ad6_of[which(rownames(ad6_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Fuerta", cex.names = 0.3)
ad6_of_Italy <- ad6_of[which(rownames(ad6_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Italy", cex.names = 0.3)
ad6_of_Tunisia <- ad6_of[which(rownames(ad6_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Tunis", cex.names = 0.3)

ad6_of_Nigeria <- ad6_of[which(rownames(ad6_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Nigera", cex.names = 0.3)
ad6_of_Ghana <- ad6_of[which(rownames(ad6_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Ghana", cex.names = 0.3)
ad6_of_Kenya <- ad6_of[which(rownames(ad6_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Kenya", cex.names = 0.3)
ad6_of_Rwanda <- ad6_of[which(rownames(ad6_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Rwanda", cex.names = 0.3)
ad6_of_RSA <- ad6_of[which(rownames(ad6_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad6_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - RSA", cex.names = 0.3)

ad6_of_StHelena <- ad6_of[which(rownames(ad6_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad6_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - StHelena", cex.names = 0.3)
ad6_of_Philippines <- ad6_of[which(rownames(ad6_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad6_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6 - Phili", cex.names = 0.3)
dev.off()

ad7 <- read.table("chr15_onlySVs/admix_in_chr15_onlySVs.7.Q")
rownames(ad7) <- fam$V1
ad7$order <- background_newordered$newnumb
ad7_or <- ad7[order(ad7$V3),]
ad7_or <- ad7_or[order(ad7_or$V2),]
ad7_or <- ad7_or[order(ad7_or$V1),]
ad7_o <- ad7[order(ad7$order),]
ad7_of <- ad7_o[,1:7]
barplot(t(as.matrix(ad7_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7", cex.names = 0.3)

tiff("all_chr15_SV_ONLY_K7.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad7_of_Fuerta <- ad7_of[which(rownames(ad7_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Fuerta", cex.names = 0.3)
ad7_of_Italy <- ad7_of[which(rownames(ad7_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Italy", cex.names = 0.3)
ad7_of_Tunisia <- ad7_of[which(rownames(ad7_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Tunis", cex.names = 0.3)

ad7_of_Nigeria <- ad7_of[which(rownames(ad7_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Nigera", cex.names = 0.3)
ad7_of_Ghana <- ad7_of[which(rownames(ad7_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Ghana", cex.names = 0.3)
ad7_of_Kenya <- ad7_of[which(rownames(ad7_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Kenya", cex.names = 0.3)
ad7_of_Rwanda <- ad7_of[which(rownames(ad7_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Rwanda", cex.names = 0.3)
ad7_of_RSA <- ad7_of[which(rownames(ad7_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad7_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - RSA", cex.names = 0.3)

ad7_of_StHelena <- ad7_of[which(rownames(ad7_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad7_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - StHelena", cex.names = 0.3)
ad7_of_Philippines <- ad7_of[which(rownames(ad7_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad7_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7 - Phili", cex.names = 0.3)
dev.off()

ad8 <- read.table("chr15_onlySVs/admix_in_chr15_onlySVs.8.Q")
rownames(ad8) <- fam$V1
ad8$order <- background_newordered$newnumb
ad8_or <- ad8[order(ad8$V1),]
ad8_o <- ad8[order(ad8$order),]
ad8_of <- ad8_o[,1:8]

barplot(t(as.matrix(ad8_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8", cex.names = 0.3)

tiff("all_chr15_SV_ONLY_K8.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad8_of_Fuerta <- ad8_of[which(rownames(ad8_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Fuerta", cex.names = 0.3)
ad8_of_Italy <- ad8_of[which(rownames(ad8_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Italy", cex.names = 0.3)
ad8_of_Tunisia <- ad8_of[which(rownames(ad8_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Tunis", cex.names = 0.3)

ad8_of_Nigeria <- ad8_of[which(rownames(ad8_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Nigera", cex.names = 0.3)
ad8_of_Ghana <- ad8_of[which(rownames(ad8_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Ghana", cex.names = 0.3)
ad8_of_Kenya <- ad8_of[which(rownames(ad8_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Kenya", cex.names = 0.3)
ad8_of_Rwanda <- ad8_of[which(rownames(ad8_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Rwanda", cex.names = 0.3)
ad8_of_RSA <- ad8_of[which(rownames(ad8_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad8_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - RSA", cex.names = 0.3)

ad8_of_StHelena <- ad8_of[which(rownames(ad8_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad8_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - StHelena", cex.names = 0.3)
ad8_of_Philippines <- ad8_of[which(rownames(ad8_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad8_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8 - Phili", cex.names = 0.3)
dev.off()

ad9 <- read.table("chr15_onlySVs/admix_in_chr15_onlySVs.9.Q")
rownames(ad9) <- fam$V1
ad9$order <- background_newordered$newnumb
ad9_or <- ad9[order(ad9$V1),]
ad9_o <- ad9[order(ad9$order),]
ad9_of <- ad9_o[,1:9]

barplot(t(as.matrix(ad9_of)), col = c("chartreuse4", "paleturquoise2", "tomato", "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9", cex.names = 0.3)

tiff("all_chr15_SV_ONLY_K9.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad9_of_Fuerta <- ad9_of[which(rownames(ad9_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Fuerta", cex.names = 0.3)
ad9_of_Italy <- ad9_of[which(rownames(ad9_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Italy", cex.names = 0.3)
ad9_of_Tunisia <- ad9_of[which(rownames(ad9_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Tunis", cex.names = 0.3)

ad9_of_Nigeria <- ad9_of[which(rownames(ad9_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Nigera", cex.names = 0.3)
ad9_of_Ghana <- ad9_of[which(rownames(ad9_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Ghana", cex.names = 0.3)
ad9_of_Kenya <- ad9_of[which(rownames(ad9_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Kenya", cex.names = 0.3)
ad9_of_Rwanda <- ad9_of[which(rownames(ad9_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Rwanda", cex.names = 0.3)
ad9_of_RSA <- ad9_of[which(rownames(ad9_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad9_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - RSA", cex.names = 0.3)

ad9_of_StHelena <- ad9_of[which(rownames(ad9_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad9_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - StHelena", cex.names = 0.3)
ad9_of_Philippines <- ad9_of[which(rownames(ad9_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad9_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9 - Phili", cex.names = 0.3)
dev.off()

ad10 <- read.table("chr15_onlySVs/admix_in_chr15_onlySVs.10.Q")
rownames(ad10) <- fam$V1
ad10$order <- background_newordered$newnumb
ad10_or <- ad10[order(ad10$V3),]
ad10_or <- ad10_or[order(ad10_or$V2),]
ad10_or <- ad10_or[order(ad10_or$V1),]
ad10_o <- ad10[order(ad10$order),]
ad10_of <- ad10_o[,1:10]

tiff("all_chr15_SV_ONLY_K10.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mfrow=c(2,5))
ad10_of_Fuerta <- ad10_of[which(rownames(ad10_of) %in% Fuerta_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Fuerta)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Fuerta", cex.names = 0.3)
ad10_of_Italy <- ad10_of[which(rownames(ad10_of) %in% Italy_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Italy)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Italy", cex.names = 0.3)
ad10_of_Tunisia <- ad10_of[which(rownames(ad10_of) %in% Tunisia_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Tunisia)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Tunis", cex.names = 0.3)

ad10_of_Nigeria <- ad10_of[which(rownames(ad10_of) %in% Nigeria_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Nigeria)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Nigera", cex.names = 0.3)
ad10_of_Ghana <- ad10_of[which(rownames(ad10_of) %in% Ghana_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Ghana)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Ghana", cex.names = 0.3)
ad10_of_Kenya <- ad10_of[which(rownames(ad10_of) %in% Kenya_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Kenya)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Kenya", cex.names = 0.3)
ad10_of_Rwanda <- ad10_of[which(rownames(ad10_of) %in% Rwanda_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Rwanda)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Rwanda", cex.names = 0.3)
ad10_of_RSA <- ad10_of[which(rownames(ad10_of) %in% RSA_indivs$indiv),]
barplot(t(as.matrix(ad10_of_RSA)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - RSA", cex.names = 0.3)

ad10_of_StHelena <- ad10_of[which(rownames(ad10_of) %in% StHelena_indivs$indiv),]
barplot(t(as.matrix(ad10_of_StHelena)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - StHelena", cex.names = 0.3)
ad10_of_Philippines <- ad10_of[which(rownames(ad10_of) %in% Philippines_indivs$indiv),]
barplot(t(as.matrix(ad10_of_Philippines)), col = c("chartreuse4", "paleturquoise2", "tomato",  "goldenrod", "darkgrey", "purple4", "hotpink", "black", "lightgreen", "yellow"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10 - Phili", cex.names = 0.3)
dev.off()


