setwd("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/plink_gwas//")
background <- read.csv("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/background/DC174_background.csv", header = T, sep = ",")
bacakground_uniq_col <- unique(background$colour)
bacakground_uniq_country <- unique(background$country)

chromosomes <- read.csv("../gemma/Dchry2.2.fa.fai", head = FALSE, sep = "\t")
rename <- read.csv("../gemma/rename_file.txt", head = FALSE, sep = " ")
rename$name <- paste(rename$V7, rename$V8)
rename$name <- gsub(" ", "", rename$name)

plot_list <- c("background", "forewing", "white")
plot_name_list <- c("background", "forewing", "white")

for(j in 1:length(plot_list)){
  filename <- paste(as.character(plot_list[j]), "_plot_emp_total.txt", sep = "")
  print(filename)
  background_points <- read.csv(filename, sep = " ")
  
  background_points <- na.omit(background_points)
  background_points$log <- -log10(background_points$P)
  background_points <- subset(background_points, background_points$log > 1)
  background_points$chrom <- paste("chr", background_points$CHR, sep = "")
  
  chromlist <- as.list(as.character(rename$name))
  b <- as.vector(chromlist)
  a <- b
  
  pdf(file=paste("./plots/", as.character(plot_name_list[j]), ".pdf", sep = ""),width=15, height=4)
  
 # par(mfrow=c(1,1),mar=c(0.01,0.05,0.05,0.01),oma=c(4,3,3,1))
  par(mfrow=c(1,1),mar=c(0.01,0.01,0.01,0.01),oma=c(4,2,2,1))
  
  l <- as.vector(chromosomes$V2)
  l <- l[1:(length(rename$V1))]
  
  #set layout parameters specifying 5 rows, one for each comarison and the widths of the chromosomes in each case
  layout(matrix(seq(1:length(rename$V1)),nrow=1,ncol=length(rename$V1),byrow=TRUE),width=c(l))
  
  altcols <- c("grey90", "grey75")
  altcols <- rep(altcols, as.integer((length(a)+1)/2))
  
  altcol_numb <- c(1,2,2,1,2,2,1,
                   2,1,2,1,2,1,2,
                   1,2,1,2,1,1,2,
                   2,1,1,2,2,1,2,
                   2,1,1,2,1,2,2,
                   1,1,2,2,1,2)
  
  print(paste("running:", filename, sep = ""))
  
  for (i in 1:length(rename$V1)){
    #extract the x/bp position for each chrom
    xFst <- background_points$BP[as.character(background_points$chrom)==as.character(a[i])]
    #extract the fst values
    yFst <- background_points$log[as.character(background_points$chrom)==as.character(a[i])]
    #then plot
    
    plot(xFst,yFst,pch=16,cex=1,ylab="-log10(p-value)",xaxt="n",axes=F,ylim=c(0,(3+max(background_points$log))),col=altcols[altcol_numb[i]])
    
    
    outliers <- subset(background_points, as.character(background_points$chrom)==as.character(a[i]))
    outliers <- subset(outliers, outliers$EMP <= as.numeric(quantile(background_points$EMP1, 0.01)[1]))
    outliers_x <- outliers$BP
    outliers_y <- outliers$log
    points(outliers_x, outliers_y, col="black", pch = 16)
    
    if (i==1){
      axis(side=2, at = c(0,20,40,60), labels = c(0,20,40,60))
    }
  }
  dev.off()
  
  #print("writing outliers")
  #outliers <- subset(background_points, background_points$log>hardcut)
  #write.csv(x = outliers, file=paste("./outlier_files/", as.character(plot_name_list[j]), ".csv", sep = ""))
  
  #max_snp <- subset(outliers, outliers$log==max(outliers$log))
  #correct_chr_name <- subset(rename, as.character(rename$name) == as.character(max_snp$chrom))
  #print(paste("the highest peak is on ", correct_chr_name$V1, sep = ""))
  #print(paste("the highest peak has a -log10(p) of ", max_snp$log, sep = ""))
  
}

#now just chr15

###background
background_points <- read.csv("background_plot_emp_total.txt", sep = " ")

background_points <- na.omit(background_points)
background_points$log <- -log10(background_points$P)
background_points <- subset(background_points, background_points$log > 1)
background_points$chrom <- paste("chr", background_points$CHR, sep = "")

background_chr15 <- subset(background_points, background_points$CHR == 17)
background_chr15pt1 <- subset(background_chr15, background_chr15$BP < 7.8e6)
background_chr15pt2 <- subset(background_chr15, background_chr15$BP > 14.8e6)
background_chr15 <- as.data.frame(rbind(background_chr15pt1, background_chr15pt2))

tiff("plots/background_chr15_new.tiff", height=4, width=18, units="in", res=300, compression="lzw")
plot(background_chr15$BP,background_chr15$log,pch=16,cex=1,
     ylab="-log10(p-value)",
     xlab="bp chr15",
     axes=T, 
     ylim=c(0,(3+max(background_chr15$log))),
     col="black",
     type = 'n')

rect(xleft = 4035981, ybottom = (2+max(background_chr15$log)), 
     xright = 5322257, ytop = (5+max(background_chr15$log)), 
     col = "#016337", border = T, lwd = 0.7)
rect(xleft = 5323805, ybottom = (2+max(background_chr15$log)),
     xright = 6220875, ytop = (5+max(background_chr15$log)), 
     col = "#016337", border = T, lwd = 0.7)
rect(xleft = 6251636, ybottom = (2+max(background_chr15$log)), 
     xright = 7829094, ytop = (5+max(background_chr15$log)), 
     col = "#a875ee", border = T, lwd = 0.7)
rect(xleft = 14803486, ybottom = (2+max(background_chr15$log)), 
     xright = 15299732, ytop = (5+max(background_chr15$log)), 
     col = "#f7614e", border = T, lwd = 0.7)

points(background_chr15$BP, background_chr15$log, pch=16, col = "grey")

outliers <- subset(background_chr15, background_chr15$EMP <= as.numeric(quantile(background_points$EMP1, 0.01)[1]))
outliers_x <- outliers$BP
outliers_y <- outliers$log
points(outliers_x, outliers_y, col="black", pch = 16)
dev.off()

#output data files for Simon
write.table(background_chr15, file = "../plink_gwas/all_background_points_chr15.csv", col.names = T, row.names = F, quote = F, sep = ",")
write.table(outliers, file = "../plink_gwas/outliers_background_points_chr15.csv", col.names = T, row.names = F, quote = F, sep = ",")

background_peaks_ordered <- background_points[order(background_points$P),]

background_peaks_ordered_top10 <- background_peaks_ordered[1:10,]
  
###forewing
background_points <- read.csv("forewing_plot_emp_total.txt", sep = " ")

background_points <- na.omit(background_points)
background_points$log <- -log10(background_points$P)
background_points <- subset(background_points, background_points$log > 1)
background_points$chrom <- paste("chr", background_points$CHR, sep = "")

background_chr15 <- subset(background_points, background_points$CHR == 17)
background_chr15pt1 <- subset(background_chr15, background_chr15$BP < 7.8e6)
background_chr15pt2 <- subset(background_chr15, background_chr15$BP > 14.8e6)
background_chr15 <- as.data.frame(rbind(background_chr15pt1, background_chr15pt2))

tiff("plots/forewing_chr15_new.tiff", height=4, width=18, units="in", res=300, compression="lzw")
plot(background_chr15$BP,background_chr15$log,pch=16,cex=1,
     ylab="-log10(p-value)",
     xlab="bp chr15",
     axes=T, 
     ylim=c(0,(3+max(background_chr15$log))),
     col="black",
     type = 'n')

rect(xleft = 4035981, ybottom = (2+max(background_chr15$log)), 
     xright = 5322257, ytop = (5+max(background_chr15$log)), 
     col = "#016337", border = T, lwd = 0.7)
rect(xleft = 5323805, ybottom = (2+max(background_chr15$log)), 
     xright = 6220875, ytop = (5+max(background_chr15$log)), 
     col = "#016337", border = T, lwd = 0.7)
rect(xleft = 6251636, ybottom = (2+max(background_chr15$log)), 
     xright = 7829094, ytop = (5+max(background_chr15$log)), 
     col = "#a875ee", border = T, lwd = 0.7)
rect(xleft = 14803486, ybottom = (2+max(background_chr15$log)), 
     xright = 15299732, ytop = (5+max(background_chr15$log)), 
     col = "#f7614e", border = T, lwd = 0.7)

points(background_chr15$BP, background_chr15$log, pch=16, col = "grey")

outliers <- subset(background_chr15, background_chr15$EMP <= as.numeric(quantile(background_points$EMP1, 0.01)[1]))
outliers_x <- outliers$BP
outliers_y <- outliers$log
points(outliers_x, outliers_y, col="black", pch = 16)
dev.off()

#output data files for Simon
write.table(background_chr15, file = "../plink_gwas/all_forewing_points_chr15.csv", col.names = T, row.names = F, quote = F, sep = ",")
write.table(outliers, file = "../plink_gwas/outliers_forewing_points_chr15.csv", col.names = T, row.names = F, quote = F, sep = ",")

forewing_peaks_ordered <- background_points[order(background_points$P),]

forewing_peaks_ordered_top10 <- forewing_peaks_ordered[1:10,]

###white
background_points <- read.csv("white_plot_emp_total.txt", sep = " ")

background_points <- na.omit(background_points)
background_points$log <- -log10(background_points$P)
background_points <- subset(background_points, background_points$log > 1)
background_points$chrom <- paste("chr", background_points$CHR, sep = "")

background_chr15 <- subset(background_points, background_points$CHR == 17)
background_chr15pt1 <- subset(background_chr15, background_chr15$BP < 7.8e6)
background_chr15pt2 <- subset(background_chr15, background_chr15$BP > 14.8e6)
background_chr15 <- as.data.frame(rbind(background_chr15pt1, background_chr15pt2))

tiff("plots/white_chr15_new.tiff", height=4, width=18, units="in", res=300, compression="lzw")
plot(background_chr15$BP,background_chr15$log,pch=16,cex=1,
     ylab="-log10(p-value)",
     xlab="bp chr15",
     axes=T, 
     ylim=c(0,(3+max(background_chr15$log))),
     col="black",
     type = 'n')

rect(xleft = 4035981, ybottom = (2+max(background_chr15$log)), 
     xright = 5322257, ytop = (4+max(background_chr15$log)), 
     col = "#016337", border = T, lwd = 0.7)
rect(xleft = 5323805, ybottom = (2+max(background_chr15$log)), 
     xright = 6220875, ytop = (4+max(background_chr15$log)), 
     col = "#016337", border = T, lwd = 0.7)
rect(xleft = 6251636, ybottom = (2+max(background_chr15$log)), 
     xright = 7829094, ytop = (4+max(background_chr15$log)), 
     col = "#a875ee", border = T, lwd = 0.7)
rect(xleft = 14803486, ybottom = (2+max(background_chr15$log)), 
     xright = 15299732, ytop = (4+max(background_chr15$log)), 
     col = "#f7614e", border = T, lwd = 0.7)

points(background_chr15$BP, background_chr15$log, pch=16, col = "grey")

outliers <- subset(background_chr15, background_chr15$EMP <= as.numeric(quantile(background_points$EMP1, 0.01)[1]))
outliers_x <- outliers$BP
outliers_y <- outliers$log
points(outliers_x, outliers_y, col="black", pch = 16)
dev.off()

#output data files for Simon
write.table(background_chr15, file = "../plink_gwas/all_white_points_chr15.csv", col.names = T, row.names = F, quote = F, sep = ",")
write.table(outliers, file = "../plink_gwas/outliers_white_points_chr15.csv", col.names = T, row.names = F, quote = F, sep = ",")

##annotation at peaks:
#FOREWING
#load annotated genes
gene_df <- read.csv("../plink_gwas/annotation/Dchry2.2.fa.masked_annotation_a001.sequences.tidy.gff3.gene.sorted.bed", sep = "\t", head = F)
gene_df_chr15 <- subset(gene_df, gene_df$V1 == "contig15.1")

gene_df_chr15$start <- c()
gene_df_chr15$end <- c()

for(i in 1:nrow(gene_df_chr15)){
  S <- gene_df_chr15$V2[i]
  E <- gene_df_chr15$V3[i]
  
  if(E>S){
    gene_df_chr15$start[i] <- S
    gene_df_chr15$end[i] <- E
  }
  
  if(E<S){
    gene_df_chr15$start[i] <- E
    gene_df_chr15$end[i] <- S
  }
}

gene_df_chr15$gene_overlap <- "empty"
gene_df_chr15$gene_overlap_name <- "empty"
gene_df_chr15$gene_overlap_log <- "empty"
gene_df_chr15$gene_overlap_order <- "empty"


for(i in 1:nrow(gene_df_chr15)){
  gstart <- gene_df_chr15$start[i]
  gend <- gene_df_chr15$end[i]
  
  set <- subset(forewing_peaks_ordered_top10, forewing_peaks_ordered_top10$BP > gstart & forewing_peaks_ordered_top10$BP < gend)
  
  if(nrow(set) > 0){
    gene_df_chr15$gene_overlap[i] <- "yes"
    gene_df_chr15$gene_overlap_name[i] <- set$BP
    gene_df_chr15$gene_overlap_log[i] <- set$log
    gene_df_chr15$gene_overlap_order[i] <- set$order
    
  }
    
}

gene_df_chr15$plot <- 1
plot(gene_df_chr15$V2, gene_df_chr15$plot, xlim = c(7100000, 7700000), frame.plot = F)
forewing_peaks_ordered_top10$plot <- 1.2
points(forewing_peaks_ordered_top10$BP, forewing_peaks_ordered_top10$plot, col = "red")

colnames(gene_df_chr15)
forewing_peaks_ordered_top10$order <- 1:10

empty_df <- gene_df_chr15[1:10,]
empty_df[1:10,] <- NA

for(i in 1:nrow(forewing_peaks_ordered_top10)){
  empty_df$V1[i] <- forewing_peaks_ordered_top10$log[i]
  empty_df$V2[i] <- NA
  empty_df$V3[i] <- NA
  empty_df$V4[i] <- forewing_peaks_ordered_top10$order[i]
  empty_df$V5[i] <- NA
  empty_df$V6[i] <- NA
  empty_df$V7[i] <- NA
  empty_df$V8[i] <- NA
  empty_df$V9[i] <- NA
  empty_df$V10[i] <- NA
  empty_df$start[i] <- forewing_peaks_ordered_top10$BP[i]
  empty_df$end[i] <- forewing_peaks_ordered_top10$BP[i]
  empty_df$plot[i] <- 1.2
}

mega_df <- rbind(gene_df_chr15, empty_df)

#BACKGROUND
#load annotated genes
gene_df <- read.csv("../plink_gwas/annotation/Dchry2.2.fa.masked_annotation_a001.sequences.tidy.gff3.gene.sorted.bed", sep = "\t", head = F)
gene_df_chr15 <- subset(gene_df, gene_df$V1 == "contig15.1")

gene_df_chr15$start <- c()
gene_df_chr15$end <- c()

for(i in 1:nrow(gene_df_chr15)){
  S <- gene_df_chr15$V2[i]
  E <- gene_df_chr15$V3[i]
  
  if(E>S){
    gene_df_chr15$start[i] <- S
    gene_df_chr15$end[i] <- E
  }
  
  if(E<S){
    gene_df_chr15$start[i] <- E
    gene_df_chr15$end[i] <- S
  }
}

gene_df_chr15$gene_overlap <- "empty"
gene_df_chr15$gene_overlap_name <- "empty"
gene_df_chr15$gene_overlap_log <- "empty"
gene_df_chr15$gene_overlap_order <- "empty"


for(i in 1:nrow(gene_df_chr15)){
  gstart <- gene_df_chr15$start[i]
  gend <- gene_df_chr15$end[i]
  
  set <- subset(background_peaks_ordered_top10, background_peaks_ordered_top10$BP > gstart & background_peaks_ordered_top10$BP < gend)
  
  if(nrow(set) > 0){
    gene_df_chr15$gene_overlap[i] <- "yes"
    gene_df_chr15$gene_overlap_name[i] <- set$BP
    gene_df_chr15$gene_overlap_log[i] <- set$log
    gene_df_chr15$gene_overlap_order[i] <- set$order
    
  }
  
}

############ SUPP FIGURE
tiff("plots/background_forewing_chr15.tiff", height=10, width=18, units="in", res=300, compression="lzw")
par(mfrow=c(2,1),mar=c(0.01,0.01,0.01,0.01),oma=c(3,2,2,1))

###background
background_points <- read.csv("background_plot_emp_total.txt", sep = " ")

background_points <- na.omit(background_points)
background_points$log <- -log10(background_points$P)
background_points <- subset(background_points, background_points$log > 1)
background_points$chrom <- paste("chr", background_points$CHR, sep = "")

background_chr15 <- subset(background_points, background_points$CHR == 17)

plot(background_chr15$BP,background_chr15$log,pch=16,cex=1,
     ylab="-log10(p-value)",
     xaxt= 'n',
     ylim=c(0,(3+max(background_chr15$log))),
     xlim=c(0, 17804635),
     col="black",
     type = 'n',
     frame.plot = FALSE)

rect(xleft = 4035981, ybottom = (2+max(background_chr15$log)), 
     xright = 5322257, ytop = (3+max(background_chr15$log)), 
     col = "#016337", border = F, lwd = 0.7)
rect(xleft = 5323805, ybottom = (2+max(background_chr15$log)),
     xright = 6220875, ytop = (3+max(background_chr15$log)), 
     col = "#016337", border = F, lwd = 0.7)
rect(xleft = 6251636, ybottom = (2+max(background_chr15$log)), 
     xright = 7829094, ytop = (3+max(background_chr15$log)), 
     col = "#a875ee", border = F, lwd = 0.7)
rect(xleft = 14803486, ybottom = (2+max(background_chr15$log)), 
     xright = 15299732, ytop = (3+max(background_chr15$log)), 
     col = "#f7614e", border = F, lwd = 0.7)

seg_df <- read.csv("../pi/50kb_core_pops/Dchry2.2_region_coordinates.tsv", sep = "\t")
seg_df <- subset(seg_df, seg_df$Region == 3)
seg_df$colour <- "empty"
for(i in 1:nrow(seg_df)){
  if(seg_df$Region[i] == 3){
    seg_df$colour[i] <- "limegreen"
  }
}

for(i in 1:nrow(seg_df)){
    a <- seg_df$Qstart[i]
    b <- (2+max(background_chr15$log))
    c <- seg_df$Qend[i]
    d <- (3+max(background_chr15$log))
    rect(a,b,c,d , col = seg_df$colour[i], border = F)
}

points(background_chr15$BP, background_chr15$log, pch=16, col = "grey")

outliers <- subset(background_chr15, background_chr15$EMP <= as.numeric(quantile(background_points$EMP1, 0.01)[1]))
outliers_x <- outliers$BP
outliers_y <- outliers$log
points(outliers_x, outliers_y, col="black", pch = 16)

v_high <- subset(background_chr15, background_chr15$log > 50)

points(v_high$BP, v_high$log, pch=16, col = "red")

###forewing
background_points <- read.csv("forewing_plot_emp_total.txt", sep = " ")

background_points <- na.omit(background_points)
background_points$log <- -log10(background_points$P)
background_points <- subset(background_points, background_points$log > 1)
background_points$chrom <- paste("chr", background_points$CHR, sep = "")

background_chr15 <- subset(background_points, background_points$CHR == 17)

plot(background_chr15$BP,background_chr15$log,pch=16,cex=1,
     ylab="-log10(p-value)",
     xlab="bp chr15",
     xlim=c(0, 17804635),
     xaxt="n", 
     ylim=c(0,(3+max(background_chr15$log))),
     col="black",
     type = 'n',
     frame.plot = FALSE)

axis(side = 1, at = c(0, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000, 11000000, 12000000, 13000000, 14000000, 15000000, 16000000, 17000000, 18000000), labels = c("0", "1", "2", "3", "4", "5", "6","7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))

points(background_chr15$BP, background_chr15$log, pch=16, col = "grey")

outliers <- subset(background_chr15, background_chr15$EMP <= as.numeric(quantile(background_points$EMP1, 0.01)[1]))
outliers_x <- outliers$BP
outliers_y <- outliers$log
points(outliers_x, outliers_y, col="black", pch = 16)

v_high <- subset(background_chr15, background_chr15$log > 70)

points(v_high$BP, v_high$log, pch=16, col = "red")

dev.off()


