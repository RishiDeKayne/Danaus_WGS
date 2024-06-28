#prepare populations for pi estimates
setwd("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/pi/50kb_core_pops/")
background <- read.csv("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/background/DC174_background.csv", header = T, sep = ",")

country_list <- c("Kenya", "RSA", "StHelena", "Nigeria", "Ghana", "Spain", "Tunisia", "Italy", "Philippines", "Rwanda")
colour_list <- c("tomato", "lightseagreen", "royalblue2", "springgreen4", "yellow3", "wheat1", "goldenrod", "hotpink", "mediumorchid1", "tomato4")

bacakground_uniq_col <- unique(background$colour)

popfile_df <- as.data.frame(cbind(background$indiv, background$country))
write.table(popfile_df, file = "popsfile.txt", col.names = F, row.names = F, quote = F, sep = "\t")

background_regions <- background
background_regions$region_colour <- "EMPTY"
background_regions$region_name <- c()
for (i in 1:nrow(background_regions)){
  col_type <- as.character(background_regions$region)[i]
  if(col_type == "darkgrey"){
    background_regions$region_colour[i] <- "black"
    background_regions$region_name[i] <- "crossed"
  }
  if(col_type == "goldenrod"){
    background_regions$region_colour[i] <- "springgreen4"
    background_regions$region_name[i] <- "med"
  }
  if(col_type == "royalblue2"){
    background_regions$region_colour[i] <- "#2143d1"
    background_regions$region_name[i] <- "orientis"
  }
  if(col_type == "springgreen4"){
    background_regions$region_colour[i] <- "#bc4754"
    background_regions$region_name[i] <- "chrysippus"
  }
  if(col_type == "tomato"){
    background_regions$region_colour[i] <- "#ffac07"
    background_regions$region_name[i] <- "klugii"
  }
  if(col_type == "tomato4"){
    background_regions$region_colour[i] <- "darkgrey"
    background_regions$region_name[i] <- "hybrid"
  }
}

popfile_df_region <- as.data.frame(cbind(background_regions$indiv, background_regions$region_name))
write.table(popfile_df_region, file = "popsfile_region.txt", col.names = F, row.names = F, quote = F, sep = "\t")

#### RUN BASH SCRIPTS TO GENERATE PI ETC.

##now plot output

##### by broad region

##and for regions: 
##now to make plot

for (i in 1:nrow(background_regions)){
  col_type <- as.character(background_regions$region)[i]
  if(col_type == "darkgrey"){
    background_regions$region_colour[i] <- "black"
    background_regions$region_name[i] <- "crossed"
  }
  if(col_type == "goldenrod"){
    background_regions$region_colour[i] <- "springgreen4"
    background_regions$region_name[i] <- "med"
  }
  if(col_type == "royalblue2"){
    background_regions$region_colour[i] <- "#2143d1"
    background_regions$region_name[i] <- "orientis"
  }
  if(col_type == "springgreen4"){
    background_regions$region_colour[i] <- "#bc4754"
    background_regions$region_name[i] <- "chrysippus"
  }
  if(col_type == "tomato"){
    background_regions$region_colour[i] <- "#ffac07"
    background_regions$region_name[i] <- "klugii"
  }
  if(col_type == "tomato4"){
    background_regions$region_colour[i] <- "darkgrey"
    background_regions$region_name[i] <- "hybrid"
  }
}

group_list <- c("WAT-Klugii", "TSW-Orientis", "NGA-Chrysippus", "Karamu", "recomb")
group_colour_list <- c("#ffac07", "#2143d1", "#bc4754", "hotpink", "black")

#get counts of individuals for each population

chromosomes <- read.csv("../../gemma/Dchry2.2.fa.fai", head = FALSE, sep = "\t")
chromosomes <- chromosomes[1:41,]
chromosomes$numb <- gsub("contig", "", chromosomes$V1)


for(j in 1:length(chromosomes$V1)){
  #get file name from contig name
  if(j == 1){
    file_name <- paste(chromosomes$V1[j], ".hap.div.regions.output.csv.gz", sep = "")
  }
  if(j > 1){
    file_name <- paste(chromosomes$V1[j], ".div.regions.output.csv.gz", sep = "")
  }
  #get contig name and length
  contig_name <- as.character(chromosomes$V1[1])
  contig_length <- as.numeric(chromosomes$V2[1])
  #read in data
  data <- read.csv(file_name)
  #create data with no Nans so we can calc plotting dimensions
  plot_parameters_df <- na.omit(data)
  #make list of all pi values to get max value
  max_val_col <- as.vector(rbind(plot_parameters_df$pi_WAT, plot_parameters_df$pi_TSW, plot_parameters_df$pi_NGA, 
                                 plot_parameters_df$pi_karamu, plot_parameters_df$pi_recomb))
  
  #plot empty plot with max value as ylim
  plot(plot_parameters_df$mid, plot_parameters_df$pi_WAT, type = 'n', ylim = c(0, max(max_val_col)),
       ylab = "pi in 100kb window", xlab = "mid bp", main = file_name)
  #add legend
  legend('topleft', 
         legend = c(group_list), 
         col = c(group_colour_list),
         lty = 1,
         lwd = 2,
         pt.cex = 1,
         cex = 0.5) 
  
  #now plot data - run a loop going country by country
  for (c in 1:length(group_list)){
    #extract the country values from data 
    country_data_vals <- data[,(5+c)]
    country_data_x <- data$mid
    #make a country df with values and midpoints
    country_data <- as.data.frame(cbind(country_data_vals, country_data_x))
    #at this point remove any missing rows with nans
    country_data <- na.omit(country_data)
    #add index column
    country_data$index <- 1:nrow(country_data)
    #now smoothing using the country-specific values
    loessMod <- loess(country_data$country_data_vals ~ index, data=country_data, span=(10000000/contig_length))
    smoothed <- predict(loessMod) 
    lines(country_data$country_data_vals, x=country_data$country_data_x, col="black", lwd = 2)
    lines(smoothed, x=country_data$country_data_x, col=group_colour_list[c], lwd = 2)
  }
}

filename <- "pi_plot_regions.pdf"
print(filename)

chromlist <- as.list(as.character(chromosomes$V1))
b <- as.vector(chromlist)
a <- b

pdf(file=filename, width=15, height=4)

par(mfrow=c(1,1),mar=c(0.01,0.01,0.01,0.01),oma=c(4,2,2,1))

l <- as.vector(chromosomes$V2)
l <- l[1:(length(chromosomes$V1))]

#set layout parameters specifying 5 rows, one for each comarison and the widths of the chromosomes in each case
layout(matrix(seq(1:length(chromosomes$V1)),nrow=1,ncol=length(chromosomes$V1),byrow=TRUE),width=c(l))

altcols <- c("grey90", "grey75")
altcols <- rep(altcols, as.integer((length(a)+1)/2))

altcol_numb <- c(1,2,2,1,2,2,1,
                 2,1,2,1,2,1,2,
                 1,2,1,2,1,1,2,
                 2,1,1,2,2,1,2,
                 2,1,1,2,1,2,2,
                 1,1,2,2,1,2)

for(j in 1:length(chromosomes$V1)){
  #get file name from contig name
  if(j == 1){
    file_name <- paste(chromosomes$V1[j], ".hap.div.regions.output.csv.gz", sep = "")
  }
  if(j > 1){
    file_name <- paste(chromosomes$V1[j], ".div.regions.output.csv.gz", sep = "")
  }  #get contig name and length
  contig_name <- as.character(chromosomes$V1[j])
  contig_length <- as.numeric(chromosomes$V2[j])
  #read in data
  data <- read.csv(file_name)
  #create data with no Nans so we can calc plotting dimensions
  plot_parameters_df <- na.omit(data)
  #make list of all pi values to get max value
  max_val_col <- as.vector(rbind(plot_parameters_df$pi_WAT, plot_parameters_df$pi_TSW, plot_parameters_df$pi_NGA, 
                                 plot_parameters_df$pi_karamu, plot_parameters_df$pi_recomb))
  
  rect_col <- altcols[altcol_numb[j]]
  #plot empty plot with ylim set by looking at individual plots
  plot(data$mid, data$pi_WAT, type = 'n', 
       ylim = c(0, 0.06), xaxt="n",axes=F, ylab = contig_name)
  rect(xleft = min(plot_parameters_df$mid), ybottom = 0, 
       xright = max(plot_parameters_df$mid), ytop = 0.06, 
       col = rect_col, border = F)
  if (j==1){
    axis(side=2)
    legend(x = 10000, y = 0.0575, 
           legend = c(group_list), 
           col = c(group_colour_list),
           lty = 1,
           lwd = 2,
           pt.cex = 1,
           cex = 0.55,
           bty = "n")}
  
  if(j==17){
    rect(xleft = 4035981, ybottom = 0, 
         xright = 5322257, ytop = 0.06, 
         col = "#016337", border = T, lwd = 0.7)
    rect(xleft = 5323805, ybottom = 0, 
         xright = 6220875, ytop = 0.06, 
         col = "#016337", border = T, lwd = 0.7)
    rect(xleft = 6251636, ybottom = 0, 
         xright = 7829094, ytop = 0.06, 
         col = "#a875ee", border = T, lwd = 0.7)
    rect(xleft = 14803486, ybottom = 0, 
         xright = 15299732, ytop = 0.06, 
         col = "#f7614e", border = T, lwd = 0.7)
  }
  
  #now plot data - run a loop going country by country
  for (c in 1:length(group_list)){
    #extract the country values from data 
    country_data_vals <- data[,(5+c)]
    country_data_x <- data$mid
    #make a country df with values and midpoints
    country_data <- as.data.frame(cbind(country_data_vals, country_data_x))
    #at this point remove any missing rows with nans
    country_data <- na.omit(country_data)
    #add index column
    country_data$index <- 1:nrow(country_data)
    #now smoothing using the country-specific values
    loessMod <- loess(country_data$country_data_vals ~ index, data=country_data, span=0.3)
    smoothed <- predict(loessMod) 
    #lines(country_data$country_data_vals, x=country_data$country_data_x, col="black", lwd = 2)
    lines(smoothed, x=country_data$country_data_x, col=group_colour_list[c], lwd = 1)
  }
  
  text(x = contig_length/2, 0.0057, chromosomes$numb[j], cex = 0.5)
}

dev.off()

#####Dxy plot
filename <- "dxy_plot_regions.pdf"
print(filename)

chromlist <- as.list(as.character(chromosomes$V1))
b <- as.vector(chromlist)
a <- b

pdf(file=filename, width=15, height=4)

par(mfrow=c(1,1),mar=c(0.01,0.01,0.01,0.01),oma=c(4,2,2,1))

l <- as.vector(chromosomes$V2)
l <- l[1:(length(chromosomes$V1))]

#set layout parameters specifying 5 rows, one for each comarison and the widths of the chromosomes in each case
layout(matrix(seq(1:length(chromosomes$V1)),nrow=1,ncol=length(chromosomes$V1),byrow=TRUE),width=c(l))

altcols <- c("grey90", "grey75")
altcols <- rep(altcols, as.integer((length(a)+1)/2))

altcol_numb <- c(1,2,2,1,2,2,1,
                 2,1,2,1,2,1,2,
                 1,2,1,2,1,1,2,
                 2,1,1,2,2,1,2,
                 2,1,1,2,1,2,2,
                 1,1,2,2,1,2)

for(j in 1:length(chromosomes$V1)){
  #get file name from contig name
  if(j == 1){
    file_name <- paste(chromosomes$V1[1], ".hap.div.regions.output.csv.gz", sep = "")
  }
  if(j > 1){
    file_name <- paste(chromosomes$V1[j], ".div.regions.output.csv.gz", sep = "")
  }  #get contig name and length
  contig_name <- as.character(chromosomes$V1[j])
  contig_length <- as.numeric(chromosomes$V2[j])
  #read in data
  data <- read.csv(file_name)
  #create data with no Nans so we can calc plotting dimensions
  plot_parameters_df <- na.omit(data)
  #make list of all pi values to get max value
  max_val_col <- as.vector(rbind(plot_parameters_df$dxy_TSW_NGA, plot_parameters_df$dxy_WAT_TSW, plot_parameters_df$dxy_WAT_NGA))
  
  rect_col <- altcols[altcol_numb[j]]
  #plot empty plot with ylim set by looking at individual plots
  plot(data$mid, data$dxy_TSW_NGA, type = 'n', 
       ylim = c(0, 0.06), xaxt="n",axes=F, ylab = contig_name)
  rect(xleft = min(plot_parameters_df$mid), ybottom = 0, 
       xright = max(plot_parameters_df$mid), ytop = 1, 
       col = rect_col, border = F)
  
  column_list <- c(15, 11, 12)
  colour_list <- c("purple4", "forestgreen", "orange")
  
  if (j==1){
    axis(side=2)
    legend(x = 10000, y = 0.0575, 
           legend = c("dxy_orientis_chrysippus", "dxy_orientis_klugii", "dxy_chrysippus_klugii"), 
           col = c(colour_list),
           lty = 1,
           lwd = 2,
           pt.cex = 1,
           cex = 0.4,
           bty = "n")}
  
  if(j==17){
    rect(xleft = 4035981, ybottom = 0, 
         xright = 5322257, ytop = 0.06, 
         col = "#016337", border = T, lwd = 0.7)
    rect(xleft = 5323805, ybottom = 0, 
         xright = 6220875, ytop = 0.06, 
         col = "#016337", border = T, lwd = 0.7)
    rect(xleft = 6251636, ybottom = 0, 
         xright = 7829094, ytop = 0.06, 
         col = "#a875ee", border = T, lwd = 0.7)
    rect(xleft = 14803486, ybottom = 0, 
         xright = 15299732, ytop = 0.06, 
         col = "#f7614e", border = T, lwd = 0.7)
  }
  
  #now plot data - run a loop going country by country
  for (c in 1:length(column_list)){
    #extract the country values from data 
    country_data_vals <- data[,(column_list[c])]
    country_data_x <- data$mid
    #make a country df with values and midpoints
    country_data <- as.data.frame(cbind(country_data_vals, country_data_x))
    #at this point remove any missing rows with nans
    country_data <- na.omit(country_data)
    #add index column
    country_data$index <- 1:nrow(country_data)
    #now smoothing using the country-specific values
    loessMod <- loess(country_data$country_data_vals ~ index, data=country_data, span=0.3)
    smoothed <- predict(loessMod) 
    #lines(country_data$country_data_vals, x=country_data$country_data_x, col="black", lwd = 2)
    lines(smoothed, x=country_data$country_data_x, col=colour_list[c], lwd = 1)
  }
  
  text(x = contig_length/2, 0.0057, chromosomes$numb[j], cex = 0.5)
}

dev.off()

#####FST plot
filename <- "fst_plot_regions.pdf"
print(filename)

chromlist <- as.list(as.character(chromosomes$V1))
b <- as.vector(chromlist)
a <- b

pdf(file=filename, width=15, height=4)

par(mfrow=c(1,1),mar=c(0.01,0.01,0.01,0.01),oma=c(4,2,2,1))

l <- as.vector(chromosomes$V2)
l <- l[1:(length(chromosomes$V1))]

#set layout parameters specifying 5 rows, one for each comarison and the widths of the chromosomes in each case
layout(matrix(seq(1:length(chromosomes$V1)),nrow=1,ncol=length(chromosomes$V1),byrow=TRUE),width=c(l))

altcols <- c("grey90", "grey75")
altcols <- rep(altcols, as.integer((length(a)+1)/2))

altcol_numb <- c(1,2,2,1,2,2,1,
                 2,1,2,1,2,1,2,
                 1,2,1,2,1,1,2,
                 2,1,1,2,2,1,2,
                 2,1,1,2,1,2,2,
                 1,1,2,2,1,2)

for(j in 1:length(chromosomes$V1)){
  #get file name from contig name
  if(j == 1){
    file_name <- paste(chromosomes$V1[1], ".hap.div.regions.output.csv.gz", sep = "")
  }
  if(j > 1){
    file_name <- paste(chromosomes$V1[j], ".div.regions.output.csv.gz", sep = "")
  }  #get contig name and length
  contig_name <- as.character(chromosomes$V1[j])
  contig_length <- as.numeric(chromosomes$V2[j])
  #read in data
  data <- read.csv(file_name)
  #create data with no Nans so we can calc plotting dimensions
  plot_parameters_df <- na.omit(data)
  #make list of all pi values to get max value
  max_val_col <- as.vector(rbind(plot_parameters_df$Fst_TSW_NGA, plot_parameters_df$Fst_WAT_TSW, plot_parameters_df$Fst_WAT_NGA))
  
  rect_col <- altcols[altcol_numb[j]]
  #plot empty plot with ylim set by looking at individual plots
  plot(data$mid, data$Fst_TSW_NGA, type = 'n', 
       ylim = c(0, 0.65), xaxt="n",axes=F, ylab = contig_name)
  rect(xleft = min(plot_parameters_df$mid), ybottom = 0, 
       xright = max(plot_parameters_df$mid), ytop = 1, 
       col = rect_col, border = F)
  
  column_list <- c(25, 21, 22)
  colour_list <- c("purple4", "forestgreen", "orange")
  
  if (j==1){
    axis(side=2)
    legend(x = 10000, y = 0.0575, 
           legend = c("fst_orientis_chrysippus", "fst_orientis_klugii", "fst_chrysippus_klugii"), 
           col = c(colour_list),
           lty = 1,
           lwd = 2,
           pt.cex = 1,
           cex = 0.4,
           bty = "n")}
  
  if(j==17){
    rect(xleft = 4035981, ybottom = 0, 
         xright = 5322257, ytop = 0.65, 
         col = "#016337", border = T, lwd = 0.7)
    rect(xleft = 5323805, ybottom = 0, 
         xright = 6220875, ytop = 0.65, 
         col = "#016337", border = T, lwd = 0.7)
    rect(xleft = 6251636, ybottom = 0, 
         xright = 7829094, ytop = 0.65, 
         col = "#a875ee", border = T, lwd = 0.7)
    rect(xleft = 14803486, ybottom = 0, 
         xright = 15299732, ytop = 0.65, 
         col = "#f7614e", border = T, lwd = 0.7)
  }
  
  #now plot data - run a loop going country by country
  for (c in 1:length(column_list)){
    #extract the country values from data 
    country_data_vals <- data[,(column_list[c])]
    country_data_x <- data$mid
    #make a country df with values and midpoints
    country_data <- as.data.frame(cbind(country_data_vals, country_data_x))
    #at this point remove any missing rows with nans
    country_data <- na.omit(country_data)
    country_data[country_data<0] <- 0
    #add index column
    country_data$index <- 1:nrow(country_data)
    #now smoothing using the country-specific values
    loessMod <- loess(country_data$country_data_vals ~ index, data=country_data, span=0.3)
    smoothed <- predict(loessMod) 
    #lines(country_data$country_data_vals, x=country_data$country_data_x, col="black", lwd = 2)
    lines(smoothed, x=country_data$country_data_x, col=colour_list[c], lwd = 1)
  }
  
  text(x = contig_length/2, 0.07, chromosomes$numb[j], cex = 0.5)
}

dev.off()

#####FST plot
filename <- "fst_plot_regions_FIG1.pdf"
print(filename)

chromlist <- as.list(as.character(chromosomes$V1))
b <- as.vector(chromlist)
a <- b

pdf(file=filename, width=15, height=4)

par(mfrow=c(1,1),mar=c(0.01,0.05,0.05,0.01),oma=c(4,3,3,1))

l <- as.vector(chromosomes$V2)
l <- l[1:(length(chromosomes$V1))]

#set layout parameters specifying 5 rows, one for each comarison and the widths of the chromosomes in each case
layout(matrix(seq(1:length(chromosomes$V1)),nrow=1,ncol=length(chromosomes$V1),byrow=TRUE),width=c(l))

altcols <- c("grey90", "grey75")
altcols <- rep(altcols, as.integer((length(a)+1)/2))

altcol_numb <- c(1,2,2,1,2,2,1,
                 2,1,2,1,2,1,2,
                 1,2,1,2,1,1,2,
                 2,1,1,2,2,1,2,
                 2,1,1,2,1,2,2,
                 1,1,2,2,1,2)

for(j in 1:length(chromosomes$V1)){
  #get file name from contig name
  if(j == 1){
    file_name <- paste(chromosomes$V1[1], ".hap.div.regions.output.csv.gz", sep = "")
  }
  if(j > 1){
    file_name <- paste(chromosomes$V1[j], ".div.regions.output.csv.gz", sep = "")
  }  #get contig name and length
  contig_name <- as.character(chromosomes$V1[j])
  contig_length <- as.numeric(chromosomes$V2[j])
  #read in data
  data <- read.csv(file_name)
  #create data with no Nans so we can calc plotting dimensions
  plot_parameters_df <- na.omit(data)
  #make list of all pi values to get max value
  max_val_col <- as.vector(rbind(plot_parameters_df$Fst_TSW_NGA, plot_parameters_df$Fst_WAT_TSW, plot_parameters_df$Fst_WAT_NGA))
  
  rect_col <- altcols[altcol_numb[j]]
  #plot empty plot with ylim set by looking at individual plots
  plot(data$mid, data$Fst_TSW_NGA, type = 'n', 
       ylim = c(0, 0.65), xaxt="n",axes=F, ylab = contig_name)
  rect(xleft = min(plot_parameters_df$mid), ybottom = 0, 
       xright = max(plot_parameters_df$mid), ytop = 1, 
       col = rect_col, border = F)
  
  column_list <- c(25, 21, 22)
  colour_list <- c("orchid3", "seagreen4", "skyblue4")
  
  if (j==1){
    axis(side=2)
    legend(x = 10000, y = 0.0575, 
           legend = c("fst_orientis_chrysippus", "fst_orientis_klugii", "fst_chrysippus_klugii"), 
           col = c(colour_list),
           lty = 1,
           lwd = 2,
           pt.cex = 1,
           cex = 0.4,
           bty = "n")}
  
  #now plot data - run a loop going country by country
  for (c in 1:length(column_list)){
    #extract the country values from data 
    country_data_vals <- data[,(column_list[c])]
    country_data_x <- data$mid
    #make a country df with values and midpoints
    country_data <- as.data.frame(cbind(country_data_vals, country_data_x))
    #at this point remove any missing rows with nans
    country_data <- na.omit(country_data)
    country_data[country_data<0] <- 0
    #add index column
    country_data$index <- 1:nrow(country_data)
    #now smoothing using the country-specific values
    loessMod <- loess(country_data$country_data_vals ~ index, data=country_data, span=0.3)
    smoothed <- predict(loessMod) 
    #lines(country_data$country_data_vals, x=country_data$country_data_x, col="black", lwd = 2)
    lines(smoothed, x=country_data$country_data_x, col=colour_list[c], lwd = 1)
  }
  
  text(x = contig_length/2, 0.07, chromosomes$numb[j], cex = 0.5)
}

dev.off()

#####FST plot
filename <- "fst_plot_regions_FIG1_NOLEGEND.pdf"
print(filename)

chromlist <- as.list(as.character(chromosomes$V1))
b <- as.vector(chromlist)
a <- b

pdf(file=filename, width=15, height=4)

#par(mfrow=c(1,1),mai=c(0.01,0,0.01,0.0),oma=c(4,3,3,1))
par(mfrow=c(1,41),mai=c(0.01,0,0,0.0),oma=c(4,3,3,1))

l <- as.vector(chromosomes$V2)
l <- l[1:(length(chromosomes$V1))]

#set layout parameters specifying 5 rows, one for each comarison and the widths of the chromosomes in each case
#layout(matrix(seq(1:length(chromosomes$V1)),nrow=1,ncol=length(chromosomes$V1),byrow=TRUE),widths=c(l))
layout(matrix(seq(1:length(chromosomes$V1)),nrow=1,ncol=length(chromosomes$V1),byrow=TRUE),width=c(l))

altcols <- c("grey90", "grey75")
altcols <- rep(altcols, as.integer((length(a)+1)/2))

altcol_numb <- c(1,2,2,1,2,2,1,
                 2,1,2,1,2,1,2,
                 1,2,1,2,1,1,2,
                 2,1,1,2,2,1,2,
                 2,1,1,2,1,2,2,
                 1,1,2,2,1,2)

for(j in 1:length(chromosomes$V1)){
  #get file name from contig name
  if(j == 1){
    file_name <- paste(chromosomes$V1[1], ".hap.div.regions.output.csv.gz", sep = "")
  }
  if(j > 1){
    file_name <- paste(chromosomes$V1[j], ".div.regions.output.csv.gz", sep = "")
  }  #get contig name and length
  contig_name <- as.character(chromosomes$V1[j])
  contig_length <- as.numeric(chromosomes$V2[j])
  #read in data
  data <- read.csv(file_name)
  #create data with no Nans so we can calc plotting dimensions
  plot_parameters_df <- na.omit(data)
  #make list of all pi values to get max value
  max_val_col <- as.vector(rbind(plot_parameters_df$Fst_TSW_NGA, plot_parameters_df$Fst_WAT_TSW, plot_parameters_df$Fst_WAT_NGA))
  
  rect_col <- altcols[altcol_numb[j]]
  #plot empty plot with ylim set by looking at individual plots
  #plot(data$mid, data$Fst_TSW_NGA, type = 'n', 
  #     ylim = c(0, 0.65), xaxt="n",axes=F, ylab = contig_name)
  #rect(xleft = min(plot_parameters_df$mid), ybottom = 0, 
  #     xright = max(plot_parameters_df$mid), ytop = 1, 
  #     col = rect_col, border = F)
  
  plot(data$mid, data$Fst_TSW_NGA, type = 'n', xlim = c(0, (l[j])),
       ylim = c(0, 0.65), xaxt="n", axes=F, ylab = contig_name)
  rect(xleft = 0, ybottom = 0, 
       xright = (l[j]), ytop = 1, 
       col = rect_col, border = F)
  
  column_list <- c(25, 21, 22)
  colour_list <- c("orchid3", "seagreen4", "skyblue4")
  
  if (j==1){
    axis(side=2)}
  
  #now plot data - run a loop going country by country
  for (c in 1:length(column_list)){
    #extract the country values from data 
    country_data_vals <- data[,(column_list[c])]
    country_data_x <- data$mid
    #make a country df with values and midpoints
    country_data <- as.data.frame(cbind(country_data_vals, country_data_x))
    #at this point remove any missing rows with nans
    country_data <- na.omit(country_data)
    country_data[country_data<0] <- 0
    #add index column
    country_data$index <- 1:nrow(country_data)
    #now smoothing using the country-specific values
    loessMod <- loess(country_data$country_data_vals ~ index, data=country_data, span=0.2)
    smoothed <- predict(loessMod) 
    #lines(country_data$country_data_vals, x=country_data$country_data_x, col="black", lwd = 2)
    lines(smoothed, x=country_data$country_data_x, col=colour_list[c], lwd = 4)
    #points(c(0), c(0), col = "black", pch = 16)
    #points(c(l[j]), c(0), col = "black", pch = 16)
    
  }
  
}

dev.off()

##FST Plot FIG1

test <- read.csv(paste(chromosomes$V1[1], ".hap.div.regions.output.csv.gz", sep = ""))
test_reduced <- as.data.frame(test[1:1,])
null_DF <- test_reduced
null_DF$numb <- "empty"

for(j in 1:length(chromosomes$V1)){
  #get file name from contig name
  if(j == 1){
    file_name <- paste(chromosomes$V1[1], ".hap.div.regions.output.csv.gz", sep = "")
  }
  if(j > 1){
    file_name <- paste(chromosomes$V1[j], ".div.regions.output.csv.gz", sep = "")
  } 
  #read in data
  data <- read.csv(file_name)
  data$numb <- as.character(j)
  null_DF <- rbind(null_DF, data)}

full_DF <- null_DF
full_DF <- subset(full_DF, full_DF$numb != "empty")
full_DF$numb <- as.numeric(full_DF$numb)
full_DF <- full_DF[order(full_DF$start),]
full_DF <- full_DF[order(full_DF$numb),]
full_DF <- na.omit(full_DF)

column_list <- c(colnames(full_DF)[25], colnames(full_DF)[21], colnames(full_DF)[22])

#chromosomes contains lengths and names
chrom_offset <- cumsum(chromosomes$V2) - chromosomes$V2
names(chrom_offset) <- chromosomes$V1

full_DF$adjusted_mid <- full_DF$mid + chrom_offset[full_DF$scaffold]

filename <- "fst_plot_regions_FIG1_NOLEGEND_REGULAR_SPACING.pdf"
print(filename)
pdf(file=filename, width=15, height=4)
par(mfrow=c(1,1),mar=c(0.01,0.05,0.05,0.01),oma=c(4,3,3,1))

plot(full_DF$adjusted_mid, full_DF$Fst_TSW_NGA, type = "n", xaxt = "n", yaxt = "n",ylim = c(0,0.6), ylab = "FST", xlab = "", frame.plot = FALSE)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6), labels = c("0.0", "0.2", "0.4", "0.6"))
#plot(full_DF$adjusted_mid, full_DF$Fst_TSW_NGA, type = "n", xaxt = "n", ylim = c(0,0.8), ylab = "FST", xlab = "")

altcols <- c("grey90", "grey75")
altcol_numb <- c(1,2,2,1,2,2,1,
                 2,1,2,1,2,1,2,
                 1,2,1,2,1,1,2,
                 2,1,1,2,2,1,2,
                 2,1,1,2,1,2,2,
                 1,1,2,2,1,2)

for(j in 1:length(chromosomes$V1)){
  contig_name <- as.character(chromosomes$V1[j])
  contig_length <- as.numeric(chromosomes$V2[j])
  full_DF_chromosome <- subset(full_DF, full_DF$scaffold == contig_name)
  full_DF_chromosome[full_DF_chromosome<0] <- 0
  
  #add index column
  full_DF_chromosome$index <- 1:nrow(full_DF_chromosome)
  
  ##LINE1 TSW = ORIENTIS + NGA = CHRYSIPPUS
  #now smoothing using the country-specific values
  loessMod1 <- loess(full_DF_chromosome$Fst_TSW_NGA ~ index, data=full_DF_chromosome, span=0.2)
  smoothed1 <- predict(loessMod1) 
  lines(smoothed1, x=full_DF_chromosome$adjusted_mid, col="orchid3", lwd = 4)
  
  ##LINE2 WAT = KLUGII + TSW = ORIENTIS
  #now smoothing using the country-specific values
  loessMod2 <- loess(full_DF_chromosome$Fst_WAT_TSW ~ index, data=full_DF_chromosome, span=0.2)
  smoothed2 <- predict(loessMod2) 
  lines(smoothed2, x=full_DF_chromosome$adjusted_mid, col="seagreen4", lwd = 4)
  
  ##LINE3 WAT = KLUGII + NGA = CHRYSIPPUS
  #now smoothing using the country-specific values
  loessMod3 <- loess(full_DF_chromosome$Fst_WAT_NGA ~ index, data=full_DF_chromosome, span=0.2)
  smoothed3 <- predict(loessMod3) 
  lines(smoothed3, x=full_DF_chromosome$adjusted_mid, col="skyblue4", lwd = 4)
  
  rect((chrom_offset[j]),-0.03,(chrom_offset[j]+contig_length),-0.01, col = altcols[altcol_numb[j]], border = F)
}
dev.off()

#########chr15 plot

chr15_df <- read.csv("contig15.1.div.regions.output.csv.gz")

pi_chr15 <- as.data.frame(cbind(chr15_df$mid, chr15_df$pi_NGA, chr15_df$pi_TSW, chr15_df$pi_WAT, chr15_df$pi_karamu, chr15_df$pi_recomb))
colnames(pi_chr15) <- c("mid", "pi_chrysippus", "pi_orientis", "pi_klugii", "pi_karamu", "pi_recomb")

dxy_chr15 <- as.data.frame(cbind(chr15_df$mid, chr15_df$dxy_TSW_NGA, chr15_df$dxy_WAT_TSW, 
                                 chr15_df$dxy_WAT_NGA, chr15_df$dxy_TSW_karamu, chr15_df$dxy_NGA_karamu,
                                 chr15_df$dxy_WAT_karamu))
colnames(dxy_chr15) <- c("mid", "dxy_OC", "dxy_OK", "dxy_CK", "dxy_OKa", "dxy_CKa", "dxy_KKa")

fst_chr15 <- as.data.frame(cbind(chr15_df$mid, chr15_df$Fst_TSW_NGA, chr15_df$Fst_WAT_TSW, 
                                 chr15_df$Fst_WAT_NGA,chr15_df$Fst_TSW_karamu, chr15_df$Fst_NGA_karamu,
                                 chr15_df$Fst_WAT_karamu))
colnames(fst_chr15) <- c("mid", "fst_OC", "fst_OK", "fst_CK", "fst_OKa", "fst_CKa", "fst_KKa")

par(mfrow=c(4,1), mai = c(0.3, 1, 0.1, 0.1))

plot(pi_chr15$mid, pi_chr15$pi_chrysippus, type = 'n', ylab = "nucleotide diversity (π)", xlab = "", xaxt = 'n', ylim = c(0,0.07))

rect(xleft = 4035981, ybottom = 0, 
     xright = 5322257, ytop = 0.08, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 125), border = T, lwd = 0.7)
rect(xleft = 5323805, ybottom = 0, 
     xright = 6220875, ytop = 0.08, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 125), border = T, lwd = 0.7)
rect(xleft = 6251636, ybottom = 0, 
     xright = 7829094, ytop = 0.08, 
     col = rgb(168,117,238, maxColorValue = 255, alpha = 125), border = T, lwd = 0.7)
rect(xleft = 14803486, ybottom = 0, 
     xright = 15299732, ytop = 0.08, 
     col = rgb(247,97,78, maxColorValue = 255, alpha = 125), border = T, lwd = 0.7)
pilist <- c(2,3,4,5,6)
pi_col_list <- c("#bc4754", "#2143d1", "#ffac07", "hotpink", "darkgrey")
for (c in 1:length(pilist)){
  #extract the country values from data 
  group_data_vals <- pi_chr15[,(pilist[c])]
  group_data_x <- pi_chr15$mid
  #make a country df with values and midpoints
  group_data <- as.data.frame(cbind(group_data_vals, group_data_x))
  lines(x=group_data$group_data_x, y=group_data_vals, col=pi_col_list[c], lwd = 2)
}
legend("topleft", 
       legend = c("π chrysippus", "π orientis", "π klugii", "π karamu", "π recomb"), 
       col = c(pi_col_list),
       lty = 1,
       lwd = 2,
       pt.cex = 1,
       cex = 0.5,
       bty = "n")

plot(dxy_chr15$mid, dxy_chr15$dxy_OC, type = 'n', ylab = "absolute pairwise divergence (Dxy)", xlab = '', xaxt = 'n', ylim = c(0,0.07))
rect(xleft = 4035981, ybottom = 0, 
     xright = 5322257, ytop = 0.08, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 125), border = T, lwd = 0.7)
rect(xleft = 5323805, ybottom = 0, 
     xright = 6220875, ytop = 0.08, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 125), border = T, lwd = 0.7)
rect(xleft = 6251636, ybottom = 0, 
     xright = 7829094, ytop = 0.08, 
     col = rgb(168,117,238, maxColorValue = 255, alpha = 125), border = T, lwd = 0.7)
rect(xleft = 14803486, ybottom = 0, 
     xright = 15299732, ytop = 0.08, 
     col = rgb(247,97,78, maxColorValue = 255, alpha = 125), border = T, lwd = 0.7)
dxylist <- c(2,3,4,5,6,7)
dxy_col_list <- c("purple", "green", "orange", "black", "darkgrey", "lightgrey")
for (c in 1:length(dxylist)){
  #extract the country values from data 
  group_data_vals <- dxy_chr15[,(dxylist[c])]
  group_data_x <- dxy_chr15$mid
  #make a country df with values and midpoints
  group_data <- as.data.frame(cbind(group_data_vals, group_data_x))
  lines(x=group_data$group_data_x, y=group_data_vals, col=dxy_col_list[c], lwd = 2)
}

legend("topleft", 
       legend = c("Dxy orientis-chrysippus", "Dxy orientis-klugii", "Dxy chrysippus-klugii",
                  "Dxy orientis-karamu", "Dxy chrysippus-karamu", "Dxy klugii-karamu"), 
       col = c(dxy_col_list),
       lty = 1,
       lwd = 2,
       pt.cex = 1,
       cex = 0.5,
       bty = "n")

plot(fst_chr15$mid, fst_chr15$fst_OC, type = 'n', ylab = "pairwise divergence (Fst)", xlab = '', xaxt = 'n', ylim = c(0,1))
rect(xleft = 4035981, ybottom = 0, 
     xright = 5322257, ytop = 1.1, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 125), border = T, lwd = 0.7)
rect(xleft = 5323805, ybottom = 0, 
     xright = 6220875, ytop = 1.1, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 125), border = T, lwd = 0.7)
rect(xleft = 6251636, ybottom = 0, 
     xright = 7829094, ytop = 1.1, 
     col = rgb(168,117,238, maxColorValue = 255, alpha = 125), border = T, lwd = 0.7)
rect(xleft = 14803486, ybottom = 0, 
     xright = 15299732, ytop = 1.1, 
     col = rgb(247,97,78, maxColorValue = 255, alpha = 125), border = T, lwd = 0.7)
fstlist <- c(2,3,4,5,6,7)
fst_col_list <- c("purple", "green", "orange", "black", "darkgrey", "lightgrey")
for (c in 1:length(dxylist)){
  #extract the country values from data 
  group_data_vals <- fst_chr15[,(fstlist[c])]
  group_data_x <- fst_chr15$mid
  #make a country df with values and midpoints
  group_data <- as.data.frame(cbind(group_data_vals, group_data_x))
  lines(x=group_data$group_data_x, y=group_data_vals, col=fst_col_list[c], lwd = 2)
}
legend("topleft", 
       legend = c("Fst orientis-chrysippus", "Fst orientis-klugii", "Fst chrysippus-klugii",
                  "Fxy orientis-karamu", "Fxy chrysippus-karamu", "Fxy klugii-karamu"), 
       col = c(fst_col_list),
       lty = 1,
       lwd = 2,
       pt.cex = 1,
       cex = 0.5,
       bty = "n")

#now do Da = dxy-(mean(pi1,pi2))
Da_df <- cbind(dxy_chr15, pi_chr15)
Da_df$da_OC <- (Da_df$dxy_OC-((pi_chr15$pi_orientis+pi_chr15$pi_chrysippus)/2))
Da_df$da_OK <- (Da_df$dxy_OK-((pi_chr15$pi_orientis+pi_chr15$pi_klugii)/2))
Da_df$da_CK <- (Da_df$dxy_CK-((pi_chr15$pi_klugii+pi_chr15$pi_chrysippus)/2))
Da_df$da_OKa <- (Da_df$dxy_OKa-((pi_chr15$pi_orientis+pi_chr15$pi_klugii)/2))
Da_df$da_CKa <- (Da_df$dxy_CKa-((pi_chr15$pi_karamu+pi_chr15$pi_chrysippus)/2))
Da_df$da_KKa <- (Da_df$dxy_KKa-((pi_chr15$pi_klugii+pi_chr15$pi_karamu)/2))

da_chr15 <- as.data.frame(cbind(Da_df$mid, Da_df$da_OC, Da_df$da_OK, 
                                Da_df$da_CK, Da_df$da_OKa, Da_df$da_CKa,
                                Da_df$da_KKa))
colnames(da_chr15) <- c("mid", "da_OC", "da_OK", "da_CK", "da_OKa", "da_CKa", "da_KKa")

plot(da_chr15$mid, da_chr15$da_OC, type = 'n', ylab = "absolute divergence (Da)", xlab = '', xaxt = 'n', ylim = c(0,0.07))
rect(xleft = 4035981, ybottom = 0, 
     xright = 5322257, ytop = 1.1, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 125), border = T, lwd = 0.7)
rect(xleft = 5323805, ybottom = 0, 
     xright = 6220875, ytop = 1.1, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 125), border = T, lwd = 0.7)
rect(xleft = 6251636, ybottom = 0, 
     xright = 7829094, ytop = 1.1, 
     col = rgb(168,117,238, maxColorValue = 255, alpha = 125), border = T, lwd = 0.7)
rect(xleft = 14803486, ybottom = 0, 
     xright = 15299732, ytop = 1.1, 
     col = rgb(247,97,78, maxColorValue = 255, alpha = 125), border = T, lwd = 0.7)
dalist <- c(2,3,4,5,6,7)
da_col_list <- c("purple", "green", "orange", "black", "darkgrey", "lightgrey")
for (c in 1:length(dalist)){
  #extract the country values from data 
  group_data_vals <- da_chr15[,(dalist[c])]
  group_data_x <- da_chr15$mid
  #make a country df with values and midpoints
  group_data <- as.data.frame(cbind(group_data_vals, group_data_x))
  lines(x=group_data$group_data_x, y=group_data_vals, col=da_col_list[c], lwd = 2)
}
legend("topleft", 
       legend = c("Da orientis-chrysippus", "Da orientis-klugii", "Da chrysippus-klugii",
                  "Da orientis-karamu", "Da chrysippus-karamu", "Da klugii-karamu"), 
       col = c(da_col_list),
       lty = 1,
       lwd = 2,
       pt.cex = 1,
       cex = 0.5,
       bty = "n")


#########chr15 part1 plot FIGURE

chr15_df <- read.csv("contig15.1.div.regions.output.csv.gz")

pi_chr15 <- as.data.frame(cbind(chr15_df$mid, chr15_df$pi_NGA, chr15_df$pi_TSW, chr15_df$pi_WAT, chr15_df$pi_karamu))
colnames(pi_chr15) <- c("mid", "pi_chrysippus", "pi_orientis", "pi_klugii", "pi_karamu")

dxy_chr15 <- as.data.frame(cbind(chr15_df$mid, chr15_df$dxy_TSW_NGA, chr15_df$dxy_WAT_TSW, 
                                 chr15_df$dxy_WAT_NGA, chr15_df$dxy_TSW_karamu, chr15_df$dxy_NGA_karamu,
                                 chr15_df$dxy_WAT_karamu))
colnames(dxy_chr15) <- c("mid", "dxy_OC", "dxy_OK", "dxy_CK", "dxy_OKa", "dxy_CKa", "dxy_KKa")

fst_chr15 <- as.data.frame(cbind(chr15_df$mid, chr15_df$Fst_TSW_NGA, chr15_df$Fst_WAT_TSW, 
                                 chr15_df$Fst_WAT_NGA,chr15_df$Fst_TSW_karamu, chr15_df$Fst_NGA_karamu,
                                 chr15_df$Fst_WAT_karamu))
colnames(fst_chr15) <- c("mid", "fst_OC", "fst_OK", "fst_CK", "fst_OKa", "fst_CKa", "fst_KKa")

par(mfrow=c(4,1), mai = c(0.3, 1, 0.1, 0.1))

# #first plot chr15 segments
# plot(pi_chr15$mid, pi_chr15$pi_chrysippus, type = 'n', ylab = "", xlab = "", xaxt = 'n', ylim = c(0,0.07), frame.plot = F)
# 
# seg_df <- read.csv("Dchry2.2_region_coordinates.tsv", sep = "\t")
# seg_df <- subset(seg_df, seg_df$Region != 3)
# seg_df$colour <- "empty"
# for(i in 1:nrow(seg_df)){
#   if(seg_df$Region[i] == "1"){
#     seg_df$colour[i] <- "darkgreen"
#   }
#   if(seg_df$Region[i] == 2){
#     seg_df$colour[i] <- "orchid"
#   }
#   if(seg_df$Region[i] == 3){
#     seg_df$colour[i] <- "limegreen"
#   }
#   if(seg_df$Region[i] == 4){
#     seg_df$colour[i] <- "tomato"
#   }
# }
# 
# for(i in 1:nrow(seg_df)){
#   sign <- seg_df$strand[i]
#   if(sign == "+"){
#     a <- seg_df$Qstart[i]
#     b <- seg_df$Qstart[i]
#     c <- (seg_df$Qend[i]-300000)
#     d <- seg_df$Qend[i]
#     e <- (seg_df$Qend[i]-300000)
#     f <- seg_df$Qstart[i]
#     polygon(x = c(a,b,c,d,e,f), y = c(0.07,0.06,0.06,0.065,0.07,0.07), col = seg_df$colour[i], border = F)
#   }
#   if(sign == "-"){
#     a <- seg_df$Qend[i]
#     b <- seg_df$Qend[i]
#     c <- (seg_df$Qstart[i]+300000)
#     d <- seg_df$Qstart[i]
#     e <- (seg_df$Qstart[i]+300000)
#     f <- seg_df$Qend[i]
#     polygon(x = c(a,b,c,d,e,f), y = c(0.07,0.06,0.06,0.065,0.07,0.07), col = seg_df$colour[i], border = F)
#   }
# }


plot(pi_chr15$mid, pi_chr15$pi_chrysippus, type = 'n', ylab = "nucleotide diversity (π)", xlab = "", xaxt = 'n', ylim = c(0,0.07), frame.plot = F)
#axis(side = 1, at = c(0, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000, 11000000, 12000000, 13000000, 14000000, 15000000, 16000000, 17000000, 18000000), labels = c("0", "1", "2", "3", "4", "5", "6","7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))
rect(xleft = 4035981, ybottom = 0, 
     xright = 5322257, ytop = 0.08, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 5323805, ybottom = 0, 
     xright = 6220875, ytop = 0.08, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 6251636, ybottom = 0, 
     xright = 7829094, ytop = 0.08, 
     col = rgb(168,117,238, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 14803486, ybottom = 0, 
     xright = 15299732, ytop = 0.08, 
     col = rgb(247,97,78, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
pilist <- c(2,3,4,5)
pi_col_list <- c("#bc4754", "#2143d1", "#ffac07", "hotpink")
for (c in 1:length(pilist)){
  #extract the country values from data 
  group_data_vals <- pi_chr15[,(pilist[c])]
  group_data_x <- pi_chr15$mid
  #make a country df with values and midpoints
  group_data <- as.data.frame(cbind(group_data_vals, group_data_x))
  lines(x=group_data$group_data_x, y=group_data_vals, col=pi_col_list[c], lwd = 2)
}
legend("topleft", 
       legend = c("π chrysippus", "π orientis", "π klugii", "π karamu"), 
       col = c(pi_col_list),
       lty = 1,
       lwd = 2,
       pt.cex = 1,
       cex = 0.5,
       bty = "n")

plot(dxy_chr15$mid, dxy_chr15$dxy_OC, type = 'n', ylab = "absolute pairwise divergence (Dxy)", xlab = '', xaxt = 'n', ylim = c(0,0.07), frame.plot = F)
#axis(side = 1, at = c(0, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000, 11000000, 12000000, 13000000, 14000000, 15000000, 16000000, 17000000, 18000000), labels = c("0", "1", "2", "3", "4", "5", "6","7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))

rect(xleft = 4035981, ybottom = 0, 
     xright = 5322257, ytop = 0.08, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 5323805, ybottom = 0, 
     xright = 6220875, ytop = 0.08, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 6251636, ybottom = 0, 
     xright = 7829094, ytop = 0.08, 
     col = rgb(168,117,238, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 14803486, ybottom = 0, 
     xright = 15299732, ytop = 0.08, 
     col = rgb(247,97,78, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
dxylist <- c(2,3,4,5,6,7)
dxy_col_list <- c("orchid3", "seagreen4", "skyblue4", "coral3", "burlywood", "black")
for (c in 1:length(dxylist)){
  #extract the country values from data 
  group_data_vals <- dxy_chr15[,(dxylist[c])]
  group_data_x <- dxy_chr15$mid
  #make a country df with values and midpoints
  group_data <- as.data.frame(cbind(group_data_vals, group_data_x))
  lines(x=group_data$group_data_x, y=group_data_vals, col=dxy_col_list[c], lwd = 2)
}

legend("topleft", 
       legend = c("Dxy orientis-chrysippus", "Dxy orientis-klugii", "Dxy chrysippus-klugii",
                  "Dxy orientis-karamu", "Dxy chrysippus-karamu", "Dxy klugii-karamu"), 
       col = c(dxy_col_list),
       lty = 1,
       lwd = 2,
       pt.cex = 1,
       cex = 0.5,
       bty = "n")

#now do Da = dxy-(mean(pi1,pi2))
Da_df <- cbind(dxy_chr15, pi_chr15)
Da_df$da_OC <- (Da_df$dxy_OC-((pi_chr15$pi_orientis+pi_chr15$pi_chrysippus)/2))
Da_df$da_OK <- (Da_df$dxy_OK-((pi_chr15$pi_orientis+pi_chr15$pi_klugii)/2))
Da_df$da_CK <- (Da_df$dxy_CK-((pi_chr15$pi_klugii+pi_chr15$pi_chrysippus)/2))
Da_df$da_OKa <- (Da_df$dxy_OKa-((pi_chr15$pi_orientis+pi_chr15$pi_klugii)/2))
Da_df$da_CKa <- (Da_df$dxy_CKa-((pi_chr15$pi_karamu+pi_chr15$pi_chrysippus)/2))
Da_df$da_KKa <- (Da_df$dxy_KKa-((pi_chr15$pi_klugii+pi_chr15$pi_karamu)/2))

da_chr15 <- as.data.frame(cbind(Da_df$mid, Da_df$da_OC, Da_df$da_OK, 
                                Da_df$da_CK, Da_df$da_OKa, Da_df$da_CKa,
                                Da_df$da_KKa))
colnames(da_chr15) <- c("mid", "da_OC", "da_OK", "da_CK", "da_OKa", "da_CKa", "da_KKa")

plot(da_chr15$mid, da_chr15$da_OC, type = 'n', ylab = "absolute divergence (Da)", xlab = '', xaxt = 'n', ylim = c(0,0.07), frame.plot = F)
#axis(side = 1, at = c(0, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000, 11000000, 12000000, 13000000, 14000000, 15000000, 16000000, 17000000, 18000000), labels = c("0", "1", "2", "3", "4", "5", "6","7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))

rect(xleft = 4035981, ybottom = 0, 
     xright = 5322257, ytop = 1.1, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 5323805, ybottom = 0, 
     xright = 6220875, ytop = 1.1, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 6251636, ybottom = 0, 
     xright = 7829094, ytop = 1.1, 
     col = rgb(168,117,238, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 14803486, ybottom = 0, 
     xright = 15299732, ytop = 1.1, 
     col = rgb(247,97,78, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
dalist <- c(2,3,4,5,6,7)
da_col_list <- c("orchid3", "seagreen4", "skyblue4", "coral3", "burlywood", "black")
for (c in 1:length(dalist)){
  #extract the country values from data 
  group_data_vals <- da_chr15[,(dalist[c])]
  group_data_x <- da_chr15$mid
  #make a country df with values and midpoints
  group_data <- as.data.frame(cbind(group_data_vals, group_data_x))
  lines(x=group_data$group_data_x, y=group_data_vals, col=da_col_list[c], lwd = 2)
}
legend("topleft", 
       legend = c("Da orientis-chrysippus", "Da orientis-klugii", "Da chrysippus-klugii",
                  "Da orientis-karamu", "Da chrysippus-karamu", "Da klugii-karamu"), 
       col = c(da_col_list),
       lty = 1,
       lwd = 2,
       pt.cex = 1,
       cex = 0.5,
       bty = "n")

plot(fst_chr15$mid, fst_chr15$fst_OC, type = 'n', ylab = "pairwise divergence (Fst)", xlab = '', xaxt = 'n', ylim = c(0,1), frame.plot = F)
axis(side = 1, at = c(0, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000, 11000000, 12000000, 13000000, 14000000, 15000000, 16000000, 17000000, 18000000), labels = c("0", "1", "2", "3", "4", "5", "6","7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))

rect(xleft = 4035981, ybottom = 0, 
     xright = 5322257, ytop = 1.1, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 5323805, ybottom = 0, 
     xright = 6220875, ytop = 1.1, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 6251636, ybottom = 0, 
     xright = 7829094, ytop = 1.1, 
     col = rgb(168,117,238, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 14803486, ybottom = 0, 
     xright = 15299732, ytop = 1.1, 
     col = rgb(247,97,78, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
fstlist <- c(2,3,4,5,6,7)
fst_col_list <- c("orchid3", "seagreen4", "skyblue4", "coral3", "burlywood", "black")

for (c in 1:length(dxylist)){
  #extract the country values from data 
  group_data_vals <- fst_chr15[,(fstlist[c])]
  group_data_x <- fst_chr15$mid
  #make a country df with values and midpoints
  group_data <- as.data.frame(cbind(group_data_vals, group_data_x))
  lines(x=group_data$group_data_x, y=group_data_vals, col=fst_col_list[c], lwd = 2)
}
legend("topleft", 
       legend = c("Fst orientis-chrysippus", "Fst orientis-klugii", "Fst chrysippus-klugii",
                  "Fxy orientis-karamu", "Fxy chrysippus-karamu", "Fxy klugii-karamu"), 
       col = c(fst_col_list),
       lty = 1,
       lwd = 2,
       pt.cex = 1,
       cex = 0.5,
       bty = "n")

#############################
#now without legends
par(mfrow=c(4,1), mai = c(0.3, 1, 0.1, 0.1))

plot(pi_chr15$mid, pi_chr15$pi_chrysippus, type = 'n', ylab = "nucleotide diversity (π)", xlab = "", xaxt = 'n', ylim = c(0,0.07), frame.plot = F)
#axis(side = 1, at = c(0, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000, 11000000, 12000000, 13000000, 14000000, 15000000, 16000000, 17000000, 18000000), labels = c("0", "1", "2", "3", "4", "5", "6","7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))
rect(xleft = 4035981, ybottom = 0, 
     xright = 5322257, ytop = 0.08, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 5323805, ybottom = 0, 
     xright = 6220875, ytop = 0.08, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 6251636, ybottom = 0, 
     xright = 7829094, ytop = 0.08, 
     col = rgb(168,117,238, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 14803486, ybottom = 0, 
     xright = 15299732, ytop = 0.08, 
     col = rgb(247,97,78, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
pilist <- c(2,3,4,5)
pi_col_list <- c("#bc4754", "#2143d1", "#ffac07", "hotpink")
for (c in 1:length(pilist)){
  #extract the country values from data 
  group_data_vals <- pi_chr15[,(pilist[c])]
  group_data_x <- pi_chr15$mid
  #make a country df with values and midpoints
  group_data <- as.data.frame(cbind(group_data_vals, group_data_x))
  lines(x=group_data$group_data_x, y=group_data_vals, col=pi_col_list[c], lwd = 2)
}

plot(dxy_chr15$mid, dxy_chr15$dxy_OC, type = 'n', ylab = "absolute pairwise divergence (Dxy)", xlab = '', xaxt = 'n', ylim = c(0,0.07), frame.plot = F)
#axis(side = 1, at = c(0, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000, 11000000, 12000000, 13000000, 14000000, 15000000, 16000000, 17000000, 18000000), labels = c("0", "1", "2", "3", "4", "5", "6","7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))

rect(xleft = 4035981, ybottom = 0, 
     xright = 5322257, ytop = 0.08, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 5323805, ybottom = 0, 
     xright = 6220875, ytop = 0.08, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 6251636, ybottom = 0, 
     xright = 7829094, ytop = 0.08, 
     col = rgb(168,117,238, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 14803486, ybottom = 0, 
     xright = 15299732, ytop = 0.08, 
     col = rgb(247,97,78, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
dxylist <- c(2,3,4,5,6,7)
dxy_col_list <- c("orchid3", "seagreen4", "skyblue4", "coral3", "burlywood", "black")
for (c in 1:length(dxylist)){
  #extract the country values from data 
  group_data_vals <- dxy_chr15[,(dxylist[c])]
  group_data_x <- dxy_chr15$mid
  #make a country df with values and midpoints
  group_data <- as.data.frame(cbind(group_data_vals, group_data_x))
  lines(x=group_data$group_data_x, y=group_data_vals, col=dxy_col_list[c], lwd = 2)
}

#now do Da = dxy-(mean(pi1,pi2))
Da_df <- cbind(dxy_chr15, pi_chr15)
Da_df$da_OC <- (Da_df$dxy_OC-((pi_chr15$pi_orientis+pi_chr15$pi_chrysippus)/2))
Da_df$da_OK <- (Da_df$dxy_OK-((pi_chr15$pi_orientis+pi_chr15$pi_klugii)/2))
Da_df$da_CK <- (Da_df$dxy_CK-((pi_chr15$pi_klugii+pi_chr15$pi_chrysippus)/2))
Da_df$da_OKa <- (Da_df$dxy_OKa-((pi_chr15$pi_orientis+pi_chr15$pi_klugii)/2))
Da_df$da_CKa <- (Da_df$dxy_CKa-((pi_chr15$pi_karamu+pi_chr15$pi_chrysippus)/2))
Da_df$da_KKa <- (Da_df$dxy_KKa-((pi_chr15$pi_klugii+pi_chr15$pi_karamu)/2))

da_chr15 <- as.data.frame(cbind(Da_df$mid, Da_df$da_OC, Da_df$da_OK, 
                                Da_df$da_CK, Da_df$da_OKa, Da_df$da_CKa,
                                Da_df$da_KKa))
colnames(da_chr15) <- c("mid", "da_OC", "da_OK", "da_CK", "da_OKa", "da_CKa", "da_KKa")

plot(da_chr15$mid, da_chr15$da_OC, type = 'n', ylab = "absolute divergence (Da)", xlab = '', xaxt = 'n', ylim = c(0,0.07), frame.plot = F)
#axis(side = 1, at = c(0, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000, 11000000, 12000000, 13000000, 14000000, 15000000, 16000000, 17000000, 18000000), labels = c("0", "1", "2", "3", "4", "5", "6","7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))

rect(xleft = 4035981, ybottom = 0, 
     xright = 5322257, ytop = 1.1, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 5323805, ybottom = 0, 
     xright = 6220875, ytop = 1.1, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 6251636, ybottom = 0, 
     xright = 7829094, ytop = 1.1, 
     col = rgb(168,117,238, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 14803486, ybottom = 0, 
     xright = 15299732, ytop = 1.1, 
     col = rgb(247,97,78, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
dalist <- c(2,3,4,5,6,7)
da_col_list <- c("orchid3", "seagreen4", "skyblue4", "coral3", "burlywood", "black")
for (c in 1:length(dalist)){
  #extract the country values from data 
  group_data_vals <- da_chr15[,(dalist[c])]
  group_data_x <- da_chr15$mid
  #make a country df with values and midpoints
  group_data <- as.data.frame(cbind(group_data_vals, group_data_x))
  lines(x=group_data$group_data_x, y=group_data_vals, col=da_col_list[c], lwd = 2)
}

plot(fst_chr15$mid, fst_chr15$fst_OC, type = 'n', ylab = "pairwise divergence (Fst)", xlab = '', xaxt = 'n', ylim = c(0,1), frame.plot = F)
axis(side = 1, at = c(0, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000, 11000000, 12000000, 13000000, 14000000, 15000000, 16000000, 17000000, 18000000), labels = c("0", "1", "2", "3", "4", "5", "6","7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))

rect(xleft = 4035981, ybottom = 0, 
     xright = 5322257, ytop = 1.1, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 5323805, ybottom = 0, 
     xright = 6220875, ytop = 1.1, 
     col = rgb(1,99,55, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 6251636, ybottom = 0, 
     xright = 7829094, ytop = 1.1, 
     col = rgb(168,117,238, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
rect(xleft = 14803486, ybottom = 0, 
     xright = 15299732, ytop = 1.1, 
     col = rgb(247,97,78, maxColorValue = 255, alpha = 40), border = F, lwd = 0.7)
fstlist <- c(2,3,4,5,6,7)
fst_col_list <- c("orchid3", "seagreen4", "skyblue4", "coral3", "burlywood", "black")

for (c in 1:length(dxylist)){
  #extract the country values from data 
  group_data_vals <- fst_chr15[,(fstlist[c])]
  group_data_x <- fst_chr15$mid
  #make a country df with values and midpoints
  group_data <- as.data.frame(cbind(group_data_vals, group_data_x))
  lines(x=group_data$group_data_x, y=group_data_vals, col=fst_col_list[c], lwd = 2)
}

