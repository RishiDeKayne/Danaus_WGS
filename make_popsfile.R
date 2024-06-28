#prepare populations for pi estimates
setwd("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/pi/")
background <- read.csv("/Users/rishide-kayne/Dropbox/RishiMac2/Danaus/DC174/background/DC174_background.csv", header = T, sep = ",")

country_list <- c("Kenya", "RSA", "StHelena", "Nigeria", "Ghana", "Spain", "Tunisia", "Italy", "Philippines", "Rwanda")
colour_list <- c("tomato", "lightseagreen", "royalblue2", "springgreen4", "yellow3", "wheat1", "goldenrod", "hotpink", "mediumorchid1", "tomato4")

bacakground_uniq_col <- unique(background$colour)

popfile_df <- as.data.frame(cbind(background$indiv, background$country))
write.table(popfile_df, file = "popsfile.txt", col.names = F, row.names = F, quote = F, sep = "\t")

popfile_df_2 <- as.data.frame(cbind(background$indiv, background$region_narrow))
write.table(popfile_df_2, file = "popsfile2.txt", col.names = F, row.names = F, quote = F, sep = "\t")


##now to make plot

#get counts of individuals for each population
indiv_counts <- c()
for(z in 1:length(country_list)){
  indiv_counts[z] <- length(which(background$country==country_list[z]))
}

#make labels to add to figure legends
country_labels <- c()
for(y in 1:length(country_list)){
  country_labels[y] <- paste(country_list[y], " n=", indiv_counts[y], sep = "")
}

chromosomes <- read.csv("../gemma/Dchry2.2.fa.fai", head = FALSE, sep = "\t")
chromosomes <- chromosomes[1:41,]
chromosomes$numb <- gsub("contig", "", chromosomes$V1)

for(j in 1:length(chromosomes$V1)){
  #get file name from contig name
  file_name <- paste(chromosomes$V1[j], ".div.output.csv.gz", sep = "")
  #get contig name and length
  contig_name <- as.character(chromosomes$V1[1])
  contig_length <- as.numeric(chromosomes$V2[1])
  #read in data
  data <- read.csv(file_name)
  #create data with no Nans so we can calc plotting dimensions
  plot_parameters_df <- na.omit(data)
  #make list of all pi values to get max value
  max_val_col <- as.vector(rbind(plot_parameters_df$pi_Kenya, plot_parameters_df$pi_RSA, plot_parameters_df$pi_StHelena, 
                                 plot_parameters_df$pi_Nigeria, plot_parameters_df$pi_Ghana, plot_parameters_df$pi_Spain, 
                                 plot_parameters_df$pi_Tunisia, plot_parameters_df$pi_Italy, plot_parameters_df$pi_Rwanda))
  
  #plot empty plot with max value as ylim
  plot(plot_parameters_df$mid, plot_parameters_df$pi_Kenya, type = 'n', ylim = c(0, max(max_val_col)),
       ylab = "pi in 100kb window", xlab = "mid bp", main = file_name)
  #add legend
  legend('topleft', 
    legend = c(country_list), 
    col = c(colour_list),
    lty = 1,
    lwd = 2,
    pt.cex = 1,
    cex = 0.5) 
  
  #now plot data - run a loop going country by country
  for (c in 1:length(country_list)){
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
    lines(smoothed, x=country_data$country_data_x, col=colour_list[c], lwd = 2)
  }
}

filename <- "pi_plot2"
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
  file_name <- paste(chromosomes$V1[j], ".div.output.csv.gz", sep = "")
  #get contig name and length
  contig_name <- as.character(chromosomes$V1[j])
  contig_length <- as.numeric(chromosomes$V2[j])
  #read in data
  data <- read.csv(file_name)
  #create data with no Nans so we can calc plotting dimensions
  plot_parameters_df <- na.omit(data)
  #make list of all pi values to get max value
  max_val_col <- as.vector(rbind(plot_parameters_df$pi_Kenya, plot_parameters_df$pi_RSA, plot_parameters_df$pi_StHelena, 
                                 plot_parameters_df$pi_Nigeria, plot_parameters_df$pi_Ghana, plot_parameters_df$pi_Spain, 
                                 plot_parameters_df$pi_Tunisia, plot_parameters_df$pi_Italy, plot_parameters_df$pi_Rwanda))

  rect_col <- altcols[altcol_numb[j]]
  #plot empty plot with ylim set by looking at individual plots
  plot(data$mid, data$pi_Kenya, type = 'n', 
       ylim = c(0, 0.06), xaxt="n",axes=F, ylab = contig_name)
  rect(xleft = min(plot_parameters_df$mid), ybottom = 0, 
       xright = max(plot_parameters_df$mid), ytop = 0.06, 
       col = rect_col, border = F)
  if (j==1){
    axis(side=2)
    legend(x = 10000, y = 0.0575, 
          legend = c(country_labels), 
          col = c(colour_list),
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
  for (c in 1:length(country_list)){
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
    loessMod <- loess(country_data$country_data_vals ~ index, data=country_data, span=0.6)
    smoothed <- predict(loessMod) 
    #lines(country_data$country_data_vals, x=country_data$country_data_x, col="black", lwd = 2)
    lines(smoothed, x=country_data$country_data_x, col=colour_list[c], lwd = 1)
  }
  
  text(x = contig_length/2, 0.0057, chromosomes$numb[j], cex = 0.5)
}

dev.off()


