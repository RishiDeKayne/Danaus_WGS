
plot_haps <- function(starts, ends, xlim, width, cols1, cols2, label="", label_cex=1){
    plot(0, cex = 0, xlim = xlim, ylim = c(0,2), ylab = "", xaxt="n", yaxt = "n", xlab = "")
    
    mtext(3,at=xlim[1],text=label, cex=label_cex, line=0, adj=0)
    
    idx <- which(starts >= xlim[1] & ends <= xlim[2])
    
    rect(starts[idx], 1, ends[idx], 1-width, col = cols1[idx], border=NA)
    rect(starts[idx], 1, ends[idx], 1+width, col = cols2[idx], border=NA)
    rect(xlim[1], 1-width, xlim[2], 1+width)
    segments(xlim[1], 1, xlim[2], 1)
    }


fill_gaps <- function(starts_and_ends, first=NULL, last=NULL){
    filled <- starts_and_ends
    n <- nrow(starts_and_ends)
    window_gaps <-  starts_and_ends[,1][-1] - starts_and_ends[,2][-n]
    filled[,2][-n] <- starts_and_ends[,2][-n] + window_gaps/2
    filled[,1][-1] <- starts_and_ends[,1][-1] - window_gaps/2
    if (is.null(first) == FALSE) filled[,1][1] <- first
    if (is.null(last) == FALSE) filled[,2][n] <- last
    filled
    }



samples <- read.table("../DC174_sample_data/DC174_region.txt", as.is = T, sep="\t")[,1]

populations <- read.table("../DC174_sample_data/DC174_region.txt", as.is = T, sep="\t")[,2]

regions <- read.table("../DC174_sample_data/DC174_region.txt", as.is = T, sep="\t")[,3]

names(populations) <- samples

hap_names <- paste(rep(samples, each=2), c("A","B"), sep="_")

hap_names_by_sample <- sapply(samples, function(sample) paste(sample, c("A","B"), sep = "_"), simplify=F)


dir="/home/simon/Research/Danaus_spiroplasma/HAPLOTYPE_ANCESTRY/Loter/"
setwd(dir)

prefix="contig15.1.double.phased.minVar2HET60BI.loter_K4_160noext.phasepaint_I20"

anc <- read.table(paste0(dir,prefix,".tsv.gz"), header=T, as.is = T)

anc[,c("start","end")] <- fill_gaps(anc[,c("start","end")], first = 1, last = 17804635)


cut_left <- 7829094
cut_right <- 14803486
gap <- cut_right-cut_left - 2e5

anc <- subset(anc, anc$end < cut_left | anc$start > cut_right)

anc$end_adjusted <- ifelse(anc$start > cut_left, anc$end - gap, anc$end)
anc$start_adjusted <- ifelse(anc$start > cut_left, anc$start - gap, anc$start)


# cols <- c(orientis="#2143d1",klugii="#ffac07",chrysippus="#bc4754","gray90")
# anc[anc==0] <- 4

cols <- c(orientis="#2143d1",klugii="#ffac07",chrysippus="#bc4754",karamu="#ff8cf1", "gray90")
anc[anc==0] <- 5


###############################################################

best_match_array <- array(dim = c(nrow(anc), length(samples), 2),
                          dimnames = list(1:nrow(anc), samples, c(1,2)))

for (sample in samples){
    best_match_array[,sample,] <- as.matrix(anc[,hap_names_by_sample[[sample]]])
    }


############ array plot with region names

left_trim = 2e6
right_trim = 17e6-gap

xlim = c(left_trim,right_trim)

regions_unique <- unique(regions)

samples_by_region <- split(samples, regions)

png(paste0(prefix, ".all.compact.tall.regions.png"), width = 2000, height = 2400, res = 100, bg="white")
par(mfcol = c(37,5), mar = c(.2,.2,1.2,.2), bty = "n", xpd=NA)

for (region_name in regions_unique){
    plot(0, cex=0, xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
    mtext(1, text=region_name, cex=1.5, line=-2)
    
    for (sample in samples_by_region[[region_name]]){
        plot_haps(anc$start_adjusted, anc$end_adjusted, xlim=c(left_trim,right_trim), width=0.8,
                cols1=cols[best_match_array[,sample,1]], cols2=cols[best_match_array[,sample,2]], label=paste(sample, populations[sample]))

        #inversions
        arrows(4035981, 2.05, 6220875, 2.05, length=0.1, code=3, lwd=1, )
        arrows(5322257, 2.05, 6220875, 2.05, length=0.1, code=3, lwd=1)
        arrows(6251636, 2.05, 7829094, 2.05, length=0.1, code=3, lwd=1)
        arrows(14803486-gap, 2.05, 15299732-gap, 2.05, length=0.1, code=3, lwd=1)
        
        }
    }


dev.off()

