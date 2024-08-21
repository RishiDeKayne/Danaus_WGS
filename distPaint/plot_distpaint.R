
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


dir="/home/simon/Research/Danaus_spiroplasma/HAPLOTYPE_ANCESTRY/distances/"
setwd(dir)

prefix <- "contig15.1.double.phased.distPaint.w200p01.phasepaint_I100"
# prefix <- "contig15.1.double.phased.distPaint.w200p01"

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

cols <- c(klugii="#ffac07",orientis="#2143d1",chrysippus="#bc4754",karamu="#ff8cf1", "gray90")

#set ancestry for missing to number 4
anc[,-(1:5)][anc[,-(1:5)]==-1] <- 4

#add 1 to each number to make it 1-based
anc[,-(1:5)] <- anc[,-(1:5)] + 1

###############################################################

best_match_array <- array(dim = c(nrow(anc), length(samples), 2),
                          dimnames = list(1:nrow(anc), samples, c(1,2)))

for (sample in samples){
    best_match_array[,sample,] <- as.matrix(anc[,hap_names_by_sample[[sample]]])
    }


##################### plot #########################

################ individual plots 

### complete
left_trim = 0
right_trim = 18e6

for (sample in samples){
    png(paste0("individual_plots/",prefix,sample,".png"), width = 2000, height = 200, res = 300, bg=NA)
    
    par(mar = c(0,0,0.5,0), bty = "n", xpd=NA)
    
    plot_haps(anc$start, anc$end, xlim=c(left_trim,right_trim), width=0.8,
              cols1=cols[best_match_array[,sample,1]], cols2=cols[best_match_array[,sample,2]],
              label="")

    #inversions
    arrows(4035981, 2.1, 6220875, 2.1, length=0.1, code=3, lwd=1, )
    arrows(5322257, 2.1, 6220875, 2.1, length=0.1, code=3, lwd=1)
    arrows(6251636, 2.1, 7829094, 2.1, length=0.1, code=3, lwd=1)
    arrows(14803486, 2.1, 15299732, 2.1, length=0.1, code=3, lwd=1)

    
    dev.off()
    }

### compact

left_trim = 2e6
right_trim = 17e6-gap


for (sample in samples){
    png(paste0("individual_plots/",prefix,sample,".compact.png"), width = 2000, height = 200, res = 300, bg=NA)
    
    par(mar = c(0,0,0.5,0), bty = "n", xpd=NA)
    
    plot_haps(anc$start_adjusted, anc$end_adjusted, xlim=c(left_trim,right_trim), width=0.8,
              cols1=cols[best_match_array[,sample,1]], cols2=cols[best_match_array[,sample,2]],
              label="")

    #inversions
    arrows(4035981, 2.1, 6220875, 2.1, length=0.1, code=3, lwd=1)
    arrows(5322257, 2.1, 6220875, 2.1, length=0.1, code=3, lwd=1)
    arrows(6251636, 2.1, 7829094, 2.1, length=0.1, code=3, lwd=1)
    arrows(14803486 - gap, 2.1, 15299732 - gap, 2.1, length=0.1, code=3, lwd=1)

    
    dev.off()
    }

#scale_bar
png(paste0("individual_plots/","compact.scale.png"), width = 2000, height = 200, res = 300, bg=NA)
par(mar = c(0,0,0.5,0), bty = "n", xpd=NA)
plot(0, cex = 0, xlim = xlim, ylim = c(0,2), ylab = "", xaxt="n", yaxt = "n", xlab = "")
arrows(9e6, 1, 10e6, 1, length=0.1, code=3, angle=90, lwd=1.5)
text(9.5e6, 0.7, labels="1Mb")
dev.off()

###array plot compact

left_trim = 2e6
right_trim = 17e6-gap


plot_samples <- samples[1:74]
png(paste0(prefix, "outer.compact.png"), width = 6000, height = 4000, res = 200, bg="white")
par(mfcol = c(25,4), mar = c(1,2,1.2,1), bty = "n", xpd=NA)

plot_samples <- samples[75:174]
png(paste0(prefix, ".inner.compact.png"), width = 6000, height = 4000, res = 200, bg="white")
par(mfcol = c(25,4), mar = c(1,2,1.2,1), bty = "n", xpd=NA)

plot_samples <- c("SM19M043", "IG20NY047", "SM21SP12")
png(paste0(prefix, ".all.compact.tall.png"), width = 1500, height = 1000, res = 300, bg="white")
par(mfcol = c(3,1), mar = c(1,.5,1.2,.5), bty = "n", xpd=NA)


plot_samples <- samples
png(paste0(prefix, ".all.compact.tall.png"), width = 2000, height = 2600, res = 100, bg="white")
par(mfcol = c(35,5), mar = c(.2,.2,1.2,.2), bty = "n", xpd=NA)




xlim = c(left_trim,right_trim)
for (sample in plot_samples){
    plot_haps(anc$start_adjusted, anc$end_adjusted, xlim=c(left_trim,right_trim), width=0.8,
              cols1=cols[best_match_array[,sample,1]], cols2=cols[best_match_array[,sample,2]], label=paste(sample, populations[sample]))

    #inversions
    arrows(4035981, 2.05, 6220875, 2.05, length=0.1, code=3, lwd=2, )
    arrows(5322257, 2.05, 6220875, 2.05, length=0.1, code=3, lwd=2)
    arrows(6251636, 2.05, 7829094, 2.05, length=0.1, code=3, lwd=2)
    arrows(14803486-gap, 2.05, 15299732-gap, 2.05, length=0.1, code=3, lwd=2)
    
    }

dev.off()


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


############# inferred genotypes

inferred_haps <- data.frame(ID=samples, hap1=rep(NA, length(samples)), hap2=rep(NA, length(samples)))
rownames(inferred_haps) = samples

for (sample in samples){
    tab1 <- table(best_match_array[1500:2000,sample,1])
    inferred_haps[sample,"hap1"] <- ifelse(max(tab1 / sum(tab1)) > .5, names(which.max(tab1 / sum(tab1))), 5)
    tab2 <- table(best_match_array[1500:2000,sample,2])
    inferred_haps[sample,"hap2"] <- ifelse(max(tab2 / sum(tab2)) > .5, names(which.max(tab2 / sum(tab2))), 5)
    }

inferred_haps$hap1 <- c("klugii", "orientis", "chrysippus", "karamu", "unknown")[as.numeric(inferred_haps$hap1)]
inferred_haps$hap2 <- c("klugii", "orientis", "chrysippus", "karamu", "unknown")[as.numeric(inferred_haps$hap2)]

inferred_haps[,c("hap1","hap2")] <- t(apply(inferred_haps[,c("hap1","hap2")], 1, sort))

#get frequencies

inferred_haps$geno = apply(inferred_haps[,c("hap1","hap2")], 1, paste, collapse="/")


inferred_haps[samples_by_region[["Rwanda (Hybrid Zone)"]], "geno"]

geno_table = table(inferred_haps[samples_by_region[["Rwanda (Hybrid Zone)"]], "geno"])

geno_freqs = geno_table/sum(geno_table)

#hap frequencies and hardy weinberg expectations
hap_table <- table(unlist(inferred_haps[samples_by_region[["Rwanda (Hybrid Zone)"]], c("hap1","hap2")]))

hap_freqs <- hap_table / sum(hap_table)

expected <- hap_freqs*hap_freqs
names(expected) <- paste(names(hap_freqs), names(hap_freqs), sep="/")

for (i in 1:4){
    for (j in (i+1):5){
        expected <- c(expected, 2*hap_freqs[i]*hap_freqs[j])
        names(expected)[length(expected)] <- paste(names(hap_freqs)[i], names(hap_freqs)[j], sep="/")
        }
    }


expected <- expected[order(expected, decreasing=T)]

#add missing for observed    
for (geno in names(expected)){
  if (!(geno %in% names(geno_freqs))){
      geno_freqs[geno] <- 0
      }
    }

geno_freqs <- geno_freqs[names(expected)]

png("Rwanda_genomes_BC_genotype_frequencies.png", width=2000, height=2000, res=200)
par(mar=c(6,10,1,1))
barplot(rbind(geno_freqs*51, expected*51), beside=T, col = c("black", "gray"), horiz=T, las=2, xlab="Genotype Frequency", xlim=c(0,15))
legend("topright", legend=c("observed", "HWE expectation"), fill=c("black", "gray"))
dev.off()


chisq.test(rbind(geno_freqs * 51, expected*51))
