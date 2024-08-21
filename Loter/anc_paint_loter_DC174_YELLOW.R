

plot_haps <- function(starts, ends, cols1, cols2, add=FALSE, xlim=NULL, y=1, width=0.8, label=NA, label_cex=1, box=TRUE, sep_col=NA, border=NA){
    if (is.null(xlim)) xlim=c(starts[1], ends[length(ends)])
    if (add==FALSE) {
        plot(0, cex = 0, xlim = xlim, ylim = c(0,2), ylab = "", xaxt="n", yaxt = "n", xlab = "")
        if (is.na(label) == FALSE) mtext(3,at=xlim[1],text=label, cex=label_cex, line=0, adj=0)
        }

    idx <- which(starts >= xlim[1] & ends <= xlim[2])
    
    rect(starts[idx], y, ends[idx], y-width, col = cols1[idx], border=border)
    rect(starts[idx], y, ends[idx], y+width, col = cols2[idx], border=border)
    if (box==TRUE) rect(xlim[1], y-width, xlim[2], y+width)
    segments(xlim[1], y, xlim[2], y, col=sep_col)
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



samples <- read.table("../DC174_sample_data/DC174_region.txt", as.is = T)[,1]

hap_names <- paste(rep(samples, each=2), c("A","B"), sep="_")

hap_names_by_sample <- sapply(samples, function(sample) paste(sample, c("A","B"), sep = "_"), simplify=F)


dir="/home/simon/Research/Danaus_spiroplasma/HAPLOTYPE_ANCESTRY/Loter/"

prefix="contig15.1.double.phased.minVar2HET60BI.loter_K3_140noext.phasepaint_I20"


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

anc_array <- array(dim = c(nrow(anc), length(samples), 2),
                          dimnames = list(1:nrow(anc), samples, c(1,2)))

for (sample in samples){
    anc_array[,sample,] <- as.matrix(anc[,hap_names_by_sample[[sample]]])
    }


### genes and gwas hits

########################### genes ############################

anno <- read.delim("Dchry2.2.fa.masked_annotation_a001.sequences.tidy.gff3.gz", sep="\t", as.is=T, comment.char="#", header=FALSE)

genes_f <- subset(anno, anno[,1] == "contig15.1" & anno[,3] == "gene" & anno[,7] == "+")
genes_r <- subset(anno, anno[,1] == "contig15.1" & anno[,3] == "gene" & anno[,7] == "-")

CDS_f <- subset(anno, anno[,1] == "contig15.1" & anno[,3] == "CDS" & anno[,7] == "+")
CDS_r <- subset(anno, anno[,1] == "contig15.1" & anno[,3] == "CDS" & anno[,7] == "-")

######################## gwas outliers ##########################

gwas <- read.csv("/home/simon/Research/Danaus_spiroplasma/gwas/all_background_points_chr15.csv", as.is=T)


top10 <- gwas$BP[order(gwas$P)[1:10]]


######################## klugii karamu orientis ##########################

plot_samples <- c("SM15W72", "SM15W74", "SM21SP06", "SM21SP09", "SM21TS05", "SM21TS33")

png(paste0(prefix, ".twisstPaint_karamu_yellow.png"), width = 4000, height = 800, res = 300, bg=NA)


par(mfcol = c(1,1), mar = c(2,2,1,1), bty = "n", xpd=FALSE)

bgcols=c(rep("light", 2), rep("dark", 4))

# png(paste0("images/",prefix, ".chrom.outer.twisstPaint.png"), width = 6000, height = 4000, res = 200, bg="white")

n=length(plot_samples)


plot(0, cex=0, bty="n", xlab="", ylab="", xlim=c(6250000,6350000), ylim=c(n+2,1), xaxt="n", yaxt="n")

# rect(6220875,0.5,6251636,n+.5)
# rect(6294547,0.5,6299719,n+.5)

for (i in 1:n){
    plot_haps(anc$start, anc$end, xlim=c(6250000,6350000),
              cols1=cols[anc_array[,plot_samples[i],1]], cols2=cols[anc_array[,plot_samples[i],2]],
              add=TRUE, y=i, width=0.3, border=NA, box=FALSE, sep_col="white")
    
    mtext(side=2, text=plot_samples[i], las=2, line=-4, adj=1, at=i)
#     points(6255000, i, pch=21, col="black", bg=ifelse(bgcols[i]=="dark", "black", NA), cex=3)
    }

axis(1, line=0, at = seq(6250000, 6350000, 10000), labels = seq(6250000, 6350000, 10000)/1e6)


# dev.off()

segments(genes_f[,4],n+1,genes_f[,5],n+1)
segments(genes_f[,4],n+1,genes_f[,4],n+1.3)
arrows(genes_f[,4],n+1.3,genes_f[,4]+1e3,n+1.3, length=0.1)
rect(CDS_f[,4],n+.8,CDS_f[,5],n+1.2, border=NA, col="black")

segments(genes_r[,4],n+1,genes_r[,5],n+1)
segments(genes_r[,5],n+1,genes_r[,5],n+1.3)
arrows(genes_r[,5],n+1.3,genes_r[,5]-1e3,n+1.3, length=0.1)
rect(CDS_r[,4],n+.8,CDS_r[,5],n+1.2, border=NA, col="black")

text(6297133, n+1.8, labels="yellow")

# legend("topleft", pch=21, pt.bg=c(NA,"black"), legend=c("light", "dark"), title="Wing background colouration", ncol=2, border=NA)
# legend("top", fill=cols[c("klugii", "orientis")],
#                              legend=c("klugii (Kenya)", "orientis (South Africa)"),
#                              title="Ancestry assignment", ncol=3, border=NA)

#gwas
points(top10, rep(n+.5, 10), pch=2, cex=1, col="forestgreen")

dev.off()


######################## chrysippus chrysippusDark orientis ##########################

plot_samples <- c("SM16N01", "SM16N05", "SM19SY01", "RV12N317", "SM21TS05", "SM21TS33")

png(paste0(prefix, ".twisstPaint_chryDark_yellow.png"), width = 4000, height = 800, res = 300, bg=NA)

par(mfcol = c(1,1), mar = c(2,2,1,1), bty = "n", xpd=FALSE)

bgcols=c(rep("light", 2), rep("dark", 4))

# png(paste0("images/",prefix, ".chrom.outer.twisstPaint.png"), width = 6000, height = 4000, res = 200, bg="white")

n=length(plot_samples)


plot(0, cex=0, bty="n", xlab="", ylab="", xlim=c(6250000,6350000), ylim=c(n+2,1), xaxt="n", yaxt="n")

# rect(6220875,0.5,6251636,n+.5)
# rect(6294547,0.5,6299719,n+.5)

for (i in 1:n){
    plot_haps(anc$start, anc$end, xlim=c(6250000,6350000),
              cols1=cols[anc_array[,plot_samples[i],1]], cols2=cols[anc_array[,plot_samples[i],2]],
              add=TRUE, y=i, width=0.3, border=NA, box=FALSE, sep_col="white")
    
    mtext(side=2, text=plot_samples[i], las=2, line=-4, adj=1, at=i)
#     points(6255000, i, pch=21, col="black", bg=ifelse(bgcols[i]=="dark", "black", NA), cex=3)
    }

axis(1, line=0, at = seq(6250000, 6350000, 10000), labels = seq(6250000, 6350000, 10000)/1e6)


# dev.off()

segments(genes_f[,4],n+1,genes_f[,5],n+1)
segments(genes_f[,4],n+1,genes_f[,4],n+1.3)
arrows(genes_f[,4],n+1.3,genes_f[,4]+1e3,n+1.3, length=0.1)
rect(CDS_f[,4],n+.8,CDS_f[,5],n+1.2, border=NA, col="black")

segments(genes_r[,4],n+1,genes_r[,5],n+1)
segments(genes_r[,5],n+1,genes_r[,5],n+1.3)
arrows(genes_r[,5],n+1.3,genes_r[,5]-1e3,n+1.3, length=0.1)
rect(CDS_r[,4],n+.8,CDS_r[,5],n+1.2, border=NA, col="black")

text(6297133, n+1.8, labels="yellow")

# legend("topleft", pch=21, pt.bg=c(NA,"black"), legend=c("light", "dark"), title="Wing background colouration", ncol=2, border=NA)
# legend("top", fill=cols[c("klugii", "orientis")],
#                              legend=c("klugii (Kenya)", "orientis (South Africa)"),
#                              title="Ancestry assignment", ncol=3, border=NA)

#gwas
points(top10, rep(n+.5, 10), pch=2, cex=1, col="forestgreen")

dev.off()


