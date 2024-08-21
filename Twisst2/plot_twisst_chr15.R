source("plot_twisst.R")

setwd("/home/simon/Dropbox/Research/Danaus_spiroplasma/HAPLOTYPE_ANCESTRY/Trees/")

############################## input files ######################################

prefix = "DC174OG.contig15.1.realigned.DP7GQ30het75min87minVar2.Twisst2_K3.contig15.1"
weights_file <- paste0(prefix,".topocounts.tsv.gz")
window_data_file <- paste0(prefix,".windowdata.tsv.gz")
topo_cols = c(orientis="#2143d1", chrysippus="#bc4754", klugii="#ffac07")

prefix = "DC174OG.contig15.1.realigned.DP7GQ30het75min87minVar2.Twisst2_2nd_kar_klu_ori.contig15.1"
weights_file <- paste0(prefix,".topocounts.tsv.gz")
window_data_file <- paste0(prefix,".windowdata.tsv.gz")
topo_cols = c(klugii="#ffac07", orientis="#2143d1", karamu="#CCCCCC")


prefix = "DC174OG.contig15.1.realigned.DP7GQ30het75min87minVar2.Twisst2_chryD_chry_ori.contig15.1"
weights_file <- paste0(prefix,".topocounts.tsv.gz")
window_data_file <- paste0(prefix,".windowdata.tsv.gz")
topo_cols = c(chrysippus="#bc4754", orientis="#2143d1", "#CCCCCC")


################################# import data ##################################

# The function import.twisst reads the weights, window data  files into a list object
# If there are multiple weights files, or a single file with different chromosomes/scaffolds/contigs
# in the window data file, these will be separated when importing.

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file, ignore_extra_columns=TRUE, min_combos=50, max_window = 50000)


                             
# make smooth weightings and plot those across chromosomes
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 50000, spacing = 1000)
# twisst_data_smooth5 <- smooth.twisst(twisst_data, span_bp = 5000, spacing = 1000)


################# whole chr plus focus on inv 2

cut_left <- 7829094
cut_right <- 14803486

left <- which(twisst_data_smooth$pos[[1]] < cut_left)
right <- which(twisst_data_smooth$pos[[1]] > cut_right)
gap <- cut_right - cut_left -2e5


png(paste0(prefix,".weights_inv2focus_overlay.png"), width=5000, height=2000, res=400, bg=NA)

if (length(twisst_data$topos) == 3) layout(matrix(c(1,2,3,4,4,4,5,5,5), nrow=3, byrow=T), height=c(2,3,3))
else layout(matrix(c(1:15,rep(16,15)), nrow=2, byrow=T), height=c(1,3))

par(mar=c(1,1,1,1), xpd=NA)

pars = list()

for (i in 1:length(twisst_data$topos)){
    pars[[i]] = par()
    plot.phylo(twisst_data$topos[[i]], type = "cladogram", edge.color=topo_cols[i],
            edge.width=2, label.offset=0.3, cex=1, rotate.tree = 0, no.margin=FALSE, x.lim=c(0,10))
    }

par(mar=c(4,2,2,1), xpd=NA, xaxt="n", yaxt="n")

left_trim = 0e6
right_trim = 18e6-gap

plot.weights(twisst_data_smooth$weights[[1]][left,], twisst_data_smooth$pos[[1]][left],lwd=0.75,stacked=FALSE, xlim=c(left_trim,right_trim),
             fill_cols = NA, line_cols=topo_cols, xlab="", ylab="")

plot.weights(twisst_data_smooth$weights[[1]][right,], twisst_data_smooth$pos[[1]][right] - gap,lwd=0.75,stacked=FALSE, xlim=c(left_trim,right_trim),
             fill_cols = NA, line_cols=topo_cols, add=T, xlab="", ylab="")

arrows(4035981, 1.05, 6220875, 1.05, length=0.1, code=3, lwd=2, col="#006437")
arrows(5322257, 1.05, 6220875, 1.05, length=0.1, code=3, lwd=2, col ="#006437")
arrows(6251636, 1.05, 7829094, 1.05, length=0.1, code=3, lwd=2, col="#a775ee")
arrows(14803486-gap, 1.05, 15299732-gap, 1.05, length=0.1, code=3, lwd=2, col="#f7614d")

rect(6.25e6, -0.02,  6.35e6, -0.08)

par(xaxt="s", yaxt="s")
axis(1, at = seq(0, cut_left, 0.5e6), labels=seq(0, cut_left/1e6, 0.5), line=1)
axis(1, at = seq(signif(cut_right, 2), 18e6, 0.5e6)-gap, labels=seq(signif(cut_right, 2), 18e6, 0.5e6)/1e6, line=1)
axis(2, at = c(0,0.5,1), line=-3)
mtext(side=2, text="Weighting", cex=.75, line=0)
mtext(side=1, text="Position (Mb)", cex=.75, line=3)

par(mar=c(4,2,2,1), xpd=NA, xaxt="n", yaxt="n")

left_trim = 6251636
right_trim = 7829094

rows <- which(twisst_data$window_data[[1]]$start >= left_trim & twisst_data$window_data[[1]]$end <= right_trim)

plot.weights(twisst_data$weights[[1]][rows,], twisst_data$window_data[[1]][rows,c("start","end")],lwd=0,stacked=TRUE, xlim=c(left_trim,right_trim),
             fill_cols = topo_cols, line_cols=NA, xlab="", ylab="")

arrows(6251636, 1.05, 7829094, 1.05, length=0.1, code=3, lwd=2, col="#a775ee")

# rect(6294547, 0,  6299719, -0.1, col="black")
rect(6.25e6, -0.02,  6.35e6, -0.08)

xaxfirst <- round(left_trim/1e6, 1)*1e6

par(xaxt="s", yaxt="s")
axis(1, at = seq(xaxfirst, right_trim, 1e5), labels=seq(xaxfirst, right_trim, 1e5)/1e6, line=1)
axis(2, at = c(0,0.5,1), line=-3)
mtext(side=2, text="Weighting", cex=.75, line=0)
mtext(side=1, text="Position (Mb)", cex=.75, line=3)

dev.off()



######## close up supergene region

cut_left <- 7829094
cut_right <- 14803486

left <- which(twisst_data$window_data[[1]]$start >= 4035981 & twisst_data$window_data[[1]]$start < cut_left)
right <- which(twisst_data$window_data[[1]]$end > cut_right & twisst_data$window_data[[1]]$end <= 15299732)
gap <- cut_right - cut_left -2e5


png(paste0(prefix,".weights_raw_supergene.png"), width=5000, height=2000, res=400, bg=NA)

if (length(twisst_data$topos) == 3) {
    layout(matrix(c(1,2,3,4,4,4), nrow=2, byrow=T), height=c(2,3))
    } else
    {layout(matrix(c(1:15,rep(16,15)), nrow=2, byrow=T), height=c(1,3))}

par(mar=c(1,1,1,1), xpd=NA, xaxt="n")

pars = list()

for (i in 1:length(twisst_data$topos)){
    pars[[i]] = par()
    plot.phylo(twisst_data$topos[[i]], type = "cladogram", edge.color=topo_cols[i],
            edge.width=2, label.offset=0.3, cex=1, rotate.tree = 0, no.margin=FALSE, x.lim=c(0,10))
    }

par(mar=c(4,4,2,1), xpd=NA, xaxt="n")

left_trim = 4e6
right_trim = 15.5e6-gap

plot.weights(twisst_data$weights[[1]][left,], twisst_data$window_data[[1]][left, c("start", "end")],lwd=0, stacked=TRUE, xlim=c(left_trim,right_trim),
             fill_cols = topo_cols, line_cols=NA, xlab="Position (MB)")

plot.weights(twisst_data$weights[[1]][right,], twisst_data$window_data[[1]][right, c("start", "end")] - gap, lwd=0, stacked=TRUE, xlim=c(left_trim,right_trim),
             fill_cols = topo_cols, line_cols=NA, add=T)

arrows(4035981, 1.05, 6220875, 1.05, length=0.1, code=3, lwd=2, col="#006437")
mtext(side=3, at=(4035981+5322257)/2, text="1.1", line=0.5)

arrows(5322257, 1.05, 6220875, 1.05, length=0.1, code=3, lwd=2, col ="#006437")
mtext(side=3, at=(5322257+6220875)/2, text="1.2", line=0.5)

arrows(6251636, 1.05, 7829094, 1.05, length=0.1, code=3, lwd=2, col="#a775ee")
mtext(side=3, at=(6251636+7829094)/2, text="2", line=0.5)

arrows(14803486-gap, 1.05, 15299732-gap, 1.05, length=0.1, code=3, lwd=2, col="#f7614d")
mtext(side=3, at=(14803486-gap+15299732-gap)/2, text="4", line=0.5)

# rect(6294547, 0,  6299719, -0.1, col="black")
# rect(6.25e6, -0.02,  6.35e6, -0.08)

par(xaxt="s")
axis(1, at = seq(left_trim, cut_left, 0.5e6), labels=seq(left_trim, cut_left, 0.5e6)/1e6, line=1)
axis(1, at = seq(signif(cut_right-gap, 2), right_trim, 0.5e6), labels=seq(signif(cut_right, 2), right_trim+gap, 0.5e6)/1e6, line=1)

dev.off()



################# whole chr

cut_left <- 7829094
cut_right <- 14803486

left <- which(twisst_data_smooth$pos[[1]] < cut_left)
right <- which(twisst_data_smooth$pos[[1]] > cut_right)
gap <- cut_right - cut_left -2e5


png(paste0(prefix,".weights_overlay.png"), width=5000, height=2000, res=400, bg=NA)

if (length(twisst_data$topos) == 3) layout(matrix(c(1,2,3,4,4,4), nrow=2, byrow=T), height=c(2,3))
else layout(matrix(c(1:15,rep(16,15)), nrow=2, byrow=T), height=c(1,3))

par(mar=c(1,1,1,1), xpd=NA)

pars = list()

for (i in 1:length(twisst_data$topos)){
    pars[[i]] = par()
    plot.phylo(twisst_data$topos[[i]], type = "cladogram", edge.color=topo_cols[i],
            edge.width=2, label.offset=0.3, cex=1, rotate.tree = 0, no.margin=FALSE, x.lim=c(0,10))
    }

par(mar=c(4,2,2,1), xpd=NA, xaxt="n", yaxt="n")

left_trim = 0e6
right_trim = 18e6-gap

plot.weights(twisst_data_smooth$weights[[1]][left,], twisst_data_smooth$pos[[1]][left],lwd=0.75,stacked=FALSE, xlim=c(left_trim,right_trim),
             fill_cols = NA, line_cols=topo_cols, xlab="", ylab="")

plot.weights(twisst_data_smooth$weights[[1]][right,], twisst_data_smooth$pos[[1]][right] - gap,lwd=0.75,stacked=FALSE, xlim=c(left_trim,right_trim),
             fill_cols = NA, line_cols=topo_cols, add=T, xlab="", ylab="")

arrows(4035981, 1.05, 6220875, 1.05, length=0.1, code=3, lwd=2, col="#006437")
arrows(5322257, 1.05, 6220875, 1.05, length=0.1, code=3, lwd=2, col ="#006437")
arrows(6251636, 1.05, 7829094, 1.05, length=0.1, code=3, lwd=2, col="#a775ee")
arrows(14803486-gap, 1.05, 15299732-gap, 1.05, length=0.1, code=3, lwd=2, col="#f7614d")

rect(6.25e6, -0.02,  6.35e6, -0.08)

par(xaxt="s", yaxt="s")
axis(1, at = seq(0, cut_left, 0.5e6), labels=seq(0, cut_left/1e6, 0.5), line=1)
axis(1, at = seq(signif(cut_right, 2), 18e6, 0.5e6)-gap, labels=seq(signif(cut_right, 2), 18e6, 0.5e6)/1e6, line=1)
axis(2, at = c(0,0.5,1), line=-3)
mtext(side=2, text="Weighting", cex=.75, line=0)
mtext(side=1, text="Position (Mb)", cex=.75, line=3)

dev.off()


