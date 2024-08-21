
library(ggplot2)
library(scatterpie)

hap_cols <- c(orientis="#2143d1",klugii="#ffac07",chrysippus="#bc4754",karamu="#ff8cf1", neoW="#7a52a6", other="gray80")

tip_locations <- read.table("DC174OG.contig15.1.realigned.DP7GQ30het75min87minVar2.regions1-4.dist.circles.locations.tsv", as.is=TRUE, row.names=1)

IDs <- read.table("Admixture_results_for_Simon/admix_in_chr15_onlySVs.fam", as.is=T)[,1]

file <- "admix_in_chr15_onlySVs.3.Q"
cluster_names <- c("chrysippus", "orientis", "klugii")

file <- "admix_in_chr15_onlySVs.4.Q"
cluster_names <- c("orientis", "karamu", "klugii", "chrysippus")

file <- "admix_in_chr15_onlySVs.5.Q"
cluster_names <- c("orientis", "klugii", "karamu", "chrysippus", "neoW")

file <- "admix_in_chr15_onlySVs.6.Q"
cluster_names <- c("klugii", "orientis", "chrysippus", "other", "karamu", "neoW")

file <- "admix_in_chr15_onlySVs.7.Q"
cluster_names <- c("klugii", "chrysippus", "karamu", "other", "other2", "neoW", "orientis")


proportions <- read.table(paste0("Admixture_results_for_Simon/", file), as.is=T)
proportions <- round(proportions, 3)
df <- as.data.frame(cbind(tip_locations[IDs,], proportions))
names(df) <- c("x", "y", cluster_names)


#### plot tips of neighbor net tree

ggplot() +
geom_scatterpie(aes(x = x, y = y, r=6), data = df, cols = cluster_names, color="black") +
scale_fill_manual("Cluster assignment", values = hap_cols[cluster_names]) +
coord_equal() + 
theme(panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill='transparent'),
      legend.box.background = element_rect(fill='transparent'),
      axis.line=element_blank(), axis.text.x=element_blank(),
      axis.text.y=element_blank(), axis.ticks=element_blank(),
      axis.title.x=element_blank(),axis.title.y=element_blank()
      )

ggsave(paste0(file, ".pdf"), bg='transparent')


plot(df$x, df$y, cex=0)
text(df$x, df$y, labels = rownames(df), cex=0.6)


### plot population groups

sample_data <- read.csv("DC174_data.csv", as.is=T)

popnames <- unique(sample_data$region_narrow)

samples_by_pop <- sapply(popnames, function(pop) sample_data[sample_data$region_narrow == pop,1], simplify=F)


for (pop in popnames){
    pop_df <- df[samples_by_pop[[pop]],]
    if (nrow(pop_df) > 16) pop_df <- pop_df[1:16,]
    n <- nrow(pop_df)
    s <- ceiling(sqrt(n))
    pop_df$x <- rep(1:s,s)[1:n]
    pop_df$y <- -1* rep(1:s,each=s)[1:n]
    ggplot() +
    geom_scatterpie(aes(x = x, y = y, r=0.3), data = pop_df, cols = cluster_names, color="black") +
    scale_fill_manual("Cluster assignment", values = hap_cols[cluster_names]) +
    coord_equal() + 
    expand_limits(x=c(1, 4), y=c(-1, -4)) +
    theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank()
        )
    ggsave(paste0(file,".", pop, ".pdf"), bg='transparent')
    }
