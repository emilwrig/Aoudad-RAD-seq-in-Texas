options(scipen=999)

library(ggplot2)
library(RColorBrewer)

popmap <- read.table("popmap_structure.txt", header=T)
popmap2 <- popmap[popmap$Region != "WildlifeCenterK",]

pca_all <- read.table("aoudad_75_mac2_10kbpthin_struc_plink_pca.eigenvec")
pca_all <- pca_all[,c(1,3,4,5,6)]
colnames(pca_all) <- c("ID", "PC1", "PC2", "PC3", "PC4")

pca_subset <- read.table("aoudad_75_mac2_10kbpthin_struc_subset_plink_pca.eigenvec")
pca_subset <- pca_subset[,c(1,3,4,5,6)]
colnames(pca_subset) <- c("ID", "PC1", "PC2", "PC3", "PC4")

plot_colors <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(n = 8, name = "Dark2"))

# plot PCA 
par(mfrow=c(1,2))
par(mar=c(4.5,4.5,1,1))
par(lwd=1)
# all
plot(pca_all[,2:3], pch=19, cex=0.1, col="white", xlab="PC1 (11.6% Variance Explained)", ylab="PC2 (3.7% Variance Explained)")
points(pca_all[pca_all[,2] > 0.4,2:3], pch=21, cex=1.2, bg=plot_colors[6])
points(pca_all[pca_all[,2] <= 0.4,2:3], pch=21, cex=1.2, bg=plot_colors[8])
plot_pca <- cbind(pca_subset, popmap2$Pop, popmap2$Region)
colnames(plot_pca) <- c("ID", "PC1", "PC2", "PC3", "PC4", "Population", "Region")
plot(plot_pca[,2:3], pch=19, cex=0.1, col="white", xlab="PC1 (4.0% Variance Explained)", ylab="PC2 (3.1% Variance Explained)")
points(plot_pca[plot_pca$Population == "Carrizo",2:3], pch=21, cex=1.2, bg= plot_colors[1])
points(plot_pca[plot_pca$Population == "VanHorn",2:3], pch=21, cex=1.2, bg= plot_colors[2])
points(plot_pca[plot_pca$Population == "SierraViejas",2:3], pch=21, cex=1.2, bg= plot_colors[2])
points(plot_pca[plot_pca$Population == "Chinatis",2:3], pch=21, cex=1.2, bg= plot_colors[3])
points(plot_pca[plot_pca$Population == "DavisMtns",2:3], pch=21, cex=1.2, bg= plot_colors[4])
points(plot_pca[plot_pca$Population == "GlassMtns",2:3], pch=21, cex=1.2, bg= plot_colors[5])
points(plot_pca[plot_pca$Population == "ElephantMtn",2:3], pch=21, cex=1.2, bg= plot_colors[6])
points(plot_pca[plot_pca$Population == "BlackGap",2:3], pch=21, cex=1.2, bg= plot_colors[7])
points(plot_pca[plot_pca$Population == "PDC",2:3], pch=24, cex=1.2, bg= plot_colors[1])
points(plot_pca[plot_pca$Population == "CCSP",2:3], pch=24, cex=1.2, bg= plot_colors[2])
points(plot_pca[plot_pca$Population == "Post",2:3], pch=22, cex=1.2, bg= plot_colors[3])
points(plot_pca[plot_pca$Population == "DolanFalls",2:3], pch=22, cex=1.2, bg= plot_colors[4])
points(plot_pca[plot_pca$Population == "KerrWMA",2:3], pch=22, cex=1.2, bg= plot_colors[5])
points(plot_pca[plot_pca$Population == "LoveCreek",2:3], pch=22, cex=1.2, bg= plot_colors[6])
points(plot_pca[plot_pca$Population == "FawcettWMA",2:3], pch=23, cex=1.2, bg= plot_colors[1])
points(plot_pca[plot_pca$Population == "FRWC_G",2:3], pch=4, cex=1.2, col= plot_colors[2])
points(plot_pca[plot_pca$Population == "NewMexico",2:3], pch=21, cex=1.2, col= plot_colors[7])
points(plot_pca[plot_pca$Population == "California",2:3], pch=21, cex=1.2, col= plot_colors[6])


