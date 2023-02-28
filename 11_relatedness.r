options(scipen=999)
library(ggplot2)
library(related)
library(RColorBrewer)

# read in genotypes  
data_file <- "aoudad_75_mac2_10kbpthin_struc.related"
x <- read.table(data_file, stringsAsFactors=F)

# all dots are missing data, replace those with zeros
for(a in 2:ncol(x)) {
	x[x[,a] == ".",a] <- 0
}

# input species name and dataset
species <- "aoudad"

# rename
total <- x

# measure relatedness of all individuals
total_relate <- coancestry(total, lynchli=2, wang=2)

# estimate the R0, R1, and KING-robust kinship stats for each comparison

# first rearrange the data matrix
x_names <- total[,1]
total2 <- total[,2:ncol(total)]
total2 <- data.frame(x_names, mapply(paste0, total2[c(T,F)], total2[c(F,T)]))

# convert all genotypes to 00, 01, 11, or missing
for(a in 2:ncol(total2)) {
	total2[total2[,a] == "00",a] <- "??"	
	total2[total2[,a] == "55",a] <- "00"
	total2[total2[,a] == "51",a] <- "01"
	total2[total2[,a] == "15",a] <- "01"
}

# determine all the pairwise comparisons desired
comps1 <- total2[1:73,1]
comps <- t(combn(comps1, 2))

output <- data.frame(ind1=as.character(comps[,1]), ind2=as.character(comps[,2]))

# calculate the A, B, C, D, E, F, G, H, I categories for each pairwise comparison
AA <- list()
BB <- list()
CC <- list()
DD <- list()
EE <- list()
FF <- list()
GG <- list()
HH <- list()
II <- list()
for(a in 1:nrow(output)) {
	if(a %% 1000 == 0) { print(a) }
	a_rep <- total2[total2[,1] == comps[a,1] | total2[,1] == comps[a,2], ]
	a_rep <- t(a_rep[,2:ncol(a_rep)])
	AA[[a]] <- length(a_rep[a_rep[,1] == "00" & a_rep[,2] =="00",1])
	BB[[a]] <- length(a_rep[a_rep[,1] == "01" & a_rep[,2] =="00",1])
	CC[[a]] <- length(a_rep[a_rep[,1] == "11" & a_rep[,2] =="00",1])
	DD[[a]] <- length(a_rep[a_rep[,1] == "00" & a_rep[,2] =="01",1])
	EE[[a]] <- length(a_rep[a_rep[,1] == "01" & a_rep[,2] =="01",1])
	FF[[a]] <- length(a_rep[a_rep[,1] == "11" & a_rep[,2] =="01",1])
	GG[[a]] <- length(a_rep[a_rep[,1] == "00" & a_rep[,2] =="11",1])
	HH[[a]] <- length(a_rep[a_rep[,1] == "01" & a_rep[,2] =="11",1])
	II[[a]] <- length(a_rep[a_rep[,1] == "11" & a_rep[,2] =="11",1])
}
output <- data.frame(ind1=as.character(output$ind1), ind2=as.character(output$ind2), A=as.numeric(unlist(AA)), B=as.numeric(unlist(BB)), C=as.numeric(unlist(CC)), D=as.numeric(unlist(DD)), E=as.numeric(unlist(EE)), F=as.numeric(unlist(FF)), G=as.numeric(unlist(GG)), H=as.numeric(unlist(HH)), I=as.numeric(unlist(II)))

# calculate the R0, R1, and KING-robust kinship stats for each pairwise comparison
R0 <- (output$C + output$G) / (output$E)
R1 <- (output$E) / (output$B + output$D + output$H + output$F + output$C + output$G)
KINGrobust <- (output$E - 2 * (output$C + output$G)) / (output$B + output$D + output$H + output$F + 2 * output$E)

# add these to output
output <- cbind(output, R0, R1, KINGrobust)

# add the wang estimates to the output matrix
add_related <- total_relate$relatedness

# setup for matching truncated names from the related output
output_names_reduced <- output[,1:2]
output_names_reduced[,1] <- substr(output_names_reduced[,1], 1, 20)
output_names_reduced[,2] <- substr(output_names_reduced[,2], 1, 20)

wang <- list()
for(a in 1:nrow(output)) {
	a_rep <- add_related[add_related$ind1.id == output_names_reduced$ind1[a] & add_related$ind2.id == output_names_reduced$ind2[a],]
	wang[[a]] <- a_rep$wang[1]
}
wang <- unlist(wang)
output <- cbind(output, wang)

# write the output
write.table(output, file=paste("output_", species, "__relatedness.txt", sep=""), sep="\t", quote=F, col.names=T, row.names=F)

# save the work
save.image(paste(species, "__relatedness.RData", sep=""))



x <- read.table("output_aoudad__relatedness.txt", header=T)

# pull out all comparisons inferred parent-offspring or full sibs
x_PO <- x[x$R0 <= 0.025 & x$KINGrobust >= 0.2, ]
x_sibs <- x[x$R0 > 0.025 & x$KINGrobust >= 0.2, ]

# write these to tables
write.table(x_PO, file="output_predicted_parent-offspring.txt", sep="\t", row.names=F, quote=F)
# no sibs relationships write.table(x_sibs, file="output_predicted_full-siblings.txt", sep="\t", row.names=F, quote=F)


# plot relatedness across all individuals
x <- read.table("output_aoudad__relatedness.txt", header=T)

plot_colors <- brewer.pal(12, "Paired") # plot colors

a_PO <- x[x$R0 <= 0.025 & x$KINGrobust >= 0.2, ]
a_sibs <- x[x$R0 > 0.025 & x$KINGrobust >= 0.2, ]
a_partly <- x[x$KINGrobust >= 0.1 & x$KINGrobust < 0.2, ]
a_not_closely_related <- x[x$KINGrobust < 0.1, ]
	
par(mar=c(4.5, 4.5, 2, 1))
par(mfrow=c(1,2))
plot(a_not_closely_related$R1, a_not_closely_related$R0, xlim=c(min(x$R1), max(x$R1)), ylim=c(min(x$R0), 1), pch=19, cex=0.4, xlab="R1", ylab="R0", col=plot_colors[7], cex.axis=1.2, cex.lab=1.4)
points(a_PO$R1, a_PO$R0, pch=19, cex=0.4, col=plot_colors[10])
points(a_sibs$R1, a_sibs$R0, pch=19, cex=0.4, col=plot_colors[2])
points(a_partly$R1, a_partly$R0, pch=19, cex=0.4, col=plot_colors[4])
plot(a_not_closely_related$R1, a_not_closely_related$KINGrobust, xlim=c(min(x$R1), max(x$R1)), ylim=c(min(x$KINGrobust), max(x$KINGrobust)), pch=19, cex=0.4, xlab="R1", ylab="KING-robust Kinship", col=plot_colors[7], cex.axis=1.2, cex.lab=1.4)
points(a_PO$R1, a_PO$KINGrobust, pch=19, cex=0.4, col=plot_colors[10])
points(a_sibs$R1, a_sibs$KINGrobust, pch=19, cex=0.4, col=plot_colors[2])
points(a_partly$R1, a_partly$KINGrobust, pch=19, cex=0.4, col=plot_colors[4])

plot(a_not_closely_related$wang, a_not_closely_related$KINGrobust, xlim=c(min(x$wang), max(x$wang)), ylim=c(min(x$KINGrobust), max(x$KINGrobust)), pch=19, cex=0.4, xlab="Relatedness (Wang 2002)", ylab="KING-robust Kinship", col=plot_colors[7], cex.axis=1.2, cex.lab=1.4)
points(a_PO$wang, a_PO$KINGrobust, pch=19, cex=0.4, col=plot_colors[10])
points(a_sibs$wang, a_sibs$KINGrobust, pch=19, cex=0.4, col=plot_colors[2])
points(a_partly$wang, a_partly$KINGrobust, pch=19, cex=0.4, col=plot_colors[4])
abline(lm(x$KINGrobust ~ x$wang))
summary(lm(x$KINGrobust ~ x$wang))

