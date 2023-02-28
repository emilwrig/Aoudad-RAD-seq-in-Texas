
# no scientific notation
options(scipen=999)

eems <- read.table("aoudad_75_mac2_10kbpthin_struc_eems.recode.vcf")
subset <- read.table("aoudad_75_mac2_10kbpthin_struc_subset.recode.vcf")
structure <- read.table("aoudad_75_mac2_10kbpthin_struc.recode.vcf")
phylo <- read.table("aoudad_75_mac3_10kbpthin_phylo.recode.vcf")

individuals <- read.table("individuals_phylo.txt")

mean_depth <- c()
sd_depth <- c()
phylo_depths <- c()
for(a in 10:ncol(phylo)) {
	a_rep <- phylo[,a]
	a_rep <- sapply(strsplit(a_rep, ":"), "[[", 2)
	a_rep[a_rep == "."] <- "0"
	a_rep <- as.numeric(a_rep)
	phylo_depths <- c(phylo_depths, a_rep)
	mean_depth <- c(mean_depth, mean(a_rep))
	sd_depth <- c(sd_depth, sd(a_rep))
}

structure_depths <- c()
for(a in 10:ncol(structure)) {
	a_rep <- structure[,a]
	a_rep <- sapply(strsplit(a_rep, ":"), "[[", 2)
	a_rep[a_rep == "."] <- "0"
	a_rep <- as.numeric(a_rep)
	structure_depths <- c(structure_depths, a_rep)
}

subset_depths <- c()
for(a in 10:ncol(subset)) {
	a_rep <- subset[,a]
	a_rep <- sapply(strsplit(a_rep, ":"), "[[", 2)
	a_rep[a_rep == "."] <- "0"
	a_rep <- as.numeric(a_rep)
	subset_depths <- c(subset_depths, a_rep)
}

eems_depths <- c()
for(a in 10:ncol(eems)) {
	a_rep <- eems[,a]
	a_rep <- sapply(strsplit(a_rep, ":"), "[[", 2)
	a_rep[a_rep == "."] <- "0"
	a_rep <- as.numeric(a_rep)
	eems_depths <- c(eems_depths, a_rep)
}


individuals <- cbind(individuals, mean_depth, sd_depth)
write.table(individuals, file="coverage_phylo_dataset.txt", sep="\t", row.names=F, quote=F)





datasets <- c("phylo", "structure", "structure_subset", "eems")
n_snps <- c(nrow(phylo), nrow(structure), nrow(subset), nrow(eems))
n_individuals <- c(ncol(phylo) - 9, ncol(structure) - 9, ncol(subset) - 9, ncol(eems) - 9)
mean_depth <- c(mean(phylo_depths), mean(structure_depths), mean(subset_depths), mean(eems_depths))
sd_depth <- c(sd(phylo_depths), sd(structure_depths), sd(subset_depths), sd(eems_depths))

datasets <- data.frame(datasets=datasets, n_snps=n_snps, n_individuals=n_individuals, mean_depth=mean_depth, sd_depth=sd_depth)

write.table(datasets, file="dataset_characteristics.txt", sep="\t", row.names=F, quote=F)



