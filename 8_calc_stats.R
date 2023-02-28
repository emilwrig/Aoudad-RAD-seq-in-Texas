
# add functions for calculations
source("reich_fst.r")

# no scientific notation
options(scipen=999)

# read in input file
input_file <- read.table("aoudad_75_mac2_10kbpthin_struc.simple.vcf", stringsAsFactors=F)
input_file <- input_file[,2:ncol(input_file)]

# subset input file 
input_file_genotypes <- input_file[,4:ncol(input_file)]

# read in populations
populations <- read.table("popmap_structure.txt", sep="\t", stringsAsFactors=F, header=T)

# define output name
output_name <- "aoudad_stats.txt"

# write output file
write(c("pop1", "pop2", "stat", "number_sites", "number_variable_sites", "calculated_stat"), ncolumns=6, file=output_name, sep="\t")


# calculate differentiation statistics for each pairwise comparison
all_combinations <- combn(unique(populations$Pop), 2)
for(a in 1:ncol(all_combinations)) {
	# define populations
	a_pop1 <- all_combinations[1,a]
	a_pop2 <- all_combinations[2,a]
	
	# subset vcf inputs
	a_input1 <- input_file_genotypes[,populations$Pop == a_pop1]
	a_input2 <- input_file_genotypes[,populations$Pop == a_pop2]
	
	
	# calculate only in comparisons where both pops have > 1 individual
	if(!is.null(ncol(a_input1)) & !is.null(ncol(a_input2))) {
		
		# remove sites that are not either invariant or bi-allelic SNPs
		a_input1 <- a_input1[nchar(input_file[,2]) == 1 & nchar(input_file[,3]) == 1, ]
		a_input2 <- a_input2[nchar(input_file[,2]) == 1 & nchar(input_file[,3]) == 1, ]
	
		differentiation(a_input1, a_input2, a_pop1, a_pop2, output_name)
	}
}

# calculate heterozygosity
for(a in 1:nrow(populations)) {
	a_rep <- input_file_genotypes[,a]
	a_rep <- a_rep[a_rep != "./."]
	output_rep <- c(populations[a,1], "-", "obs_het", length(a_rep), length(a_rep[a_rep == "0/1"]), length(a_rep[a_rep == "0/1"]) / length(a_rep))
	write(output_rep, ncolumns=6, file=output_name, sep="\t", append=T)
}

# calculate missingness
for(a in 1:nrow(populations)) {
	a_rep <- input_file_genotypes[,a]
	a_rep <- a_rep[a_rep == "./."]
	output_rep <- c(populations[a,1], "-", "missingness", length(a_rep), "-", length(a_rep) / nrow(input_file_genotypes))
	write(output_rep, ncolumns=6, file=output_name, sep="\t", append=T)
}