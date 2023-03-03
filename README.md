# aoudad_radseq_in_tx
Use these scripts in numerical order to analyze the raw read files in the BioProject PRJNA940156. 

METHODs
1_clean_align.sh -> takes raw reads, removes the EcoRI restriction enzyme, aligns to reference genome
2_stacks.sh -> takes resulting bam files for input into Stacks and produces files for the populations script
3_populations75.sh -> uses the output from Stacks to produce vcf file for further filtering
4_aoudad_filtering.sh -> filter data for all analyses.  Use Plink to convert files for admixture, PCA, and EEMS.  Use Plink PCA analyses.  Use admixture for population structure analyses.  Also convert files for relatedness, raxml, and treemix analyses.
5a_EEMS_matrix.R -> setting the parameters for the EEMS script, calculates genetic distance matrix for the output to be used in 5b_eems.sh
5b_eems.sh -> running EEMS 
5c_plot_eems.R -> plots results of the EEMS in R and used shapefiles to plot results on a map of Texas
6_run_treemix.sh -> script to run treemix analysis. Specified a maximum of 5 migration edges.  100 bootstraps for 200 SNP bootrap blocks.
7_raxml.sh -> running raxml  
8_calc_stats.R -> calculate Reich's Fst, missingness, and observed heterozygosity
9_ibd.R -> using output Fst values and geographic distance to run a regression and making a plot 
10_plot_structure.R -> use the plink outputs of the PCA and plot
11_relatedness.r -> use the filtered vcf files to estimate relatedness values using the different estimation methods (i.e., Wang, King-robust, etc.)
12_coverage_missing_snps_summaries.R -> determine summary stats for each dataset (i.e., number of SNPs, mean coverage depth, , etc.)
