interactive -p quanah

workdir=/lustre/scratch/jmanthey/17_aoudad

cd ${workdir}


###########################
#### prepping and filtering files
###########################


# copy and rename vcf file
cp /lustre/scratch/emilwrig/aoudad/final_analyses/04_populations/populations_75/populations.snps.vcf .
mv populations.snps.vcf aoudad_75.snps.vcf

# get vcf header from STACKS output (to keep track)
grep "#C" aoudad_75.snps.vcf > header.txt


# run vcftools to filter the data
# all samples including outgroup for phylogenetics (75% complete matrix)
vcftools --vcf ${workdir}/aoudad_75.snps.vcf --keep individuals_phylo.txt --max-missing 0.75 --minGQ 20 --minDP 8 --max-meanDP 100 --min-alleles 2 --max-alleles 2 --mac 3 --thin 10000 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/aoudad_75_mac3_10kbpthin_phylo

# all samples excluding outgroup for structure (75% complete matrix)
vcftools --vcf ${workdir}/aoudad_75.snps.vcf --keep individuals_structure.txt --max-missing 0.75 --minGQ 20 --minDP 8 --max-meanDP 100 --min-alleles 2 --max-alleles 2 --mac 2 --thin 10000 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/aoudad_75_mac2_10kbpthin_struc

# all samples excluding outgroup and three divergent individuals for structure (75% complete matrix)
vcftools --vcf ${workdir}/aoudad_75.snps.vcf --keep individuals_structure_subset.txt --max-missing 0.75 --minGQ 20 --minDP 8 --max-meanDP 100 --min-alleles 2 --max-alleles 2 --mac 2 --thin 10000 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/aoudad_75_mac2_10kbpthin_struc_subset

# only Texas for EEMS (75% complete matrix)
vcftools --vcf ${workdir}/aoudad_75.snps.vcf --keep individuals_structure_eems.txt --max-missing 0.75 --minGQ 20 --minDP 8 --max-meanDP 100 --min-alleles 2 --max-alleles 2 --mac 2 --thin 10000 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/aoudad_75_mac2_10kbpthin_struc_eems

# make chromosome map for the vcfs
grep -v "#" ${workdir}/aoudad_75_mac2_10kbpthin_struc.recode.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt
grep -v "#" ${workdir}/aoudad_75_mac2_10kbpthin_struc_eems.recode.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map2.txt


# bgzip and index the vcf files
bgzip -c ${workdir}/aoudad_75_mac2_10kbpthin_struc.recode.vcf > ${workdir}/aoudad_75_mac2_10kbpthin_struc.recode.vcf.gz

bgzip -c ${workdir}/aoudad_75_mac3_10kbpthin_phylo.recode.vcf > ${workdir}/aoudad_75_mac3_10kbpthin_phylo.recode.vcf.gz

bgzip -c ${workdir}/aoudad_75_mac2_10kbpthin_struc_subset.recode.vcf > ${workdir}/aoudad_75_mac2_10kbpthin_struc_subset.recode.vcf.gz

bgzip -c ${workdir}/aoudad_75_mac2_10kbpthin_struc_eems.recode.vcf > ${workdir}/aoudad_75_mac2_10kbpthin_struc_eems.recode.vcf.gz

tabix -p vcf ${workdir}/aoudad_75_mac2_10kbpthin_struc.recode.vcf.gz

tabix -p vcf ${workdir}/aoudad_75_mac3_10kbpthin_phylo.recode.vcf.gz

tabix -p vcf ${workdir}/aoudad_75_mac2_10kbpthin_struc_subset.recode.vcf.gz

tabix -p vcf ${workdir}/aoudad_75_mac2_10kbpthin_struc_eems.recode.vcf.gz


# run bcftools to simplify the vcftools output for the 10kbp spacing for each dataset
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n ' ${workdir}/aoudad_75_mac2_10kbpthin_struc.recode.vcf.gz > ${workdir}/aoudad_75_mac2_10kbpthin_struc.simple.vcf

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n ' ${workdir}/aoudad_75_mac3_10kbpthin_phylo.recode.vcf.gz > ${workdir}/aoudad_75_mac3_10kbpthin_phylo.simple.vcf

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n ' ${workdir}/aoudad_75_mac2_10kbpthin_struc_subset.recode.vcf.gz > ${workdir}/aoudad_75_mac2_10kbpthin_struc_subset.simple.vcf

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n ' ${workdir}/aoudad_75_mac2_10kbpthin_struc_eems.recode.vcf.gz > ${workdir}/aoudad_75_mac2_10kbpthin_struc_eems.simple.vcf


# run vcftools for the vcfs (intended for admixture and PCA and EEMS) with plink output
vcftools --vcf ${workdir}/aoudad_75_mac2_10kbpthin_struc.recode.vcf --plink --chrom-map chrom_map.txt --out ${workdir}/aoudad_75_mac2_10kbpthin_struc

vcftools --vcf ${workdir}/aoudad_75_mac2_10kbpthin_struc_subset.recode.vcf --plink --chrom-map chrom_map.txt --out ${workdir}/aoudad_75_mac2_10kbpthin_struc_subset

vcftools --vcf ${workdir}/aoudad_75_mac2_10kbpthin_struc_eems.recode.vcf --plink --chrom-map chrom_map2.txt --out ${workdir}/aoudad_75_mac2_10kbpthin_eems


# convert with plink for each of the datasets that we will estimate admixture and pca and eems
# all 
plink --file ${workdir}/aoudad_75_mac2_10kbpthin_struc --recode12 --allow-extra-chr --out ${workdir}/aoudad_75_mac2_10kbpthin_struc_plink
# w/o 3 divergent individuals
plink --file ${workdir}/aoudad_75_mac2_10kbpthin_struc_subset --recode12 --allow-extra-chr --out ${workdir}/aoudad_75_mac2_10kbpthin_struc_subset_plink
# Texas EEMS
plink --file ${workdir}/aoudad_75_mac2_10kbpthin_eems --recode12 --allow-extra-chr --out ${workdir}/aoudad_75_mac2_10kbpthin_eems_plink




###########################
#### estimating structure with PCA and ADMIXTURE
###########################



# run each dataset for pca
# all
plink --file ${workdir}/aoudad_75_mac2_10kbpthin_struc_plink --pca --allow-extra-chr --out ${workdir}/aoudad_75_mac2_10kbpthin_struc_plink_pca
# w/o 3 divergent individuals
plink --file ${workdir}/aoudad_75_mac2_10kbpthin_struc_subset_plink --pca --allow-extra-chr --out ${workdir}/aoudad_75_mac2_10kbpthin_struc_subset_plink_pca

# run each dataset for admixture
# all
for K in 1 2 3 4; do admixture --cv ${workdir}/aoudad_75_mac2_10kbpthin_struc_plink.ped $K  | tee log_10kbpthin_aoudad_75_all_${K}.out; done
# w/o 3 divergent individuals
for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv ${workdir}/aoudad_75_mac2_10kbpthin_struc_subset_plink.ped $K  | tee log_10kbpthin_aoudad_75_subset_${K}.out; done



###########################
#### converting files for other analyses
###########################

# load R
module load intel R

R

source("vcf_convert.r")

# convert to related package style format for relatedness estimates
vcf_to_related("aoudad_75_mac2_10kbpthin_struc.simple.vcf", "aoudad_75_mac2_10kbpthin_struc.related", "individuals_structure.txt")

# convert to aligned fasta for raxml
create_fasta_from_vcf("aoudad_75_mac3_10kbpthin_phylo.simple.vcf", "individuals_phylo.txt", "aoudad_75_mac3_10kbpthin_phylo.fasta", 1)

# convert to treemix format
vcf_to_treemix("aoudad_75_mac3_10kbpthin_phylo.simple.vcf", "aoudad.treemix", "popmap_treemix.txt")

q()

gzip aoudad.treemix




