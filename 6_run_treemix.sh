# run with up to 5 migration edges:
src/treemix -i input/aoudad.treemix.gz -root Outgroup -o output/aoudad.treemix
src/treemix -i input/aoudad.treemix.gz -m 1 -g output/aoudad.treemix.vertices.gz output/aoudad.treemix.edges.gz -o output/aoudad_m1.treemix
src/treemix -i input/aoudad.treemix.gz -m 1 -g output/aoudad_m1.treemix.vertices.gz output/aoudad_m1.treemix.edges.gz -o output/aoudad_m2.treemix
src/treemix -i input/aoudad.treemix.gz -m 1 -g output/aoudad_m2.treemix.vertices.gz output/aoudad_m2.treemix.edges.gz -o output/aoudad_m3.treemix
src/treemix -i input/aoudad.treemix.gz -m 1 -g output/aoudad_m3.treemix.vertices.gz output/aoudad_m3.treemix.edges.gz -o output/aoudad_m4.treemix
src/treemix -i input/aoudad.treemix.gz -m 1 -g output/aoudad_m4.treemix.vertices.gz output/aoudad_m4.treemix.edges.gz -o output/aoudad_m5.treemix


#bootstraps over 200s snps 
for i in {1..100}; do
    src/treemix -i input/aoudad.treemix.gz -bootstrap -k 200 -o output/$i.treemix
done;

# unzip the tree files
for i in $( ls *treeout.gz ); do
    gzip -d $i
done;

# in R:
#summarize bootstraps
x <- list.files(pattern="*treeout")
for(a in 1:length(x)) {
	if (a==1) {
		output <- scan(x[a], what="character")[1]
	} else {
		output <- c(output, scan(x[a], what="character")[1])
	}
}
write(output, file="aoudad_bootstraps.trees", ncolumns=1)

# in bash 
# summarize bootstraps
sumtrees.py --output=aoudad_treemix_summed.tre --min-clade-freq=0.01 aoudad_bootstraps.trees