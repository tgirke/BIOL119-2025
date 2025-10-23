#GitHub Repository Contains Needed Files

#If setting up workspace, need to run the following in Posit Cloud:
#install.packages("rentrez")
#install.packages("seqinr") #for intermediate file creation
#install.packages("ape")
#install.packages("pheatmap")
#
#install.packages("BiocManager")
######################################################################
### [the steps below can take some time - perhaps several minutes] ###
######################################################################
#BiocManager::install("org.Hs.eg.db") #for gene symbol mapping
#BiocManager::install("msa")

result_ID_file = "data/phyloAnalysis/sequence.gi"
FASTA_file = "data/phyloAnalysis/rentrez-P450_UniProtKB_SwissProt-gene_symbol.fasta"

library(msa)
library(ape)

# Perform Multiple Sequence Alignment (using `msa`)
# NOTE : Starting Analysis from Gene-Symbol FASTA File
P450_SeqObj = readAAStringSet(FASTA_file)
msa_ClustalOmega = msa(P450_SeqObj, method = "ClustalOmega")

#SKIP - Visualize Alignment Near Invariant Cysteine
#...this is not affected by Clustering Method

#SKIP -  Visualize Distance Matrix (using `pheatmap`)
#...this is not affected by Clustering Method

# Create Phylogenetic Tree with Bootstrap Support (using `ape`)
ape_tree_file = "./data/phyloAnalysis/rentrez-msa_ClustalOmega-P452_UniProtKB_SwissProt.ape_tree-UPGMA.bootstrap-gene_symbol.png" #for OUTPUT, not input
bootstrap_n = 100 #should really use at least 1000, but we are reducing time and resources

upgma = function(x) hclust(x, "average")
#use `as.phylo`, as decribed here : https://cran.r-project.org/web/packages/ape/refman/ape.html#as.phylo
tree_function = function(x) as.phylo(upgma(dist.aa(x)))
ape_bootstrap_wrapper = function(msa.ape, boot_function, B){
	msa.dist = dist.aa(msa.ape, pairwise.deletion = FALSE, scaled = FALSE)
	msa.tree = tree_function(msa.ape)

	msa.boot_obj = boot.phylo(msa.tree, msa.ape,
									boot_function, B = B, trees = TRUE)
	msa.tree_clade = prop.clades(msa.tree, msa.boot_obj$trees, rooted = FALSE)
	msa.tree_boot = prop.clades(msa.tree, msa.boot_obj$trees)

	layout(1)
	par(mar = rep(2, 4))
	plot(msa.tree, main = "")
	#drawSupportOnEdges(round(100 * msa.tree_boot / B))
	nodelabels(round(100* msa.tree_clade / B))
	#legend("bottomleft", legend = c("Bipartitions", "Clades"), pch = 22,
	#	   pt.bg = c("green", "lightblue"), pt.cex = 2.5)
	legend("topleft", legend = c("Clades"), pch = 22,
		   pt.bg = c("lightblue"), pt.cex = 2.5)
}#end def ape_bootstrap_wrapper

msa_ClustalOmega.ape = msaConvert(msa_ClustalOmega, type="ape::AAbin") #conversion of R object for next step
msa_ClustalOmega.dist = dist.aa(msa_ClustalOmega.ape, pairwise.deletion = FALSE, scaled = FALSE)
png(ape_tree_file, height=400, width=800)
ape_bootstrap_wrapper(msa_ClustalOmega.ape, tree_function, bootstrap_n)
dev.off()