#Upload file following these instructions
#2 Files to Upload:
#--> sequence.gi
#--> [this file]

#Need to run the following in Posit Cloud:
#install.packages("rentrez")
#install.packages("seqinr") #for intermediate file creation
#install.packages("ape")
#
#install.packages("BiocManager")
#####################################################################
### [the steps below can take some time - perhaps several minutes] ###
#####################################################################
#BiocManager::install("org.Hs.eg.db") #for gene symbol mapping
#BiocManager::install("msa")
#
#source("[renamed version of this file].R")
#
#[I can also create an R Markdown file, to help with step-wise execution]
#
###################################################################
### Otherwise, this WORKED and completed ALMOST INSTANTANEOUSLY ###
###################################################################

result_ID_file = "sequence.gi"
FASTA_file = "rentrez-P450_UniProtKB_SwissProt-gene_symbol.fasta"

library(rentrez)
library(org.Hs.eg.db)
library(seqinr)#to create intermediate FASTA file
library(msa)
library(ape)

#######################################
### Add Code for `rentrez` Analysis ###
#######################################

#Consider the following resource : https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html

proteinIDs = readLines(result_ID_file)
all_recs = entrez_fetch(db="protein", id=proteinIDs, rettype="fasta")
linked_gene_IDs = entrez_link(dbfrom="protein", id=proteinIDs, db="gene")
write(all_recs, file="TEMP.fasta")

gene_symbols = mapIds( x = org.Hs.eg.db, keys = linked_gene_IDs$links$protein_gene,
					column = "SYMBOL", keytype = "ENTREZID")

temp_SeqObj = readAAStringSet("TEMP.fasta")
temp_SeqObj.seq = as.character(temp_SeqObj)
temp_SeqObj.seq = list()
for (i in 1:length(temp_SeqObj)){
	temp_SeqObj.seq[i]=as.character(temp_SeqObj)[i]
}
write.fasta(temp_SeqObj.seq, gene_symbols, file.out=FASTA_file)
unlink("TEMP.fasta") # Delete temporary file
#after creating the FASTA file, it can be downloaded with the following instructions : https://forum.posit.co/t/how-does-one-download-files-from-rstudio-cloud-onto-desktop/52132

#########################################
###    Alignment steps from `msa`     ###
### Choose One Method : Clustal Omega ###
#########################################
tree_function = function(x) nj(dist.aa(x))
bootstrap_n = 100 #should really use at least 1000, but we are reducing time and resources

P450_SeqObj = readAAStringSet(FASTA_file)

msa_ClustalOmega = msa(P450_SeqObj, method = "ClustalOmega")
msa_ClustalOmega.ape = msaConvert(msa_ClustalOmega, type="ape::AAbin")
msa_ClustalOmega.dist = dist.aa(msa_ClustalOmega.ape, pairwise.deletion = FALSE, scaled = FALSE)
msa_ClustalOmega.tree = tree_function(msa_ClustalOmega.ape)
msa_ClustalOmega.boot_obj = boot.phylo(msa_ClustalOmega.tree, msa_ClustalOmega.ape,
								tree_function, B = bootstrap_n, trees = TRUE)
msa_ClustalOmega.tree_clade = prop.clades(msa_ClustalOmega.tree, msa_ClustalOmega.boot_obj$trees, rooted = FALSE)
msa_ClustalOmega.tree_boot = prop.clades(msa_ClustalOmega.tree, msa_ClustalOmega.boot_obj$trees)
ape_tree_file = "rentrez-msa_ClustalOmega-P450_UniProtKB_SwissProt.ape_tree.bootstrap-gene_symbol.png"
png(ape_tree_file, height=1200, width=2000)
layout(1)
par(mar = rep(2, 4))
plot(msa_ClustalOmega.tree, main = "")
drawSupportOnEdges(round(100 * msa_ClustalOmega.tree_boot / bootstrap_n))
nodelabels(round(100* msa_ClustalOmega.tree_clade / bootstrap_n))
legend("bottomleft", legend = c("Bipartitions", "Clades"), pch = 22,
       pt.bg = c("green", "lightblue"), pt.cex = 2.5)
dev.off()