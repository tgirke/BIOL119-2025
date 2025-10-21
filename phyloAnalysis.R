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

library(rentrez)
library(Biostrings)
library(org.Hs.eg.db)
library(seqinr)#to create intermediate FASTA file
library(msa)
library(ape)
library(pheatmap)
source("phyloAnalysis_Fct.R") # as long as this is saved in the same location

# Define Gene IDs (using `rentrez`)
linked_gene_IDs = geneID_wrapper(result_ID_file)


# Define Gene Symbols (using `org.Hs.eg.db`)
gene_symbols = mapIds( x = org.Hs.eg.db, keys = linked_gene_IDs,
					column = "SYMBOL", keytype = "ENTREZID")
FASTA_symbol_wrapper(result_ID_file, gene_symbols, FASTA_file)


# Perform Multiple Sequence Alignment (using `msa`)
P450_SeqObj = readAAStringSet(FASTA_file)
msa_ClustalOmega = msa(P450_SeqObj, method = "ClustalOmega")

#Visualize Alignment Near Invariant Cysteine
msa_basename = "rentrez-msa_ClustalOmega-P452_UniProtKB_SwissProt.msaPrettyPrint.pdf" #please note the folder that is saved for GitHub visualization
msa_pdf = paste("./data/phyloAnalysis",msa_basename,sep="/") #please note the folder that is saved for GitHub visualization
msaPrettyPrint(msa_ClustalOmega, file = msa_basename,
			y=c(540, 640), output="pdf", showNames="left", 
			logoColors="rasmol", showLogo="top",
			shadingMode="identical", shadingModeArg =c(50,90), shadingColors = "grays",
			askForOverwrite=FALSE)
file_copy(msa_basename, msa_pdf) #make sure that the PDF is saved in the intended location

msa_ClustalOmega.ape = msaConvert(msa_ClustalOmega, type="ape::AAbin") #conversion of R object for next step

# Visualize Distance Matrix (using `pheatmap`)
msa_ClustalOmega.dist = dist.aa(msa_ClustalOmega.ape, pairwise.deletion = FALSE, scaled = FALSE)
pheatmap(as.matrix(msa_ClustalOmega.dist))

# Create Phylogenetic Tree with Bootstrap Support (using `ape`)
ape_tree_file = "rentrez-msa_ClustalOmega-P452_UniProtKB_SwissProt.ape_tree.bootstrap-gene_symbol.png" #for OUTPUT, not input
tree_function = function(x) nj(dist.aa(x))
bootstrap_n = 100 #should really use at least 1000, but we are reducing time and resources
ape_bootstrap_wrapper(msa_ClustalOmega.ape, tree_function, bootstrap_n)