geneID_wrapper = function(input_file){
	proteinIDs = readLines(result_ID_file)

	linked_gene_IDs.unlisted = unlist(entrez_link(dbfrom="protein", id=proteinIDs, db="gene", by_id = TRUE))
	mappable_index = seq(1,length(linked_gene_IDs.unlisted),2)
	linked_gene_IDs = unlist(linked_gene_IDs.unlisted[mappable_index])
	
	return(linked_gene_IDs)
}#end def geneID_wrapper

FASTA_symbol_wrapper = function(starting_GI_file, symbol_mapping, output_file){
	proteinIDs = readLines(starting_GI_file)
	all_recs = entrez_fetch(db="protein", id=proteinIDs, rettype="fasta")
	write(all_recs, file="TEMP.fasta") #used to assist with conversion between R objects

	temp_SeqObj = readAAStringSet("TEMP.fasta")
	temp_SeqObj.seq = as.character(temp_SeqObj)
	temp_SeqObj.seq = list()
	for (i in 1:length(temp_SeqObj)){
		temp_SeqObj.seq[i]=as.character(temp_SeqObj)[i]
	}
	write.fasta(temp_SeqObj.seq, symbol_mapping, file.out=output_file)
	unlink("TEMP.fasta") # Delete temporary file
}#end def FASTA_symbol_wrapper

file_copy = function(file_source, file_destination){
	command = paste("cp", file_source, file_destination, sep=" ")
	system(command)
}#end def file_copy

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
	drawSupportOnEdges(round(100 * msa.tree_boot / B))
	nodelabels(round(100* msa.tree_clade / B))
	legend("bottomleft", legend = c("Bipartitions", "Clades"), pch = 22,
		   pt.bg = c("green", "lightblue"), pt.cex = 2.5)
}#end def ape_bootstrap_wrapper

ape_bootstrap_wrapper.clade_only = function(msa.ape, boot_function, B){
	msa.dist = dist.aa(msa.ape, pairwise.deletion = FALSE, scaled = FALSE)
	msa.tree = tree_function(msa.ape)

	msa.boot_obj = boot.phylo(msa.tree, msa.ape,
									boot_function, B = B, trees = TRUE)
	msa.tree_clade = prop.clades(msa.tree, msa.boot_obj$trees, rooted = FALSE)
	msa.tree_boot = prop.clades(msa.tree, msa.boot_obj$trees)

	layout(1)
	par(mar = rep(2, 4))
	plot(msa.tree, main = "")
	nodelabels(round(100* msa.tree_clade / B))
	legend("bottomleft", legend = c( "Clades"), pch = 22,
		   pt.bg = c( "lightblue"), pt.cex = 2.5)
}#end def ape_bootstrap_wrapper.clade_only

##################
## T Girke Code ##
##################
## UPGMA Tree Building from MSA. Dist ma computed internally.
upgma_bootstrap_aa <- function(aln_mat, B = 100, model = "JTT", seed = 1) {
  build <- function(mat) {
    mat <- toupper(mat)        # lower -> upper
    mat[mat == "-"] <- "?"     # treat gaps as missing (phyDat accepts "?")
    phy <- phyDat(mat, type = "AA")
    as.phylo(hclust(dist.ml(phy, model = model), method = "average"))  # UPGMA
  }
  tr <- build(aln_mat)
  set.seed(seed)
  bs <- boot.phylo(tr, x = aln_mat, FUN = build, B = B)
  tr$node.label <- round(100 * bs / B)
  tr
}
## Usage:
# msa_bios <- readAAMultipleAlignment("alignment.fasta", format="fasta")
# m <- as.matrix(unmasked(msa_bios))                               # each row = sequence
# tree_upgma <- upgma_bootstrap_aa(m)

## Plot UPGMA Tree with Bootstrap Values
plot_tree_bootstrap <- function(tree, tip_cex=.9, bs_cex=.7,
                                 cutoff=50, main="UPGMA with bootstrap (%)",
                                 legend_outside=TRUE) {
  if (is.null(tree$node.label)) stop("tree$node.label is empty.")
  bs  <- as.numeric(tree$node.label)
  lab <- ifelse(bs >= cutoff, as.character(bs), "")
  col <- ifelse(bs < cutoff, "grey70", ifelse(bs < 90, "#3182BD", "#08306B"))

  # add right margin if we want the legend outside
  par(xpd = NA, mar = c(2.5, 1, 2, if (legend_outside) 10 else 2))
  tr <- ladderize(tree)
  plot(tr, cex = tip_cex, no.margin = TRUE); title(main)

  ints <- (Ntip(tr) + 1):(Ntip(tr) + length(tr$node.label))
  nodelabels(text = lab, node = ints, frame = "n",
             adj = c(1.1, -0.25), cex = bs_cex, col = col)
  add.scale.bar()

  if (legend_outside) {
    legend("right", inset = c(-0.16, 0), xpd = NA, bty = "n",
           legend = c(paste0("<", cutoff, "%"), "≥70%", "≥90%"),
           pch = 16, pt.cex = 1,
           col = c("grey70", "#3182BD", "#08306B"),
           title = "Bootstrap")
  } else {
    legend("topright", bty = "n",
           legend = c(paste0("<", cutoff, "%"), "≥70%", "≥90%"),
           pch = 16, pt.cex = 1,
           col = c("grey70", "#3182BD", "#08306B"),
           title = "Bootstrap")
  }
}
# Usage:
# plot_tree_bootstrap(tree_upgma, cutoff = 70, main = "")


## NJ from amino acid alignment (matrix object) with bootstrapping
nj_bootstrap_aa <- function(aln_mat, B = 100, model = "JTT", seed = 1) {
  build <- function(mat) {
    mat <- toupper(mat); mat[mat == "-"] <- "?"     # gaps -> missing, standard
    phy <- phyDat(mat, type = "AA")
    D   <- dist.ml(phy, model = model)
    nj(D)                                            # Neighbor-Joining
  }
  tr <- build(aln_mat)
  set.seed(seed)
  bs <- boot.phylo(tr, x = aln_mat, FUN = build, B = B)
  tr$node.label <- round(100 * bs / B)               # attach bootstrap %
  tr
}
# Usage:
# msa_bios <- readAAMultipleAlignment("alignment.fasta", format="fasta")
# m  <- as.matrix(unmasked(msa_bios))                               # each row = sequence
# tree_nj <- nj_bootstrap_aa(m, B = 500, model = "JTT")
# Midpoint rooting 
# tree_nj_rooted <- phangorn::midpoint(tree_nj)
# Outgroup rooting
# tree_nj_rooted <- ape::root(tree_nj, outgroup = "YourOutgroupName", resolve.root = TRUE)
# Plot tree with same function as UPGMA tree
# plot_tree_bootstrap(tree_nj_rooted, cutoff = 70, main = "")



