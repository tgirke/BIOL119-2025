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
