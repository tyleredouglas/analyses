#takes log2_drop.R output file and gtf file, outputs list of genes within log2drop window

get_genes <- function(log2_file, gtf, output_name, dest) {
	hits <- read.table(log2_file, header = TRUE) 
	dmel <- read.table(gtf, header = FALSE)

	#initialize dataframe to store candidate genes in
	genes <- tibble(gene = character(), peak_location = character(), logp = numeric())

	for (i in 1:nrow(hits)) {
  
	  hit_row <- hits[i,]
	  hit_loc <- as.character(hit_row[1])
	  hit_chr <- as.character(hit_row[3])
	  hit_min <- as.numeric(hit_row[4])
	  hit_max <- as.numeric(hit_row[5])
	  hit_p <- as.numeric(hit_row[2])
  
	  genes_to_add <- data.frame()

	  #subset gtf by log2drop window
	  within_window <- subset(dmel, V1 == hit_chr & V4 > hit_min & V4 < hit_max)

	  #get unique genes (column V10) from gtf
	  genes_to_add <- data.frame(unique(as.character(within_window$V10)))
  
	  if (nrow(genes_to_add) > 0) {

		#retain hit location to ID genes from a shared peak (hitchhiking)
	    genes_to_add$peak_location <- as.character(hit_loc)
	    genes_to_add$logp <- as.numeric(hit_p)
	    genes <- rbind(genes, genes_to_add)
    
 	 }	
}		

	genes <- rename(genes, gene_name = unique.as.character.within_window.V10..)
	genes_filename = paste(dest, "/", output_name, sep = "")
	write.table(genes, genes_filename, row.names=FALSE, sep="\t", quote=FALSE)
}      
