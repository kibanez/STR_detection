require(dplyr)

data = read.csv("/Users/kibanez/Documents/GEL_STR/RpanelApp_dmontaner/panelapp_genes_and_panels.tsv", sep = "\t", header= T)


l_genes = read.csv("/Users/kibanez/Documents/GEL_STR/RpanelApp_dmontaner/list_trinucleotide-repeat-expansions.tsv", sep = '\t', header = F)
l_genes = as.character(l_genes$V1)


for (i in 1:length(l_genes)){

	df_gene = filter(data, GeneSymbol == l_genes[i])

	if (i == 1){

		df = df_gene

	}else{

		df = rbind(df, df_gene)

	}
	

}

write.table(df, '/Users/kibanez/Documents/GEL_STR/RpanelApp_dmontaner/list_panels_including_trinucleotide-repeat-expansions.tsv', sep = '\t', col.names = T, row.names = F, quote = F)


# And also, we are going to generate/create a table containing, in the 1st column the panel and in the 2nd column the list of genes within the NRE that is included in the panel


panel_names = unique(df$Panel_Name)

for (i in 1:length(panel_names)){

	gene_names = as.character(filter(df, Panel_Name == panel_names[i])$GeneSymbol)

	gene_names_string = paste(shQuote(gene_names), collapse=",")

	if (i == 1){

		df_panels = data.frame(panel = as.character(panel_names[i]), genes = gene_names_string)

	}else{

		df_panels = rbind(df_panels, cbind(panel = as.character(panel_names[i]), genes = gene_names_string))
	}

}

write.table(df_panels, '/Users/kibanez/Documents/GEL_STR/RpanelApp_dmontaner/table_list_panel_its_genes_in_NRE.tsv', quote = F, col.names=T, row.names=F, sep = '\t')