################################################################################################################################################################################################################################################
# 2017.02.03
# csv's provided by David Montaner (/Users/kibanez/Documents/GEL_STR/labKey_Catalog_2016.12.20/rd_pilot_labkey2catalog_2016_12_20/*.csv)

# panels.csv contains `gelFamilyId`, `gelId` (to extract) and `panel`. We will extract all those participants (`gelFamilyId` and `gelId`) with `ID` and `HA` as panel assigned
# samples.csv contains `gelID` and `plateKey`. With the list of `gelId` retrieved from panels.csv, we'll extract all those `plateKey` corresponding to `gelId`
# registration.csv contains `gelId` and `gelFamilyId`, in order to extract all those `gelID` that are members within a family given a `gelFamilyId` 


# we are not interested only in the participants, but also in the rest of their familiars. With the `gelFamilyId` we will extract all the members of the family as well
# the difference comparing to the `ID` and `HA` is that here we are going to pass to the script a table or list containing a the name of the panel of names, and we are going to 
# generate/create an excel or tsv for each of them, containing the probands/patients + familiars that have been assigned that panel of genes

################################################################################################################################################################################################################################################
require(dplyr)

setwd("/Users/kibanez/Documents/GEL_STR/labKey_Catalog_2016.12.20/rd_pilot_labkey2catalog_2016_12_20/")

df_panels = read.csv('panels.csv', sep=',', header=T)

df_samples = read.csv('samples.csv', sep=',', header=T)

df_registration = read.csv('registration.csv', sep=',', header=T)



l_panel_genes = read.csv("/Users/kibanez/Documents/GEL_STR/RpanelApp_dmontaner/list_panels_including_trinucleotide-repeat-expansions.tsv", sep = '\t', header = T)
l_panel_genes = as.character(l_panel_genes$Panel_Name)

# we want to keep only the uniqeu ones, remove the repetitives
l_panel_genes = unique(l_panel_genes)


for (i in 1:length(l_panel_genes)){

	panel_name = l_panel_genes[i]
	df_aux = filter(df_panels, panel == panel_name)

	# 1) We first extract the `gelFamilyId` and `gelId` 
	l_gel_family_id = as.character(df_aux$gelFamilyId)

	# 2) Extract all the `gelId`-s included in each `gelFamilyId`
	l_gel_id = c()
	l_affected = c()

	for (j in 1:length(l_gel_family_id)){

		index = which(as.character(df_registration$gelFamilyId) %in% l_gel_family_id[j])

		l_gel_id = c(l_gel_id, as.character(df_registration$gelId[index]))
		l_affected = c(l_affected, as.character(df_registration$disease_status[index]))

	} # for 

	# 3) Extract from df_samples the LP numbers
	l_lp_numbers = as.character(df_samples$plateKey[pmatch(l_gel_id, df_samples$gelID)])

	# 4) we incorporate now the path info as well
	data = read.csv('/Users/kibanez/Documents/GEL_STR/labKey_Catalog_2016.12.20/RD_pilot_sample_metrics.txt',sep='\t', header = T)

	l_path_id = c()

	for (j in 1:length(l_lp_numbers)){
		
		index = which(data$iid == l_lp_numbers[j])

		if (length(index) == 0){

			l_path_id = c(l_path_id, '.')

		}
		else if (length(index) > 1){

			# we include in the analysis the most recent one, the latest one
			index2 = index[length(index)]
			l_path_id = c(l_path_id, as.character(data$path[index2]))

		}else{

			l_path_id = c(l_path_id, as.character(data$path[index]))

		}
	}# for

	df_output = cbind(lp = l_lp_numbers, path = l_path_id, status = l_affected)

	output_tsv = paste(panel_name, 'participants_and_familiars_enriched.tsv' , sep= '_')

	write.table(df_output, output_tsv , sep='\t', quote = F, row.names = F, col.names = T)

	if (i == 1){

		df_all = df_output
	}else if (length(l_gel_family_id) > 0) {

		df_all = rbind(df_all, df_output)
	}

} # for i


# special cases
# IN panelApp even if the panel name is one, sometimes in Catalog the name in which this panel appears, is one of the name appeared in `relevant disorders`
# For instance, for parkinson's disease; 
# Relevant disorders: Complex Parkinsonism (includes pallido-pyramidal syndromes), Early onset and familial Parkinson's Disease 

df_aux = df_panels[grep("Parkinson", df_panels$panel),]
df_aux = df_panels[grep("sudden death", df_panels$panel),]
df_aux = filter(df_panels, panel == 'Infantile pseudo-obstruction')

df_aux = filter(df_panels, panel == 'Lactic acidosis')
df_aux = rbind(df_aux, filter(df_panels, panel == 'All recognised syndromes and those with suggestive features'))

df_aux = filter(df_panels, panel == 'Epilepsy plus other features')
df_aux = filter(df_panels, panel == "Early onset familial premature ovarian failure")


df_aux = filter(df_panels, panel == 'Vein of Galen malformation')
df_aux = filter(df_panels, panel == 'Cerebral vascular malformations')
df_aux = rbind(df_aux, filter(df_panels, panel == 'Cerebral arteriovenous malformations'))
df_aux = rbind(df_aux, filter(df_panels, panel == 'Moyamoya disease'))

 
df_aux = filter(df_panels, panel == 'Bilateral microtia')
#df_aux = filter(df_panels, panel == 'Ear malformations with hearing impairment')
df_aux = rbind(df_aux, filter(df_panels, panel == 'Ear malformations'))
df_aux = rbind(df_aux, filter(df_panels, panel == 'Familial hemifacial microsomia'))




	l_gel_family_id = as.character(df_aux$gelFamilyId)

	# 2) Extract all the `gelId`-s included in each `gelFamilyId`
	l_gel_id = c()
	l_affected = c()

	for (j in 1:length(l_gel_family_id)){

		index = which(as.character(df_registration$gelFamilyId) %in% l_gel_family_id[j])

		l_gel_id = c(l_gel_id, as.character(df_registration$gelId[index]))
		l_affected = c(l_affected, as.character(df_registration$disease_status[index]))

	} # for 

	# 3) Extract from df_samples the LP numbers
	l_lp_numbers = as.character(df_samples$plateKey[pmatch(l_gel_id, df_samples$gelID)])

	# 4) we incorporate now the path info as well
	data = read.csv('/Users/kibanez/Documents/GEL_STR/labKey_Catalog_2016.12.20/RD_pilot_sample_metrics.txt',sep='\t', header = T)

	l_path_id = c()

	for (j in 1:length(l_lp_numbers)){
		
		index = which(data$iid == l_lp_numbers[j])

		if (length(index) == 0){

			l_path_id = c(l_path_id, '.')

		}
		else if (length(index) > 1){

			# we include in the analysis the most recent one, the latest one
			index2 = index[length(index)]
			l_path_id = c(l_path_id, as.character(data$path[index2]))

		}else{

			l_path_id = c(l_path_id, as.character(data$path[index]))

		}
	}# for

	df_output = cbind(lp = l_lp_numbers, path = l_path_id, status = l_affected)

	output_tsv = paste(panel_name, 'participants_and_familiars_enriched.tsv' , sep= '_')

	write.table(df_output, output_tsv , sep='\t', quote = F, row.names = F, col.names = T)

	if (i == 1){

		df_all = df_output
	}else if (length(l_gel_family_id) > 0) {

		df_all = rbind(df_all, df_output)
	}


