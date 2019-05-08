################################################################################################################################################################################################################################################
# 2017.01.17
# csv's provided by David Montaner (/Users/kibanez/Documents/GEL_STR/labKey_Catalog_2016.12.20/rd_pilot_labkey2catalog_2016_12_20/*.csv)

# panels.csv contains `gelFamilyId`, `gelId` (to extract) and `panel`. We will extract all those participants (`gelFamilyId` and `gelId`) with `ID` and `HA` as panel assigned
# samples.csv contains `gelID` and `plateKey`. With the list of `gelId` retrieved from panels.csv, we'll extract all those `plateKey` corresponding to `gelId`
# registration.csv contains `gelId` and `gelFamilyId`, in order to extract all those `gelID` that are members within a family given a `gelFamilyId` 

# we are not interested only in the participants, but also in the rest of their familiars. With the `gelFamilyId` we will extract all the members of the family as well

################################################################################################################################################################################################################################################

setwd("/Users/kibanez/Documents/GEL_STR/labKey_Catalog_2016.12.20/rd_pilot_labkey2catalog_2016_12_20/")

df_panels = read.csv('panels.csv', sep=',', header=T)

df_samples = read.csv('samples.csv', sep=',', header=T)

df_registration = read.csv('registration.csv', sep=',', header=T)

## 1) We first extract the `gelFamilyId` and `gelId` for those participants that have been assigned `HA` or `ID`

# ID 
index_ID = which(df_panels$panel == "Intellectual disability")

l_gel_family_id_ID = as.character(df_panels$gelFamilyId[index_ID])

# HA
index_HA = which(df_panels$panel == "Hereditary ataxia")

l_gel_family_id_HA = as.character(df_panels$gelFamilyId[index_HA])


## 2) Extract all the `gelId`-s included in each `gelFamilyId`

# ID
#index_gel_id_ID = which(l_gel_family_id_ID %in% as.character(df_registration$gelFamilyId))

l_gel_id_ID = c()
l_affected_ID = c()

for (i in 1:length(l_gel_family_id_ID)){

	index = which(as.character(df_registration$gelFamilyId) %in% l_gel_family_id_ID[i])

	l_gel_id_ID = c(l_gel_id_ID, as.character(df_registration$gelId[index]))
	l_affected_ID = c(l_affected_ID, as.character(df_registration$disease_status[index]))

}


# HA
#index_gel_id_HA = which(l_gel_family_id_HA %in% as.character(df_registration$gelFamilyId))

l_gel_id_HA = c()
l_affected_HA = c()

for (i in 1:length(l_gel_family_id_HA)){

	index = which(as.character(df_registration$gelFamilyId) %in% l_gel_family_id_HA[i])

	l_gel_id_HA = c(l_gel_id_HA, as.character(df_registration$gelId[index]))
	l_affected_HA = c(l_affected_HA, as.character(df_registration$disease_status[index]))

}


## 3) Extract from df_samples the LP numbers

l_lp_numbers_ID = as.character(df_samples$plateKey[pmatch(l_gel_id_ID, df_samples$gelID)])


l_lp_numbers_HA = as.character(df_samples$plateKey[pmatch(l_gel_id_HA, df_samples$gelID)])


write.table(l_lp_numbers_ID, '../list_plateKey_participants_ID.tsv', quote = F, row.names = F, col.names = F)

write.table(l_lp_numbers_HA, '../list_plateKey_participants_HA.tsv', quote = F, row.names = F, col.names = F)


df_status_ID = data.frame(lp_id = l_lp_numbers_ID, status = l_affected_ID)
df_status_HA = data.frame(lp_id = l_lp_numbers_HA, status = l_affected_HA)


write.table(df_status_ID, '../list_plateKey_status_participants_ID.tsv', sep='\t', quote = F, row.names = F, col.names = T)
write.table(df_status_HA, '../list_plateKey_status_participants_HA.tsv', sep='\t', quote = F, row.names = F, col.names = T)


# we incorporate now the path info as well


setwd("/Users/kibanez/Documents/GEL_STR/labKey_Catalog_2016.12.20")


# in list_id and list_ha we have the platekeys or LP identifiers for the samples/participants to include in the `HA` and `ID` cohorts

data = read.csv('RD_pilot_sample_metrics.txt',sep='\t', header = T)

# ID

l_path_id = c()

for (i in 1:length(l_lp_numbers_ID)){
	
	index = which(data$iid == l_lp_numbers_ID[i])

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
}


# HA

l_path_ha = c()

for (i in 1:length(l_lp_numbers_HA)){
	
	index = which(data$iid == l_lp_numbers_HA[i])

	if (length(index) == 0){

		l_path_ha = c(l_path_ha, '.')

	}else if (length(index) > 1){

		# we include in the analysis the most recent one, the latest one
		index2 = index[length(index)]
		l_path_ha = c(l_path_ha, as.character(data$path[index2]))

	}else{

		l_path_ha = c(l_path_ha, as.character(data$path[index]))

	}
}

df_id = cbind(lp = l_lp_numbers_ID, path = l_path_id, status = l_affected_ID)

df_ha = cbind(lp = l_lp_numbers_HA, path = l_path_ha, status = l_affected_HA)

## we write in a simple plain tsv file, lp platekeys and their corresponding paths

write.table(df_id,'participants_and_familiars_intellectual_disability_january2017_enriched.tsv', quote=F, col.names=T, row.names=F, sep='\t')

write.table(df_ha,'participants_and_familiars_hereditary_ataxia_january2017_enriched.tsv', quote=F, col.names=T, row.names=F, sep='\t')

