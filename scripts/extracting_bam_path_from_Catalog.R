# we have the list of all the participants for which have assigned a panel `HA` or `ID`
# we retrieve from Catalog the bam path for them


# (and also some family-level data derived from the participant-level data)
# R CMD BATCH --no-save --no-restore /Users/kibanez/git/families_for_interpretation/merge_RD_pilot_participant_data.R 


load("/Users/kibanez/Documents/GEL_STR/labKey_Catalog_2016.12.20/RD_pilot_participant_data.RData")



# run `R CMD BATCH --no-save --no-restore /Users/kibanez/git/families_for_interpretation/merge_RD_pilot_participant_data.R ` the at the cluster

> pwd
> /genomes/scratch/kgarikano/labkKey-Catalog



`/home/ksmith/plink1_9/current`, e.g. `/home/ksmith/plink1_9/current/GEL_RD_pilot_sex_queries.csv`. 


`current` is just a symbolic link that points a folder in the same directory that is named after a date, e.g. it currently points to `/home/ksmith/plink1_9/2016_12_08/GEL_RD_pilot_sex_queries.csv`

That’s for the pilot - the corresponding files for the main programme live in or are copied to `/home/ksmith/plink1_9/GEL_RD_MP/current/`


The workflow still does generate a variable called `path` in object `pids` that stores the path to the delivery directory. The original source of the path is the file `~/Documents/sample_metrics/merge/current/RD_pilot_sample_metrics.txt`.
That is not in `/home/ksmith/plink1_9/current` and doesn’t particularly belong there - I’ll attach it here



==> al final no he accedido a Catalog, el fichero `RD_pilot_sample_metrics.txt` contiene todos esos paths y las muestras


####


setwd("/Users/kibanez/Documents/GEL_STR/labKey_Catalog_2016.12.20")

list_id = read.table('./list_plateKey_participants_ID.tsv',stringsAsFactors=F)$V1

list_ha = read.table('./list_plateKey_participants_HA.tsv',stringsAsFactors=F)$V1


list_id = list_id[!is.na(list_id)]

list_ha = list_ha[!is.na(list_ha)]


# in list_id and list_ha we have the platekeys or LP identifiers for the samples/participants to include in the `HA` and `ID` cohorts

data = read.csv('RD_pilot_sample_metrics.txt',sep='\t', header = T)

# ID

l_path_id = c()

for (i in 1:length(list_id)){
	
	index = which(data$iid == list_id[i])

	if (length(index) > 1){

		# we include in the analysis the most recent one, the latest one
		index2 = index[length(index)]
		l_path_id = c(l_path_id, as.character(data$path[index2]))

	}else{

		l_path_id = c(l_path_id, as.character(data$path[index]))

	}
}


# HA

l_path_ha = c()

for (i in 1:length(list_ha)){
	
	index = which(data$iid == list_ha[i])

	if (length(index) > 1){

		# we include in the analysis the most recent one, the latest one
		index2 = index[length(index)]
		l_path_ha = c(l_path_ha, as.character(data$path[index2]))

	}else{

		l_path_ha = c(l_path_ha, as.character(data$path[index]))

	}
}

df_id = cbind(lp = list_id, path = l_path_id)

df_ha = cbind(lp = list_ha, path = l_path_ha)

## we write in a simple plain tsv file, lp platekeys and their corresponding paths

write.table(df_id,'participants_and_familiars_intellectual_disability_january2017.tsv', quote=F, col.names=T, row.names=F, sep='\t')

write.table(df_ha,'participants_and_familiars_hereditary_ataxia_january2017.tsv', quote=F, col.names=T, row.names=F, sep='\t')



