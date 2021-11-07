library(dplyr)
library(tibble)
library(rjson)
library(argparse)

options(stringsAsFactors = F)
options(header = T)

parser <- ArgumentParser(description="Associate ATC to a compound based on chemb ID match")
parser$add_argument("--compounds_json", type = "character", help = "")
parser$add_argument("--atc_compounds_json", type = "character", help = "")
parser$add_argument("--atc_csv", type = "character", help = "")
parser$add_argument("--out", type = "character", help = "")

args <- parser$parse_args()
compounds_json <- args$compounds_json
atc_compounds_json <- args$atc_compounds_json
atc_csv <- args$atc_csv
out <- args$out

###############
#### input ####
###############

# load atc table
atc_table <- read.csv(atc_csv)

# load json files
compounds_to_annotate <- fromJSON(file = compounds_json)
compounds_atc <- fromJSON(file = atc_compounds_json)

#################
### functions ###
#################

# convert json to dfs
convert_json_to_df <- function(list_json){
    
    n_entry <- length(list_json)
    names_feat <- names(list_json[[1]])
    df <- matrix(NA, ncol = length(names_feat), nrow = n_entry)
    colnames(df) <- names_feat
    df <- as.data.frame(df)

    for(id_entry in 1:n_entry){

        tmp <- list_json[[id_entry]]
        for(col_name in names_feat){
            if(length(tmp[[col_name]]) == 0){
                value_assing <- NA
            }else if(length(tmp[[col_name]]) == 1) {
               value_assing <- as.character(tmp[[col_name]])
            }else{
                value_assing <- paste0(as.character(tmp[[col_name]]),collapse = ',')
            }
            df[id_entry,col_name] <- value_assing
        }
    }
    # remove chembl ID NAs
    df <- df %>% filter(!is.na(molecule_chembl_id))

    return(df)
}


##########################

df_new <- convert_json_to_df(compounds_to_annotate)
df_atc <- convert_json_to_df(compounds_atc)

# match by chembl id
match_by_chembl <- inner_join(x = df_new, y = df_atc, by = 'molecule_chembl_id')
match_by_chembl <- match_by_chembl[!duplicated(match_by_chembl$compound.x),]

# attach ATC code(s) to compounds list 
df_new_with_atc <- data.frame(compound = match_by_chembl$compound.x, 
                                 chembl_id = match_by_chembl$molecule_chembl_id, 
                                 atc_code = NA, n_atc_code = NA,  
                                 atc_compound = match_by_chembl$compound.y)

for(idx_compound in 1:nrow(df_new_with_atc)){

    compound_name <- df_new_with_atc$atc_compound[idx_compound]
    atc_search <- atc_table$ATC[atc_table$substance %in% compound_name]
    n_out <- length(atc_search)

    df_new_with_atc$n_atc_code[idx_compound] <- n_out

    if(n_out>1){
        atc_search <- paste0(atc_search, collapse = '-')
    }
    df_new_with_atc$atc_code[idx_compound] <- atc_search
}
# save:
write.table(x = df_new_with_atc, file = sprintf("%s.csv",out), sep = ',', col.names = T, 
            row.names = F, quote = T)
