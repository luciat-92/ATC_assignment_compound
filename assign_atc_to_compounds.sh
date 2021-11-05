#!/bin/bash

#SBATCH --job-name=atc
#SBATCH --partition=cpuq
#SBATCH --time=62:00:00
#SBATCH --nodes=1
#SBATCH --output=atc_%j.out
#SBATCH --error=atc_%j.err
#SBATCH --mem=10G

module load R

source /home/lucia.trastulla/miniconda3/etc/profile.d/conda.sh
conda activate py38

cd /home/lucia.trastulla/workspace/ATC_assignment_compound/

# extract ATC compounds chemblIDs
python src/extract_chembl_id.py --compound_list_csv data/covid_compound_names.csv --output_file_name data/compounds_covid_chemblid.json
# extract merged compounds lists (Aurora) chemblIDs
python src/extract_chembl_id.py --compound_list_csv data/atc_table_who_20211101_filtered.csv --output_file_name data/compounds_atc_who_chemblid.json

# combine lists
Rscript src/assign_atc_to_compounds_by_chembl.R \
	--compounds_json
	--atc_compounds_json
	--atc_csv
	--out


