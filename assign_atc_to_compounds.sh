#!/bin/bash

#SBATCH --job-name=atc
#SBATCH --partition=cpuq
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --output=atc_%j.out
#SBATCH --error=atc_%j.err
#SBATCH --mem=10G

module load R

source /home/lucia.trastulla/miniconda3/etc/profile.d/conda.sh
conda activate py38

cd /home/lucia.trastulla/workspace/ATC_assignment_compound/

# extract merged compounds lists (Aurora) chemblIDs
file_new=/group/iorio/lucia/cMAP_shared/list_compounds_cMAP.csv
python src/extract_chembl_id.py --compound_list_csv ${file_new} --output_file_name data/compounds_cMAP_chemblid

# extract ATC compounds chemblIDs
file_atc=/home/lucia.trastulla/datasets/WHO_ATC_compounds/atc_table_who_20211101_filtered.csv
python src/extract_chembl_id.py --compound_list_csv ${file_atc} --output_file_name data/compounds_atc_who_chemblid

echo "chembl ID extracted"

# combine lists
Rscript src/assign_atc_to_compounds_by_chembl.R \
	--compounds_json data/compounds_cMAP_chemblid.json \
	--atc_compounds_json data/compounds_atc_who_chemblid.json \
	--atc_csv ${file_atc} \
	--out data/compounds_cMAP_with_atc

echo "ATC assigned to compound list"

cp data/compounds_cMAP_with_atc.csv /group/iorio/lucia/cMAP_shared/
