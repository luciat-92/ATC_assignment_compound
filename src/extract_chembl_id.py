# Data manipulation
import pandas as pd
import numpy as np
import os
import json
import math
import time
import argparse

# chembl client
from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.utils import utils

parser = argparse.ArgumentParser(description='Convert a list of compounds to dictiornaries with chembl ID')
parser.add_argument('--compound_list_csv', type = str, help='compound list in csv format, must contain "substance" column')
parser.add_argument('--output_file_name', type = str, help='path and filename to save json output')

args = parser.parse_args()
compound_list_csv = args.compound_list_csv
output_file_name = args.output_file_name

print("input file:" + compound_list_csv)
print("output file:" + output_file_name + ".json")

# import input drug lists
input_table = pd.read_csv(compound_list_csv, header=0)
compound_names = input_table['substance'].unique().tolist()

# maximum number of downloaded chembl ids to remove compound
MAX_N_META = 10000

print('input data loaded')

####################
#### functions #####
####################

# extract chembl id for a compound name
def get_list_compound_meta_first(compounds_list, MAX_N_META):
    
    molecule = new_client.molecule
    
    res = []
    count = 0
    
    for compound_name in compounds_list:
        
        print(count)
        
        chembl_search_compound = molecule.search(compound_name)
        
        try:
            len_search = len(chembl_search_compound)
            if len_search == 0 or len_search > MAX_N_META:
                res.append([])
            else:
                best = pd.DataFrame.from_dict([chembl_search_compound[0]])
                res.append(best.to_dict(orient='records')[0])
        except:
            res.append('error for ' + compound_name)
            print(compound_name)
        time.sleep(1)
        count += 1

    return res

# search twice chembl id if skipped first time
def get_chemblid(compounds_list, MAX_N_META):
    
    res = get_list_compound_meta_first(compounds_list, MAX_N_META)
    
    # check if some element are missing (only once), if so run again search
    missed_compounds_list = []
    idx_missed_compounds_list = []

    for idx_compound in range(len(compounds_list)):
        if res[idx_compound] == str('error for ' + compounds_list[idx_compound]):
            missed_compounds_list.append(compounds_list[idx_compound])
            idx_missed_compounds_list.append(idx_compound)

    if len(missed_compounds_list) > 0:
        print('second search for ' +  str(len(missed_compounds_list)) + ' missed compounds')
        res_first_missed = get_list_compound_meta_first(missed_compounds_list, MAX_N_META)
        
        for old_index, new_res in zip(idx_missed_compounds_list, res_first_missed): 
            res[old_index] = new_res
    
    return res

# create dictionaries with compound name and chembl id:
def create_dict_compound_chemblid(compounds_list, meta_search_output):
    compound_matched = []
    for idx_compound in range(len(compounds_list)):
        #print(idx_compound)
        if len(meta_search_output[idx_compound]) == 0 or meta_search_output[idx_compound] == str('error for ' + compounds_list[idx_compound]):
            compound_matched.append({'compound': compounds_list[idx_compound],
            'molecule_chembl_id': [], 
            'atc_chembl': []})
        else:
            compound_matched.append({'compound': compounds_list[idx_compound],
            'molecule_chembl_id': meta_search_output[idx_compound]['molecule_chembl_id'], 
            'atc_chembl': meta_search_output[idx_compound]['atc_classifications']})
    return compound_matched

#############################

## get chemblid 
# for mantra
chembl_compounds = get_chemblid(compound_names, MAX_N_META=MAX_N_META)
print('all chembl ids from compounds extracted')

# create dictionary
dictionary_compound_with_chembl = create_dict_compound_chemblid(compound_names, chembl_compounds)

# save (json format)
output_file_name = output_file_name + ".json"
out_file = open(output_file_name, "w") 
json.dump(dictionary_compound_with_chembl, out_file, indent = 6)   
out_file.close() 
