import pandas as pd
import numpy as np
import os 


dir = ["HS2_RE2", "HS3RE_sub1", "HS4", "HS5", "HS6", "HS7", "HS8", "HS9", "HS10", "HS11"]
file_names = [[f"ProTurn_Output_hs{i}_imp{j+1}.txt" for j in range(10) ] for i in range(2, 12)]
hl = 'hl.out'
sample_unique_proteins_dmi = {}
total_unique_proteins = []
out_dir = "output"
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

for d in dir: 
    #get unique proteins from the raw data
    raw = pd.read_csv(f"{d}/{hl}", sep="\t")
    #filter raw to include R2 greater than 0.8 or SS lesss than 0.05, and DP > 4 
    filtered_raw = raw[((raw['R2'] > 0.8) | (raw['SS'] < 0.05)) & (raw['DP'] > 4)]
    filtered_raw_proteins = filtered_raw['Uniprot'].unique()
    #array above to csv 
    filtered_raw_proteins = pd.DataFrame(filtered_raw_proteins, columns=['Uniprot'])
    filtered_raw_proteins.to_csv(f"{out_dir}/raw_pro_{d}.csv", sep=",", index=False)  
    
    print(f"{d} no-impute : {len(filtered_raw_proteins)}")
    
#get the unique proteins from the imputed data 
for imps, d in zip(file_names, dir): 
    imputed_proteins = []
    for f in imps:
        proturn_out = pd.read_csv(f"{d}/{f}", sep="\t") 
        #same filtering as raw but DP >= 2
        filtered_out = proturn_out[((proturn_out['R2'] > 0.8) | (proturn_out['SS'] < 0.05)) & (proturn_out['DP'] >= 2)]
        imputed_proteins.extend(filtered_out['Uniprot'].unique())
    print(f"{d} impute : {len(set(imputed_proteins))}")
    imputed_proteins = pd.DataFrame(set(imputed_proteins), columns=['Uniprot'])
    imputed_proteins.to_csv(f"{out_dir}/out_pro_{d}.csv", sep=",", index=False)
    