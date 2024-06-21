import pandas as pd
import requests
import json

#Get Marker DB data
file_path_markers = "marker_proteins.tsv"
markers = pd.read_csv(file_path_markers, sep='\t', header=None)
protein_biomarkers = markers[1].unique()
print(f'We obtain {len(protein_biomarkers)} unique protein biomarkers curated from Marker DB.')


# Get UniProt ID for a given protein name
def get_uniprot_id(protein_name):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={protein_name}&format=list&fields=accession"
    response = requests.get(url)
    if response.status_code == 200:
        result = response.text.strip().split('\n')
        return result[0] if result else "Not Found"
    else:
        return "Error"

# Call Fxn above to map protein names to UniProt IDs
# Uncomment these lines if it has not been run before
# protein_to_uniprot = {}
# for protein in protein_biomarkers:
#     uniprot_id = get_uniprot_id(protein)
#     protein_to_uniprot[protein] = uniprot_id
#write the protein to uniprot dictionary to a json file 
# with open('MarkerDB_Proteins_with_Uniprot_Mapping.json', 'w') as f:
#     json.dump(protein_to_uniprot, f)

#open the json file and load the protein to uniprot dictionary
with open('MarkerDB_Proteins_with_Uniprot_Mapping.json', 'r') as f:
    protein_markers = json.load(f)

#import the pre and post imputation proteins from each human sample 
dir = "output/"
samples = ["HS2_RE2", "HS3RE_sub1", "HS4", "HS5", "HS6", "HS7", "HS8", "HS9", "HS10", "HS11"]
file_suffix = [f"_pro_{s}.csv" for s in samples]
file_prefix = ["out", "raw"] # out is imputed protein list and raw is the non imputed protein list
imputed = [f"{dir}{file_prefix[0]}{file_suffix[i]}" for i in range(10)]
non_imputed = [f"{dir}{file_prefix[1]}{file_suffix[i]}" for i in range(10)]

# read in the protein list from each human sample
# calculate the intersection of the imputed/non imputed list and the values in the protein_markers dictionary
for i in range(10):
    pre = pd.read_csv(non_imputed[i])
    post = pd.read_csv(imputed[i])
    pre_proteins = pre['Uniprot'].values
    post_proteins = post['Uniprot'].values
    pre_proteins = [uniprot for uniprot in pre_proteins if uniprot in protein_markers.values()]
    post_proteins = [uniprot for uniprot in post_proteins if uniprot in protein_markers.values()]
    with open("Imputed_vs_Non_Imputed_Protein_Markers.csv", 'a') as f:
        f.write(f"{samples[i]}, wo_Imputation, {len(pre_proteins)}, {pre_proteins}\n")
        f.write(f"{samples[i]}, w_Imputation, {len(post_proteins)}, {post_proteins}\n")  
