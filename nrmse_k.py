import os
import pandas as pd
import numpy as np
import argparse

class ImputationAnalysis:
    def __init__(self, imputation):
        self.imputation = imputation
        self.filepath_common = 'proturn_output_'
        self.strains = ['aj'] #, 'balbc', 'c57', 'cej', 'dba', 'fvb']
        self.exp = ['ctrl', 'iso']
        self.true_data_files = 'hl.out'
        self.mask = [1, 2, 3, 4, 5]
        self.output_dir = 'nrmse_k'
        self.folders = [self.filepath_common + s + '_' + e + '/' for s in self.strains for e in self.exp]

        os.makedirs(self.output_dir, exist_ok=True)
    
    def calculate_nrmse(self, imputed_data, true_data):
        #filter by well fit peptides 
        imputed_data = imputed_data.loc[(imputed_data['R2'] >= 0.8) | (imputed_data['SS'] <= 0.05)]
        #imputed_data = imputed_data.loc[(imputed_data['k'] >= 0) & (imputed_data['k'] <= 100)] 
        true_data = true_data.loc[(true_data['R2'] >= 0.8)] 

        sse = 0 
        nf = 0 
        for prot in imputed_data['Uniprot'].unique():
                if self.imputation == 'imp': 
                    #imputed data for a specific protein
                    imp_prot = imputed_data.loc[imputed_data['Uniprot'] == prot]#.groupby(['ID', 'Uniprot']).agg({'k': 'mean'}).reset_index()
                    #true data for a specific protein
                    t = true_data.loc[true_data['Uniprot'] == prot]
                    # filter t so it only contains rows where ID is in imp['ID']
                    t = t.loc[t['ID'].isin(imp_prot['ID'])]
                    error = np.median(t['k']) - np.median(imp_prot['k'])
                    sse += error**2
                    nf += np.median(t['k'])**2
                else: 
                    #imputed data for a specific protein
                    imp_prot = imputed_data.loc[imputed_data['Uniprot'] == prot]
                    #true data for a specific protein
                    t = true_data.loc[true_data['Uniprot'] == prot]
                    # filter t so it only contains rows where ID is in imp['ID']
                    t = t.loc[t['ID'].isin(imp_prot['ID'])]
                    error = np.median(t['k']) - np.median(imp_prot['k'])
                    sse += error**2
                    nf += np.median(t['k'])**2
        
        return sse, nf 

    def run_analysis(self):
        for m in self.mask:
            sse_mask = []
            nf_mask = []
            for f in self.folders:
                try:
                    # Read in file
                    imputed_data = pd.read_csv(f'{f}ProTurn_{self.imputation}_{m}.txt', sep='\t')
                    true_data = pd.read_csv(f'{f}hl.out', sep='\t')
                except FileNotFoundError:
                    print(f"File not found: {f}")
                    continue

                sse,nf = self.calculate_nrmse(imputed_data, true_data) 

                sse_mask.append(sse)
                nf_mask.append(nf)

            nrmse = sum(sse_mask) / sum(nf_mask)
            print(f'{nrmse} for {m}')

            with open(f'{self.output_dir}/nrmse.csv', 'a') as f:
                f.write(f'{nrmse},{m},{self.imputation}\n')

def main():
    parser = argparse.ArgumentParser(description='Run imputation analysis.')
    parser.add_argument('imputation', type=str, choices=['imp', 'si', 'knn'], help='Type of imputation method: {imp, si, knn}')
    args = parser.parse_args()

    analysis = ImputationAnalysis(args.imputation)
    analysis.run_analysis()

if __name__ == "__main__":
    main()
