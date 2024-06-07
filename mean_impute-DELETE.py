import pandas as pd 
import numpy as np
import os
import sys
from sklearn.impute import SimpleImputer
import argparse

class ImputeData:
    def __init__(self, n_masked):
        self.n_masked = n_masked
        self.filepath_common = 'proturn_output_'
        self.strains = ['aj', 'balbc', 'c57', 'cej', 'dba', 'fvb']
        self.exp_type = ['ctrl', 'iso']
        self.mask = f'hl-data_masked_{n_masked}.csv'
        self.true = 'hl-data_filtered.csv'
        self.imputer_mean = SimpleImputer(strategy='mean')

    def get_file_paths(self):
        file_paths_masked = [f'{self.filepath_common}{s}_{e}/{self.mask}' for s in self.strains for e in self.exp_type]
        return file_paths_masked

    def impute_data(self, file_paths_masked):
        for mask in file_paths_masked:
            try:
                data_mask = pd.read_csv(mask)
                data_mask = data_mask[['ID', 't', 'A0']]
            except FileNotFoundError:
                print(f"File not found: {mask}")
                continue
            except Exception as e:
                print(f"An error occurred while reading {mask}: {e}")
                continue
            
            print(f'Imputing for {mask}...\n')
            
            for ID in data_mask['ID'].unique():
                idx = data_mask.loc[data_mask['ID'] == ID].index
                peptide_dat = data_mask.loc[idx, ['t', 'A0']]
                imputed = self.imputer_mean.fit_transform(peptide_dat[['A0']])[:,0]
                
                imputed_df = pd.DataFrame(imputed, columns=['A0']) 
                imputed_df['ID'] = [ID] * len(imputed)
                imputed_df['t'] = peptide_dat['t'].values
                
                imputed_file_path = f'{mask.split('/')[0]}/mean_si_{self.n_masked}.csv'
                if not os.path.isfile(imputed_file_path):
                    imputed_df.to_csv(imputed_file_path, index=False, mode='a', header=True)
                else:
                    imputed_df.to_csv(imputed_file_path, index=False, mode='a', header=False)

def main(n_masked):
    imputer = ImputeData(n_masked)
    file_paths_masked = imputer.get_file_paths()
    imputer.impute_data(file_paths_masked)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run mean imputation.')
    parser.add_argument('n_masked', type=int, choices=[1, 2, 3, 4], help='Number of masked values')
    args = parser.parse_args()
    main(args.n_masked)
