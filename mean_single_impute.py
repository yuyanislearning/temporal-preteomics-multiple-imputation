import pandas as pd 
import numpy as np
import os
import sys
from sklearn.impute import SimpleImputer
import argparse

class ImputeAndCalculateRMSE:
    def __init__(self, n_masked, exp_type):
        self.n_masked = n_masked
        self.filepath_common = 'proturn_output_'
        self.strains = ['aj', 'balbc', 'c57', 'cej', 'dba', 'fvb']
        self.exp_type = exp_type
        self.mask = f'hl-data_masked_{n_masked}.csv'
        self.true = 'hl-data_filtered.csv'
        self.imputer_mean = SimpleImputer(strategy='mean')
        self.output_dir = 'rmse'
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def get_file_paths(self):
        file_paths_masked = [f'{self.filepath_common}{s}_{self.exp_type}/{self.mask}' for s in self.strains]
        file_paths_true = [f'{self.filepath_common}{s}_{self.exp_type}/{self.true}' for s in self.strains]
        return file_paths_masked, file_paths_true

    def calculate_rmse(self, file_paths_masked, file_paths_true):
        n_ctrl = []
        mse_sum_strain = []

        for mask, un_masked in zip(file_paths_masked, file_paths_true):
            try:
                data_mask = pd.read_csv(mask)
                #filter the dataframe to only include columns ID, t and A0
                data_mask = data_mask[['ID', 't', 'A0']]
        
                data_true = pd.read_csv(un_masked)
                #filter the dataframe to only include columns ID, t and A0
                data_true = data_true[['ID', 't', 'A0']]
                
            except FileNotFoundError:
                print(f"File not found: {mask} or {un_masked}")
                continue
            except Exception as e:
                print(f"An error occurred while reading {mask} or {un_masked}: {e}")
                continue
            
            n_ctrl.append(len(data_mask['A0']))
            mse_k = []
            print(f'imputing for {mask}...\n')
            for ID in data_mask['ID'].unique():
                idx = data_mask.loc[data_mask['ID'] == ID].index
                imputed = self.imputer_mean.fit_transform(data_mask[['A0']])[:,0]
                true = data_true['A0'].values 
                mse = np.sum((true - imputed) ** 2)
                mse_k.append(mse)

            mse_sum_strain.append(sum(mse_k))

        mse_sum_ctrl = sum(mse_sum_strain)
        n_ctrl_total = sum(n_ctrl) - 1 # -1 to account for the degrees of freedom
        rmse_ctrl = np.sqrt(mse_sum_ctrl / n_ctrl_total)
        return rmse_ctrl

    def save_rmse(self, rmse_value):
        with open(f'./{self.output_dir}/rmse_{self.exp_type}.csv', 'a') as f:
            f.write(f'{rmse_value}, {self.n_masked}, mean_single_impute\n') 

def main(n_masked, exp_type):
    imputer = ImputeAndCalculateRMSE(n_masked, exp_type)
    file_paths_masked, file_paths_true = imputer.get_file_paths()
    rmse_value = imputer.calculate_rmse(file_paths_masked, file_paths_true)
    imputer.save_rmse(rmse_value)
    print(f"RMSE for n_masked={n_masked}, exp_type={exp_type}: {rmse_value}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run mean imputation and calculate RMSE.')
    parser.add_argument('n_masked', type=int, choices=[1, 2, 3, 4], help='Number of masked values')
    parser.add_argument('exp_type', type=str, choices=['ctrl', 'iso'], help='Experiment type')
    args = parser.parse_args()
    main(args.n_masked, args.exp_type)
