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
        self.imputed_dir = 'mean_imputation_masked_data'
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        if not os.path.exists(self.imputed_dir):
            os.makedirs(self.imputed_dir)

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
                #get t and A0 vales for the idx values
                peptide_dat = data_mask.loc[idx, ['t', 'A0']]
                imputed = self.imputer_mean.fit_transform(peptide_dat[['A0']])[:,0]
                
                #write the imputed array to a csv file (appending to the csv each loop)
                # add corresponding columns for ID and time {0,1,3,5,7,10,14} as well 
                imputed_df = pd.DataFrame(imputed, columns=['A0']) 
                imputed_df['ID'] = [ID] * len(imputed)
                imputed_df['t'] = peptide_dat['t'].values
                
                # Check if file exists to determine if header should be included
                imputed_file_path = f'./{self.imputed_dir}/{mask.split('/')[0]}_mean_imputed_{self.n_masked}.csv'
                if not os.path.isfile(imputed_file_path):
                    imputed_df.to_csv(imputed_file_path, index=False, mode='a', header=True)
                else:
                    imputed_df.to_csv(imputed_file_path, index=False, mode='a', header=False)
                
                #calculate squared error for the imputed values for each peptide, these will be summed up later            
                true = data_true.loc[idx, ['A0']].values 
                mse = np.sum((true - imputed) ** 2)
                mse_k.append(mse)
            #sum all the squared errors for each peptide in the strain 
            mse_sum_strain.append(sum(mse_k))
            
        #sum all the squared errors for each strain in the experiment group 
        mse_sum_ctrl = sum(mse_sum_strain)
        n_ctrl_total = sum(n_ctrl) - 1 # -1 to account for the degrees of freedom
        #calculate the root mean squared error for the experiment group
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
