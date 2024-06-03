#import necessary modules 
import pandas as pd 
import numpy as np
import os
from sklearn.impute import KNNImputer, SimpleImputer

# Call KNN imputer to impute missing values
imputer_knn = KNNImputer(n_neighbors=3)

# Number of items masked 
n_masked = 3

### Systematic way to get filepaths for the hl-data.out files (which contain A0 values) in each folder
filepath_common = 'proturn_output_'
strains = ['aj', 'balbc', 'c57', 'cej', 'dba', 'fvb']
exp = ['ctrl', 'iso'] 
mask = 'hl-data_masked.csv'
true = 'hl-data_filtered.csv'

# List filepaths to the ctrl data
ctrl_mask = [f'{filepath_common}{s}_{exp[0]}/{mask}' for s in strains]
ctrl_true = [f'{filepath_common}{s}_{exp[0]}/{true}' for s in strains]

# KNN imputer for ctrl data
# Initialize lists to accumulate results
n_masked_ctrl = []
mse_sum_strain = []

# Loop through each strain in control 
for mask, un_masked in zip(ctrl_mask, ctrl_true):
    # Read in the file 
    # data has ID, t, and A0 columns
    data_mask = pd.read_csv(mask)
    data_true = pd.read_csv(un_masked)
    
    # Calculate number of missing values in the mouse strain's dataset and add it to n_masked_ctrl 
    n_masked_ctrl.append(data_mask['A0'].isnull().sum())
    
    # Initialize list to store MSE for each peptide in the strain's dataset
    mse_k = []

    # Impute values for each peptide in the strain's dataset
    for ID in data_mask['ID'].unique():
        # Get the indices of the A0 values for the current ID
        idx = data_mask.loc[data_mask['ID'] == ID].index
        # KNN impute
        imputed = imputer_knn.fit_transform(data_mask.loc[idx, 'A0'].values.reshape(-1, 1)).flatten()
        true = data_true.loc[idx, 'A0'].values.flatten()
        
        # Calculate mean squared error for the peptide and append to mse_k
        mse = np.mean((true - imputed) ** 2)
        mse_k.append(mse)

    # Sum the mse values for each peptide stored in mse_k
    mse_sum_strain.append(sum(mse_k))

# Sum up mse values for each strain in ctrl data to get a total mse value for ctrl data
mse_sum_ctrl = sum(mse_sum_strain)

# Combined RMSE is square root of total mse divided by the total number of missing values
n_masked_ctrl_total = sum(n_masked_ctrl)
rmse_ctrl = np.sqrt(mse_sum_ctrl / n_masked_ctrl_total)
#write rmse_ctrl to a csv file with the rmse value and n_masked value as items
#ensure the file is able to be appended to
#make an rmse directory 
if not os.path.exists('rmse'):
    os.makedirs('rmse')
with open ('./rmse/rmse_ctrl.csv', 'a') as f:
    f.write(f'{rmse_ctrl}, {n_masked}, KNN\n') 

