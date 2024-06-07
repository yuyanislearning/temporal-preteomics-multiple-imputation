import pandas as pd 
import numpy as np

# file paths to the true A0 data, imputed A0 data, masked data 
imputation = 'mean_si'
filepath_common = 'proturn_output_'
strains = ['aj', 'balbc', 'c57', 'cej', 'dba', 'fvb']
exp = ['ctrl', 'iso'] 
true_data_files = 'hl-data_filtered.csv'
mask = [1, 2, 3, 4] 
masked_data_files = 'hl-data_masked_'
imputed_data_files = f'{imputation}'

folders = [filepath_common + s + '_' + e + '/' for s in strains for e in exp]

# calculate sum of squared errors
def calculate_sse(true, imputed):
    sse = np.sum((true - imputed) ** 2)
    return sse

for m in mask: 
    sse_mask = []
    #observed values sum of squares
    ovss = []
    for f in folders: 
        # calculate number of missing values
        masked_data = pd.read_csv(f + masked_data_files + str(m) + '.csv')
        #find the index where values are nan
        masked_indices = masked_data.loc[masked_data['A0'].isnull(), 'A0'].index 
        
        # get true data
        true_data = pd.read_csv(f + true_data_files)
        #get observed values sum of squares
        ovss.append(np.sum(true_data.loc[masked_indices, 'A0'] ** 2))
        
        # get imputed data
        imputed_data = pd.read_csv(f + imputed_data_files + "_" + str(m) + '.csv')
        imputed_copy = imputed_data.copy()
        
        true = np.array(true_data['A0'].values)
        imputed = np.array(imputed_copy['A0'].values)
        
        # Calculate sum of squared errors 
        sse_n = calculate_sse(true, imputed)
        sse_mask.append(sse_n)

    #Print nrmse for n masked datasets 
    print(f'{np.sum(sse_mask)/sum(ovss)} for {m}')
