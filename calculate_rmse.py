import pandas as pd 
import numpy as np

# calculate sum of squared errors
def calculate_sse(true, imputed):
    sse = np.sum((true - imputed) ** 2)
    return sse

# file paths to the true A0 data, imputed A0 data, masked data 
imputation = 'dmi'
filepath_common = 'proturn_output_'
strains = ['aj', 'balbc', 'c57', 'cej', 'dba', 'fvb']
exp = ['ctrl', 'iso'] 
true_data_files = 'hl-data_filtered.csv'
mask = [1, 2, 3, 4] 
masked_data_files = 'hl-data_masked_'
imputed_data_files = f'{imputation}'
id_paths = 'IDs_'

folders = [filepath_common + s + '_' + e + '/' for s in strains for e in exp]

for m in mask: 
    sse_mask = []
    #observed values sum of squares
    ovss = []
    for f in folders: 
        # calculate number of missing values
        masked_data = pd.read_csv(f + masked_data_files + str(m) + '.csv')
        #find the index where values are nan
        masked_indices = masked_data.loc[masked_data['A0'].isnull(), 'A0'].index 
        
        # get ids and imputed data
        ids = pd.read_csv(f + id_paths + str(m) + ".csv")
        imputed_data = pd.read_csv(f + imputed_data_files + "_" + str(m) + '.csv')
        
        # get true data
        true_data = pd.read_csv(f + true_data_files)
        #get observed values sum of squares
        ovss.append(np.sum(true_data.loc[masked_indices, 'A0'] ** 2))
        
        #pivot true data to match imputed data format 
        true_data = true_data.pivot(index='ID', columns='t', values='A0').reset_index()
    
        # Need to calculate SSE for each of the n imputed matrices 
        sse_n_total = []
        
        for i in range(10): 
            imputed_subset = imputed_data.copy()
            imputed_subset = imputed_subset[imputed_subset['.imp'] == i + 1]
            
            # Add id column to imputed data
            imputed_subset.loc[:, 'x'] = ids['x'].values
            
            # print(imputed_subset)
            # print(true_data)
            true = np.array(true_data.drop(columns=['ID']))
            imputed = np.array(imputed_subset.drop(columns=['x', 'N', 'a', '.id', '.imp', 'Unnamed: 0']))
            
            # Calculate sum of squared errors for imputed matrix i 
            sse_n = calculate_sse(true, imputed)
            sse_n_total.append(sse_n)


        # Take average of sse_n_total
        sse_n_total = np.mean(sse_n_total)
        sse_mask.append(sse_n_total)
    # Print sum of squared errors for n masked datasets 
    print(f'{np.sqrt(sum(sse_mask))/sum(ovss)} for {m}')
