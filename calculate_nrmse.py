import pandas as pd 
import numpy as np

# calculate sum of squared errors
def calculate_sse(true, imputed):
    sse = np.sum((true - imputed) ** 2)
    return sse

# file paths to the true A0 data, imputed A0 data, masked data 
imputation = 'dmi_np'
n_imputations = 1
filepath_common = 'proturn_output_'
strains = ['aj', 'balbc', 'c57', 'cej', 'dba', 'fvb']
exp = ['ctrl', 'iso'] 
true_data_files = 'hl-data_filtered.csv'
mask = [1, 2, 3, 4, 5] 
masked_data_files = 'hl-data_masked_'
imputed_data_files = f'{imputation}'
id_paths = 'IDs_'
output_dir = 'rmse'
#create output directory
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

folders = [filepath_common + s + '_' + e + '/' for s in strains for e in exp]

for m in mask: 
    sse_mask = []
    #observed values sum of squares
    ovss = []
    for f in folders: 
        # get masked data 
        masked_data = pd.read_csv(f + masked_data_files + str(m) + '.csv')
        #find the index where values are nan
        masked_indices = masked_data.loc[masked_data['A0'].isnull(), 'A0'].index 
        
        # get imputed data
        imputed_data = pd.read_csv(f + imputed_data_files + "_" + str(m) + '.csv')
        
        # get true data
        true_data = pd.read_csv(f + true_data_files)
        #get observed values sum of squares
        ovss.append(np.sum(true_data.loc[masked_indices, 'A0'] ** 2))
        
        #pivot true data to match imputed data format 
        true_data = true_data.pivot(index='ID', columns='t', values='A0').reset_index()
    
        # Need to calculate SSE for each of the n imputed matrices 
        sse_n_total = []
        
        for i in range(n_imputations): 
            imputed_subset = imputed_data.copy()
            
            if imputation == 'dmi_np':
                imputed_subset = imputed_subset[imputed_subset['.imp'] == i+1]
            
            true = np.array(true_data.drop(columns=['ID']))
            #imputed = np.array(imputed_subset.drop(columns=['N', 'a', '.id', '.imp', 'Unnamed: 0']))
            if imputation == 'dmi_np':
                imputed = np.array(imputed_subset.drop(columns=['.id', '.imp', 'Unnamed: 0']))
            if imputation == 'si_mean':
                imputed = np.array(imputed_subset.drop(columns=['Unnamed: 0']))
      
            # Calculate sum of squared errors for imputed matrix i 
            sse_n = calculate_sse(true, imputed)
            sse_n_total.append(sse_n)

        # Take average of sse_n_total
        sse_n_total = np.mean(sse_n_total)
        sse_mask.append(sse_n_total)
        
    # Print nrmse for n masked datasets 
    print(f'{np.sum(sse_mask)/sum(ovss)} for {m}')
    
    #save nrmse to file in output_dir, save the nrmse, #masked values, and imputation method
    with open(f'{output_dir}/nrmse_{imputation}.txt', 'a') as f:
        f.write(f'{np.sum(sse_mask)/sum(ovss)} for {m}\n')
        f.close()
