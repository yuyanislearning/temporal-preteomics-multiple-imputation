import pandas as pd 
import numpy as np
import os

### Systematic way to get filepaths for the hl-data.out files (which contain A0 values) in each folder
filepath_common = 'proturn_output_'
strains = ['aj', 'balbc', 'c57', 'cej', 'dba', 'fvb']
exp = ['ctrl', 'iso'] 
file = 'hl-data.out'

# List filepaths to the ctrl data
ctrl_paths = [f'{filepath_common}{s}_{exp[0]}/{file}' for s in strains]
# List filepaths to the iso data
iso_paths = [f'{filepath_common}{s}_{exp[1]}/{file}' for s in strains]

### Function to process and mask data
def process_and_mask_data(file_paths, exp_type, mask_counts):
    for file in file_paths: 
        try:
            # Read in the file m
            data = pd.read_csv(file, sep='\t')
        except FileNotFoundError:
            print(f"File not found: {file}")
            continue
        except Exception as e:
            print(f"An error occurred while reading {file}: {e}")
            continue
        
        # hl-data.out files have ID, t, and A0 columns
        # Get the IDs of items with 7 items present
        non_missing_ids = data['ID'].value_counts().loc[lambda x: x == 7].index

        # Filter data to only include items with 7 items present
        data = data[data['ID'].isin(non_missing_ids)]
        # Save the filtered data to a new file
        filtered_file_path = f'{file[:-4]}_filtered.csv'
        try:
            data.to_csv(filtered_file_path, index=False)
        except Exception as e:
            print(f"An error occurred while writing {filtered_file_path}: {e}")
            continue

        # For each number of missing values in mask_counts
        for n_mask in mask_counts:
            masked_data = data.copy()
            # For each ID in the dataset, randomly mask n_mask A0 values
            for ID in masked_data['ID'].unique():
                # Get the indices of the A0 values for the current ID
                idx = masked_data.loc[masked_data['ID'] == ID].index
                # Randomly select n_mask indices to mask
                mask_idx = np.random.choice(idx, n_mask, replace=False)
                # Mask the A0 values at the selected indices
                masked_data.loc[mask_idx, 'A0'] = np.nan

            # Save the masked data to a new file
            masked_file_path = f'{file[:-4]}_masked_{n_mask}.csv'
            try:
                masked_data.to_csv(masked_file_path, index=False)
            except Exception as e:
                print(f"An error occurred while writing {masked_file_path}: {e}")

### Define the number of missing values to mask
mask_counts = [1, 2, 3, 4]

### Process CTRL Data
process_and_mask_data(ctrl_paths, 'ctrl', mask_counts)

### Process ISO Data
process_and_mask_data(iso_paths, 'iso', mask_counts)
