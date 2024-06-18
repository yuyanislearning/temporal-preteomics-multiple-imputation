import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
df = pd.read_csv('nrmse/nrmse.csv')

# Set up the plot
plt.figure(figsize=(10, 6))
#set y -axis for plot from 0 to 0.35 
plt.ylim(0, 0.35)

#change the imputation_method column names in df to be more descriptive for the plot
df['imputation_method'] = df['imputation_method'].replace({'dmi_np': 'Data Multiple Imputation', 
                                                           'si_mean': 'Single Imputation with Mean',
                                                           'knn_30': 'KNN imputation with 30 neighbors'})

# Plot each imputation method
methods = df['imputation_method'].unique()
for method in methods:
    subset = df[df['imputation_method'] == method]
    plt.scatter(subset['n_masked'], subset['NRMSE'], label=method)
    plt.plot(subset['n_masked'], subset['NRMSE'])

# Add labels and title
plt.xlabel('Number of Masked Time Points for Each Peptide (n_masked)')
plt.ylabel('NRMSE')
#plt.title('NRMSE vs Number of Masked Values for Different Imputation Methods')
plt.legend(title='Imputation Method')
#save the plot as a pdf 
plt.savefig('nrmse/nrmse_plot.pdf')

# Show the plot
plt.show()


