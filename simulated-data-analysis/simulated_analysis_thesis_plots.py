#!/bin/env python3

import pandas as pd
import glob
import matplotlib as plt
from plotnine import *
from plotnine.data import *
from scipy.stats import t
import scipy
import numpy as np

############################# Read in Bayes & GBLUP Metrics ##########################################
######################################################################################################
dataframes = []

for ancestry in {"CEU", "CHB", "YRI", "TSI"}:
    for bayes in {"bayes_r", "gblup"}:
        for h2 in {"0.4", "0.8"}:
            for ncs in {"100", "10k"}:
                for p in {"2"}:
                    bayes_directory = f"/home/duffme/Output_Error_Files/target_comparison_linear_phenotype_neural_networks/{ancestry}"
                    # Define the file pattern
                    file_pattern = f"{bayes}_metrics_h2{h2}_ncs{ncs}_pvalue{p}.txt"
                    # Use glob to get a list of file paths that match the pattern
                    file_paths = glob.glob(f"{bayes_directory}/**/{file_pattern}", recursive=True)
                    # Loop through each file path and read it as a DataFrame
                    for file_path in file_paths:
                        df = pd.read_csv(file_path)
                        df['method'] = f'{bayes}'
                        df['ancestry'] = f'{ancestry}'
                        dataframes.append(df)
                    
combined_bayes = pd.concat(dataframes, ignore_index=True)

combined_bayes.head()

count_per_group = combined_bayes.groupby(['h2', 'num_causal_snps', 'ancestry', 'method']).size().reset_index(name='row_count')

print(count_per_group)

################################ Read in MLP Metrics #################################################
######################################################################################################

dataframes = []
                    
# for sim in range(1,110):
#     for h2 in {"0.4", "0.8"}:
#         for ncs in {"100", "10k"}:
#             for p in {"4"}:
#                 nn_directory = f"/home/duffme/Output_Error_Files/target_MLP_linear_phenotype_neural_networks"
#                 # Define the file pattern
#                 file_pattern = f"MLP_linear_pheno_sim{sim}_ncs{ncs}_h2{h2}_pvalue_{p}_bs64_e"
#                 # Use glob to get a list of file paths that match the pattern
#                 file_paths = glob.glob(f"{nn_directory}/**/{file_pattern}"+ "*", recursive=True)
#                 target_snps_paths = [path for path in file_paths if "target_snps" in path]
#                 # Loop through each file path and read it as a DataFrame
#                 for file_path in target_snps_paths:
#                     # Read the CSV file into a DataFrame
#                     df = pd.read_csv(file_path)
#                     df['sim']=f'{sim}'
#                     df['method'] = 'MLP'
#                     # Append the DataFrame to the list
#                     dataframes.append(df)


for ancestry in {"CEU", "CHB", "YRI", "TSI"}:                  
  for sim in range(1,110):
    for h2 in {"0.4", "0.8"}:
      for ncs in {"100", "10k"}:
        for p in {"4"}:
          nn_directory = f"/scratch/duffme_gensim/Output_Error_Files/target_MLP_linear_phenotype_neural_networks/{ancestry}/sim_{sim}"
          # Define the file pattern
          file_pattern = f"MLP_2_linear_pheno_sim{sim}_ncs{ncs}_h2{h2}_pvalue_{p}_bs64_e"
          # Use glob to get a list of file paths that match the pattern
          file_paths = glob.glob(f"{nn_directory}/**/{file_pattern}"+ "*", recursive=True)
          target_snps_paths = [path for path in file_paths if "target_snps" in path]
          # Loop through each file path and read it as a DataFrame
          for file_path in target_snps_paths:
            # Read the CSV file into a DataFrame
            df = pd.read_csv(file_path)
            df['sim']=f'{sim}'
            df['method'] = 'MLP'
            # Append the DataFrame to the list
            dataframes.append(df)

combined_MLP_temp = pd.concat(dataframes, ignore_index=True)

combined_MLP_temp = combined_MLP_temp.dropna()

grouped = combined_MLP_temp.groupby(['sim', 'h2', 'num_causal_snps', 'Ancestry', 'p-value'])

combined_MLP = grouped.apply(lambda x: x.sort_values(by='Validation_Correlation', ascending=False).head(1))

combined_MLP.reset_index(drop=True, inplace=True)

count_per_group = combined_MLP.groupby(['h2', 'num_causal_snps', 'Ancestry', 'p-value']).size().reset_index(name='row_count')

print(count_per_group)

################################### Read in TL Metrics ###############################################
######################################################################################################

dataframes = []

for ancestry in {"CEU", "CHB", "YRI", "TSI"}:                    
    for sim in range(1,110):
        for h2 in {"0.4", "0.8"}:
            for ncs in {"100", "10k"}:
                for p in {"4"}:
                    nn_directory = f"/home/duffme/Output_Error_Files/target_TL_MLP_linear_phenotype_neural_networks/{ancestry}/sim_{sim}/"
                    # Define the file pattern
                    file_pattern_1 = f"{ancestry}_TL_MLP_linear_pheno_sim{sim}_ncs{ncs}_h2{h2}_pvalue_{p}_bs64"
                    # Use glob to get a list of file paths that match the pattern
                    file_paths_1 = glob.glob(f"{nn_directory}/**/{file_pattern_1}"+ "*", recursive=True)
                    # Loop through each file path and read it as a DataFrame
                    if file_paths_1:
                        for file_path in file_paths_1:
                            # Read the CSV file into a DataFrame
                            df=pd.read_csv(file_path)
                            df['sim']=f'{sim}'
                            df['method']="TL_MLP_" + df['Approach']
                            # Append the DataFrame to the list
                            dataframes.append(df)
                        
combined_TL_temp = pd.concat(dataframes, ignore_index=True)

grouped = combined_TL_temp.groupby(['sim', 'h2', 'num_causal_snps', 'Ancestry', 'p-value', 'Approach'])

combined_TL = grouped.apply(lambda x: x.sort_values(by='Validation_Correlation', ascending=False).head(1))

combined_TL.reset_index(drop=True, inplace=True)

count_per_group = combined_TL.groupby(['h2', 'num_causal_snps', 'Ancestry', 'p-value', 'Approach']).size().reset_index(name='row_count')

print(count_per_group)

# dataframes = []
#                         
# for ancestry in {"CEU", "CHB", "YRI", "TSI"}:                    
#     for sim in range(1,110):
#         for h2 in {"0.4", "0.8"}:
#             for ncs in {"100", "10k"}:
#                 for p in {"4", "8"}:
#                     nn_directory = f"/home/duffme/Output_Error_Files/target_TL_MLP_linear_phenotype_neural_networks/{ancestry}/sim_{sim}/"
#                     # Define the file pattern
#                     file_pattern_2 = f"{ancestry}_TL_2_MLP_linear_pheno_sim{sim}_ncs{ncs}_h2{h2}_pvalue_{p}_bs16000"
#                     # Use glob to get a list of file paths that match the pattern
#                     file_paths_2 = glob.glob(f"{nn_directory}/**/{file_pattern_2}"+ "*", recursive=True)
#                     # Loop through each file path and read it as a DataFrame
#                     if file_paths_2:
#                         file_path_2=file_paths_2[0]
#                         # Read the CSV file into a DataFrame
#                         df = pd.read_csv(file_path_2)
#                         df['sim']=f'{sim}'
#                         # Append the DataFrame to the list
#                         dataframes.append(df)
# 
# combined_TL_2 = pd.concat(dataframes, ignore_index=True)
# 
# combined_TL_2.head()
# 
# count_per_group = combined_TL_2.groupby(['h2', 'num_causal_snps', 'Ancestry', 'p-value']).size().reset_index(name='row_count')
# 
# print(count_per_group)
# 
# merged_df = pd.merge(combined_TL_1, combined_TL_2, on=['Ancestry', 'Approach', 'h2', 'num_causal_snps', 'Ancestry', 'p-value', 'sim'], how='outer',suffixes=('_1', '_2'))
# 
# merged_df.head()
# 
# filtered_rows = merged_df[(merged_df['Validation_Correlation_1'] < 0) & (merged_df['Approach'] == "Fine_Tune")]
# 
# print(filtered_rows)
# 
# # Iterate through each row and select columns with '_1' suffix if val_correlation_1 > val_correlation_2
# selected_rows = []
# 
# for index, row in merged_df.iterrows():
#     if row['Validation_Correlation_1'] > row['Validation_Correlation_2']:
#         selected_columns = ['Ancestry', 'Approach', 'Train_Correlation_1', 'Validation_Correlation_1', 'Test_Correlation_1', 'Train_MSE_1', 'Validation_MSE_1', 'Test_MSE_1', 'Train_R2_1', 'Validation_R2_1', 'Test_R2_1', 'time_1', 'h2', 'num_causal_snps', 'p-value', 'sim']
#     else:
#         selected_columns = ['Ancestry', 'Approach', 'Train_Correlation_2', 'Validation_Correlation_2', 'Test_Correlation_2', 'Train_MSE_2', 'Validation_MSE_2', 'Test_MSE_2', 'Train_R2_2', 'Validation_R2_2', 'Test_R2_2', 'time_2', 'h2', 'num_causal_snps', 'p-value', 'sim']
#     selected_row = row.loc[selected_columns]
#     row_dataframe = pd.DataFrame(selected_row).transpose()
#     row_dataframe.columns = ['Ancestry', 'method', 'Train_Correlation', 'Validation_Correlation', 'Test_Correlation', 'Train_MSE', 'Validation_MSE', 'Test_MSE', 'Train_R2', 'Validation_R2', 'Test_R2', 'time', 'h2', 'num_causal_snps',  'p-value', 'sim']
#     selected_rows.append(row_dataframe)
# 
# combined_TL_best = pd.concat(selected_rows, ignore_index=True)
# print(combined_TL_best.columns)

combined_TL_best = combined_TL[['Ancestry', 'method', 'Train_Correlation', 'Validation_Correlation', 'Test_Correlation', 'Train_MSE', 'Validation_MSE', 'Test_MSE', 'Train_R2', 'Validation_R2', 'Test_R2', 'h2', 'num_causal_snps',  'p-value', 'time', 'sim']]
combined_TL_best['Ancestry'] = combined_TL_best['Ancestry'].astype('category')
combined_TL_best['h2'] = combined_TL_best['h2'].astype('category')
combined_TL_best['num_causal_snps'] = combined_TL_best['num_causal_snps'].astype('category')
combined_TL_best['p-value'] = combined_TL_best['p-value'].astype('category')
# count_per_group = combined_TL_best.groupby(['h2', 'num_causal_snps', 'Ancestry', 'p-value']).size().reset_index(name='row_count')
# 
# print(count_per_group)

################################### Read in CNN Metrics ###############################################
######################################################################################################

dataframes = []

for ancestry in {"CEU", "CHB", "YRI", "TSI"}:                    
    for sim in range(1,20):
        for h2 in {"0.4", "0.8"}:
            for ncs in {"100", "10k"}:
                for p in {"4"}:
                    nn_directory = f"/scratch/duffme_gensim/Output_Error_Files/target_CNN_linear_phenotype_neural_networks/{ancestry}/sim_{sim}/"
                    # Define the file pattern
                    file_pattern_1 = f"{ancestry}_CNN_linear_pheno_sim{sim}_ncs{ncs}_h2{h2}_pvalue_{p}_bs64"
                    # Use glob to get a list of file paths that match the pattern
                    file_paths_1 = glob.glob(f"{nn_directory}/**/{file_pattern_1}"+ "*nn_results.csv", recursive=True)
                    # Loop through each file path and read it as a DataFrame
                    if file_paths_1:
                        for file_path in file_paths_1:
                            # Read the CSV file into a DataFrame
                            df=pd.read_csv(file_path)
                            df['sim']=f'{sim}'
                            df['method']="CNN"
                            df['Ancestry']=f'{ancestry}'
                            # Append the DataFrame to the list
                            dataframes.append(df)
                        
combined_CNN = pd.concat(dataframes, ignore_index=True)

combined_CNN.head()

count_per_group = combined_CNN.groupby(['h2', 'num_causal_snps', 'Ancestry', 'method']).size().reset_index(name='row_count')

print(count_per_group)


################################### Read in CNN TL Metrics ###############################################
######################################################################################################

dataframes = []

for ancestry in {"CEU", "CHB", "YRI", "TSI"}:                    
    for sim in range(1,20):
        for h2 in {"0.4", "0.8"}:
            for ncs in {"100", "10k"}:
                for p in {"4"}:
                    nn_directory = f"/scratch/duffme_gensim/Output_Error_Files/target_TL_CNN_linear_phenotype_neural_networks/{ancestry}/sim_{sim}/"
                    # Define the file pattern
                    file_pattern_1 = f"{ancestry}_TL_CNN_linear_pheno_sim{sim}_ncs{ncs}_h2{h2}_pvalue_{p}_bs64"
                    # Use glob to get a list of file paths that match the pattern
                    file_paths_1 = glob.glob(f"{nn_directory}/**/{file_pattern_1}"+ "*_nn_results.csv", recursive=True)
                    # Loop through each file path and read it as a DataFrame
                    if file_paths_1:
                        for file_path in file_paths_1:
                            # Read the CSV file into a DataFrame
                            df=pd.read_csv(file_path)
                            df['sim']=f'{sim}'
                            df['method']="TL_CNN_"+ df['Approach']
                            # Append the DataFrame to the list
                            dataframes.append(df)
                        
combined_TL_CNN = pd.concat(dataframes, ignore_index=True)

combined_TL_CNN.head()

count_per_group = combined_TL_CNN.groupby(['h2', 'num_causal_snps', 'Ancestry', 'method']).size().reset_index(name='row_count')

print(count_per_group)
   
#combined_TL_CNN.loc[(combined_TL_CNN["Ancestry"] == "YRI") & (combined_TL_CNN["h2"] == 0.4) & (combined_TL_CNN["num_causal_snps"] == "10k")]


################################### Combine Metric Dataframes ###############################################
############################################################################################################

combined_bayes.columns = ["Train_MSE", "Train_Correlation", "Train_R2", "Test_MSE", "Test_Correlation", "Test_R2", "sim", "h2", "num_causal_snps", "p-value", "time", "method", "Ancestry"]

subset_combined_bayes = combined_bayes[['Train_Correlation', 'Test_Correlation',
       'Train_MSE','Test_MSE', 'Train_R2',
       'Test_R2', 'sim', 'h2', 'num_causal_snps', 'p-value', 'method', 'Ancestry']]
       
subset_combined_MLP = combined_MLP[['Ancestry', 'Train_Correlation', 'Test_Correlation',
       'Train_MSE','Test_MSE', 'Train_R2',
       'Test_R2', 'h2', 'num_causal_snps', 'p-value', 'sim', 'method']]
       
subset_combined_TL_best = combined_TL_best[['Ancestry', 'method', 'Train_Correlation', 'Test_Correlation',
       'Train_MSE','Test_MSE', 'Train_R2',
       'Test_R2', 'h2', 'num_causal_snps', 'p-value', 'sim']]
       
subset_combined_CNN = combined_CNN[['Ancestry', 'method', 'Train_Correlation', 'Test_Correlation',
       'Train_MSE','Test_MSE', 'Train_R2',
       'Test_R2', 'h2', 'num_causal_snps', 'p-value', 'sim']]
       
subset_combined_TL_CNN = combined_TL_CNN[['Ancestry', 'method', 'Train_Correlation', 'Test_Correlation',
       'Train_MSE','Test_MSE', 'Train_R2',
       'Test_R2', 'h2', 'num_causal_snps', 'p-value', 'sim']]

# combined_bayes['Train_MSE'] = pd.to_numeric(combined_bayes['Train_MSE'])
# combined_bayes['Train_Correlation'] = pd.to_numeric(combined_bayes['Train_Correlation'])
# combined_bayes['Train_R2'] = pd.to_numeric(combined_bayes['Train_R2'])
# combined_bayes['Test_MSE'] = pd.to_numeric(combined_bayes['Test_MSE'])
# combined_bayes['Test_Correlation'] = pd.to_numeric(combined_bayes['Test_Correlation'])
# combined_bayes['Test_R2'] = pd.to_numeric(combined_bayes['Test_R2'])
# combined_bayes['h2'] = pd.to_numeric(combined_bayes['h2'])
# combined_bayes['p-value'] = pd.to_numeric(combined_bayes['p-value'])

linear_target_sets = pd.concat([combined_bayes, subset_combined_MLP, subset_combined_TL_best, subset_combined_CNN, subset_combined_TL_CNN])


#################################################### READ IN NON-LINEAR PHENOTYPE RESULTS ###############################################################################################
#########################################################################################################################################################################################
#########################################################################################################################################################################################


############################# Read in Bayes & GBLUP Metrics ##########################################
######################################################################################################
dataframes = []

for ancestry in {"CEU", "CHB", "YRI", "TSI"}:
    for bayes in {"bayes_r", "gblup"}:
        for h2 in {"0.05", "0.1"}:
            for ncs in {"1k", "10k"}:
                for p in {"2"}:
                    bayes_directory = f"/home/duffme/Output_Error_Files/target_comparison_linear_phenotype_neural_networks/{ancestry}"
                    # Define the file pattern
                    file_pattern = f"non_linear_{bayes}_metrics_h2{h2}_ncs{ncs}_pvalue{p}.txt"
                    # Use glob to get a list of file paths that match the pattern
                    file_paths = glob.glob(f"{bayes_directory}/**/{file_pattern}", recursive=True)
                    # Loop through each file path and read it as a DataFrame
                    for file_path in file_paths:
                        df = pd.read_csv(file_path)
                        df['method'] = f'{bayes}'
                        df['ancestry'] = f'{ancestry}'
                        dataframes.append(df)
                    
combined_bayes = pd.concat(dataframes, ignore_index=True)

combined_bayes.head()

count_per_group = combined_bayes.groupby(['h2', 'num_causal_snps', 'ancestry', 'method']).size().reset_index(name='row_count')

print(count_per_group)

################################ Read in MLP Metrics #################################################
######################################################################################################

dataframes = []

for ancestry in {"CEU", "CHB", "YRI", "TSI"}:                 
  for sim in range(1,110):
    for h2 in {"0.05", "0.1"}:
      for ncs in {"1k", "10k"}:
        for p in {"4"}:
          nn_directory = f"/scratch/duffme_gensim/Output_Error_Files/target_MLP_non_linear_phenotype_neural_networks/{ancestry}"
          # Define the file pattern
          file_pattern = f"MLP_non_linear_pheno_sim{sim}_ncs{ncs}_h2{h2}_pvalue_{p}_bs64_e"
          # Use glob to get a list of file paths that match the pattern
          file_paths = glob.glob(f"{nn_directory}/**/{file_pattern}"+ "*", recursive=True)
          # Loop through each file path and read it as a DataFrame
          if file_paths:
            file_path=file_paths[0]
            # Read the CSV file into a DataFrame
            df = pd.read_csv(file_path)
            df['sim']=f'{sim}'
            df['method'] = 'MLP'
            # Append the DataFrame to the list
            dataframes.append(df)

combined_MLP = pd.concat(dataframes, ignore_index=True)

combined_MLP.head()

combined_MLP=combined_MLP.dropna()

count_per_group = combined_MLP.groupby(['h2', 'num_causal_snps', 'Ancestry', 'p-value']).size().reset_index(name='row_count')

print(count_per_group)

################################### Read in TL Metrics ###############################################
######################################################################################################

dataframes = []

for ancestry in {"CEU", "CHB", "YRI", "TSI"}:                    
    for sim in range(1,110):
        for h2 in {"0.05", "0.1"}:
            for ncs in {"1k", "10k"}:
                for p in {"4"}:
                    nn_directory = f"/scratch/duffme_gensim/Output_Error_Files/target_TL_MLP_non_linear_phenotype_neural_networks/{ancestry}/sim_{sim}/"
                    # Define the file pattern
                    file_pattern_1 = f"{ancestry}_TL_MLP_non_linear_pheno_sim{sim}_ncs{ncs}_h2{h2}_pvalue_{p}_bs64"
                    # Use glob to get a list of file paths that match the pattern
                    file_paths_1 = glob.glob(f"{nn_directory}/**/{file_pattern_1}"+ "*", recursive=True)
                    # Loop through each file path and read it as a DataFrame
                    if file_paths_1:
                        file_path_1=file_paths_1[0]
                        # Read the CSV file into a DataFrame
                        df=pd.read_csv(file_path_1)
                        df['sim']=f'{sim}'
                        df['method']="TL_MLP_"+ df['Approach']
                        # Append the DataFrame to the list
                        dataframes.append(df)
                        
combined_TL_1 = pd.concat(dataframes, ignore_index=True)

combined_TL_1.head()

count_per_group = combined_TL_1.groupby(['h2', 'num_causal_snps', 'Ancestry', 'p-value']).size().reset_index(name='row_count')

print(count_per_group)


# dataframes = []
#                         
# for ancestry in {"CEU", "CHB", "YRI", "TSI"}:                    
#     for sim in range(1,110):
#         for h2 in {"0.05", "0.1"}:
#             for ncs in {"1k", "10k"}:
#                 for p in {"4"}:
#                     nn_directory = f"/home/duffme/Output_Error_Files/target_TL_MLP_non_linear_phenotype_neural_networks/{ancestry}/sim_{sim}/"
#                     # Define the file pattern
#                     file_pattern_2 = f"{ancestry}_TL_2_MLP_non_linear_pheno_sim{sim}_ncs{ncs}_h2{h2}_pvalue_{p}_bs16000"
#                     # Use glob to get a list of file paths that match the pattern
#                     file_paths_2 = glob.glob(f"{nn_directory}/**/{file_pattern_2}"+ "*", recursive=True)
#                     # Loop through each file path and read it as a DataFrame
#                     if file_paths_2:
#                         file_path_2=file_paths_2[0]
#                         # Read the CSV file into a DataFrame
#                         df = pd.read_csv(file_path_2)
#                         df['sim']=f'{sim}'
#                         # Append the DataFrame to the list
#                         dataframes.append(df)
# 
# combined_TL_2 = pd.concat(dataframes, ignore_index=True)
# 
# combined_TL_2.head()
# 
# count_per_group = combined_TL_2.groupby(['h2', 'num_causal_snps', 'Ancestry', 'p-value']).size().reset_index(name='row_count')
# 
# print(count_per_group)
# 
# merged_df = pd.merge(combined_TL_1, combined_TL_2, on=['Ancestry', 'Approach', 'h2', 'num_causal_snps', 'Ancestry', 'p-value', 'sim'], how='outer',suffixes=('_1', '_2'))
# 
# merged_df.head()
# 
# filtered_rows = merged_df[(merged_df['Validation_Correlation_1'] < 0) & (merged_df['Approach'] == "Fine_Tune")]
# 
# print(filtered_rows)
# 
# # Iterate through each row and select columns with '_1' suffix if val_correlation_1 > val_correlation_2
# selected_rows = []
# 
# for index, row in merged_df.iterrows():
#     if row['Validation_Correlation_1'] > row['Validation_Correlation_2']:
#         selected_columns = ['Ancestry', 'Approach', 'Train_Correlation_1', 'Validation_Correlation_1', 'Test_Correlation_1', 'Train_MSE_1', 'Validation_MSE_1', 'Test_MSE_1', 'Train_R2_1', 'Validation_R2_1', 'Test_R2_1', 'time_1', 'h2', 'num_causal_snps', 'p-value', 'sim']
#     else:
#         selected_columns = ['Ancestry', 'Approach', 'Train_Correlation_2', 'Validation_Correlation_2', 'Test_Correlation_2', 'Train_MSE_2', 'Validation_MSE_2', 'Test_MSE_2', 'Train_R2_2', 'Validation_R2_2', 'Test_R2_2', 'time_2', 'h2', 'num_causal_snps', 'p-value', 'sim']
#     selected_row = row.loc[selected_columns]
#     row_dataframe = pd.DataFrame(selected_row).transpose()
#     row_dataframe.columns = ['Ancestry', 'method', 'Train_Correlation', 'Validation_Correlation', 'Test_Correlation', 'Train_MSE', 'Validation_MSE', 'Test_MSE', 'Train_R2', 'Validation_R2', 'Test_R2', 'time', 'h2', 'num_causal_snps',  'p-value', 'sim']
#     selected_rows.append(row_dataframe)

combined_TL_best = combined_TL_1[['Ancestry', 'method', 'Train_Correlation', 'Validation_Correlation', 'Test_Correlation', 'Train_MSE', 'Validation_MSE', 'Test_MSE', 'Train_R2', 'Validation_R2', 'Test_R2', 'h2', 'num_causal_snps',  'p-value', 'time', 'sim']]
combined_TL_best['Ancestry'] = combined_TL_best['Ancestry'].astype('category')
combined_TL_best['h2'] = combined_TL_best['h2'].astype('category')
combined_TL_best['num_causal_snps'] = combined_TL_best['num_causal_snps'].astype('category')
combined_TL_best['p-value'] = combined_TL_best['p-value'].astype('category')
count_per_group = combined_TL_best.groupby(['h2', 'num_causal_snps', 'Ancestry', 'p-value']).size().reset_index(name='row_count')

print(count_per_group)
   
################################### Read in CNN Metrics ###############################################
######################################################################################################

dataframes = []

for ancestry in {"CEU", "CHB", "YRI", "TSI"}:                    
    for sim in range(1,20):
        for h2 in {"0.05", "0.1"}:
            for ncs in {"1k", "10k"}:
                for p in {"4"}:
                    nn_directory = f"/scratch/duffme_gensim/Output_Error_Files/target_CNN_non_linear_phenotype_neural_networks/{ancestry}/sim_{sim}/"
                    # Define the file pattern
                    file_pattern_1 = f"{ancestry}_CNN_linear_pheno_sim{sim}_ncs{ncs}_h2{h2}_pvalue_{p}_bs64"
                    # Use glob to get a list of file paths that match the pattern
                    file_paths_1 = glob.glob(f"{nn_directory}/**/{file_pattern_1}"+ "*", recursive=True)
                    # Loop through each file path and read it as a DataFrame
                    if file_paths_1:
                        for file_path in file_paths_1:
                            # Read the CSV file into a DataFrame
                            df=pd.read_csv(file_path)
                            df['sim']=f'{sim}'
                            df['method']="CNN"
                            df['Ancestry']=f'{ancestry}'
                            # Append the DataFrame to the list
                            dataframes.append(df)
                        
combined_CNN = pd.concat(dataframes, ignore_index=True)

combined_CNN.head()

count_per_group = combined_CNN.groupby(['h2', 'num_causal_snps', 'Ancestry', 'method']).size().reset_index(name='row_count')

print(count_per_group)


################################### Read in CNN TL Metrics ###############################################
######################################################################################################

dataframes = []

for ancestry in {"CEU", "CHB", "YRI", "TSI"}:                    
    for sim in range(1,20):
        for h2 in {"0.05", "0.1"}:
            for ncs in {"1k", "10k"}:
                for p in {"4"}:
                    nn_directory = f"/scratch/duffme_gensim/Output_Error_Files/target_TL_CNN_non_linear_phenotype_neural_networks/{ancestry}/sim_{sim}/"
                    # Define the file pattern
                    file_pattern_1 = f"{ancestry}_TL_CNN_non_linear_pheno_sim{sim}_ncs{ncs}_h2{h2}_pvalue_{p}_bs64"
                    # Use glob to get a list of file paths that match the pattern
                    file_paths_1 = glob.glob(f"{nn_directory}/**/{file_pattern_1}"+ "*", recursive=True)
                    # Loop through each file path and read it as a DataFrame
                    if file_paths_1:
                        for file_path in file_paths_1:
                            # Read the CSV file into a DataFrame
                            df=pd.read_csv(file_path)
                            df['sim']=f'{sim}'
                            df['method']="TL_CNN_"+ df['Approach']
                            # Append the DataFrame to the list
                            dataframes.append(df)
                        
combined_TL_temp = pd.concat(dataframes, ignore_index=True)

grouped = combined_TL_temp.groupby(['sim', 'h2', 'num_causal_snps', 'Ancestry', 'p-value', 'Approach'])

combined_TL_CNN = grouped.apply(lambda x: x.sort_values(by='Validation_Correlation', ascending=False).head(1))

combined_TL_CNN.reset_index(drop=True, inplace=True)

count_per_group = combined_TL_CNN.groupby(['h2', 'num_causal_snps', 'Ancestry', 'p-value', 'Approach']).size().reset_index(name='row_count')

print(count_per_group)

################################### Combine Metric Dataframes ###############################################
############################################################################################################

combined_bayes.columns = ["Train_MSE", "Train_Correlation", "Train_R2", "Test_MSE", "Test_Correlation", "Test_R2", "sim", "h2", "num_causal_snps", "p-value", "time", "method", "Ancestry"]

subset_combined_bayes = combined_bayes[['Train_Correlation', 'Test_Correlation',
       'Train_MSE','Test_MSE', 'Train_R2',
       'Test_R2', 'sim', 'h2', 'num_causal_snps', 'p-value', 'method', 'Ancestry']]
       
subset_combined_MLP = combined_MLP[['Ancestry', 'Train_Correlation', 'Test_Correlation',
       'Train_MSE','Test_MSE', 'Train_R2',
       'Test_R2', 'h2', 'num_causal_snps', 'p-value', 'sim', 'method']]
       
subset_combined_TL_best = combined_TL_best[['Ancestry', 'method', 'Train_Correlation', 'Test_Correlation',
       'Train_MSE','Test_MSE', 'Train_R2',
       'Test_R2', 'h2', 'num_causal_snps', 'p-value', 'sim']]
       
subset_combined_CNN = combined_CNN[['Ancestry', 'method', 'Train_Correlation', 'Test_Correlation',
       'Train_MSE','Test_MSE', 'Train_R2',
       'Test_R2', 'h2', 'num_causal_snps', 'p-value', 'sim']]
       
subset_combined_TL_CNN = combined_TL_CNN[['Ancestry', 'method', 'Train_Correlation', 'Test_Correlation',
       'Train_MSE','Test_MSE', 'Train_R2',
       'Test_R2', 'h2', 'num_causal_snps', 'p-value', 'sim']]

# combined_bayes['Train_MSE'] = pd.to_numeric(combined_bayes['Train_MSE'])
# combined_bayes['Train_Correlation'] = pd.to_numeric(combined_bayes['Train_Correlation'])
# combined_bayes['Train_R2'] = pd.to_numeric(combined_bayes['Train_R2'])
# combined_bayes['Test_MSE'] = pd.to_numeric(combined_bayes['Test_MSE'])
# combined_bayes['Test_Correlation'] = pd.to_numeric(combined_bayes['Test_Correlation'])
# combined_bayes['Test_R2'] = pd.to_numeric(combined_bayes['Test_R2'])
# combined_bayes['h2'] = pd.to_numeric(combined_bayes['h2'])
# combined_bayes['p-value'] = pd.to_numeric(combined_bayes['p-value'])

non_linear_target_sets = pd.concat([combined_bayes, subset_combined_MLP, subset_combined_TL_best, subset_combined_CNN, subset_combined_TL_CNN])


###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
target_sets = pd.concat([linear_target_sets, non_linear_target_sets])

target_sets['num_causal_snps'] = target_sets['num_causal_snps'].astype('category')
target_sets['Ancestry'] = target_sets['Ancestry'].astype('category')
target_sets['method'] = target_sets['method'].astype('category')
target_sets['Train_Correlation'] = pd.to_numeric(target_sets['Train_Correlation'])
target_sets['Test_Correlation'] = pd.to_numeric(target_sets['Test_Correlation'])
target_sets['p-value'] = target_sets['p-value'].astype('category')
target_sets['h2'] = target_sets['h2'].astype('category')
target_sets['Test_R2'] = pd.to_numeric(target_sets['Test_R2'])
target_sets['Train_R2'] = pd.to_numeric(target_sets['Train_R2'])
target_sets['Train_MSE'] = pd.to_numeric(target_sets['Train_MSE'])
target_sets['Test_MSE'] = pd.to_numeric(target_sets['Test_MSE'])

target_sets['p-value'] = target_sets['p-value'].cat.rename_categories({
    '0.0001': '1e-04'
})

target_sets['method'] = target_sets['method'].replace({
    'bayes_r': 'BayesR',
    'gblup': 'GBLUP'
})

print(target_sets.columns)

unique_categories = target_sets['num_causal_snps'].cat.categories
print("num-casual-snps:")
print(unique_categories)

unique_categories = target_sets['p-value'].cat.categories
print("p-value:")
print(unique_categories)

unique_categories = target_sets['h2'].cat.categories
print("h2:")
print(unique_categories)

print(target_sets['Test_Correlation'].describe())
print(target_sets['Test_R2'].describe())

##################################################################################### PLOTS #############################################################################################
#########################################################################################################################################################################################
#########################################################################################################################################################################################

####################### Test Correlation Means with Error bounds -- bar plot ########################
# Function to calculate mean and confidence interval
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h, n

target_sets_2 = target_sets.dropna(subset=['Test_Correlation'])

# Calculate mean and confidence intervals for each group
summary_df = target_sets_2.groupby(['method', 'Ancestry', 'num_causal_snps', 'h2']).agg(
    Test_Correlation_Mean=('Test_Correlation', lambda x: mean_confidence_interval(x)[0]),
    lower=('Test_Correlation', lambda x: mean_confidence_interval(x)[1]),
    upper=('Test_Correlation', lambda x: mean_confidence_interval(x)[2]),
    count=('Test_Correlation', lambda x: mean_confidence_interval(x)[3])
).reset_index()

############################################### CREATE COUNT OF COMPLETED REPLICATES ####################################################
#########################################################################################################################################
#########################################################################################################################################
counts_table = summary_df.groupby(['Ancestry', 'num_causal_snps', 'h2', 'method']).size().reset_index(name='Count')

# Display the table
print(counts_table)



############################################### PLOT 1: Subset 100 causal SNPs ##########################################################
#########################################################################################################################################
#########################################################################################################################################

#summary_df_subset = summary_df[(summary_df['num_causal_snps']=="10k")&(summary_df['Ancestry']!="CHB")&(summary_df['Ancestry']!="TSI")&(summary_df['method']!="TL_CNN_Freeze")&(summary_df['method']!="TL_MLP_Freeze")&(summary_df['method']!="gblup")]

summary_df_subset = summary_df[(summary_df['num_causal_snps']==100) &
                               (summary_df['Ancestry']!="CHB") &
                               (summary_df['Ancestry']!="TSI") &
                               (summary_df['method']!="CNN") &
                               (summary_df['method']!="TL_CNN_Fine_Tune") &
                               (summary_df['method']!="TL_CNN_Freeze") &
                               (summary_df['method']!="TL_CNN_Fine_Tune") &
                               (summary_df['method']!="TL_MLP_Freeze") &
                               (summary_df['method']!="GBLUP") &
                               (summary_df['h2']!=0.05) &
                               (summary_df['h2']!=0.1)
                               ]

summary_df_subset['Ancestry'] = summary_df_subset['Ancestry'].cat.remove_unused_categories()
summary_df_subset['h2'] = summary_df_subset['h2'].cat.remove_unused_categories()
summary_df_subset['method'] = summary_df_subset['method'].cat.remove_unused_categories()

summary_df_subset['method'] = pd.Categorical(summary_df_subset['method'], 
                                             categories=['TL_MLP_Fine_Tune', 'MLP', 'BayesR'],
                                             ordered=True)

summary_df_subset['method'] = summary_df_subset['method'].cat.rename_categories({
    'TL_MLP_Fine_Tune': 'TL_MLP_FT',
    'MLP': 'MLP',
    'BayesR': 'Bayes R'
})

plotted_ancestries = summary_df_subset['Ancestry'].unique()

color_mapping = {
    'CEU': '#E69F00',  # Orange
    'YRI': '#56B4E9',  # Sky Blue
    'CHB': '#009E73',  # Bluish Green
    'TSI': '#D55E00'   # Vermilion
}

# Ensure that only plotted ancestries are included in the color mapping
color_mapping_filtered = {key: color_mapping[key] for key in plotted_ancestries if key in color_mapping}
#Additive Effects Phenotypes: 

# Define custom labeller for facet grids
summary_df_subset['pheno_format'] = summary_df_subset['num_causal_snps'].astype(str) + "/" + summary_df_subset['h2'].astype(str) + "/" + "0"

# Create the ggplot object
gg = (
    ggplot(summary_df_subset)
    + aes(x="method", y="Test_Correlation_Mean", fill="Ancestry", group="Ancestry")
    + facet_grid(". ~ pheno_format")
    + labs(
        x="Method",
        y="Test Correlation",
        title="Test Correlation Mean (with 95% error bounds)"
    )
    + theme(
        text=element_text(size=12),
        axis_text_x=element_text(size=10),
        axis_text_y=element_text(size=10),
        axis_title_x=element_text(size=12),
        axis_title_y=element_text(size=12),
        legend_text=element_text(size=10),
        legend_title=element_text(size=12),
        plot_title=element_text(size=12),
        strip_text=element_text(size=12)
    )
    + geom_bar(stat="identity", position=position_dodge(width=0.8))
    + geom_errorbar(aes(ymin='lower', ymax='upper'), width=0.2, position=position_dodge(width=0.8))
    + scale_fill_manual(values=color_mapping_filtered)
)

# Save the plot
gg.save(f'/scratch/duffme_gensim/Plots/linear_phenotype_target_correlation_mean_error_bounds_bar_plot_snapshot_ncs_100.png', width=8, height=5, dpi=600)

############################################### PLOT 2: Effect Size Distribution ##########################################################
#########################################################################################################################################
#########################################################################################################################################
dataframes = []

for ancestry in {'GBR'}:   
    for h2 in {'0.4', '0.8'}:
        for ncs in {100, '10k'}:
            effect_df = pd.read_csv(f"/scratch/duffme_gensim/Simulations/{ancestry}/sim_1/subset_1_{ancestry}_{ncs}_{h2}_phenotypes.par", sep = "\t")
            effect_df['ncs'] = ncs
            effect_df['h2'] = h2
            effect_df['Ancestry']=ancestry
            dataframes.append(effect_df)

combined_effect_df = pd.concat(dataframes, ignore_index=True)


plotted_ancestries = combined_effect_df['Ancestry'].unique()

color_mapping = {
    'CEU': '#E69F00',  # Orange
    'YRI': '#56B4E9',  # Sky Blue
    'CHB': '#009E73',  # Bluish Green
    'TSI': '#D55E00',  # Vermilion
    'GBR': '#CC79A7'   # Reddish Purple
}

# Ensure that only plotted ancestries are included in the color mapping
color_mapping_filtered = {key: color_mapping[key] for key in plotted_ancestries if key in color_mapping}
#Additive Effects Phenotypes: 
# Plotting
combined_effect_df_subset=combined_effect_df[combined_effect_df['Ancestry']=='GBR']

combined_effect_df_subset['pheno_format'] = combined_effect_df_subset['ncs'].astype(str) + "/" + combined_effect_df_subset['h2'].astype(str) + "/" + "0"

gg = (
    ggplot(combined_effect_df_subset)
    + aes(x="Effect", color="Ancestry", fill="Ancestry", group="Ancestry")
    + facet_grid("pheno_format ~ ", scales='free')
    + labs(
        x="Effect Size",
        y="Count",
        title="Phenotype Effect Size Distribution"
    )
    + theme(
        text=element_text(size=12),
        axis_text_x=element_text(size=10),
        axis_text_y=element_text(size=8),
        axis_title_x=element_text(size=12),
        axis_title_y=element_text(size=12),
        legend_text=element_text(size=10),
        legend_title=element_text(size=12),
        plot_title=element_text(size=12),
        strip_text=element_text(size=12),
        strip_text_y=element_text(size=12),  # Ensure facet labels are visible
        strip_text_x=element_text(size=12)   # Ensure facet labels are visible
    )
    + geom_histogram(position=position_dodge(width=0.8), bins=100)  # Adjust 'bins' as needed
    + scale_color_manual(values=color_mapping_filtered)
    + scale_fill_manual(values=color_mapping_filtered)  # Adding fill scale to match color scale
)

# Save the plot
gg.save('/scratch/duffme_gensim/Plots/linear_effect_vs_maf.png', width=8, height=5, dpi=600)

############################################### PLOT 3: Subset 10k causal SNPs ##########################################################
#########################################################################################################################################
#########################################################################################################################################

summary_df_subset = summary_df[(summary_df['num_causal_snps']=="10k") &
                               (summary_df['Ancestry']!="CHB") &
                               (summary_df['Ancestry']!="TSI") &
                               (summary_df['method']!="CNN") &
                               (summary_df['method']!="TL_CNN_Fine_Tune") &
                               (summary_df['method']!="TL_CNN_Freeze") &
                               (summary_df['method']!="TL_CNN_Fine_Tune") &
                               (summary_df['method']!="TL_MLP_Freeze") &
                               (summary_df['method']!="GBLUP") &
                               (summary_df['h2']!=0.05) &
                               (summary_df['h2']!=0.1)
                               ]

summary_df_subset['Ancestry'] = summary_df_subset['Ancestry'].cat.remove_unused_categories()
summary_df_subset['h2'] = summary_df_subset['h2'].cat.remove_unused_categories()
summary_df_subset['method'] = summary_df_subset['method'].cat.remove_unused_categories()

summary_df_subset['method'] = pd.Categorical(summary_df_subset['method'], 
                                             categories=['TL_MLP_Fine_Tune', 'MLP', 'BayesR'],
                                             ordered=True)

summary_df_subset['method'] = summary_df_subset['method'].cat.rename_categories({
    'MLP': 'MLP',
    'TL_MLP_Fine_Tune': 'TL_MLP_FT',
    'BayesR': 'Bayes R'
})

plotted_ancestries = summary_df_subset['Ancestry'].unique()

color_mapping = {
    'CEU': '#E69F00',  # Orange
    'YRI': '#56B4E9',  # Sky Blue
    'CHB': '#009E73',  # Bluish Green
    'TSI': '#D55E00'   # Vermilion
}

# Ensure that only plotted ancestries are included in the color mapping
color_mapping_filtered = {key: color_mapping[key] for key in plotted_ancestries if key in color_mapping}

summary_df_subset['pheno_format'] = summary_df_subset['num_causal_snps'].astype(str) + "/" + summary_df_subset['h2'].astype(str) + "/" + "0"

#Additive Effects Phenotypes: 
# Plotting
gg = (
    ggplot(summary_df_subset)
    + aes(x="method", y="Test_Correlation_Mean", fill="Ancestry", group="Ancestry")
    + facet_grid(". ~ pheno_format")
    + labs(
        x="Method",
        y="Test Correlation",
        title="Test Correlation Mean (with 95% error bounds)"
    )
    + theme(
        text=element_text(size=12),
        axis_text_x=element_text(size=10),
        axis_text_y=element_text(size=10),
        axis_title_x=element_text(size=12),
        axis_title_y=element_text(size=12),
        legend_text=element_text(size=10),
        legend_title=element_text(size=12),
        plot_title=element_text(size=12),
        strip_text=element_text(size=12)
    )
    + geom_bar(stat="identity", position=position_dodge(width=0.8))
    + geom_errorbar(aes(ymin='lower', ymax='upper'), width=0.2, position=position_dodge(width=0.8))
    + scale_fill_manual(values=color_mapping_filtered)
)

gg.save(f'/scratch/duffme_gensim/Plots/linear_phenotype_target_correlation_mean_error_bounds_bar_plot_snapshot_ncs_10k.png', width=8, height=5, dpi=600)

############################################### PLOT 4: Subset additive & interaction effect phenotypes causal SNPs ##########################################################
##############################################################################################################################################################################
##############################################################################################################################################################################

summary_df_subset = summary_df[
                               (summary_df['num_causal_snps']!=100) &
                               (summary_df['Ancestry']!="CHB") &
                               (summary_df['Ancestry']!="TSI") &
                               (summary_df['method']!="CNN") &
                               (summary_df['method']!="TL_CNN_Fine_Tune") &
                               (summary_df['method']!="TL_CNN_Freeze") &
                               (summary_df['method']!="TL_CNN_Fine_Tune") &
                               (summary_df['method']!="TL_MLP_Freeze") &
                               (summary_df['method']!="GBLUP") &
                               (summary_df['h2']!=0.4) &
                               (summary_df['h2']!=0.8)
                               ]

summary_df_subset['Ancestry'] = summary_df_subset['Ancestry'].cat.remove_unused_categories()
summary_df_subset['h2'] = summary_df_subset['h2'].cat.remove_unused_categories()
summary_df_subset['method'] = summary_df_subset['method'].cat.remove_unused_categories()

summary_df_subset['method'] = pd.Categorical(summary_df_subset['method'], 
                                             categories=['TL_MLP_Fine_Tune', 'MLP', 'BayesR'],
                                             ordered=True)

summary_df_subset['method'] = summary_df_subset['method'].cat.rename_categories({
    'MLP': 'MLP',
    'TL_MLP_Fine_Tune': 'TL_MLP_FT',
    'BayesR': 'Bayes R'
})

summary_df_subset['h2'] = summary_df_subset['h2'].cat.reorder_categories(
    [0.05, 0.1], ordered=True)

summary_df_subset['num_causal_snps'] = summary_df_subset['num_causal_snps'].cat.remove_unused_categories()
summary_df_subset['num_causal_snps'] = summary_df_subset['num_causal_snps'].cat.reorder_categories(
    ['1k', '10k'], ordered=True)

plotted_ancestries = summary_df_subset['Ancestry'].unique()

color_mapping = {
    'CEU': '#E69F00',  # Orange
    'YRI': '#56B4E9',  # Sky Blue
    'CHB': '#009E73',  # Bluish Green
    'TSI': '#D55E00'   # Vermilion
}

# Ensure that only plotted ancestries are included in the color mapping
color_mapping_filtered = {key: color_mapping[key] for key in plotted_ancestries if key in color_mapping}

# Additive heritability column 
summary_df_subset['add_h2'] = None

# Apply conditions for assignment
summary_df_subset.loc[summary_df_subset['h2'] == 0.05, 'add_h2'] = 0.475
summary_df_subset.loc[summary_df_subset['h2'] == 0.1, 'add_h2'] = 0.45

# Interaction heritability column 
summary_df_subset['int_h2'] = None

# Apply conditions for assignment
summary_df_subset.loc[summary_df_subset['h2'] == 0.05, 'int_h2'] = 0.025
summary_df_subset.loc[summary_df_subset['h2'] == 0.1, 'int_h2'] = 0.05

# Create new variable to name the facets how I want 
summary_df_subset['pheno_format'] = summary_df_subset['num_causal_snps'].astype(str) + "/" + summary_df_subset['add_h2'].astype(str) + "/" + summary_df_subset['int_h2'].astype(str)

summary_df_subset['pheno_format'] = pd.Categorical(summary_df_subset['pheno_format'], 
                                             categories=['1k/0.45/0.05', '1k/0.475/0.025', '10k/0.45/0.05', '10k/0.475/0.025'],
                                             ordered=True)

#Additive Effects Phenotypes: 
# Plotting
gg = (
    ggplot(summary_df_subset)
    + aes(x="method", y="Test_Correlation_Mean", fill="Ancestry", group="Ancestry")
    + facet_grid(". ~ pheno_format")
    + labs(
        x="Method",
        y="Test Correlation",
        title="Test Correlation Mean (with 95% error bounds)"
    )
    + theme(
        text=element_text(size=12),
        axis_text_x=element_text(size=10, angle=45),
        axis_text_y=element_text(size=10),
        axis_title_x=element_text(size=12),
        axis_title_y=element_text(size=12),
        legend_text=element_text(size=10),
        legend_title=element_text(size=12),
        plot_title=element_text(size=12),
        strip_text=element_text(size=12)
    )
    + geom_bar(stat="identity", position=position_dodge(width=0.8))
    + geom_errorbar(aes(ymin='lower', ymax='upper'), width=0.2, position=position_dodge(width=0.8))
    + scale_fill_manual(values=color_mapping_filtered)
)

gg.save(f'/scratch/duffme_gensim/Plots/non_linear_phenotypes_target_correlation_mean_error_bounds_bar_plot_snapshot.png', width=8, height=5, dpi=600)

############################################### PLOT 5: Test Correlation Relative Improvement ################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
target_sets_subset = target_sets[(target_sets["method"]!="TL_CNN_Freeze") & (target_sets["method"]!="TL_MLP_Freeze") & (target_sets["num_causal_snps"]!=100)]
target_sets_subset['sim'] = target_sets_subset['sim'].astype(str)

#def rel_improve_by_bayes_r(group):
#    # Find the Test_Correlation for BayesR
#    bayes_r_value = group.loc[group['method'] == 'bayes_r', 'Test_Correlation']
#    if not bayes_r_value.empty:
#        # If BayesR is present, normalize by its Test_Correlation
#        bayes_r_value = bayes_r_value.iloc[0]
#        group['Relative_Improvement_Test_Correlation_baseline_bayesR'] = group['Test_Correlation'] / bayes_r_value
#    else:
#        # If BayesR is not present, set Normalized_Test_Correlation to NaN or 1
#        group['Relative_Improvement_Test_Correlation_baseline_bayesR'] = pd.NA  # or group['Test_Correlation']
#    return group

def rel_improve_by_best_method(group):
    # Find the maximum Test_Correlation value in the group
    subset = group[(group['method']!="TL_MLP_Fine_Tune") & (group['method']!="TL_CNN_Fine_Tune")]
    best_value = subset['Test_Correlation'].max()
    # Normalize each method's Test_Correlation by the best value
    rel_subset = group[(group['method']=="TL_MLP_Fine_Tune") | 
                       (group['method']=="TL_CNN_Fine_Tune")]
    rel_subset['Relative_Improvement_Test_Correlation_baseline_best'] = rel_subset['Test_Correlation'] / best_value
    return rel_subset

rel_improve_target_sets = target_sets_subset.groupby(["sim", "Ancestry", "num_causal_snps", "h2"]).apply(rel_improve_by_best_method)
rel_improve_target_sets.reset_index(drop=True, inplace=True)
rel_improve_target_sets['method'] = rel_improve_target_sets['method'].cat.remove_unused_categories()

rel_improve_target_sets['method'] = rel_improve_target_sets['method'].cat.rename_categories({
    'TL_CNN_Fine_Tune': 'TL_CNN_FT',
    'TL_MLP_Fine_Tune': 'TL_MLP_FT'
})

rel_improve_target_sets['method'] = rel_improve_target_sets['method'].cat.reorder_categories(
    ['TL_MLP_FT', 'TL_CNN_FT'], ordered=True)

rel_improve_target_sets['h2'] = rel_improve_target_sets['h2'].cat.reorder_categories(
    [0.4, 0.05, 0.1, 0.8], ordered=True)

rel_improve_target_sets['num_causal_snps'] = rel_improve_target_sets['num_causal_snps'].cat.remove_unused_categories()
rel_improve_target_sets['num_causal_snps'] = rel_improve_target_sets['num_causal_snps'].cat.reorder_categories(
    ['1k', '10k'], ordered=True)

plotted_ancestries = rel_improve_target_sets['Ancestry'].unique()

color_mapping = {
    'CEU': '#E69F00',  # Orange
    'YRI': '#56B4E9',  # Sky Blue
    'CHB': '#009E73',  # Bluish Green
    'TSI': '#D55E00'   # Vermilion
}

# Ensure that only plotted ancestries are included in the color mapping
color_mapping_filtered = {key: color_mapping[key] for key in plotted_ancestries if key in color_mapping}

# Additive heritability column 
rel_improve_target_sets['add_h2'] = None

# Apply conditions for assignment
rel_improve_target_sets.loc[rel_improve_target_sets['h2'] == 0.05, 'add_h2'] = 0.475
rel_improve_target_sets.loc[rel_improve_target_sets['h2'] == 0.1, 'add_h2'] = 0.45
rel_improve_target_sets.loc[rel_improve_target_sets['h2'] == 0.4, 'add_h2'] = 0.4
rel_improve_target_sets.loc[rel_improve_target_sets['h2'] == 0.8, 'add_h2'] = 0.8

# Interaction heritability column 
rel_improve_target_sets['int_h2'] = None

# Apply conditions for assignment
rel_improve_target_sets.loc[rel_improve_target_sets['h2'] == 0.05, 'int_h2'] = 0.025
rel_improve_target_sets.loc[rel_improve_target_sets['h2'] == 0.1, 'int_h2'] = 0.05
rel_improve_target_sets.loc[rel_improve_target_sets['h2'] == 0.4, 'int_h2'] = 0
rel_improve_target_sets.loc[rel_improve_target_sets['h2'] == 0.8, 'int_h2'] = 0

# Create new variable to name the facets how I want 
rel_improve_target_sets['pheno_format'] = rel_improve_target_sets['num_causal_snps'].astype(str) + "/" + rel_improve_target_sets['add_h2'].astype(str) + "/" + rel_improve_target_sets['int_h2'].astype(str)

rel_improve_target_sets['pheno_format'] = pd.Categorical(rel_improve_target_sets['pheno_format'], 
                                             categories=['1k/0.45/0.05', '1k/0.475/0.025', '10k/0.4/0', '10k/0.45/0.05', '10k/0.475/0.025', '10k/0.8/0'],
                                             ordered=True)


gg = (
  ggplot(rel_improve_target_sets)
  + aes(x="method", y="Relative_Improvement_Test_Correlation_baseline_best", color="Ancestry")
  + facet_grid(". ~ pheno_format")
  + labs(
    x="Method",
    y="Relative Improvement of Test Correlation",
    title=f"Test Correlation Relative Improvement (Baseline: Comparison Method with Best Performance)",
    )
  + theme(
        text=element_text(size=11),
        axis_text_x=element_text(size=8, angle=45),
        axis_text_y=element_text(size=8),
        axis_title_x=element_text(size=10),
        axis_title_y=element_text(size=10),
        legend_text=element_text(size=8),
        legend_title=element_text(size=10),
        plot_title=element_text(size=10),
        strip_text=element_text(size=8)
    )
  + geom_boxplot()
  + geom_hline(yintercept=1, linetype='dashed', color='grey')
  + scale_color_manual(values=color_mapping_filtered)
    )
    
gg.save(f'/scratch/duffme_gensim/Plots/snapshot_target_check_test_corr_rel_improve.png', width=8, height=5, dpi=600)

############################################### PLOT 6: Test correlation standard deviation ################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################

summary_statistics = target_sets.groupby(['h2', 'num_causal_snps', 'Ancestry', 'p-value', 'method'])['Test_Correlation'].describe().reset_index()

summary_df_subset=summary_statistics[(summary_statistics['method']!="TL_MLP_Freeze") &
                                     (summary_statistics['method']!="TL_CNN_Freeze") &
                                     (summary_statistics['num_causal_snps']!=100) &
                                     (summary_statistics['h2']!=0.05) & 
                                     (summary_statistics['h2']!=0.1)
]

plotted_ancestries = summary_df_subset['Ancestry'].unique()

color_mapping = {
    'CEU': '#E69F00',  # Orange
    'YRI': '#56B4E9',  # Sky Blue
    'CHB': '#009E73',  # Bluish Green
    'TSI': '#D55E00'   # Vermilion
}

# Ensure that only plotted ancestries are included in the color mapping
color_mapping_filtered = {key: color_mapping[key] for key in plotted_ancestries if key in color_mapping}

summary_df_subset['method'] = summary_df_subset['method'].cat.rename_categories({
    'TL_CNN_Fine_Tune': 'TL_CNN_FT',
    'TL_MLP_Fine_Tune': 'TL_MLP_FT',
    'BayesR': 'Bayes R'
})

# Additive summary_df_subset column 
summary_df_subset['add_h2'] = None

# Apply conditions for assignment
summary_df_subset.loc[summary_df_subset['h2'] == 0.05, 'add_h2'] = 0.475
summary_df_subset.loc[summary_df_subset['h2'] == 0.1, 'add_h2'] = 0.45
summary_df_subset.loc[summary_df_subset['h2'] == 0.4, 'add_h2'] = 0.4
summary_df_subset.loc[summary_df_subset['h2'] == 0.8, 'add_h2'] = 0.8

# Interaction heritability column 
summary_df_subset['int_h2'] = None

# Apply conditions for assignment
summary_df_subset.loc[summary_df_subset['h2'] == 0.05, 'int_h2'] = 0.025
summary_df_subset.loc[summary_df_subset['h2'] == 0.1, 'int_h2'] = 0.05
summary_df_subset.loc[summary_df_subset['h2'] == 0.4, 'int_h2'] = 0
summary_df_subset.loc[summary_df_subset['h2'] == 0.8, 'int_h2'] = 0

# Create new variable to name the facets how I want 
summary_df_subset['pheno_format'] = summary_df_subset['num_causal_snps'].astype(str) + "/" + summary_df_subset['add_h2'].astype(str) + "/" + summary_df_subset['int_h2'].astype(str)

summary_df_subset['pheno_format'] = pd.Categorical(summary_df_subset['pheno_format'], 
                                             categories=['1k/0.45/0.05', '1k/0.475/0.025', '10k/0.4/0', '10k/0.45/0.05', '10k/0.475/0.025', '10k/0.8/0'],
                                             ordered=True)

summary_df_subset['method'] = pd.Categorical(summary_df_subset['method'], 
                                             categories=['MLP', 'CNN', 'TL_MLP_FT', 'TL_CNN_FT',  'Bayes R', 'GBLUP'],
                                             ordered=True)

gg = (
    ggplot(summary_df_subset)
    + aes(x="method", y="std", fill="Ancestry")
    + facet_grid(". ~ pheno_format")
    + labs(
        x="Method",
        y="Test Correlation Standard Deviation",
        title="Test Correlation Standard Deviation",
    )
    + theme(
        text=element_text(size=11),
        axis_text_x=element_text(size=8, angle=45),
        axis_text_y=element_text(size=8),
        axis_title_x=element_text(size=10),
        axis_title_y=element_text(size=10),
        legend_text=element_text(size=8),
        legend_title=element_text(size=10),
        plot_title=element_text(size=10),
        strip_text=element_text(size=10),
    )
    + geom_col(position=position_dodge(width=0.7))
    + scale_fill_manual(values=color_mapping_filtered)
)
    
gg.save(f'/scratch/duffme_gensim/Plots/std_correlation_target_check.png', width=8, height=5, dpi=600)

############################################### PLOT 7: MLP vs CNN ######################################################################
#########################################################################################################################################
#########################################################################################################################################

summary_df_subset = summary_df[
                               (summary_df['Ancestry']!="CHB") &
                               (summary_df['Ancestry']!="TSI") &
                               (summary_df['method']!="BayesR") &
                               (summary_df['method']!="TL_CNN_Freeze") &
                               (summary_df['method']!="TL_MLP_Freeze") &
                               (summary_df['method']!="GBLUP") &
                               (summary_df['h2']!=0.4) &
                               (summary_df['h2']!=0.8) &
                               (summary_df['num_causal_snps']!=100)
                               ]

summary_df_subset['Ancestry'] = summary_df_subset['Ancestry'].cat.remove_unused_categories()
summary_df_subset['h2'] = summary_df_subset['h2'].cat.remove_unused_categories()
summary_df_subset['num_causal_snps'] = summary_df_subset['num_causal_snps'].cat.remove_unused_categories()
summary_df_subset['method'] = summary_df_subset['method'].cat.remove_unused_categories()

summary_df_subset['method'] = summary_df_subset['method'].cat.rename_categories({
    'MLP': 'MLP',
    'CNN': 'CNN',
    'TL_MLP_Fine_Tune': 'TL_MLP_FT',
    'TL_CNN_Fine_Tune': 'TL_CNN_FT'
})

plotted_ancestries = summary_df_subset['Ancestry'].unique()

color_mapping = {
    'CEU': '#E69F00',  # Orange
    'YRI': '#56B4E9',  # Sky Blue
    'CHB': '#009E73',  # Bluish Green
    'TSI': '#D55E00'   # Vermilion
}

# Ensure that only plotted ancestries are included in the color mapping
color_mapping_filtered = {key: color_mapping[key] for key in plotted_ancestries if key in color_mapping}

summary_df_subset['method'] = pd.Categorical(summary_df_subset['method'], 
                                             categories=['MLP', 'CNN', 'TL_MLP_FT', 'TL_CNN_FT'],
                                             ordered=True)

# Additive summary_df_subset column 
summary_df_subset['add_h2'] = None

# Apply conditions for assignment
summary_df_subset.loc[summary_df_subset['h2'] == 0.05, 'add_h2'] = 0.475
summary_df_subset.loc[summary_df_subset['h2'] == 0.1, 'add_h2'] = 0.45
summary_df_subset.loc[summary_df_subset['h2'] == 0.4, 'add_h2'] = 0.4
summary_df_subset.loc[summary_df_subset['h2'] == 0.8, 'add_h2'] = 0.8

# Interaction heritability column 
summary_df_subset['int_h2'] = None

# Apply conditions for assignment
summary_df_subset.loc[summary_df_subset['h2'] == 0.05, 'int_h2'] = 0.025
summary_df_subset.loc[summary_df_subset['h2'] == 0.1, 'int_h2'] = 0.05
summary_df_subset.loc[summary_df_subset['h2'] == 0.4, 'int_h2'] = 0
summary_df_subset.loc[summary_df_subset['h2'] == 0.8, 'int_h2'] = 0

# Create new variable to name the facets how I want 
summary_df_subset['pheno_format'] = summary_df_subset['num_causal_snps'].astype(str) + "/" + summary_df_subset['add_h2'].astype(str) + "/" + summary_df_subset['int_h2'].astype(str)

summary_df_subset['pheno_format'] = pd.Categorical(summary_df_subset['pheno_format'], 
                                             categories=['1k/0.45/0.05', '1k/0.475/0.025', '10k/0.4/0', '10k/0.45/0.05', '10k/0.475/0.025', '10k/0.8/0'],
                                             ordered=True)

#Additive Effects Phenotypes: 
# Plotting
gg = (
    ggplot(summary_df_subset)
    + aes(x="method", y="Test_Correlation_Mean", fill="Ancestry", group="Ancestry")
    + facet_grid(". ~ pheno_format")
    + labs(
        x="Method",
        y="Test Correlation",
        title="Test Correlation Mean (with 95% error bounds)"
    )
    + theme(
        text=element_text(size=12),
        axis_text_x=element_text(size=10, angle=45),
        axis_text_y=element_text(size=10),
        axis_title_x=element_text(size=12),
        axis_title_y=element_text(size=12),
        legend_text=element_text(size=10),
        legend_title=element_text(size=12),
        plot_title=element_text(size=12),
        strip_text=element_text(size=12)
    )
    + geom_bar(stat="identity", position=position_dodge(width=0.8))
    + geom_errorbar(aes(ymin='lower', ymax='upper'), width=0.2, position=position_dodge(width=0.8))
    + scale_fill_manual(values=color_mapping_filtered)
)

gg.save(f'/scratch/duffme_gensim/Plots/cnn_vs_mlp.png', width=8, height=5, dpi=600)

############################################### PLOT 8: Freeze vs. Fine Tune ############################################################
#########################################################################################################################################
#########################################################################################################################################

summary_df_subset = summary_df[
                               (summary_df['Ancestry']!="CHB") &
                               (summary_df['Ancestry']!="TSI") &
                               (summary_df['method']!="BayesR") &
                               (summary_df['method']!="MLP") &
                               (summary_df['method']!="CNN") &
                               (summary_df['method']!="GBLUP") &
                               (summary_df['method']!="TL_CNN_Fine_Tune") &
                               (summary_df['method']!="TL_CNN_Freeze") &
                               (summary_df['h2']!=0.1) &
                               (summary_df['h2']!=0.05) &
                               (summary_df['num_causal_snps']!="1k")
                               ]

summary_df_subset['Ancestry'] = summary_df_subset['Ancestry'].cat.remove_unused_categories()
summary_df_subset['h2'] = summary_df_subset['h2'].cat.remove_unused_categories()
summary_df_subset['num_causal_snps'] = summary_df_subset['num_causal_snps'].cat.remove_unused_categories()
summary_df_subset['method'] = summary_df_subset['method'].cat.remove_unused_categories()

summary_df_subset['method'] = summary_df_subset['method'].cat.rename_categories({
    'TL_MLP_Fine_Tune': 'TL_MLP_FT',
    'TL_MLP_Freeze': 'TL_MLP_Freeze'
})

plotted_ancestries = summary_df_subset['Ancestry'].unique()

color_mapping = {
    'CEU': '#E69F00',  # Orange
    'YRI': '#56B4E9',  # Sky Blue
    'CHB': '#009E73',  # Bluish Green
    'TSI': '#D55E00'   # Vermilion
}

# Ensure that only plotted ancestries are included in the color mapping
color_mapping_filtered = {key: color_mapping[key] for key in plotted_ancestries if key in color_mapping}


summary_df_subset['method'] = pd.Categorical(summary_df_subset['method'], 
                                             categories=['TL_MLP_Freeze', 'TL_MLP_FT'],
                                             ordered=True)

# Additive summary_df_subset column 
summary_df_subset['add_h2'] = None

# Apply conditions for assignment
summary_df_subset.loc[summary_df_subset['h2'] == 0.05, 'add_h2'] = 0.475
summary_df_subset.loc[summary_df_subset['h2'] == 0.1, 'add_h2'] = 0.45
summary_df_subset.loc[summary_df_subset['h2'] == 0.4, 'add_h2'] = 0.4
summary_df_subset.loc[summary_df_subset['h2'] == 0.8, 'add_h2'] = 0.8

# Interaction heritability column 
summary_df_subset['int_h2'] = None

# Apply conditions for assignment
summary_df_subset.loc[summary_df_subset['h2'] == 0.05, 'int_h2'] = 0.025
summary_df_subset.loc[summary_df_subset['h2'] == 0.1, 'int_h2'] = 0.05
summary_df_subset.loc[summary_df_subset['h2'] == 0.4, 'int_h2'] = 0
summary_df_subset.loc[summary_df_subset['h2'] == 0.8, 'int_h2'] = 0

# Create new variable to name the facets how I want 
summary_df_subset['pheno_format'] = summary_df_subset['num_causal_snps'].astype(str) + "/" + summary_df_subset['add_h2'].astype(str) + "/" + summary_df_subset['int_h2'].astype(str)

summary_df_subset['pheno_format'] = pd.Categorical(summary_df_subset['pheno_format'], 
                                             categories=['100/0.4/0', '100/0.8/0', '10k/0.4/0', '10k/0.45/0.05', '10k/0.475/0.025', '10k/0.8/0'],
                                             ordered=True)

#Additive Effects Phenotypes: 
# Plotting
gg = (
    ggplot(summary_df_subset)
    + aes(x="method", y="Test_Correlation_Mean", fill="Ancestry", group="Ancestry")
    + facet_grid(". ~ pheno_format")
    + labs(
        x="Method",
        y="Test Correlation",
        title="Test Correlation Mean (with 95% error bounds)"
    )
    + theme(
        text=element_text(size=12),
        axis_text_x=element_text(size=10, angle=45),
        axis_text_y=element_text(size=10),
        axis_title_x=element_text(size=12),
        axis_title_y=element_text(size=12),
        legend_text=element_text(size=10),
        legend_title=element_text(size=12),
        plot_title=element_text(size=12),
        strip_text=element_text(size=12)
    )
    + geom_bar(stat="identity", position=position_dodge(width=0.8))
    + geom_errorbar(aes(ymin='lower', ymax='upper'), width=0.2, position=position_dodge(width=0.8))
    + scale_fill_manual(values=color_mapping_filtered)
)

gg.save(f'/scratch/duffme_gensim/Plots/freeze_vs_fine_tune.png', width=8, height=5, dpi=600)

############################################### PLOT 8: Freeze vs. Fine Tune ############################################################
#########################################################################################################################################
#########################################################################################################################################
target_sets_subset = target_sets[
                               (target_sets['method']!="MLP") &
                               (target_sets['method']!="CNN") &
                               (target_sets['method']!="GBLUP") &
                               (target_sets['method']!="TL_CNN_Fine_Tune") &
                               (target_sets['method']!="TL_CNN_Freeze") &
                               (target_sets['method']!="TL_MLP_Freeze") &
                               (target_sets['h2']!=0.1) &
                               (target_sets['h2']!=0.05) &
                               (target_sets['h2']!=0.4) &
                               (target_sets['num_causal_snps']!=100) &
                               (target_sets['num_causal_snps']!="1k")
                               ]

target_sets_subset['Ancestry'] = target_sets_subset['Ancestry'].cat.remove_unused_categories()
target_sets_subset['h2'] = target_sets_subset['h2'].cat.remove_unused_categories()
target_sets_subset['num_causal_snps'] = target_sets_subset['num_causal_snps'].cat.remove_unused_categories()
target_sets_subset['method'] = target_sets_subset['method'].cat.remove_unused_categories()

target_sets_subset['method'] = target_sets_subset['method'].cat.rename_categories({
    'TL_MLP_Fine_Tune': 'TL_MLP_FT',
    'BayesR': 'Bayes R'
})

plotted_ancestries = target_sets_subset['Ancestry'].unique()

color_mapping = {
    'CEU': '#E69F00',  # Orange
    'YRI': '#56B4E9',  # Sky Blue
    'CHB': '#009E73',  # Bluish Green
    'TSI': '#D55E00'   # Vermilion
}

# Ensure that only plotted ancestries are included in the color mapping
color_mapping_filtered = {key: color_mapping[key] for key in plotted_ancestries if key in color_mapping}

# Additive summary_df_subset column 
target_sets_subset['add_h2'] = None

# Apply conditions for assignment
target_sets_subset.loc[target_sets_subset['h2'] == 0.05, 'add_h2'] = 0.475
target_sets_subset.loc[target_sets_subset['h2'] == 0.1, 'add_h2'] = 0.45
target_sets_subset.loc[target_sets_subset['h2'] == 0.4, 'add_h2'] = 0.4
target_sets_subset.loc[target_sets_subset['h2'] == 0.8, 'add_h2'] = 0.8

# Interaction heritability column 
target_sets_subset['int_h2'] = None

# Apply conditions for assignment
target_sets_subset.loc[target_sets_subset['h2'] == 0.05, 'int_h2'] = 0.025
target_sets_subset.loc[target_sets_subset['h2'] == 0.1, 'int_h2'] = 0.05
target_sets_subset.loc[target_sets_subset['h2'] == 0.4, 'int_h2'] = 0
target_sets_subset.loc[target_sets_subset['h2'] == 0.8, 'int_h2'] = 0

# Create new variable to name the facets how I want 
target_sets_subset['pheno_format'] = target_sets_subset['num_causal_snps'].astype(str) + "/" + target_sets_subset['add_h2'].astype(str) + "/" + target_sets_subset['int_h2'].astype(str)

target_sets_subset['pheno_format'] = pd.Categorical(target_sets_subset['pheno_format'], 
                                             categories=['100/0.4/0', '100/0.8/0', '10k/0.4/0', '10k/0.45/0.05', '10k/0.475/0.025', '10k/0.8/0'],
                                             ordered=True)

#Additive Effects Phenotypes: 
# Plotting
gg = (
    ggplot(target_sets_subset)
    + aes(x="method", y="Test_Correlation", color="Ancestry")
    + facet_grid(". ~ pheno_format")
    + labs(
        x="Method",
        y="Test Correlation",
        title="Ancestry Comparison Snapshot"
    )
    + theme(
        text=element_text(size=12),
        axis_text_x=element_text(size=10, angle=45),
        axis_text_y=element_text(size=10),
        axis_title_x=element_text(size=12),
        axis_title_y=element_text(size=12),
        legend_text=element_text(size=10),
        legend_title=element_text(size=12),
        plot_title=element_text(size=12),
        strip_text=element_text(size=12)
    )
    + geom_boxplot()
    + scale_color_manual(values=color_mapping_filtered)
)

gg.save(f'/scratch/duffme_gensim/Plots/ancestry_comparison.png', width=8, height=5, dpi=600)

############################################### Appendix Plots 1-2: Additive Effects Only  ############################################################
########################################################################################################################################################
########################################################################################################################################################
target_sets_subset = target_sets[
                               (target_sets['h2']!=0.1) &
                               (target_sets['h2']!=0.05) 
                               ]

target_sets_subset['Ancestry'] = target_sets_subset['Ancestry'].cat.remove_unused_categories()
target_sets_subset['h2'] = target_sets_subset['h2'].cat.remove_unused_categories()
target_sets_subset['num_causal_snps'] = target_sets_subset['num_causal_snps'].cat.remove_unused_categories()
target_sets_subset['method'] = target_sets_subset['method'].cat.remove_unused_categories()

target_sets_subset['method'] = target_sets_subset['method'].cat.rename_categories({
    'TL_MLP_Fine_Tune': 'TL_MLP_FT',
    'TL_CNN_Fine_Tune': 'TL_CNN_FT',
    'TL_MLP_Freeze': 'TL_MLP_FR',
    'TL_CNN_Freeze': 'TL_CNN_FR',
    'BayesR': 'Bayes R'
})

plotted_ancestries = target_sets_subset['Ancestry'].unique()

color_mapping = {
    'CEU': '#E69F00',  # Orange
    'YRI': '#56B4E9',  # Sky Blue
    'CHB': '#009E73',  # Bluish Green
    'TSI': '#D55E00'   # Vermilion
}

# Ensure that only plotted ancestries are included in the color mapping
color_mapping_filtered = {key: color_mapping[key] for key in plotted_ancestries if key in color_mapping}

# Additive summary_df_subset column 
target_sets_subset['add_h2'] = None

# Apply conditions for assignment
target_sets_subset.loc[target_sets_subset['h2'] == 0.05, 'add_h2'] = 0.475
target_sets_subset.loc[target_sets_subset['h2'] == 0.1, 'add_h2'] = 0.45
target_sets_subset.loc[target_sets_subset['h2'] == 0.4, 'add_h2'] = 0.4
target_sets_subset.loc[target_sets_subset['h2'] == 0.8, 'add_h2'] = 0.8

# Interaction heritability column 
target_sets_subset['int_h2'] = None

# Apply conditions for assignment
target_sets_subset.loc[target_sets_subset['h2'] == 0.05, 'int_h2'] = 0.025
target_sets_subset.loc[target_sets_subset['h2'] == 0.1, 'int_h2'] = 0.05
target_sets_subset.loc[target_sets_subset['h2'] == 0.4, 'int_h2'] = 0
target_sets_subset.loc[target_sets_subset['h2'] == 0.8, 'int_h2'] = 0

# Create new variable to name the facets how I want 
target_sets_subset['pheno_format'] = target_sets_subset['num_causal_snps'].astype(str) + "/" + target_sets_subset['add_h2'].astype(str) + "/" + target_sets_subset['int_h2'].astype(str)

target_sets_subset['pheno_format'] = pd.Categorical(target_sets_subset['pheno_format'], 
                                             categories=['100/0.4/0', '100/0.8/0', '10k/0.4/0', '10k/0.45/0.05', '10k/0.475/0.025', '10k/0.8/0'],
                                             ordered=True)

target_sets_subset['method'] = pd.Categorical(target_sets_subset['method'], 
                                             categories=['TL_MLP_FR', 'TL_MLP_FT', 'TL_CNN_FR', 'TL_CNN_FT', 'MLP', 'CNN', 'Bayes R', 'GBLUP'],
                                             ordered=True)

gg = (
  ggplot(target_sets_subset)
  + aes(x="method", y="Test_Correlation", color="Ancestry")
  + facet_grid("pheno_format ~ .")
  + labs(
    x="Method",
    y="Test Correlation",
    title="Additive Effects Only Phenotypes: Test Correlation",
    )
  + theme(
        text=element_text(size=12),
        axis_text_x=element_text(size=10, angle=20, hjust = 1),
        axis_text_y=element_text(size=10),
        axis_title_x=element_text(size=12),
        axis_title_y=element_text(size=12),
        legend_text=element_text(size=10),
        legend_title=element_text(size=12),
        plot_title=element_text(size=12),
        strip_text=element_text(size=12)
    )
    + geom_boxplot(outlier_size=1) 
    + scale_color_manual(values=color_mapping_filtered)
    )

    
gg.save(f'/scratch/duffme_gensim/Plots/linear_phenotype_target_correlation.png', width=8, height=11, dpi=600)

gg = (
  ggplot(target_sets_subset)
  + aes(x="method", y="Test_Correlation", color="Ancestry", group="Ancestry")
  + facet_grid("pheno_format ~ .")
  + labs(
    x="Method",
    y="Test Correlation",
    title="Additive Effects Only Phenotypes: Test Correlation Mean (with 95% Confidence Intervals)",
    )
  + theme(
        text=element_text(size=12),
        axis_text_x=element_text(size=10, angle=20, hjust = 1),
        axis_text_y=element_text(size=10),
        axis_title_x=element_text(size=12),
        axis_title_y=element_text(size=12),
        legend_text=element_text(size=10),
        legend_title=element_text(size=12),
        plot_title=element_text(size=11),
        strip_text=element_text(size=12)
    )
    + stat_summary(fun_data="mean_cl_normal", geom="errorbar", width=1, position=position_dodge(width=0.5)) # Add error bars for mean
    + stat_summary(fun_data="mean_cl_normal", geom="point", position=position_dodge(width=0.5), size=0.5) # Add mean points
    + scale_color_manual(values=color_mapping_filtered)
    )

gg.save(f'/scratch/duffme_gensim/Plots/linear_phenotype_target_correlation_mean_error_bounds.png', width=8, height=11, dpi=600)

############################################### Appendix Plots 3: Additive Effects Only STD  ############################################################
########################################################################################################################################################
########################################################################################################################################################
summary_statistics = target_sets.groupby(['h2', 'num_causal_snps', 'Ancestry', 'p-value', 'method'])['Test_Correlation'].describe().reset_index()

summary_df_subset=summary_statistics[
                                     (summary_statistics['h2']!=0.05) & 
                                     (summary_statistics['h2']!=0.1)
]


summary_df_subset['h2'] = summary_df_subset['h2'].cat.remove_unused_categories()
plotted_ancestries = summary_df_subset['Ancestry'].unique()

color_mapping = {
    'CEU': '#E69F00',  # Orange
    'YRI': '#56B4E9',  # Sky Blue
    'CHB': '#009E73',  # Bluish Green
    'TSI': '#D55E00'   # Vermilion
}

# Ensure that only plotted ancestries are included in the color mapping
color_mapping_filtered = {key: color_mapping[key] for key in plotted_ancestries if key in color_mapping}

summary_df_subset['method'] = summary_df_subset['method'].cat.rename_categories({
    'TL_MLP_Fine_Tune': 'TL_MLP_FT',
    'TL_CNN_Fine_Tune': 'TL_CNN_FT',
    'TL_MLP_Freeze': 'TL_MLP_FR',
    'TL_CNN_Freeze': 'TL_CNN_FR',
    'BayesR': 'Bayes R'
})

# Additive summary_df_subset column 
summary_df_subset['add_h2'] = None

# Apply conditions for assignment
summary_df_subset.loc[summary_df_subset['h2'] == 0.05, 'add_h2'] = 0.475
summary_df_subset.loc[summary_df_subset['h2'] == 0.1, 'add_h2'] = 0.45
summary_df_subset.loc[summary_df_subset['h2'] == 0.4, 'add_h2'] = 0.4
summary_df_subset.loc[summary_df_subset['h2'] == 0.8, 'add_h2'] = 0.8

# Interaction heritability column 
summary_df_subset['int_h2'] = None

# Apply conditions for assignment
summary_df_subset.loc[summary_df_subset['h2'] == 0.05, 'int_h2'] = 0.025
summary_df_subset.loc[summary_df_subset['h2'] == 0.1, 'int_h2'] = 0.05
summary_df_subset.loc[summary_df_subset['h2'] == 0.4, 'int_h2'] = 0
summary_df_subset.loc[summary_df_subset['h2'] == 0.8, 'int_h2'] = 0

# Create new variable to name the facets how I want 
summary_df_subset['pheno_format'] = summary_df_subset['num_causal_snps'].astype(str) + "/" + summary_df_subset['add_h2'].astype(str) + "/" + summary_df_subset['int_h2'].astype(str)

summary_df_subset['pheno_format'] = pd.Categorical(summary_df_subset['pheno_format'], 
                                             categories=['100/0.4/0', '100/0.8/0', '10k/0.4/0', '10k/0.45/0.05', '10k/0.475/0.025', '10k/0.8/0'],
                                             ordered=True)

summary_df_subset['method'] = pd.Categorical(summary_df_subset['method'], 
                                             categories=['TL_MLP_FR', 'TL_MLP_FT', 'TL_CNN_FR', 'TL_CNN_FT', 'MLP', 'CNN', 'Bayes R', 'GBLUP'],
                                             ordered=True)

gg = (
    ggplot(summary_df_subset)
    + aes(x="method", y="std", fill="Ancestry")
    + facet_grid("pheno_format ~ .")
    + labs(
        x="Method",
        y="Test Correlation Standard Deviation",
        title="Additive Effects Only Phenotypes: Test Correlation Standard Deviation",
    )
  + theme(
        text=element_text(size=12),
        axis_text_x=element_text(size=10, angle=20, hjust = 1),
        axis_text_y=element_text(size=10),
        axis_title_x=element_text(size=12),
        axis_title_y=element_text(size=12),
        legend_text=element_text(size=10),
        legend_title=element_text(size=12),
        plot_title=element_text(size=11),
        strip_text=element_text(size=12)
    )
    + geom_col(position=position_dodge(width=0.7))
    + scale_fill_manual(values=color_mapping_filtered)
)
    
gg.save(f'/scratch/duffme_gensim/Plots/linear_phenotype_target_correlation_std.png', width=8, height=11, dpi=600)

############################################### Appendix Plot 4-6: Additive & Interaction  #########################################################################
####################################################################################################################################################################
####################################################################################################################################################################

target_sets_subset = target_sets[
                               (target_sets['h2']!=0.4) &
                               (target_sets['h2']!=0.8) 
                               ]

target_sets_subset['Ancestry'] = target_sets_subset['Ancestry'].cat.remove_unused_categories()
target_sets_subset['h2'] = target_sets_subset['h2'].cat.remove_unused_categories()
target_sets_subset['num_causal_snps'] = target_sets_subset['num_causal_snps'].cat.remove_unused_categories()
target_sets_subset['method'] = target_sets_subset['method'].cat.remove_unused_categories()

target_sets_subset['method'] = target_sets_subset['method'].cat.rename_categories({
    'TL_MLP_Fine_Tune': 'TL_MLP_FT',
    'TL_CNN_Fine_Tune': 'TL_CNN_FT',
    'TL_MLP_Freeze': 'TL_MLP_FR',
    'TL_CNN_Freeze': 'TL_CNN_FR',
    'BayesR': 'Bayes R'
})

plotted_ancestries = target_sets_subset['Ancestry'].unique()

color_mapping = {
    'CEU': '#E69F00',  # Orange
    'YRI': '#56B4E9',  # Sky Blue
    'CHB': '#009E73',  # Bluish Green
    'TSI': '#D55E00'   # Vermilion
}

# Ensure that only plotted ancestries are included in the color mapping
color_mapping_filtered = {key: color_mapping[key] for key in plotted_ancestries if key in color_mapping}

# Additive summary_df_subset column 
target_sets_subset['add_h2'] = None

# Apply conditions for assignment
target_sets_subset.loc[target_sets_subset['h2'] == 0.05, 'add_h2'] = 0.475
target_sets_subset.loc[target_sets_subset['h2'] == 0.1, 'add_h2'] = 0.45
target_sets_subset.loc[target_sets_subset['h2'] == 0.4, 'add_h2'] = 0.4
target_sets_subset.loc[target_sets_subset['h2'] == 0.8, 'add_h2'] = 0.8

# Interaction heritability column 
target_sets_subset['int_h2'] = None

# Apply conditions for assignment
target_sets_subset.loc[target_sets_subset['h2'] == 0.05, 'int_h2'] = 0.025
target_sets_subset.loc[target_sets_subset['h2'] == 0.1, 'int_h2'] = 0.05
target_sets_subset.loc[target_sets_subset['h2'] == 0.4, 'int_h2'] = 0
target_sets_subset.loc[target_sets_subset['h2'] == 0.8, 'int_h2'] = 0

# Create new variable to name the facets how I want 
target_sets_subset['pheno_format'] = target_sets_subset['num_causal_snps'].astype(str) + "/" + target_sets_subset['add_h2'].astype(str) + "/" + target_sets_subset['int_h2'].astype(str)

target_sets_subset['pheno_format'] = pd.Categorical(target_sets_subset['pheno_format'], 
                                             categories=['1k/0.45/0.05', '1k/0.475/0.025', '10k/0.45/0.05', '10k/0.475/0.025'],
                                             ordered=True)

target_sets_subset['method'] = pd.Categorical(target_sets_subset['method'], 
                                             categories=['TL_MLP_FR', 'TL_MLP_FT', 'TL_CNN_FR', 'TL_CNN_FT', 'MLP', 'CNN', 'Bayes R', 'GBLUP'],
                                             ordered=True)

gg = (
  ggplot(target_sets_subset)
  + aes(x="method", y="Test_Correlation", color="Ancestry")
  + facet_grid("pheno_format ~ .")
  + labs(
    x="Method",
    y="Test Correlation",
    title="Additive and Interaction Effects Phenotypes: Test Correlation"
    )
  + theme(
        text=element_text(size=12),
        axis_text_x=element_text(size=10, angle=20, hjust = 1),
        axis_text_y=element_text(size=10),
        axis_title_x=element_text(size=12),
        axis_title_y=element_text(size=12),
        legend_text=element_text(size=10),
        legend_title=element_text(size=12),
        plot_title=element_text(size=12),
        strip_text=element_text(size=12)
    )
    + geom_boxplot(outlier_size=1) 
    + scale_color_manual(values=color_mapping_filtered)
    )

    
gg.save(f'/scratch/duffme_gensim/Plots/non_linear_phenotype_target_correlation.png', width=8, height=11, dpi=600)

gg = (
  ggplot(target_sets_subset)
  + aes(x="method", y="Test_Correlation", color="Ancestry", group="Ancestry")
  + facet_grid("pheno_format ~ .")
  + labs(
    x="Method",
    y="Test Correlation",
    title="Additive and Interaction Effects Phenotypes: ",
    subtitle="Test Correlation Mean (with 95% Confidence Intervals)"
    )
  + theme(
        text=element_text(size=12),
        axis_text_x=element_text(size=10, angle=20, hjust = 1),
        axis_text_y=element_text(size=10),
        axis_title_x=element_text(size=12),
        axis_title_y=element_text(size=12),
        legend_text=element_text(size=10),
        legend_title=element_text(size=12),
        plot_title=element_text(size=11),
        plot_subtitle = element_text(size = 10),
        strip_text=element_text(size=12)
    )
    + stat_summary(fun_data="mean_cl_normal", geom="errorbar", width=1, position=position_dodge(width=0.5)) # Add error bars for mean
    + stat_summary(fun_data="mean_cl_normal", geom="point", position=position_dodge(width=0.5), size=0.5) # Add mean points
    + scale_color_manual(values=color_mapping_filtered)
    )

gg.save(f'/scratch/duffme_gensim/Plots/non_linear_phenotype_target_correlation_mean_error_bounds.png', width=8, height=11, dpi=600)

############################################### Appendix Plot 6: Additive and Interaction STD  ############################################################
###########################################################################################################################################################
###########################################################################################################################################################
summary_statistics = target_sets.groupby(['h2', 'num_causal_snps', 'Ancestry', 'p-value', 'method'])['Test_Correlation'].describe().reset_index()

summary_df_subset=summary_statistics[
                                     (summary_statistics['h2']!=0.4) & 
                                     (summary_statistics['h2']!=0.8) & 
                                     (summary_statistics['num_causal_snps']!=100)
]

summary_df_subset['h2'] = summary_df_subset['h2'].cat.remove_unused_categories()
plotted_ancestries = summary_df_subset['Ancestry'].unique()

color_mapping = {
    'CEU': '#E69F00',  # Orange
    'YRI': '#56B4E9',  # Sky Blue
    'CHB': '#009E73',  # Bluish Green
    'TSI': '#D55E00'   # Vermilion
}

# Ensure that only plotted ancestries are included in the color mapping
color_mapping_filtered = {key: color_mapping[key] for key in plotted_ancestries if key in color_mapping}

summary_df_subset['method'] = summary_df_subset['method'].cat.rename_categories({
    'TL_MLP_Fine_Tune': 'TL_MLP_FT',
    'TL_CNN_Fine_Tune': 'TL_CNN_FT',
    'TL_MLP_Freeze': 'TL_MLP_FR',
    'TL_CNN_Freeze': 'TL_CNN_FR',
    'BayesR': 'Bayes R'
})

# Additive summary_df_subset column 
summary_df_subset['add_h2'] = None

# Apply conditions for assignment
summary_df_subset.loc[summary_df_subset['h2'] == 0.05, 'add_h2'] = 0.475
summary_df_subset.loc[summary_df_subset['h2'] == 0.1, 'add_h2'] = 0.45
summary_df_subset.loc[summary_df_subset['h2'] == 0.4, 'add_h2'] = 0.4
summary_df_subset.loc[summary_df_subset['h2'] == 0.8, 'add_h2'] = 0.8

# Interaction heritability column 
summary_df_subset['int_h2'] = None

# Apply conditions for assignment
summary_df_subset.loc[summary_df_subset['h2'] == 0.05, 'int_h2'] = 0.025
summary_df_subset.loc[summary_df_subset['h2'] == 0.1, 'int_h2'] = 0.05
summary_df_subset.loc[summary_df_subset['h2'] == 0.4, 'int_h2'] = 0
summary_df_subset.loc[summary_df_subset['h2'] == 0.8, 'int_h2'] = 0

# Create new variable to name the facets how I want 
summary_df_subset['pheno_format'] = summary_df_subset['num_causal_snps'].astype(str) + "/" + summary_df_subset['add_h2'].astype(str) + "/" + summary_df_subset['int_h2'].astype(str)

summary_df_subset['pheno_format'] = pd.Categorical(summary_df_subset['pheno_format'], 
                                             categories=['1k/0.45/0.05', '1k/0.475/0.025', '10k/0.45/0.05', '10k/0.475/0.025'],
                                             ordered=True)

summary_df_subset['method'] = pd.Categorical(summary_df_subset['method'], 
                                             categories=['TL_MLP_FR', 'TL_MLP_FT', 'TL_CNN_FR', 'TL_CNN_FT', 'MLP', 'CNN', 'Bayes R', 'GBLUP'],
                                             ordered=True)


gg = (
    ggplot(summary_df_subset)
    + aes(x="method", y="std", fill="Ancestry")
    + facet_grid("pheno_format ~ .")
    + labs(
        x="Method",
        y="Test Correlation Standard Deviation",
        title="Additive and Interaction Effects Phenotypes: Standard Deviation",
    )
  + theme(
        text=element_text(size=12),
        axis_text_x=element_text(size=10, angle=20, hjust = 1),
        axis_text_y=element_text(size=10),
        axis_title_x=element_text(size=12),
        axis_title_y=element_text(size=12),
        legend_text=element_text(size=10),
        legend_title=element_text(size=12),
        plot_title=element_text(size=11),
        strip_text=element_text(size=12)
    )
    + geom_col(position=position_dodge(width=0.7))
    + scale_fill_manual(values=color_mapping_filtered)
)
    
gg.save(f'/scratch/duffme_gensim/Plots/non_linear_phenotype_target_correlation_std.png', width=8, height=11, dpi=600)
