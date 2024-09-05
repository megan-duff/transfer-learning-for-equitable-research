
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Define arguments
dataframes=[]
for ancestry in {'african', 'chinese', 'irish', 'british_rural'}:
    for pheno_name, pheno_code in [['standing_height', '50'], ['bmi', '21001'], ['diastolic_blood_pressure', '4079']]:

        ################################ Phenotype Data ###############################################
        # Read in phenotype files
        print(f"{ancestry}: {pheno_code}")
        train_pheno_path = f'/mnt/project/Data/phenotypes/{ancestry}_corrected_predicted_{pheno_name}_train.txt'
        val_pheno_path = f'/mnt/project/Data/phenotypes/{ancestry}_corrected_predicted_{pheno_name}_val.txt'
        test_pheno_path = f'/mnt/project/Data/phenotypes/{ancestry}_corrected_predicted_{pheno_name}_test.txt'

        if os.path.exists(train_pheno_path):
            train_pheno_file = pd.read_csv(train_pheno_path, sep='\s+')

        if os.path.exists(val_pheno_path):
            val_pheno_file = pd.read_csv(val_pheno_path, sep='\s+')

        if os.path.exists(test_pheno_path):
            test_pheno_file = pd.read_csv(test_pheno_path, sep='\s+')

        # Grab phenotype value
        train_pheno = train_pheno_file.iloc[:, 1].astype(float)
        val_pheno = val_pheno_file.iloc[:, 1].astype(float)
        test_pheno = test_pheno_file.iloc[:, 1].astype(float)
        Y_train = train_pheno.values
        Y_val = val_pheno.values
        Y_test = test_pheno.values

        ################################################ Bayes R & GBLUP scores ################################################
        ########################################################################################################################
        train_bayes_r_path = f'/mnt/project/Model_Output/BayesR/{ancestry}_pheno_{pheno_code}_bayes_r_clump_snps_1e-2/train_pred_{ancestry}_pheno_{pheno_code}_bayes_r_snps_clump_snps_1e-2.profile'
        val_bayes_r_path = f'/mnt/project/Model_Output/BayesR/{ancestry}_pheno_{pheno_code}_bayes_r_clump_snps_1e-2/val_pred_{ancestry}_pheno_{pheno_code}_bayes_r_snps_clump_snps_1e-2.profile'
        test_bayes_r_path = f'/mnt/project/Model_Output/BayesR/{ancestry}_pheno_{pheno_code}_bayes_r_clump_snps_1e-2/test_pred_{ancestry}_pheno_{pheno_code}_bayes_r_snps_clump_snps_1e-2.profile'

        train_gblup_path = f'/mnt/project/Model_Output/GBLUP/{ancestry}_pheno_{pheno_code}_gblup_clump_snps_1e-2/train_pred_{ancestry}_pheno_{pheno_code}_gblup_clump_snps_1e-2.profile'
        val_gblup_path = f'/mnt/project/Model_Output/GBLUP/{ancestry}_pheno_{pheno_code}_gblup_clump_snps_1e-2/val_pred_{ancestry}_pheno_{pheno_code}_gblup_clump_snps_1e-2.profile'
        test_gblup_path = f'/mnt/project/Model_Output/GBLUP/{ancestry}_pheno_{pheno_code}_gblup_clump_snps_1e-2/test_pred_{ancestry}_pheno_{pheno_code}_gblup_clump_snps_1e-2.profile'

        # Initialize variables
        train_bayes_r_scores = None
        val_bayes_r_scores = None
        test_bayes_r_scores = None

        train_gblup_scores = None
        val_gblup_scores = None
        test_gblup_scores = None

        # Read in the files if they exist
        if os.path.exists(train_bayes_r_path):
            train_bayes_r_scores = pd.read_csv(train_bayes_r_path, sep='\s+')

        if os.path.exists(val_bayes_r_path):
            val_bayes_r_scores = pd.read_csv(val_bayes_r_path, sep='\s+')

        if os.path.exists(test_bayes_r_path):
            test_bayes_r_scores = pd.read_csv(test_bayes_r_path, sep='\s+')

        if os.path.exists(train_gblup_path):
            train_gblup_scores = pd.read_csv(train_gblup_path, sep='\s+')

        if os.path.exists(val_gblup_path):
            val_gblup_scores = pd.read_csv(val_gblup_path, sep='\s+')

        if os.path.exists(test_gblup_path):
            test_gblup_scores = pd.read_csv(test_gblup_path, sep='\s+')


        ################################################ Order scores to match phenotypes ################################################
        ##################################################################################################################################

        if test_gblup_scores is not None:

            ############ Order Bayes R scores
            merged_train_bayes_r = train_pheno_file.merge(train_bayes_r_scores, on="IID")
            merged_train_bayes_r_scores = merged_train_bayes_r[['IID', 'SCORESUM']]

            merged_val_bayes_r = val_pheno_file.merge(val_bayes_r_scores, on="IID")
            merged_val_bayes_r_scores = merged_val_bayes_r[['IID', 'SCORESUM']]

            merged_test_bayes_r = test_pheno_file.merge(test_bayes_r_scores, on="IID")
            merged_test_bayes_r_scores = merged_test_bayes_r[['IID', 'SCORESUM']]

            ############ Order GBLUP scores
            merged_train_gblup = train_pheno_file.merge(train_gblup_scores, on="IID")
            merged_train_gblup_scores = merged_train_gblup[['IID', 'SCORESUM']]

            merged_val_gblup = val_pheno_file.merge(val_gblup_scores, on="IID")
            merged_val_gblup_scores = merged_val_gblup[['IID', 'SCORESUM']]

            merged_test_gblup = test_pheno_file.merge(test_gblup_scores, on="IID")
            merged_test_gblup_scores = merged_test_gblup[['IID', 'SCORESUM']]

            #### Grab scores
            bayes_r_train_score = merged_train_bayes_r_scores['SCORESUM'].values
            bayes_r_val_score = merged_val_bayes_r_scores['SCORESUM'].values
            bayes_r_test_score = merged_test_bayes_r_scores['SCORESUM'].values

            gblup_train_score = merged_train_gblup_scores['SCORESUM'].values
            gblup_val_score = merged_val_gblup_scores['SCORESUM'].values
            gblup_test_score = merged_test_gblup_scores['SCORESUM'].values
            
            print(merged_val_bayes_r.columns)
            print(merged_val_gblup.columns)

            Y_train_bayes = merged_train_bayes_r.iloc[:,1].values
            Y_val_bayes = merged_val_bayes_r.iloc[:,1].values
            Y_test_bayes = merged_test_bayes_r.iloc[:,1].values

            Y_train_gblup = merged_train_gblup.iloc[:,1].values
            Y_val_gblup = merged_val_gblup.iloc[:,1].values
            Y_test_gblup = merged_test_gblup.iloc[:,1].values
            
            
            #################### Correlation / MSE #############################
            # Compute MSE for train/val/test sets
            bayes_r_MSE_trn = np.mean((bayes_r_train_score - Y_train_bayes) ** 2)
            bayes_r_MSE_val = np.mean((bayes_r_val_score - Y_val_bayes) ** 2)
            bayes_r_MSE_tst = np.mean((bayes_r_test_score - Y_test_bayes) ** 2)

            gblup_MSE_trn = np.mean((gblup_train_score - Y_train_gblup) ** 2)
            gblup_MSE_val = np.mean((gblup_val_score - Y_val_gblup) ** 2)
            gblup_MSE_tst = np.mean((gblup_test_score - Y_test_gblup) ** 2)

            # Compute correlation for train/val/test sets
            bayes_r_COR_trn = np.corrcoef(bayes_r_train_score, Y_train_bayes)[0, 1]
            bayes_r_COR_val = np.corrcoef(bayes_r_val_score, Y_val_bayes)[0, 1]
            bayes_r_COR_tst = np.corrcoef(bayes_r_test_score, Y_test_bayes)[0, 1]

            gblup_COR_trn = np.corrcoef(gblup_train_score, Y_train_gblup)[0, 1]
            gblup_COR_val = np.corrcoef(gblup_val_score, Y_val_gblup)[0, 1]
            gblup_COR_tst = np.corrcoef(gblup_test_score, Y_test_gblup)[0, 1]

            # Compute R2 for train/val/test sets
            bayes_r_R2_trn = bayes_r_COR_trn**2
            bayes_r_R2_val = bayes_r_COR_val**2
            bayes_r_R2_tst = bayes_r_COR_tst**2

            gblup_R2_trn = gblup_COR_trn**2
            gblup_R2_val = gblup_COR_val**2
            gblup_R2_tst = gblup_COR_tst**2

            data = {
            'Ancestry': ancestry,
            'pheno_code': pheno_code,
            'pheno_name': pheno_name,
            'Method': ['Bayes_R', 'GBLUP'],
            'MSE_train': [bayes_r_MSE_trn, gblup_MSE_trn],
            'MSE_val': [bayes_r_MSE_val, gblup_MSE_val],
            'MSE_test': [bayes_r_MSE_tst, gblup_MSE_tst],
            'COR_train': [bayes_r_COR_trn, gblup_COR_trn],
            'COR_val': [bayes_r_COR_val, gblup_COR_val],
            'COR_test': [bayes_r_COR_tst, gblup_COR_tst],
            'R2_train': [bayes_r_R2_trn, gblup_R2_trn],
            'R2_val': [bayes_r_R2_val, gblup_R2_val],
            'R2_test': [bayes_r_R2_tst, gblup_R2_tst]
            }

            df = pd.DataFrame(data, index=['BayesR', 'GBLUP'])  
            dataframes.append(df)
        else:
            continue


combined_bayes_gblup = pd.concat(dataframes, ignore_index=True)

dataframes=[]
for ancestry in {'african', 'chinese', 'irish', 'british_rural'}:
    for pheno_name, pheno_code in [['standing_height', '50'], ['bmi', '21001'], ['diastolic_blood_pressure', '4079']]:
        nn_path = f'/mnt/project/Model_Output/Transfer_Learning_Output/{ancestry}_MLP_{pheno_code}_nn_results.csv'
        nn_df = pd.read_csv(nn_path, sep='\,')
        dataframes.append(nn_df)
        

combined_TL = pd.concat(dataframes, ignore_index=True)

# Write 'metrics_data' to a text file with tab-separated values
metrics_file_path = f"bayes_r_gblup_metrics.txt"

print("Writing metrics file...")
combined_bayes_gblup.to_csv(metrics_file_path, sep=",", index=False, header=False, quoting=False)

df = pd.DataFrame(combined_bayes_gblup)

# Set the style of the visualization
sns.set(style="whitegrid")

# Melt the DataFrame for easier plotting with seaborn
df_melted = df.melt(id_vars=['Ancestry', 'pheno_code', 'pheno_name', 'Model'],
                    value_vars=['MSE_train', 'MSE_val', 'MSE_test', 
                                'COR_train', 'COR_val', 'COR_test', 
                                'R2_train', 'R2_val', 'R2_test'],
                    var_name='Metric', value_name='Value')

# Create a FacetGrid to show the metrics for each model
g = sns.catplot(x='Metric', y='Value', hue='Model', data=df_melted, kind='bar', height=6, aspect=2)
g.set_axis_labels("Metric", "Value")
g.set_titles("Comparison of Models by Different Metrics")
plt.xticks(rotation=45)

# Display the plot
plt.show()



    # Print results
    print("Bayes R")
    print("-----------------------")
    print(f"Train correlation: {bayes_r_COR_trn}")
    print(f"Val correlation: {bayes_r_COR_val}")
    print(f"Test correlation: {bayes_r_COR_tst}")

    print(f"Train R2: {bayes_r_R2_trn}")
    print(f"Val R2: {bayes_r_R2_val}")
    print(f"Test R2: {bayes_r_R2_tst}")

    print("GBLUP")
    print("-----------------------")
    print(f"Train correlation: {gblup_COR_trn}")
    print(f"Val correlation: {gblup_COR_val}")
    print(f"Test correlation: {gblup_COR_tst}")

    print(f"Train R2: {gblup_R2_trn}")
    print(f"Val R2: {gblup_R2_val}")
    print(f"Test R2: {gblup_R2_tst}")

###################################### Create data frame to save aggregated results ######################################

metrics_data = {
    "method": ["bayes_r", "gblup"],
    "train_MSE": [bayes_r_MSE_trn, gblup_MSE_trn],
    "train_COR": [bayes_r_COR_trn, gblup_COR_trn],
    "train_R2": [bayes_r_R2_trn, gblup_R2_trn],
    "test_MSE": [bayes_r_MSE_tst, gblup_MSE_tst],
    "test_COR": [bayes_r_COR_tst, gblup_COR_tst],
    "test_R2": [bayes_r_R2_tst, gblup_R2_tst]
}

metrics_df = pd.DataFrame(metrics_data)

# Write 'metrics_data' to a text file with tab-separated values
metrics_file_path = f"{ancestry}_pheno_{pheno_code}_bayes_r_gblup_metrics.txt"

print("Writing metrics file...")
metrics_df.to_csv(metrics_file_path, sep=",", index=False, header=False, quoting=False)