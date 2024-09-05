import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

dataframes=[]
for ancestry in {'african', 'chinese', 'irish', 'british_rural'}:
    for pheno_name, pheno_code in [['standing_height', '50'], ['bmi', '21001'], ['diastolic_blood_pressure', '4079']]:
        nn_path = f'/mnt/project/Model_Output/Transfer_Learning_Output/attempt_1/{ancestry}_MLP_{pheno_code}_nn_results.csv'
        nn_df = pd.read_csv(nn_path, sep=',')
        nn_df['Ancestry']=ancestry
        nn_df['Pheno_Name']=pheno_name
        nn_df['Pheno_Code']=pheno_code
        dataframes.append(nn_df)
        
combined_TL = pd.concat(dataframes, ignore_index=True)


dataframes=[]
for ancestry in {'african', 'chinese', 'irish', 'british_rural'}:
    for pheno_name, pheno_code in [['standing_height', '50'], ['bmi', '21001'], ['diastolic_blood_pressure', '4079']]:
        nn_path = f'/mnt/project/Model_Output/Transfer_Learning_Output/p_1e-2/{ancestry}_MLP_{pheno_code}_nn_results.csv'
        nn_df = pd.read_csv(nn_path, sep=',')
        nn_df['Ancestry']=ancestry
        nn_df['Pheno_Name']=pheno_name
        nn_df['Pheno_Code']=pheno_code
        dataframes.append(nn_df)
        
combined_TL_1eneg2 = pd.concat(dataframes, ignore_index=True)

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

        if test_gblup_scores is not None and train_bayes_r_path is not None:

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

for nn_type in ['MLP', 'CNN']:
    for ancestry in ['chinese', 'african' ,'irish', 'british_rural']:
        for pheno_name, pheno_code in [['standing_height', '50'], ['diastolic_blood_pressure', '4079']]:
            nn_path = f'/mnt/project/Model_Output/Target_{nn_type}/{ancestry}_MLP_{pheno_code}_nn_results.csv'
            if os.path.exists(nn_path):
                nn_df = pd.read_csv(nn_path, sep=',')
                nn_df['Pheno_Name']=pheno_name
                nn_df['Method']=f"Target_{nn_type}"
                nn_df['Ancestry']=ancestry
                dataframes.append(nn_df)
            else:
                continue
        
combined_target_nn = pd.concat(dataframes, ignore_index=True)

combined_bayes_gblup.columns=['Ancestry','Pheno_Code',
       'Pheno_Name', 'Method', 'Train_MSE', 'Validation_MSE', 'Test_MSE','Train_Correlation', 'Validation_Correlation',
       'Test_Correlation', 'Train_R2', 'Validation_R2', 'Test_R2']

combined_bayes_gblup_subset=combined_bayes_gblup[['Ancestry',
       'Pheno_Name', 'Method', 'Train_MSE', 'Validation_MSE', 'Test_MSE','Train_Correlation', 'Validation_Correlation',
       'Test_Correlation', 'Train_R2', 'Validation_R2', 'Test_R2']]

combined_bayes_gblup_subset['p_value'] = np.nan

combined_TL_subset=combined_TL[['Ancestry',
       'Pheno_Name', 'Method', 'Train_MSE', 'Validation_MSE', 'Test_MSE','Train_Correlation', 'Validation_Correlation',
       'Test_Correlation', 'Train_R2', 'Validation_R2', 'Test_R2']]

combined_TL_subset['p_value'] = "1e-4"

combined_target_nn_subset=combined_target_nn[['Ancestry',
       'Pheno_Name', 'Method', 'Train_MSE', 'Validation_MSE', 'Test_MSE','Train_Correlation', 'Validation_Correlation',
       'Test_Correlation', 'Train_R2', 'Validation_R2', 'Test_R2']]

combined_target_nn_subset['p_value'] = np.nan

combined_TL_1eneg2=combined_TL_1eneg2[['Ancestry',
       'Pheno_Name', 'Method', 'Train_MSE', 'Validation_MSE', 'Test_MSE','Train_Correlation', 'Validation_Correlation',
       'Test_Correlation', 'Train_R2', 'Validation_R2', 'Test_R2']]

combined_TL_1eneg2['p_value'] = "1e-2"

combined_sets = pd.concat([combined_TL_subset, combined_TL_1eneg2, combined_bayes_gblup_subset, combined_target_nn_subset])

import matplotlib.pyplot as plt

# Initialize the figure and axes
fig, axs = plt.subplots(1, len(combined_sets['Pheno_Name'].unique()), figsize=(18, 6), sharey=True)

# Get unique phenotypes
phenotypes = combined_sets['Pheno_Name'].unique()

# Loop over each phenotype
for i, phenotype in enumerate(phenotypes):
    # Filter data for the current phenotype
    phenotype_data = combined_sets[combined_sets['Pheno_Name'] == phenotype]
    
    # Create scatter plot for each method
    for ancestry, color, label in zip(['irish', 'african', 'chinese', 'british_rural'],
                                      ['blue', 'red', 'green', 'orange'],
                                      ['Irish', 'African', 'Chinese', 'British Rural']):
        axs[i].scatter(phenotype_data.loc[phenotype_data['Ancestry'] == ancestry, 'Method'],
                       phenotype_data.loc[phenotype_data['Ancestry'] == ancestry, 'Test_Correlation'],
                       label=label, marker='o', color=color)
    
    # Set titles and labels
    axs[i].set_title(f'Phenotype: {phenotype}')
    axs[i].set_xlabel('Method')
    axs[i].set_ylabel('Test Correlation')
    
    # Add legend
    axs[i].legend()

# Adjust layout
fig.suptitle('Transfer Learning Performance', fontsize=16)
fig.tight_layout()

# Save the plot
plt.savefig('TL_performance_check_scatter.png', dpi=600)

# Show the plot (optional)
plt.show()




dataframes=[]
for pheno_code in {'50', '4079', '21001'}:
    source = pd.read_csv(f"/mnt/project/Genotype_Data/british_urban/pheno_set_files/british_urban_pheno_{pheno_code}_val.bim", sep="\t", header=None)
    source_rs = source.iloc[:,1]
    for ancestry in {'african', 'chinese', 'irish', 'british_rural'}:
        target_path=f"/mnt/project/Model_Output/BayesR/{ancestry}_pheno_{pheno_code}_bayes_r_clump_snps_1e-2/{ancestry}_pheno_{pheno_code}_bayes_r_snps_clump_snps_1e-2.snpRes"
        if os.path.exists(target_path):
            target = pd.read_csv(target_path, sep='\s+')
            matching_count = target[target.iloc[:,1].isin(source.iloc[:,1])]
            percent = matching_count.shape[0]/target.shape[0]
            df = pd.DataFrame({
            'Pheno_Code': [pheno_code],
            'Ancestry': [ancestry],
            'Percent': [percent*100]
            })
            dataframes.append(df)


dataframes=[]
for pheno_code in {'50', '4079', '21001'}:
    source_list = pd.read_csv(f"source_{pheno_code}_input_snps.txt", header=None)
    for ancestry in {'chinese', 'african', 'irish', 'british_rural', 'british_urban'}:
        frq = pd.read_csv(f"{ancestry}_whole_genome_genotype_calls_train.frq", sep="\s+")
        subset_frq = frq[frq['SNP'].isin(source_list.iloc[:,0])]
        subset_frq['Ancestry']=ancestry
        subset_frq['Pheno_Code']=pheno_code
        dataframes.append(subset_frq)

maf = pd.concat(dataframes, ignore_index=True)   

import seaborn as sns

# Create a FacetGrid for separate violin plots with box plots overlaid by Pheno_Code
g = sns.FacetGrid(maf, row='Pheno_Code', hue='Ancestry', height=4, aspect=2, palette='Set3')

# Map the histogram plot
g.map(sns.histplot, 'MAF', bins=100, alpha=0.6)
# Add vertical line at MAF = 0.01
for ax in g.axes.ravel():
    ax.axvline(x=0.01, color='r', linestyle='--', linewidth=1)
# Add titles and labels
g.set_axis_labels('MAF', 'Count')
g.set_titles(row_template='{row_name}')
g.add_legend(title='Ancestry')

plt.subplots_adjust(top=0.9)
g.fig.suptitle('Histogram of MAF Faceted by Pheno_Code and Colored by Ancestry')

plt.show()