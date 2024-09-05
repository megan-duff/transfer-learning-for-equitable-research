import tensorflow as tf
import pandas_plink as ppl
from pandas_plink import read_plink1_bin
from tensorflow import keras
from keras import backend as K
from tensorflow.keras import layers, regularizers, Model, models
from tensorflow.keras.layers import Activation, Dense, InputLayer
from tensorflow.keras.models import load_model, clone_model
import keras_tuner as kt
import pandas as pd
from sklearn.preprocessing import StandardScaler
import dask.array as da
import numpy as np
import math
from keras_tuner.tuners import RandomSearch
from keras_tuner.tuners import Hyperband
from tensorflow.keras.layers import LayerNormalization
import os, psutil
import datetime
from sklearn.model_selection import train_test_split
from tensorflow.keras.utils import plot_model
from tensorflow.keras import metrics
import pickle
import sys
import signal
import time
import matplotlib.pyplot as plt
import csv
import glob
from keras.optimizers import Adam, Adamax

############################################ DEFINE FUNCTIONS ##################################################################
  
# Calculate correlation 
def calculate_correlation_mse(y_true, y_pred):
    correlation = np.corrcoef(y_true, y_pred, rowvar=False)[0, 1]
    mse = np.mean(np.square(y_true - y_pred))
    return correlation, mse

############################################ BEGIN TARGET NN SCRIPT ##################################################################

## Define variables for neural network 

# 21001 -- BMI
# 4079 -- Diastolic Blood Pressure
# 50 -- Standing height  

num_batch_size = 64
num_epochs = 200
ancestry=sys.argv[1]
pheno_name=sys.argv[2]
pheno_code=sys.argv[3]

if pheno_name=="standing_height":
    pheno_name_adjust="Predicted_Height"
elif pheno_name=="bmi":
    pheno_name_adjust="BMI"
else:
    pheno_name_adjust=pheno_name

train_std_matrix_file = f"/mnt/project/Genotype_Data/{ancestry}/target_model_files/{ancestry}_whole_genome_maf_unrelated_pheno_{pheno_code}_1neg2_train.npy"
val_std_matrix_file = f"/mnt/project/Genotype_Data/{ancestry}/target_model_files/{ancestry}_whole_genome_maf_unrelated_pheno_{pheno_code}_1neg2_val.npy"
test_std_matrix_file = f"/mnt/project/Genotype_Data/{ancestry}/target_model_files/{ancestry}_whole_genome_maf_unrelated_pheno_{pheno_code}_1neg2_test.npy"
train_pheno_file = f"/mnt/project/Data/phenotypes/{ancestry}_corrected_predicted_{pheno_name}_train.txt"
val_pheno_file = f"/mnt/project/Data/phenotypes/{ancestry}_corrected_predicted_{pheno_name}_val.txt"
test_pheno_file = f"/mnt/project/Data/phenotypes/{ancestry}_corrected_predicted_{pheno_name}_test.txt"

timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

p = psutil.Process(os.getpid())

start_time = datetime.datetime.now()

print("Start time:" + str(datetime.datetime.now()), flush=True)

## Read in standardized genotype files for train/val/test sets

scaled_train_genotypes = np.load(train_std_matrix_file)
print("Finished reading in train set", flush = True)
scaled_val_genotypes = np.load(val_std_matrix_file)
print("Finished reading in validation set", flush = True)
scaled_test_genotypes = np.load(test_std_matrix_file)
print("Finished reading in test set", flush = True)

## Read in standardized phenotype files for train/val/test sets

scaled_train_phenotypes = pd.read_csv(train_pheno_file, sep = " ")
scaled_val_phenotypes = pd.read_csv(val_pheno_file, sep = " ")
scaled_test_phenotypes = pd.read_csv(test_pheno_file, sep = " ")

# Define pheno_name_adjust (this is for 4079 script specifically)
pheno_name_adjust="BMI"

## Ensure corrected phenotypes are in the correct order 
train_fam_file=pd.read_csv(f"/mnt/project/Genotype_Data/{ancestry}/target_model_files/{ancestry}_whole_genome_maf_unrelated_clump_snps_1e-2_train.fam", header=None, sep=" ")
train_fam_file.columns=['NA1', 'IID', 'NA2' ,'NA3', 'NA4', 'NA5']
merged_train_data = train_fam_file.merge(scaled_train_phenotypes, on="IID")
merged_train_phenotype_data = merged_train_data[['IID', f'Corrected_{pheno_name}']]

val_fam_file=pd.read_csv(f"/mnt/project/Genotype_Data/{ancestry}/target_model_files/{ancestry}_whole_genome_maf_unrelated_clump_snps_1e-2_val.fam", sep=" ", header=None)
val_fam_file.columns=['NA1', 'IID', 'NA2' ,'NA3', 'NA4', 'NA5']
merged_val_data = val_fam_file.merge(scaled_val_phenotypes, on="IID")
merged_val_phenotype_data = merged_val_data[['IID', f'Corrected_{pheno_name_adjust}']]

test_fam_file=pd.read_csv(f"/mnt/project/Genotype_Data/{ancestry}/target_model_files/{ancestry}_whole_genome_maf_unrelated_clump_snps_1e-2_test.fam", sep=" ", header=None)
test_fam_file.columns=['NA1', 'IID', 'NA2' ,'NA3', 'NA4', 'NA5']
merged_test_data = test_fam_file.merge(scaled_test_phenotypes, on="IID")
merged_test_phenotype_data = merged_test_data[['IID', f'Corrected_{pheno_name_adjust}']]

print("Time after data set reading:" + str(datetime.datetime.now()), flush=True)

# Remove individuals in genotype matrix not in phenotype file
matching_rows = train_fam_file[train_fam_file['IID'].isin(merged_train_phenotype_data['IID'])].index.tolist()
print(f"Number of Samples : {len(matching_rows)}", flush = True)
# Subset the PLINK files based on the location of the top k SNPs list
scaled_train_genotypes_subset = scaled_train_genotypes[matching_rows,:]
train_fam_file=train_fam_file.iloc[matching_rows,:]

matching_rows = val_fam_file[val_fam_file['IID'].isin(merged_val_phenotype_data['IID'])].index.tolist()
print(f"Number of Samples : {len(matching_rows)}", flush = True)
# Subset the PLINK files based on the location of the top k SNPs list
scaled_val_genotypes_subset = scaled_val_genotypes[matching_rows,:]
val_fam_file=val_fam_file.iloc[matching_rows,:]

matching_rows = test_fam_file[test_fam_file['IID'].isin(merged_test_phenotype_data['IID'])].index.tolist()
print(f"Number of Samples : {len(matching_rows)}", flush = True)
# Subset the PLINK files based on the location of the top k SNPs list
scaled_test_genotypes_subset = scaled_test_genotypes[matching_rows,:]
test_fam_file=test_fam_file.iloc[matching_rows,:]

# Rename the variables for neural network training! 
X_train = scaled_train_genotypes_subset
Y_train = pd.DataFrame(merged_train_phenotype_data.iloc[:,1].values.reshape(-1, 1))
X_val = scaled_val_genotypes_subset
Y_val = pd.DataFrame(merged_val_phenotype_data.iloc[:,1].values.reshape(-1, 1))
X_test = scaled_test_genotypes_subset
Y_test = pd.DataFrame(merged_test_phenotype_data.iloc[:,1].values.reshape(-1, 1))

# Check shape of variables 
print("X_train shape:", X_train.shape, flush=True)
print("Y_train shape:", Y_train.shape, flush=True)
print("X_val shape:", X_val.shape, flush=True)
print("Y_val shape:", Y_val.shape, flush=True)
print("X_test shape:", X_test.shape, flush=True)
print("Y_test shape:", Y_test.shape, flush=True)

# Save shape of X_train to use in NN model to specifiy the start number of nodes
num_input_nodes = X_train.shape[1]

# Ensure the X and Y dataset are of the same type
Y_train = tf.cast(Y_train, tf.float32)
Y_val = tf.cast(Y_val, tf.float32)
Y_test = tf.cast(Y_test, tf.float32)

# Convert nd.arrays to tensorflow dataset
train_dataset = tf.data.Dataset.from_tensor_slices((X_train, Y_train))
val_dataset = tf.data.Dataset.from_tensor_slices((X_val, Y_val))
test_dataset = tf.data.Dataset.from_tensor_slices((X_test, Y_test))

# Shuffle (only the training set) and batch the tensorflow datasets
SHUFFLE_BUFFER_SIZE = 100

train_dataset = train_dataset.shuffle(SHUFFLE_BUFFER_SIZE).batch(num_batch_size)
val_dataset = val_dataset.batch(num_batch_size)
test_dataset = test_dataset.batch(num_batch_size)

################################# Before begining - ensure genotypes and phenotypes were in the correct order #######################################
if np.array_equal(train_fam_file['IID'].values, merged_train_phenotype_data.iloc[:, 0].values):
    print("Train set is in correct order!")
else:
    print("Train Set Error!!")

if np.array_equal(val_fam_file['IID'].values, merged_val_phenotype_data.iloc[:, 0].values):
    print("Val set is in correct order!")
else:
    print("Val Set Error!!")

if np.array_equal(test_fam_file['IID'].values, merged_test_phenotype_data.iloc[:, 0].values):
    print("Test set is in correct order!")
else:
    print("Test Set Error!!")

################################ Neural Network Optimization ###############################################
############################################################################################################

print("############################### Neural Network Optimization ###############################")

### Define the model building function

def build_model_with_L2_and_dropout(hp):
  #Define model to be of sequential nature
  model = tf.keras.Sequential()
  # Create input layer with appropriate number of nodes
  model.add(InputLayer((num_input_nodes,)))
  # Tune the number of layers and nodes in the hidden layers
  # Within each layer tune choice of activation function and strength of L2 regularization 
  for i in range(hp.Int('num_layers', min_value=1, max_value=5, step=1)):
    model.add(layers.Dropout(hp.Choice("drop_rate", [0.2, 0.3, 0.5])))
    model.add(tf.keras.layers.Dense(units=hp.Choice('num_of_nodes', [16, 32, 64, 128]),
                                    activation=hp.Choice('activation_function', ["relu", "elu", "tanh", "sigmoid", "linear"]), 
                                    kernel_regularizer=keras.regularizers.l2(hp.Choice('l2', [0.01, 0.15, 0.225, 0.3]))))
    # Add dropout layer - tune the drop out rate 
    model.add(tf.keras.layers.BatchNormalization())
  #Add the output layer      
  model.add(layers.Dropout(hp.Choice("drop_rate", [0.2, 0.3, 0.5])))
  model.add(tf.keras.layers.Dense(1, activation='linear'))
  #Compile the model with appropriate loss function and optimizer
  optimizer = hp.Choice('optim', ['adam', 'adamax'])
  # Tune the learning rate
  learning_rate = 1e-4
  if optimizer == 'adam':
        optimizer = Adam(learning_rate=learning_rate)
  else:
        optimizer = Adamax(learning_rate=learning_rate)
  model.compile(optimizer=optimizer,
                  loss='mean_squared_error',
                  metrics=['mean_squared_error', tf.keras.metrics.RootMeanSquaredError()])
  return model

### Define the tuner
# Choose hyperband bc it is more efficient than randomsearch 
tuner = Hyperband(
  hypermodel=build_model_with_L2_and_dropout,
  objective='val_mean_squared_error',
  max_epochs=100,
  hyperband_iterations=1,
  project_name=f'{ancestry}_MLP_{pheno_code}',
  overwrite=True)

# made patience 10% of number of epochs  
stop_early = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=num_epochs*0.10)

print("Start hyperparameter search!", flush=True)

print("Memory Usage (before hp opt):", p.memory_info().rss/1024/1024, "MB")

### Search for the best hyperparameters
tuner.search(train_dataset, validation_data=val_dataset, batch_size=num_batch_size, verbose = 2, callbacks=[stop_early])

print("Memory Usage (after hp opt):", p.memory_info().rss/1024/1024, "MB")

print("Retrieving hyperparameters...")
# Get the best hyperparameters
tf.keras.backend.clear_session()
best_hps=tuner.get_best_hyperparameters(num_trials=1)[0]

print("The hyperparameter search is complete.")
# Print the hyperparameters and their values
print(best_hps.get_config()['values'])

print("Getting the best model")
# get the best model
best_model = tuner.get_best_models(1)[0]

# Saving the model object
print("Save best model parameters...")
save_best_param=open(f'{ancestry}_MLP_{pheno_code}_params.obj', 'wb')
pickle.dump(best_hps, save_best_param)

print("Building best model...")
# Building best model 
best_model_params = tuner.hypermodel.build(best_hps)

################################## Neural Network Training ###############################################
##########################################################################################################

print("################################### Neural Network Training ########################################")
print("Fitting the best model...")
# Early Stopping callback
early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=num_epochs*0.20, verbose=0, restore_best_weights=True)

## Note, to find the optimal epochs I am setting the value of epochs large and then using early stopping to determine optimal number for training 

# Train the best model with early stopping to find the optimal number of epochs for training
history = best_model_params.fit(train_dataset, validation_data=val_dataset, epochs=num_epochs, batch_size=num_batch_size, verbose=2, callbacks=[early_stopping])

print("Plotting results...")
# Extract loss and validation loss from the history object
train_loss = history.history['loss']
val_loss = history.history['val_loss']
# Get the number of epochs
num_epochs = len(train_loss)
# Plot training and validation loss over epochs
plt.figure(figsize=(10, 6))
plt.plot(range(1, num_epochs + 1), train_loss, label='Training Loss')
plt.plot(range(1, num_epochs + 1), val_loss, label='Validation Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.title('Training and Validation Loss over Epochs')
plt.suptitle(f'Phenotype: {pheno_code}')
plt.legend()
plt.grid(True)
# Save the plot as an image file
plt.savefig(f'{ancestry}_MLP_{pheno_code}_training_loss_plot.png')
# Close the plot to free up resources (optional)
plt.close()

# Saving best model
print("Save best model...")
best_model_params.save(f'{ancestry}_MLP_{pheno_code}_best_model_trained.keras')

# View NN architecture of best model 
best_model_params.summary()

################################ Neural Network Testing ####################################################
############################################################################################################

print("###################################### Neural Network Testing #####################################")

# Get predictions for X_train, X_val, and X_test and reshape to 1-D array
Y_train_pred_NN = best_model_params.predict(X_train).reshape(-1)
Y_val_pred_NN = best_model_params.predict(X_val).reshape(-1)
Y_test_pred_NN = best_model_params.predict(X_test).reshape(-1)

# Calculate correlation and MSE for training set
train_correlation_NN, train_mse_NN = calculate_correlation_mse(Y_train, Y_train_pred_NN)
train_r_squared = train_correlation_NN**2
# Calculate correlation and MSE for validation set
val_correlation_NN, val_mse_NN = calculate_correlation_mse(Y_val, Y_val_pred_NN)
val_r_squared = val_correlation_NN**2
# Calculate correlation and MSE for test set
test_correlation_NN, test_mse_NN = calculate_correlation_mse(Y_test, Y_test_pred_NN)
test_r_squared = test_correlation_NN**2

# Evaluate models
#test_evaluate = best_model_params.evaluate(X_test, return_dict=True)
#print("Test evaluate: ", flush=True)
#print(test_evaluate)

# Record end time 
end_time = datetime.datetime.now()

# Calculate the time elapsed
time_elapsed = end_time - start_time

################################ Save & Print Results ####################################################
##########################################################################################################

# Define the data to be saved
data = [
    ["Train_Correlation", "Validation_Correlation", "Test_Correlation", "Train_MSE", "Validation_MSE", "Test_MSE", "Train_R2", "Validation_R2", "Test_R2", "Time"],
    [train_correlation_NN, val_correlation_NN, test_correlation_NN, train_mse_NN, val_mse_NN, test_mse_NN, train_r_squared, val_r_squared, test_r_squared, time_elapsed]
]

# Specify the CSV file name
csv_file = f"{ancestry}_MLP_{pheno_code}_nn_results.csv"

# Write the data to the CSV file
with open(csv_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(data)

print(f"Results saved to {csv_file}")
#
print("FOR NEURAL NETWORK:")
print(f"Training Set - Correlation: {train_correlation_NN}, MSE: {train_mse_NN}, R2: {train_r_squared}")
print(f"Validation Set - Correlation: {val_correlation_NN}, MSE: {val_mse_NN}, R2: {val_r_squared}")
print(f"Test Set - Correlation: {test_correlation_NN}, MSE: {test_mse_NN}, R2: {test_r_squared}")

# Print the total time elapsed in seconds
print("Total time elapsed (seconds):", time_elapsed.total_seconds())