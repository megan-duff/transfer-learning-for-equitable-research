#!/bin/env python3

import tensorflow as tf
import pandas_plink as ppl
from pandas_plink import read_plink1_bin
from tensorflow import keras
from keras import backend as K
from tensorflow.keras import layers, regularizers
from tensorflow.keras.layers import experimental, InputLayer, Flatten, Conv1D
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
from functions import *
import csv
from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession


############################## Define variables ##############################
##############################################################################
nn_type="CNN"
#training set's bed/bim/fam file prefix
train_file_prefix=""
#validation set's bed/bim/fam file prefix
val_file_prefix=""
#test set's bed/bim/fam file prefix
test_file_prefix=""
# GWAS results file produced from plink
snp_file=""
#p-value threshold to subset SNPs based on their GWAS p-values
p_value_threshold=1e-04
#training set's allele frequency files as produced by plink
train_af_file = f"/scratch/duffme_gensim/Simulations/GBR/sim_1/train_WG_af.afreq"
#training set's phenotype file
train_pheno_file_name=""
#validation set's phenotype file
val_pheno_file_name=""
#test set's phenotype file
test_pheno_file_name=""
#directory to save neural network logs 
nn_directory_path=""
#project name to save neural network logs
nn_project_name=""
#number of epochs to be used in neural network training 
num_epochs=200
#number of batch size to be used in neural network training 
num_batch_size=64
#name of file to save model parameters
save_model_parameter_file=""
#name of file to save model after hp optimization 
save_model_file=""
#name of file to save loss plot (.png file)
loss_plot_figure_file=""
#name of file to save the trained best model (.keras file)
save_trained_best_model_file=""
#name of file to save metrics based on the model (.txt file)
csv_file = ""


################################ Start time tracker ########################################################
timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

start_time = datetime.datetime.now()

###################################################### Start Script ######################################################
timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

start_time = datetime.datetime.now()

########################## Read in PLINK files ##########################################
print("Time before data set reading:" + str(datetime.datetime.now()), flush=True)

(train_bim, train_fam, train_bed) = ppl.read_plink(train_file_prefix)
print("Finished reading in train sets")
(val_bim, val_fam, val_bed) = ppl.read_plink(val_file_prefix)
print("Finished reading in validation sets")
(test_bim, test_fam, test_bed) = ppl.read_plink(test_file_prefix)
print("Finished reading in test sets", flush = True)

print("Time after data set reading:" + str(datetime.datetime.now()), flush=True)

print("################### Load in GWAS analysis files to get top SNPs ##############################")

(
    filtered_train_bim,
    filtered_train_bed,
    filtered_train_fam,
    filtered_val_bim,
    filtered_val_bed,
    filtered_val_fam,
    filtered_test_bim,
    filtered_test_bed,
    filtered_test_fam,
    num_train_snps,
    num_train_samples,
    num_val_snps,
    num_val_samples,
    num_test_snps,
    num_test_samples,
    matching_rows
) = subset_bed_bim_fam(snp_file, p_value_threshold, train_bim, train_bed, train_fam, val_bim, val_bed, val_fam, test_bim, test_bed, test_fam)

print("matching_rows length: " + str(len(matching_rows)), flush=True)

print("################################ To read in training genotypes ######################################")

train_genotypes_modified = get_genotypes(filtered_train_bed)

# Print the shape of the genotype matrix
print("Train genotype shape:" + str(train_genotypes_modified.shape))

print("############################# To read in validation genotypes #######################################")
print("Start validation genotypes...")

val_genotypes_modified = get_genotypes(filtered_val_bed, 10)
         
# Print the shape of the genotype matrix
print("Validation genotype shape:" + str(val_genotypes_modified.shape))

print("################################ To read in test genotypes #######################################")

print("Start test genotypes...")

test_genotypes_modified = get_genotypes(filtered_test_bed, 10)

# Print the shape of the genotype matrix
print("Test genotype shape:" + str(test_genotypes_modified.shape))

####################################### Standardize Data #######################################
##################################################################################################

scaled_train_genotypes, scaled_val_genotypes, scaled_test_genotypes = standardize_geno(train_af_file, matching_rows, train_genotypes_modified, val_genotypes_modified, test_genotypes_modified)

################################ Phenotype Data ###############################################
##################################################################################################

scaled_train_pheno, scaled_val_pheno, scaled_test_pheno = read_and_standardize_pheno(train_pheno_file_name, val_pheno_file_name, test_pheno_file_name)

# Load the dataset and split into training and validation sets
X_train = scaled_train_genotypes
Y_train = scaled_train_pheno
X_val = scaled_val_genotypes
Y_val = scaled_val_pheno
X_test = scaled_test_genotypes
Y_test = scaled_test_pheno

# Check shape of variables 
print("X_train shape:", X_train.shape)
print("Y_train shape:", Y_train.shape)
print("X_val shape:", X_val.shape)
print("Y_val shape:", Y_val.shape)
print("X_test shape:", X_test.shape)
print("Y_test shape:", Y_test.shape)

# Reshape the input data to have the appropriate shape
X_train_reshaped = np.reshape(X_train, (X_train.shape[0], X_train.shape[1], 1))
Y_train_array = Y_train.to_numpy()
Y_train_reshaped = np.reshape(Y_train_array, (Y_train_array.shape[0], 1))

# Reshape the validation data as well
X_val_reshaped = np.reshape(X_val, (X_val.shape[0], X_val.shape[1], 1))
Y_val_array = Y_val.to_numpy()
Y_val_reshaped = np.reshape(Y_val_array, (Y_val_array.shape[0], 1))

# Reshape the validation data as well
X_test_reshaped = np.reshape(X_test, (X_test.shape[0], X_test.shape[1], 1))
Y_test_array = Y_test.to_numpy()
Y_test_reshaped = np.reshape(Y_test_array, (Y_test_array.shape[0], 1))

# Convert to tensors for memory efficiency
X_train_tensor = tf.convert_to_tensor(X_train_reshaped, dtype=tf.float32)
X_val_tensor = tf.convert_to_tensor(X_val_reshaped, dtype=tf.float32)
X_test_tensor = tf.convert_to_tensor(X_test_reshaped, dtype=tf.float32)

Y_train_tensor = tf.convert_to_tensor(Y_train_reshaped, dtype=tf.float32)
Y_val_tensor = tf.convert_to_tensor(Y_val_reshaped, dtype=tf.float32)
Y_test_tensor = tf.convert_to_tensor(Y_test_reshaped, dtype=tf.float32)

print("############################### Neural Network Optimization ###############################")

### Define the model building function

def build_CNN(hp):
    # Define model to be of sequential nature
    model = tf.keras.Sequential()
    # Create input layer with appropriate number of nodes
    model.add(InputLayer(input_shape=(X_train_tensor.shape[1], X_train_tensor.shape[2])))
    # Add convolutional and batch normalization layers
    for i in range(2):
        model.add(tf.keras.layers.Conv1D(filters=hp.Choice("filters_" + str(i), [5, 8, 16]),
                                         kernel_size=hp.Choice("conv_kernel_size_" + str(i), [4, 5, 8]),
                                         activation=hp.Choice('conv_activation_function', ["relu", "tanh", "linear"]),
                                         padding=hp.Choice("conv_padding_" + str(i), ["same"]),
                                         strides=hp.Choice("conv_stride_" + str(i), [1, 2, 3])))
        model.add(tf.keras.layers.BatchNormalization())
        model.add(tf.keras.layers.MaxPooling1D(pool_size=hp.Choice("pool_size_" + str(i), [2, 3, 5]),
                                                padding="same",
                                                strides=hp.Choice("pool_stride_" + str(i), [2, 3, 5])))
        model.add(layers.Dropout(hp.Float("conv_drop_rate", min_value=0.05, max_value=0.55, step=0.10)))
    model.add(Flatten())
    # Tune the number of layers and nodes in the hidden layers
    # Within each layer tune choice of activation function and strength of L2 regularization
    for i in range(hp.Int('fc_num_layers', min_value=1, max_value=3, step=1)):
        model.add(tf.keras.layers.Dense(units=hp.Int('fc_num_of_nodes', min_value=16, max_value=208, step=32),
                                        activation=hp.Choice('fc_activation_function', ["relu", "tanh", "linear"]),
                                        kernel_regularizer=keras.regularizers.l2(hp.Choice('fc_l2', [0.01, 0.15, 0.225, 0.3]))))
        model.add(tf.keras.layers.BatchNormalization())
        # Tune if dropout layer should be added and if so tune the drop out rate
        model.add(layers.Dropout(hp.Float("fc_drop_rate", min_value=0.05, max_value=0.55, step=0.10)))
    # Add the output layer
    model.add(tf.keras.layers.Dense(1, activation='linear'))
    # Compile the model with appropriate loss function and optimizer
    model.compile(optimizer=hp.Choice('optim', ['adam', 'adamax', 'sgd']),
                  loss='mean_squared_error',
                  metrics=['mean_squared_error', tf.keras.metrics.RootMeanSquaredError()])
    # Tune the learning rate
    learning_rate = hp.Choice('learning_rate', values=[1e-3, 1e-4])
    K.set_value(model.optimizer.learning_rate, learning_rate)
    return model

### Define the tuner
# Choose hyperband bc it is more efficient than randomsearch 
tuner = Hyperband(
  hypermodel=build_CNN,
  objective='val_mean_squared_error',
  max_epochs=10,
  hyperband_iterations=1,
  directory=nn_directory_path,
  project_name=nn_project_name,
  overwrite=True)

# made patience 10% of number of epochs
stop_early = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=num_epochs*0.10)

### Search for the best hyperparameters
tuner.search(X_train_tensor, Y_train_tensor, validation_data=(X_val_tensor, Y_val_tensor), batch_size=num_batch_size, epochs=num_epochs, verbose = 2, callbacks=[stop_early])

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
save_best_param=open(save_model_parameter_file, 'wb')
pickle.dump(best_hps, save_best_param)

print("Save best model...")
save_best_model=open(save_model_file, 'wb')
pickle.dump(best_model, save_best_model)

print("Building best model...")
# Building best model 
best_model_params = tuner.hypermodel.build(best_hps)

print("################################### Neural Network Training ########################################")
print("Fitting the best model...")
# Early Stopping callback
early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=.10*num_epochs, verbose=0, restore_best_weights=True)

## Note, to find the optimal epochs I am setting the value of epochs large and then using early stopping to determine optimal number for training 

# Train the best model with early stopping to find the optimal number of epochs for training
history = best_model_params.fit(X_train_tensor, Y_train_tensor, validation_data=(X_val_tensor, Y_val_tensor), batch_size=num_batch_size, epochs=num_epochs, verbose=0, callbacks=[early_stopping])

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
plt.legend()
plt.grid(True)
# Save the plot as an image file
plt.savefig(loss_plot_figure_file)
# Close the plot to free up resources (optional)
plt.close()

# Saving best model
print("Save best model...")
best_model_params.save(save_trained_best_model_file)

# View NN architecture of best model 
best_model_params.summary()

print("###################################### Neural Network Testing #####################################")

# Get predictions for X_train, X_val, and X_test and reshape to 1-D array
Y_train_pred_NN = best_model_params.predict(X_train_tensor).reshape(-1)
Y_val_pred_NN = best_model_params.predict(X_val_tensor).reshape(-1)
Y_test_pred_NN = best_model_params.predict(X_test_tensor).reshape(-1)

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
test_evaluate = best_model_params.evaluate(X_test_tensor, Y_test_tensor, return_dict=True)
print("Test evaluate: ", flush=True)
print(test_evaluate)

############################## SAVE RESULTS ##############################

# Define the data to be saved
data = [
    ["Train_Correlation", "Validation_Correlation", "Test_Correlation", "Train_MSE", "Validation_MSE", "Test_MSE", "Train_R2", "Validation_R2", "Test_R2"],
    [train_correlation_NN, val_correlation_NN, test_correlation_NN, train_mse_NN, val_mse_NN, test_mse_NN, train_r_squared, val_r_squared, test_r_squared]
]

# Write the data to the CSV file
with open(csv_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(data)

print(f"Results saved to {csv_file}")

############################## PRINT RESULTS ##############################
print("FOR NEURAL NETWORK:")
print(f"Training Set - Correlation: {train_correlation_NN}, MSE: {train_mse_NN}, R2: {train_r_squared}")
print(f"Validation Set - Correlation: {val_correlation_NN}, MSE: {val_mse_NN}, R2: {val_r_squared}")
print(f"Test Set - Correlation: {test_correlation_NN}, MSE: {test_mse_NN}, R2: {test_r_squared}")

# Record the end time
end_time = datetime.datetime.now()

# Calculate the time elapsed
time_elapsed = end_time - start_time

# Print the total time elapsed in seconds
print("Total time elapsed (seconds):", time_elapsed.total_seconds())
