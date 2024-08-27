#!/bin/env python3

import tensorflow as tf
import pandas_plink as ppl
from pandas_plink import read_plink1_bin
from tensorflow import keras
from keras import backend as K
from tensorflow.keras import layers, regularizers, Model, models
from tensorflow.keras.layers import experimental, Activation, Dense
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
from functions import *
import csv
import glob

############################## Define variables ##############################
##############################################################################
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
train_af_file = f""
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
#source model's best parameter object path (.obj)
source_model_best_param_file_path = ""
#source model's best pre-trained model path (.keras)
pre_trained_file_pattern = ""
#name of file to save loss plot (.png file) for TL Freeze Model
TL_freeze_loss_plot_figure_file=""
#name of file to save loss plot (.png file) for TL Fine_Tune Model
TL_fine_tune_loss_plot_figure_file=""
#name of file to save the trained TL freeze best model (.keras file)
save_trained_best_model_TL_freeze_file=""
#name of file to save the trained TL fine tune best model (.keras file)
save_trained_best_model_TL_fine_tune_file=""
#name of file to save metrics based on the models (.txt file)
csv_file = ""

# Set variables to track time
timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

############################## Read in target sets data ########################################################
################################################################################################################

print("Time before data set reading:" + str(datetime.datetime.now()))

# Read in PLINK files
(train_bim, train_fam, train_bed) = ppl.read_plink(train_file_prefix)
print("Finished reading in train sets")

(val_bim, val_fam, val_bed) = ppl.read_plink(val_file_prefix)
print("Finished reading in validation sets")

(test_bim, test_fam, test_bed) = ppl.read_plink(test_file_prefix)
print("Finished reading in test sets")

print("Time after data set reading:" + str(datetime.datetime.now()))


########################### Load in GWAS analysis files to get top SNPs & subset files #######################
##############################################################################################################

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

################################ To read in training genotypes #######################################
######################################################################################################

train_genotypes_modified = get_genotypes(filtered_train_bed, 16)

# Print the shape of the genotype matrix
print("Train genotype shape:" + str(train_genotypes_modified.shape))

############################# To read in validation genotypes #######################################
#####################################################################################################

print("Start validation genotypes...")

val_genotypes_modified = get_genotypes(filtered_val_bed, 2)
         
# Print the shape of the genotype matrix
print("Validation genotype shape:" + str(val_genotypes_modified.shape))

################################ To read in test genotypes ####################################### 
##################################################################################################

print("Start test genotypes...")

test_genotypes_modified = get_genotypes(filtered_test_bed, 2)
         
#Print the shape of the genotype matrix
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

# Convert to tensors for memory efficiency
X_train_tensor = tf.convert_to_tensor(X_train, dtype=tf.float32)
X_val_tensor = tf.convert_to_tensor(X_val, dtype=tf.float32)
X_test_tensor = tf.convert_to_tensor(X_test, dtype=tf.float32)

Y_train_tensor = tf.convert_to_tensor(Y_train, dtype=tf.float32)
Y_val_tensor = tf.convert_to_tensor(Y_val, dtype=tf.float32)
Y_test_tensor = tf.convert_to_tensor(Y_test, dtype=tf.float32)


#################################### Transfer Learning ##########################################
#################################################################################################

# Find the matching file(s)
matching_files = glob.glob(source_model_best_param_file_path)

file_path = matching_files[0]

# Open the file in binary read mode
with open(file_path, 'rb') as file:
    # Load the object from the file
    best_hps = pickle.load(file)

pretrained_optimizer = best_hps.get_config()['values']['optim']

#################################### Freeze Approach ####################################
freeze_start_time = datetime.datetime.now()

# Find the matching file(s)
matching_files = glob.glob(pre_trained_file_pattern)

model_path = matching_files[0]

freeze_transfer_model_history, freeze_transfer_model = freeze_transfer_learning_train(
    model_path,
    layer_count=0,
    pretrained_optimizer=pretrained_optimizer,
    X_train=X_train_tensor,
    Y_train=Y_train_tensor,
    X_val=X_val_tensor,
    Y_val=Y_val_tensor,
    X_test=X_test_tensor,
    Y_test=Y_test_tensor,
    num_batch_size=num_batch_size,
    num_epochs=num_epochs
)

freeze_end_time = datetime.datetime.now()
#################################### Fine-tune Approach ####################################

ft_start_time = datetime.datetime.now()

# Grab number of layers
layer_count = len(freeze_transfer_model.layers)

best_correlation=0 
best_layer=0

for layer_value in range(2,layer_count):
    test_correlation = fine_tuning_transfer_learning_comparison(
      freeze_transfer_model=freeze_transfer_model,
      layer_count=layer_value,
      pretrained_optimizer=pretrained_optimizer,
      X_train=X_train_tensor,
      Y_train=Y_train_tensor,
      X_val=X_val_tensor,
      Y_val=Y_val_tensor,
      X_test=X_test_tensor,
      Y_test=Y_test_tensor,
      num_batch_size=num_batch_size,
      num_epochs=num_epochs
    )
    print(f"Test correlation: {test_correlation}", flush=True)
    if test_correlation > best_correlation:
      best_correlation=test_correlation
      best_layer=layer_value

print(f"Best Layer for Fine Tuning: {best_layer}", flush=True)

best_fine_tune_model_history, best_fine_tune_model = fine_tuning_transfer_learning_train(
    freeze_transfer_model=freeze_transfer_model,
    layer_count=best_layer,
    pretrained_optimizer=pretrained_optimizer,
    X_train=X_train_tensor,
    Y_train=Y_train_tensor,
    X_val=X_val_tensor,
    Y_val=Y_val_tensor,
    X_test=X_test_tensor,
    Y_test=Y_test_tensor,
    num_batch_size=num_batch_size,
    num_epochs=num_epochs
)

ft_end_time = datetime.datetime.now()
############################ Plot Results ###################################

############################ Freeze #########################################
# Plot training vs. validation loss by epoch
print("Plotting results...")
# Extract loss and validation loss from the history object
train_loss = freeze_transfer_model_history.history['loss']
val_loss = freeze_transfer_model_history.history['val_loss']
# Get the number of epochs
num_epochs = len(train_loss)
# Plot training and validation loss over epochs
plt.figure(figsize=(10, 6))
plt.plot(range(1, num_epochs + 1), train_loss, label='Training Loss')
plt.plot(range(1, num_epochs + 1), val_loss, label='Validation Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.title('TL Freeze Training and Validation Loss over Epochs')
plt.legend()
plt.grid(True)
# Save the plot as an image file
plt.savefig(TL_freeze_loss_plot_figure_file)
# Close the plot to free up resources (optional)
plt.close()

###################### Saving best model ############################
print("Save best transfer learning freeze model...")
freeze_transfer_model.save(save_trained_best_model_TL_freeze_file)

############################ Fine Tune #########################################
# Plot training vs. validation loss by epoch
print("Plotting results...")
# Extract loss and validation loss from the history object
train_loss = best_fine_tune_model_history.history['loss']
val_loss = best_fine_tune_model_history.history['val_loss']
# Get the number of epochs
num_epochs = len(train_loss)
# Plot training and validation loss over epochs
plt.figure(figsize=(10, 6))
plt.plot(range(1, num_epochs + 1), train_loss, label='Training Loss')
plt.plot(range(1, num_epochs + 1), val_loss, label='Validation Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.title('TL Fine Tune Training and Validation Loss over Epochs')
plt.legend()
plt.grid(True)
# Save the plot as an image file
plt.savefig(TL_fine_tune_loss_plot_figure_file)
# Close the plot to free up resources (optional)
plt.close()

###################### Saving best model ############################
print("Save best transfer learning fine tune model...")
best_fine_tune_model.save(save_trained_best_model_TL_fine_tune_file)

###################################### Neural Network Testing #####################################
###################################################################################################

# Get predictions for X_train, X_val, and X_test and reshape to 1-D array
Y_train_pred_freeze_TL = freeze_transfer_model.predict(X_train_tensor).reshape(-1)
Y_val_pred_freeze_TL = freeze_transfer_model.predict(X_val_tensor).reshape(-1)
Y_test_pred_freeze_TL = freeze_transfer_model.predict(X_test_tensor).reshape(-1)

Y_train_pred_finetune_TL = best_fine_tune_model.predict(X_train_tensor).reshape(-1)
Y_val_pred_finetune_TL = best_fine_tune_model.predict(X_val_tensor).reshape(-1)
Y_test_pred_finetune_TL = best_fine_tune_model.predict(X_test_tensor).reshape(-1)

# Calculate correlation and MSE for training set
train_correlation_freeze_TL, train_mse_freeze_TL = calculate_correlation_mse(Y_train, Y_train_pred_freeze_TL)
train_r_squared_freeze_TL = train_correlation_freeze_TL**2

train_correlation_finetune_TL, train_mse_finetune_TL = calculate_correlation_mse(Y_train, Y_train_pred_finetune_TL)
train_r_squared_finetune_TL = train_r_squared_finetune_TL**2

# Calculate correlation and MSE for validation set
val_correlation_freeze_TL, val_mse_freeze_TL = calculate_correlation_mse(Y_val, Y_val_pred_freeze_TL)
val_r_squared_freeze_TL = val_correlation_freeze_TL**2

val_correlation_finetune_TL, val_mse_finetune_TL = calculate_correlation_mse(Y_val, Y_val_pred_finetune_TL)
val_r_squared_finetune_TL = val_correlation_finetune_TL**2

# Calculate correlation and MSE for test set
test_correlation_freeze_TL, test_mse_freeze_TL = calculate_correlation_mse(Y_test, Y_test_pred_freeze_TL)
test_r_squared_freeze_TL = test_correlation_freeze_TL**2

test_correlation_finetune_TL, test_mse_finetune_TL = calculate_correlation_mse(Y_test, Y_test_pred_finetune_TL)
test_r_squared_finetune_TL = test_correlation_finetune_TL**2

# Evaluate models 
freeze_test_evaluate = freeze_transfer_model.evaluate(X_test_tensor, Y_test_tensor, return_dict=True)
print("Freeze Test Evaluations:")
print(freeze_test_evaluate)

fine_tune_evaluate = best_fine_tune_model.evaluate(X_test_tensor, Y_test_tensor, return_dict=True)
print("Fine Tune Test Evaluations:")
print(fine_tune_evaluate)

############################## SAVE RESULTS ##############################
freeze_time = freeze_end_time - freeze_start_time
ft_time = ft_end_time - ft_start_time

# Define the data to be saved
data = [
    ["Approach","Train_Correlation", "Validation_Correlation", "Test_Correlation", "Train_MSE", "Validation_MSE", "Test_MSE", "Train_R2", "Validation_R2", "Test_R2", "p_value_threshold", "time"],
    ["Freeze", train_correlation_freeze_TL, val_correlation_freeze_TL, test_correlation_freeze_TL, train_mse_freeze_TL, val_mse_freeze_TL, test_mse_freeze_TL, train_r_squared_freeze_TL, val_r_squared_freeze_TL, test_r_squared_freeze_TL, p_value_threshold, freeze_time.total_seconds()],
    ["Fine_Tune", train_correlation_finetune_TL, val_correlation_finetune_TL, test_correlation_finetune_TL, train_mse_finetune_TL, val_mse_finetune_TL, test_mse_finetune_TL, train_r_squared_finetune_TL, val_r_squared_finetune_TL, test_r_squared_finetune_TL, p_value_threshold, ft_time.total_seconds()]
]

# Write the data to the CSV file
with open(csv_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(data)

print(f"Results saved to {csv_file}")

############################## PRINT RESULTS ##############################

print("FOR FREEZE APPROACH:")
print(f"Training Set - Correlation: {train_correlation_freeze_TL}, MSE: {train_mse_freeze_TL}, R2: {train_r_squared_freeze_TL}")
print(f"Validation Set - Correlation: {val_correlation_freeze_TL}, MSE: {val_mse_freeze_TL}, R2: {val_r_squared_freeze_TL}")
print(f"Test Set - Correlation: {test_correlation_freeze_TL}, MSE: {test_mse_freeze_TL}, R2: {test_r_squared_freeze_TL}")

print("FOR FINE-TUNE APPROACH:")
print(f"Training Set - Correlation: {train_correlation_finetune_TL}, MSE: {train_mse_finetune_TL}, R2: {train_r_squared_finetune_TL}")
print(f"Validation Set - Correlation: {val_correlation_finetune_TL}, MSE: {val_mse_finetune_TL}, R2: {val_r_squared_finetune_TL}")
print(f"Test Set - Correlation: {test_correlation_finetune_TL}, MSE: {test_mse_finetune_TL}, R2: {test_r_squared_finetune_TL}")
