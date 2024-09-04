#!/bin/env python3

import tensorflow as tf
import pandas_plink as ppl
from pandas_plink import read_plink1_bin
from tensorflow import keras
from keras import backend as K
from tensorflow.keras import layers, regularizers
from tensorflow.keras.layers import experimental, InputLayer, Dropout
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
import random

def set_seed(seed=None):
    if seed is None:
        seed = int(time.time())  # If seed is set to None (aka not set by the user) a random seed is set based on the current time
    np.random.seed(seed)
    tf.random.set_seed(seed)
    random.seed(seed)
    print(f"Seed set to: {seed}")
    return seed

def get_timestamp():
    return datetime.datetime.now().strftime("%Y%m%d%H%M%S")

def print_current_time(label):
    print(f"{label}: {str(datetime.datetime.now())}", flush=True)

def read_plink_files(train_file_prefix, val_file_prefix, test_file_prefix):
    print_current_time("Time before data set reading")
    
    train_bim, train_fam, train_bed = ppl.read_plink(train_file_prefix)
    print("Finished reading in train sets")
    
    val_bim, val_fam, val_bed = ppl.read_plink(val_file_prefix)
    print("Finished reading in validation sets")
    
    test_bim, test_fam, test_bed = ppl.read_plink(test_file_prefix)
    print("Finished reading in test sets", flush=True)
    
    print_current_time("Time after data set reading")
    
    return (train_bim, train_fam, train_bed), (val_bim, val_fam, val_bed), (test_bim, test_fam, test_bed)

def load_gwas_and_subset(snp_file, p_value_threshold, train_data, val_data, test_data):
    train_bim, train_fam, train_bed = train_data
    val_bim, val_fam, val_bed = val_data
    test_bim, test_fam, test_bed = test_data
    
    print("################### Load in GWAS analysis files to get top SNPs ##############################")
    
    return subset_bed_bim_fam(
        snp_file, p_value_threshold, train_bim, train_bed, train_fam, 
        val_bim, val_bed, val_fam, test_bim, test_bed, test_fam
    )

def read_genotypes(filtered_bed, label):
    print(f"################################ To read in {label} genotypes ######################################")
    print(f"Start {label} genotypes...")
    
    genotypes = get_genotypes(filtered_bed, 10)
    
    print(f"{label.capitalize()} genotype shape:" + str(genotypes.shape))
    
    return genotypes

def standardize_and_prepare_data(train_af_file, matching_rows, train_genotypes, val_genotypes, test_genotypes):
    print("-- Standardize Gentoype Data")
    
    return standardize_geno(train_af_file, matching_rows, train_genotypes, val_genotypes, test_genotypes)

def read_and_standardize_phenotypes(train_pheno_file_name, val_pheno_file_name, test_pheno_file_name):
    print("-- Standardize Phenotype Data")
    
    return read_and_standardize_pheno(train_pheno_file_name, val_pheno_file_name, test_pheno_file_name)

def prepare_tensors(X_train, Y_train, X_val, Y_val, X_test, Y_test):
    X_train_tensor = tf.convert_to_tensor(X_train, dtype=tf.float32)
    X_val_tensor = tf.convert_to_tensor(X_val, dtype=tf.float32)
    X_test_tensor = tf.convert_to_tensor(X_test, dtype=tf.float32)
    
    Y_train_tensor = tf.convert_to_tensor(Y_train, dtype=tf.float32)
    Y_val_tensor = tf.convert_to_tensor(Y_val, dtype=tf.float32)
    Y_test_tensor = tf.convert_to_tensor(Y_test, dtype=tf.float32)
    
    return X_train_tensor, Y_train_tensor, X_val_tensor, Y_val_tensor, X_test_tensor, Y_test_tensor

def load_best_hyperparameters(source_model_best_param_file_path):
    matching_files = glob.glob(source_model_best_param_file_path)
    file_path = matching_files[0]

    with open(file_path, 'rb') as file:
        best_hps = pickle.load(file)
    
    return best_hps.get_config()['values']['optim']

def load_pretrained_model(pre_trained_file_pattern):
    matching_files = glob.glob(pre_trained_file_pattern)
    return matching_files[0]

def plot_training_validation_loss(history, title, save_path):
    train_loss = history.history['loss']
    val_loss = history.history['val_loss']
    num_epochs = len(train_loss)

    plt.figure(figsize=(10, 6))
    plt.plot(range(1, num_epochs + 1), train_loss, label='Training Loss')
    plt.plot(range(1, num_epochs + 1), val_loss, label='Validation Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(save_path)
    plt.close()

def save_model(model, save_path):
    model.save(save_path)

def calculate_performance_metrics(model, X_train, X_val, X_test, Y_train, Y_val, Y_test):
    Y_train_pred = model.predict(X_train).reshape(-1)
    Y_val_pred = model.predict(X_val).reshape(-1)
    Y_test_pred = model.predict(X_test).reshape(-1)

    train_correlation, train_mse = calculate_correlation_mse(Y_train, Y_train_pred)
    train_r_squared = train_correlation ** 2

    val_correlation, val_mse = calculate_correlation_mse(Y_val, Y_val_pred)
    val_r_squared = val_correlation ** 2

    test_correlation, test_mse = calculate_correlation_mse(Y_test, Y_test_pred)
    test_r_squared = test_correlation ** 2

    return {
        "train": (train_correlation, train_mse, train_r_squared),
        "val": (val_correlation, val_mse, val_r_squared),
        "test": (test_correlation, test_mse, test_r_squared)
    }

def save_results_to_csv(data, csv_file):
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(data)

def print_results(approach, metrics):
    print(f"FOR {approach.upper()} APPROACH:")
    print(f"Training Set - Correlation: {metrics['train'][0]}, MSE: {metrics['train'][1]}, R2: {metrics['train'][2]}")
    print(f"Validation Set - Correlation: {metrics['val'][0]}, MSE: {metrics['val'][1]}, R2: {metrics['val'][2]}")
    print(f"Test Set - Correlation: {metrics['test'][0]}, MSE: {metrics['test'][1]}, R2: {metrics['test'][2]}")
