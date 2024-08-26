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

def build_model_with_L2_and_dropout(hp, input_shape):
    """Builds a Keras model with hyperparameters tuned for L2 regularization and dropout."""
    
    model = tf.keras.Sequential()
    model.add(InputLayer(input_shape=(input_shape,)))
    
    for i in range(hp.Int('num_layers', min_value=1, max_value=5, step=1)):
        model.add(tf.keras.layers.Dense(
            units=hp.Int('num_of_nodes', min_value=16, max_value=512, step=32),
            activation=hp.Choice('activation_function', ["relu", "leaky_relu", "elu", "tanh", "sigmoid", "hard_sigmoid", "softplus", "linear"]),
            kernel_regularizer=regularizers.l2(hp.Choice('l2', [0.01, 0.15, 0.225, 0.3]))
        ))
        if hp.Boolean("dropout"):
            model.add(Dropout(hp.Float("drop_rate", min_value=0.05, max_value=0.55, step=0.10)))
    
    model.add(tf.keras.layers.Dense(1, activation='linear'))
    
    model.compile(
        optimizer=hp.Choice('optim', ['adam', 'adamax']),
        loss='mean_squared_error',
        metrics=['mean_squared_error', tf.keras.metrics.RootMeanSquaredError()]
    )
    
    learning_rate = hp.Choice('learning_rate', values=[1e-2, 1e-3, 1e-4])
    tf.keras.backend.set_value(model.optimizer.learning_rate, learning_rate)
    
    return model

def setup_tuner(num_epochs, nn_directory_path, nn_project_name, input_shape):
    """Sets up the Hyperband tuner for hyperparameter optimization."""
    
    return Hyperband(
        hypermodel=lambda hp: build_model_with_L2_and_dropout(hp, input_shape),
        objective='val_mean_squared_error',
        max_epochs=num_epochs,
        hyperband_iterations=1,
        directory=nn_directory_path,
        project_name=nn_project_name
    )

def setup_early_stopping(num_epochs):
    """Sets up early stopping callback."""
    
    return EarlyStopping(monitor='val_loss', patience=int(num_epochs * 0.10))

def monitor_memory_usage(label):
    """Prints memory usage with a given label."""
    
    p = psutil.Process(os.getpid())
    print(f"Memory Usage ({label}): {p.memory_info().rss / 1024 / 1024:.2f} MB")

def search_best_hyperparameters(tuner, X_train_tensor, Y_train_tensor, X_val_tensor, Y_val_tensor, num_batch_size, num_epochs):
    """Searches for the best hyperparameters using the tuner."""
    
    stop_early = setup_early_stopping(num_epochs)
    monitor_memory_usage("before hp opt")
    
    tuner.search(
        X_train_tensor, Y_train_tensor,
        validation_data=(X_val_tensor, Y_val_tensor),
        batch_size=num_batch_size,
        epochs=num_epochs,
        verbose=1,
        callbacks=[stop_early]
    )
    
    monitor_memory_usage("after hp opt")

def save_best_hyperparameters_and_model(tuner, save_model_parameter_file, save_model_file):
    """Saves the best hyperparameters and model."""
    
    print("Retrieving hyperparameters...")
    tf.keras.backend.clear_session()
    best_hps = tuner.get_best_hyperparameters(num_trials=1)[0]
    
    print("The hyperparameter search is complete.")
    print(best_hps.get_config()['values'])
    
    print("Getting the best model")
    best_model = tuner.get_best_models(1)[0]
    
    print("Saving best model parameters...")
    with open(save_model_parameter_file, 'wb') as f:
        pickle.dump(best_hps, f)
    
    print("Saving best model...")
    with open(save_model_file, 'wb') as f:
        pickle.dump(best_model, f)
    
    print("Building best model...")
    return tuner.hypermodel.build(best_hps)

def setup_early_stopping(num_epochs):
    """Sets up early stopping callback with given patience."""
    return tf.keras.callbacks.EarlyStopping(
        monitor='val_loss',
        patience=int(0.10 * num_epochs),
        verbose=0,
        restore_best_weights=True
    )

def train_best_model(best_model_params, X_train_tensor, Y_train_tensor, X_val_tensor, Y_val_tensor, num_epochs, num_batch_size, early_stopping):
    """Trains the best model with early stopping to find the optimal number of epochs."""
    return best_model_params.fit(
        X_train_tensor,
        Y_train_tensor,
        validation_data=(X_val_tensor, Y_val_tensor),
        batch_size=num_batch_size,
        epochs=num_epochs,
        verbose=0,
        callbacks=[early_stopping]
    )

def plot_training_history(history, loss_plot_figure_file):
    """Plots training and validation loss over epochs and saves the plot to a file."""
    train_loss = history.history['loss']
    val_loss = history.history['val_loss']
    num_epochs = len(train_loss)
    
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, num_epochs + 1), train_loss, label='Training Loss')
    plt.plot(range(1, num_epochs + 1), val_loss, label='Validation Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('Training and Validation Loss over Epochs')
    plt.legend()
    plt.grid(True)
    plt.savefig(loss_plot_figure_file)
    plt.close()

def save_best_model(best_model_params, save_trained_best_model_file):
    """Saves the best trained model to a file."""
    best_model_params.save(save_trained_best_model_file)

def display_model_summary(best_model_params):
    """Prints the summary of the best model."""
    best_model_params.summary()


def get_predictions(model, X_train_tensor, X_val_tensor, X_test_tensor):
    """Get model predictions and reshape them to 1-D array."""
    Y_train_pred = model.predict(X_train_tensor).reshape(-1)
    Y_val_pred = model.predict(X_val_tensor).reshape(-1)
    Y_test_pred = model.predict(X_test_tensor).reshape(-1)
    return Y_train_pred, Y_val_pred, Y_test_pred

def calculate_metrics(Y_true, Y_pred):
    """Calculate correlation and MSE between true and predicted values."""
    correlation = np.corrcoef(Y_true, Y_pred)[0, 1]
    mse = np.mean((Y_true - Y_pred) ** 2)
    return correlation, mse

def evaluate_model(model, X_test_tensor, Y_test_tensor):
    """Evaluate the model on the test dataset."""
    return model.evaluate(X_test_tensor, Y_test_tensor, return_dict=True)

def save_results_to_csv(csv_file, train_correlation, val_correlation, test_correlation,
                        train_mse, val_mse, test_mse, train_r2, val_r2, test_r2):
    """Save evaluation results to a CSV file."""
    data = [
        ["Train_Correlation", "Validation_Correlation", "Test_Correlation", "Train_MSE", "Validation_MSE", "Test_MSE", "Train_R2", "Validation_R2", "Test_R2"],
        [train_correlation, val_correlation, test_correlation, train_mse, val_mse, test_mse, train_r2, val_r2, test_r2]
    ]
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(data)
    print(f"Results saved to {csv_file}")

def print_results(train_correlation, val_correlation, test_correlation,
                  train_mse, val_mse, test_mse, train_r2, val_r2, test_r2):
    """Print the evaluation results."""
    print("FOR NEURAL NETWORK:")
    print(f"Training Set - Correlation: {train_correlation}, MSE: {train_mse}, R2: {train_r2}")
    print(f"Validation Set - Correlation: {val_correlation}, MSE: {val_mse}, R2: {val_r2}")
    print(f"Test Set - Correlation: {test_correlation}, MSE: {test_mse}, R2: {test_r2}")
