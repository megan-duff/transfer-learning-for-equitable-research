#!/bin/env python
import numpy as np
import pandas as pd
import pandas_plink as ppl
import tensorflow as tf
from tensorflow.keras.optimizers import Adam, Adamax, SGD
from tensorflow.keras import layers, regularizers, Model, models
from tensorflow.keras.layers import experimental, Activation, Dense
from tensorflow.keras.models import load_model, clone_model
from tensorflow.keras import metrics

def get_genotypes(filtered_bed, num_chunks=100):
    """
    Reads and processes genotype data in chunks to manage memory usage.

    Args:
        filtered_bed (dask.array): Dask array of genotype data from PLINK bed file.
        num_chunks (int): Number of chunks to split the data into for processing. Default is 100.

    Returns:
        np.ndarray: A modified genotype matrix where alleles are flipped (0 <-> 2).
    """
    # Initialize variables
    num_snps, num_samples = filtered_bed.shape
    samples_per_chunk = num_samples // num_chunks
    
    # Initialize the list to store chunked data
    sample_wide_genotypes = []
    #print("Start to read in data")
    
    # Process genotype data in chunks
    for i in range(num_chunks):
        #print(f"Start chunk {i}")
        
        # Calculate the start and end indices for the current chunk
        start_index = i * samples_per_chunk
        end_index = (i + 1) * samples_per_chunk
        
        # Subset the data for the current chunk
        #print("Starting to read genotypes...")
        chunk_compute = filtered_bed[:, start_index:end_index].compute()
        #print(f"Got genotypes of chunk {i}")
        
        # Transpose the chunked genotypes to match the format (sample x SNPs)
        chunk_genotypes = np.transpose(chunk_compute)
        #print(f"Transposed chunk {i}")
        
        # Append the chunked data to the list
        sample_wide_genotypes.append(chunk_genotypes)
        #print(f"Appended chunk {i}")
        
        # Remove variables from memory to free up space
        del chunk_genotypes
        del chunk_compute
        #print("-------------------")
    
    # Combine all chunks and flip alleles
    genotypes = np.vstack(sample_wide_genotypes)
    genotypes_modified = np.where(genotypes == 0, 2, np.where(genotypes == 2, 0, genotypes))
    return genotypes_modified

def calculate_correlation_mse(y_true, y_pred):
    """
    Calculates the Pearson correlation coefficient and mean squared error (MSE) between true and predicted values.

    Args:
        y_true (np.ndarray): True phenotype values.
        y_pred (np.ndarray): Predicted phenotype values.

    Returns:
        tuple: Correlation coefficient and MSE.
    """
    # Calculate the Pearson correlation coefficient
    correlation = np.corrcoef(y_true, y_pred, rowvar=False)[0, 1]
    
    # Calculate the mean squared error
    mse = np.mean(np.square(y_true - y_pred))
    return correlation, mse

def subset_bed_bim_fam(snp_file, p_value_threshold, train_bim, train_bed, train_fam, val_bim, val_bed, val_fam, test_bim, test_bed, test_fam):
    """
    Subsets the PLINK files (bed, bim, fam) to a subset of SNPs based on GWAS p-value thresholds.

    Args:
        snp_file (str): Path to the file containing SNPs and their associated p-values.
        p_value_threshold (float): P-value threshold for SNP inclusion.
        train_bim, train_bed, train_fam, val_bim, val_bed, val_fam, test_bim, test_bed, test_fam (pd.DataFrame or dask.array): 
            PLINK bim, bed, and fam files for training, validation, and testing datasets.

    Returns:
        tuple: Subsets of bim, bed, and fam files for training, validation, and testing, 
               along with the number of SNPs and samples.
    """
    print("Reading in top SNPs found from GWAS analysis", flush=True)
    print("SNP File:" + snp_file, flush=True)

    # Read in SNPs filtered by p-value threshold
    dtype_dict = {0: str, 1: str, 7: str, 8: str, 9: str, 10: str, 11: str}
    snp_filter_dataset = pd.read_csv(snp_file, header=None, delimiter="\t", dtype=dtype_dict)
    
    # Drop the first row (with header) and remove rows with NaN values
    snp_filter_dataset = snp_filter_dataset[snp_filter_dataset.iloc[:, 11] != 'P']
    snp_filter_dataset = snp_filter_dataset.dropna()
    
    # Convert p-values to float type
    snp_filter_dataset.iloc[:, 11] = snp_filter_dataset.iloc[:, 11].astype(float)
    
    # Filter to SNPs with linear effects (ADD) and apply p-value threshold
    linear_rows = snp_filter_dataset[snp_filter_dataset.iloc[:, 6] == 'ADD']
    cut_off_rows = linear_rows[linear_rows.iloc[:, 11] < p_value_threshold]
    
    # Extract rs IDs of SNPs passing the p-value threshold
    snp_filter_rs_ids = cut_off_rows.iloc[:, 2].tolist()
    
    # Subset the bim file based on the selected SNPs
    matching_rows = train_bim[train_bim['snp'].isin(snp_filter_rs_ids)].index.tolist()
    print(f"Number of top SNPs identified in data set : {len(matching_rows)}", flush=True)
    
    # Subset the PLINK files based on the location of the selected SNPs
    filtered_train_bim = train_bim.iloc[matching_rows, :]
    filtered_train_bed = train_bed[matching_rows, :]
    filtered_train_fam = train_fam 
    filtered_val_bim = val_bim.iloc[matching_rows, :]
    filtered_val_bed = val_bed[matching_rows, :]
    filtered_val_fam = val_fam 
    filtered_test_bim = test_bim.iloc[matching_rows, :]
    filtered_test_bed = test_bed[matching_rows, :]
    filtered_test_fam = test_fam
    
    # Determine the size of the genotype matrix
    num_train_snps, num_train_samples = filtered_train_bed.shape
    num_val_snps, num_val_samples = filtered_val_bed.shape
    num_test_snps, num_test_samples = filtered_test_bed.shape
    return (
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
    )

def standardize_geno(train_af_file, matching_rows, train_genotypes_modified, val_genotypes_modified, test_genotypes_modified):    
    # Read in the training sets allele frequencies
    af_df = pd.read_csv(train_af_file, delimiter='\t')
    af_df_subset = af_df.iloc[matching_rows,:] 
    train_freq_vector = np.array(af_df_subset.iloc[:, 4])
    scale_factor = np.array(np.sqrt(train_freq_vector * (1 - train_freq_vector) * 2))
    # Create a mask for variants with MAF of 0 or 1 
    mask = (scale_factor == 0)
    # Standardize train genotype counts
    train_numerator = (train_genotypes_modified - 2 * train_freq_vector)
    scaled_train_genotypes = np.divide(train_numerator, scale_factor, out=np.zeros_like(train_numerator), where=~mask)
    scaled_train_genotypes[:,mask] = train_genotypes_modified[:,mask]
    # Standardize validation genotype counts
    val_numerator = (val_genotypes_modified - 2 * train_freq_vector)
    scaled_val_genotypes = np.divide(val_numerator, scale_factor, out=np.zeros_like(val_numerator), where=~mask)
    scaled_val_genotypes[:,mask] = val_genotypes_modified[:,mask]
    # Standardize validation genotype counts
    test_numerator = (test_genotypes_modified - 2 * train_freq_vector)
    scaled_test_genotypes = np.divide(test_numerator, scale_factor, out=np.zeros_like(test_numerator), where=~mask)
    scaled_test_genotypes[:,mask] = test_genotypes_modified[:,mask]
    return scaled_train_genotypes, scaled_val_genotypes, scaled_test_genotypes

def read_and_standardize_pheno(train_pheno_file_name, val_pheno_file_name, test_pheno_file_name):  
    # Read in phenotype files 
    # I had to add the appropriate header for GBR train files to perfrom GWAS 
    # that is why the format is different
    train_pheno_file = pd.read_csv(train_pheno_file_name, 
          delimiter='\t', 
          header=None)
    val_pheno_file = pd.read_csv(val_pheno_file_name, 
    delimiter='\t', 
    header=None)
    test_pheno_file = pd.read_csv(test_pheno_file_name, 
    delimiter='\t', 
    header=None)
    # Grab phenotype value
    train_pheno = train_pheno_file.iloc[:,2]
    val_pheno = val_pheno_file.iloc[:,2]
    test_pheno = test_pheno_file.iloc[:,2]
    #Standardize phenotype files
    scaled_train_pheno = (train_pheno - train_pheno.mean()) / train_pheno.std()
    scaled_val_pheno = (val_pheno - val_pheno.mean()) / val_pheno.std()
    scaled_test_pheno = (test_pheno - test_pheno.mean()) / test_pheno.std()
    return scaled_train_pheno, scaled_val_pheno, scaled_test_pheno


def freeze_transfer_learning_train(pre_trained_file, layer_count, pretrained_optimizer, X_train, Y_train, X_val, Y_val, X_test, Y_test, num_batch_size, num_epochs):
  # Load the pre-trained model
  pretrained_model = load_model(pre_trained_file, compile=False)
  # Remove prediction layer 
  pretrained_model.pop()
  pretrained_model.trainable = False
  # Create a prediction layer 
  predictions = Dense(1, activation='linear')(pretrained_model.output) 
  # Assemble transfer learning model 
  freeze_transfer_model = Model(inputs=[pretrained_model.input], outputs=[predictions])
  # Check architecture of transfer learning model 
  freeze_transfer_model.summary()  
  # Compile the transfer learning model with your optimizer, loss function, and metrics
  freeze_transfer_model.compile(optimizer=pretrained_optimizer, loss='mean_squared_error', metrics=['mean_squared_error', tf.keras.metrics.RootMeanSquaredError()])
    # Train TL Freeze Model 
  early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=.1*num_epochs, min_delta=0.001, verbose=1, restore_best_weights=True)
  freeze_transfer_model_history = freeze_transfer_model.fit(X_train, Y_train, validation_data=(X_val, Y_val), batch_size=num_batch_size, epochs=num_epochs, verbose=2, callbacks=[early_stopping])
  return freeze_transfer_model_history, freeze_transfer_model


def fine_tuning_transfer_learning_comparison(freeze_transfer_model, layer_count, pretrained_optimizer, X_train, Y_train, X_val, Y_val, X_test, Y_test, num_batch_size, num_epochs, learning_rate=1e-3):
  #learning_rate = freeze_transfer_model.optimizer.lr.numpy()
  print(f"Original learning rate: {freeze_transfer_model.optimizer.lr.numpy()}", flush=True)
  fine_tune_transfer_model = tf.keras.models.clone_model(freeze_transfer_model)
  for freeze_layer, fine_tune_layer in zip(freeze_transfer_model.layers, fine_tune_transfer_model.layers):
    fine_tune_layer.set_weights(freeze_layer.get_weights())
  fine_tune_transfer_model.summary()
  fine_tune_transfer_model.trainable=True
  # Define amount of fine-tuning
  for layer in fine_tune_transfer_model.layers[:-layer_count]:
        layer.trainable = False
  # Check architecture of transfer learning model 
  for layer in fine_tune_transfer_model.layers:
        print(f"Layer Name: {layer.name}, Trainable: {layer.trainable}")
  # Compile the transfer learning model with your optimizer, loss function, and metrics
  optimizers = {
    'adamax': Adamax,
    'adam': Adam,
    'sgd': SGD
  }
  print(pretrained_optimizer)
  selected_optimizer = optimizers.get(pretrained_optimizer.lower())
  print(selected_optimizer)
  optimizer_with_lr = selected_optimizer(learning_rate=learning_rate)
  fine_tune_transfer_model.compile(optimizer=optimizer_with_lr, loss='mean_squared_error', metrics=['mean_squared_error', tf.keras.metrics.RootMeanSquaredError()])
    # Train TL Freeze Model 
  early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=.10*num_epochs, min_delta=0.001, verbose=1, restore_best_weights=True)
  fine_tune_transfer_model_history = fine_tune_transfer_model.fit(X_train, Y_train, validation_data=(X_val, Y_val), batch_size=num_batch_size, epochs=num_epochs, verbose=2, callbacks=[early_stopping])
  Y_val_pred_NN = np.array(fine_tune_transfer_model.predict(X_val)).reshape(-1)
  val_correlation, val_mse = calculate_correlation_mse(Y_val, Y_val_pred_NN)
  return val_correlation

def fine_tuning_transfer_learning_train(freeze_transfer_model, layer_count, pretrained_optimizer, X_train, Y_train, X_val, Y_val, X_test, Y_test, num_batch_size, num_epochs, learning_rate=1e-3):
  #learning_rate = freeze_transfer_model.optimizer.lr.numpy() 
  fine_tune_transfer_model = tf.keras.models.clone_model(freeze_transfer_model)
  for freeze_layer, fine_tune_layer in zip(freeze_transfer_model.layers, fine_tune_transfer_model.layers):
    fine_tune_layer.set_weights(freeze_layer.get_weights())
  fine_tune_transfer_model.summary()
  fine_tune_transfer_model.trainable=True
  # Define amount of fine-tuning
  for layer in fine_tune_transfer_model.layers[:-layer_count]:
        layer.trainable = False
  # Check architecture of transfer learning model 
  for layer in fine_tune_transfer_model.layers:
        print(f"Layer Name: {layer.name}, Trainable: {layer.trainable}")
  # Compile the transfer learning model with your optimizer, loss function, and metrics
  optimizers = {
    'adamax': Adamax,
    'adam': Adam,
    'sgd': SGD
  }
  new_learning_rate=learning_rate
  selected_optimizer = optimizers.get(pretrained_optimizer.lower())
  optimizer_with_lr = selected_optimizer(learning_rate=new_learning_rate)
  fine_tune_transfer_model.compile(optimizer=optimizer_with_lr, loss='mean_squared_error', metrics=['mean_squared_error', tf.keras.metrics.RootMeanSquaredError()])
  # Train TL Freeze Model
  early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=.10*num_epochs, min_delta=0.001, verbose=1, restore_best_weights=True)
  fine_tune_transfer_model_history = fine_tune_transfer_model.fit(X_train, Y_train, validation_data=(X_val, Y_val), batch_size=num_batch_size, epochs=num_epochs, verbose=2, callbacks=[early_stopping])
  return fine_tune_transfer_model_history, fine_tune_transfer_model

def standardize_geno(train_af_file, matching_rows, train_genotypes_modified, val_genotypes_modified, test_genotypes_modified):
    """
    Standardizes genotype data using allele frequencies from the training set.

    Parameters:
    -----------
    train_af_file : str
        File path to the allele frequency data of the training set.
    matching_rows : list
        List of row indices to match variants between genotype data and allele frequencies.
    train_genotypes_modified : np.ndarray
        Modified genotype data for the training set.
    val_genotypes_modified : np.ndarray
        Modified genotype data for the validation set.
    test_genotypes_modified : np.ndarray
        Modified genotype data for the test set.

    Returns:
    --------
    scaled_train_genotypes : np.ndarray
        Standardized training genotype data.
    scaled_val_genotypes : np.ndarray
        Standardized validation genotype data.
    scaled_test_genotypes : np.ndarray
        Standardized test genotype data.
    """
    # Read in the training set's allele frequencies
    af_df = pd.read_csv(train_af_file, delimiter='\t')
    
    # Select the rows corresponding to matching variants
    af_df_subset = af_df.iloc[matching_rows, :]
    
    # Extract the frequency vector for the selected variants
    train_freq_vector = np.array(af_df_subset.iloc[:, 4])
    
    # Calculate the scale factor for standardization
    scale_factor = np.array(np.sqrt(train_freq_vector * (1 - train_freq_vector) * 2))
    
    # Create a mask for variants with minor allele frequency (MAF) of 0 or 1
    mask = (scale_factor == 0)
    
    # Standardize train genotype counts
    train_numerator = (train_genotypes_modified - 2 * train_freq_vector)
    scaled_train_genotypes = np.divide(train_numerator, scale_factor, out=np.zeros_like(train_numerator), where=~mask)
    scaled_train_genotypes[:, mask] = train_genotypes_modified[:, mask]
    
    # Standardize validation genotype counts
    val_numerator = (val_genotypes_modified - 2 * train_freq_vector)
    scaled_val_genotypes = np.divide(val_numerator, scale_factor, out=np.zeros_like(val_numerator), where=~mask)
    scaled_val_genotypes[:, mask] = val_genotypes_modified[:, mask]
    
    # Standardize test genotype counts
    test_numerator = (test_genotypes_modified - 2 * train_freq_vector)
    scaled_test_genotypes = np.divide(test_numerator, scale_factor, out=np.zeros_like(test_numerator), where=~mask)
    scaled_test_genotypes[:, mask] = test_genotypes_modified[:, mask]
    
    return scaled_train_genotypes, scaled_val_genotypes, scaled_test_genotypes


def read_and_standardize_pheno(train_pheno_file_name, val_pheno_file_name, test_pheno_file_name):
    """
    Reads and standardizes phenotype data.

    Parameters:
    -----------
    train_pheno_file_name : str
        File path to the training phenotype data.
    val_pheno_file_name : str
        File path to the validation phenotype data.
    test_pheno_file_name : str
        File path to the test phenotype data.

    Returns:
    --------
    scaled_train_pheno : pd.Series
        Standardized training phenotype data.
    scaled_val_pheno : pd.Series
        Standardized validation phenotype data.
    scaled_test_pheno : pd.Series
        Standardized test phenotype data.
    """
    # Read in phenotype files
    train_pheno_file = pd.read_csv(train_pheno_file_name, delimiter='\t', header=None)
    val_pheno_file = pd.read_csv(val_pheno_file_name, delimiter='\t', header=None)
    test_pheno_file = pd.read_csv(test_pheno_file_name, delimiter='\t', header=None)
    
    # Extract the phenotype values
    train_pheno = train_pheno_file.iloc[:, 2]
    val_pheno = val_pheno_file.iloc[:, 2]
    test_pheno = test_pheno_file.iloc[:, 2]
    
    # Standardize phenotype data
    scaled_train_pheno = (train_pheno - train_pheno.mean()) / train_pheno.std()
    scaled_val_pheno = (val_pheno - val_pheno.mean()) / val_pheno.std()
    scaled_test_pheno = (test_pheno - test_pheno.mean()) / test_pheno.std()
    
    return scaled_train_pheno, scaled_val_pheno, scaled_test_pheno


def freeze_transfer_learning_train(pre_trained_file, layer_count, pretrained_optimizer, X_train, Y_train, X_val, Y_val, X_test, Y_test, num_batch_size, num_epochs):
    """
    Trains a transfer learning model using a pre-trained model with frozen layers.

    Parameters:
    -----------
    pre_trained_file : str
        File path to the pre-trained model.
    layer_count : int
        Number of layers to freeze in the pre-trained model.
    pretrained_optimizer : str
        Name of the optimizer to use for training.
    X_train : np.ndarray
        Training data features.
    Y_train : np.ndarray
        Training data labels.
    X_val : np.ndarray
        Validation data features.
    Y_val : np.ndarray
        Validation data labels.
    X_test : np.ndarray
        Test data features.
    Y_test : np.ndarray
        Test data labels.
    num_batch_size : int
        Batch size for training.
    num_epochs : int
        Number of epochs to train the model.

    Returns:
    --------
    freeze_transfer_model_history : tf.keras.callbacks.History
        Training history of the transfer learning model.
    freeze_transfer_model : tf.keras.Model
        The trained transfer learning model.
    """
    # Load the pre-trained model
    pretrained_model = load_model(pre_trained_file, compile=False)
    
    # Remove the prediction layer from the pre-trained model
    pretrained_model.pop()
    
    # Freeze the layers of the pre-trained model
    pretrained_model.trainable = False
    
    # Create a new prediction layer
    predictions = Dense(1, activation='linear')(pretrained_model.output)
    
    # Assemble the transfer learning model
    freeze_transfer_model = Model(inputs=[pretrained_model.input], outputs=[predictions])
    
    # Check the architecture of the transfer learning model
    freeze_transfer_model.summary()
    
    # Compile the transfer learning model with the specified optimizer, loss function, and metrics
    freeze_transfer_model.compile(optimizer=pretrained_optimizer, loss='mean_squared_error', metrics=['mean_squared_error', tf.keras.metrics.RootMeanSquaredError()])
    
    # Train the transfer learning model
    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=int(0.1 * num_epochs), min_delta=0.001, verbose=1, restore_best_weights=True)
    freeze_transfer_model_history = freeze_transfer_model.fit(X_train, Y_train, validation_data=(X_val, Y_val), batch_size=num_batch_size, epochs=num_epochs, verbose=2, callbacks=[early_stopping])
    
    return freeze_transfer_model_history, freeze_transfer_model

def fine_tuning_transfer_learning_comparison(freeze_transfer_model, layer_count, pretrained_optimizer, X_train, Y_train, X_val, Y_val, X_test, Y_test, num_batch_size, num_epochs, learning_rate=1e-3):
    """
    Fine-tunes a pre-trained transfer learning model and compares its performance.

    Parameters:
    -----------
    freeze_transfer_model : tf.keras.Model
        The frozen transfer learning model to be fine-tuned.
    layer_count : int
        Number of layers to fine-tune in the pre-trained model.
    pretrained_optimizer : str
        Name of the optimizer used in the original frozen model.
    X_train : np.ndarray
        Training data features.
    Y_train : np.ndarray
        Training data labels.
    X_val : np.ndarray
        Validation data features.
    Y_val : np.ndarray
        Validation data labels.
    X_test : np.ndarray
        Test data features.
    Y_test : np.ndarray
        Test data labels.
    num_batch_size : int
        Batch size for training.
    num_epochs : int
        Number of epochs to train the model.
    learning_rate : float, optional
        Learning rate for fine-tuning (default is 1e-3).

    Returns:
    --------
    val_correlation : float
        Correlation between the predicted and actual values on the validation set.
    """
    # Print the original learning rate of the frozen model
    print(f"Original learning rate: {freeze_transfer_model.optimizer.lr.numpy()}", flush=True)
    
    # Clone the frozen model to create a new fine-tuning model
    fine_tune_transfer_model = tf.keras.models.clone_model(freeze_transfer_model)
    
    # Transfer weights from the frozen model to the fine-tuning model
    for freeze_layer, fine_tune_layer in zip(freeze_transfer_model.layers, fine_tune_transfer_model.layers):
        fine_tune_layer.set_weights(freeze_layer.get_weights())
    
    # Print the architecture of the fine-tuning model
    fine_tune_transfer_model.summary()
    
    # Set the entire model as trainable
    fine_tune_transfer_model.trainable = True
    
    # Freeze layers except for the last 'layer_count' layers
    for layer in fine_tune_transfer_model.layers[:-layer_count]:
        layer.trainable = False
    
    # Print the trainability of each layer
    for layer in fine_tune_transfer_model.layers:
        print(f"Layer Name: {layer.name}, Trainable: {layer.trainable}")
    
    # Select the optimizer based on the input string
    optimizers = {
        'adamax': Adamax,
        'adam': Adam,
        'sgd': SGD
    }
    selected_optimizer = optimizers.get(pretrained_optimizer.lower())
    
    # Compile the fine-tuning model with the new learning rate
    optimizer_with_lr = selected_optimizer(learning_rate=learning_rate)
    fine_tune_transfer_model.compile(optimizer=optimizer_with_lr, loss='mean_squared_error', metrics=['mean_squared_error', tf.keras.metrics.RootMeanSquaredError()])
    
    # Train the fine-tuning model with early stopping
    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=int(0.10 * num_epochs), min_delta=0.001, verbose=1, restore_best_weights=True)
    fine_tune_transfer_model_history = fine_tune_transfer_model.fit(X_train, Y_train, validation_data=(X_val, Y_val), batch_size=num_batch_size, epochs=num_epochs, verbose=2, callbacks=[early_stopping])
    
    # Predict on the validation set and calculate correlation and mean squared error
    Y_val_pred_NN = np.array(fine_tune_transfer_model.predict(X_val)).reshape(-1)
    val_correlation, val_mse = calculate_correlation_mse(Y_val, Y_val_pred_NN)
    
    return val_correlation


def fine_tuning_transfer_learning_train(freeze_transfer_model, layer_count, pretrained_optimizer, X_train, Y_train, X_val, Y_val, X_test, Y_test, num_batch_size, num_epochs, learning_rate=1e-3):
    """
    Fine-tunes a pre-trained transfer learning model and trains it on new data.

    Parameters:
    -----------
    freeze_transfer_model : tf.keras.Model
        The frozen transfer learning model to be fine-tuned.
    layer_count : int
        Number of layers to fine-tune in the pre-trained model.
    pretrained_optimizer : str
        Name of the optimizer used in the original frozen model.
    X_train : np.ndarray
        Training data features.
    Y_train : np.ndarray
        Training data labels.
    X_val : np.ndarray
        Validation data features.
    Y_val : np.ndarray
        Validation data labels.
    X_test : np.ndarray
        Test data features.
    Y_test : np.ndarray
        Test data labels.
    num_batch_size : int
        Batch size for training.
    num_epochs : int
        Number of epochs to train the model.
    learning_rate : float, optional
        Learning rate for fine-tuning (default is 1e-3).

    Returns:
    --------
    fine_tune_transfer_model_history : tf.keras.callbacks.History
        Training history of the fine-tuned transfer learning model.
    fine_tune_transfer_model : tf.keras.Model
        The trained fine-tuned transfer learning model.
    """
    # Clone the frozen model to create a new fine-tuning model
    fine_tune_transfer_model = tf.keras.models.clone_model(freeze_transfer_model)
    
    # Transfer weights from the frozen model to the fine-tuning model
    for freeze_layer, fine_tune_layer in zip(freeze_transfer_model.layers, fine_tune_transfer_model.layers):
        fine_tune_layer.set_weights(freeze_layer.get_weights())
    
    # Print the architecture of the fine-tuning model
    fine_tune_transfer_model.summary()
    
    # Set the entire model as trainable
    fine_tune_transfer_model.trainable = True
    
    # Freeze layers except for the last 'layer_count' layers
    for layer in fine_tune_transfer_model.layers[:-layer_count]:
        layer.trainable = False
    
    # Print the trainability of each layer
    for layer in fine_tune_transfer_model.layers:
        print(f"Layer Name: {layer.name}, Trainable: {layer.trainable}")
    
    # Select the optimizer based on the input string
    optimizers = {
        'adamax': Adamax,
        'adam': Adam,
        'sgd': SGD
    }
    
    # Compile the fine-tuning model with the new learning rate
    optimizer_with_lr = selected_optimizer(learning_rate=learning_rate)
    fine_tune_transfer_model.compile(optimizer=optimizer_with_lr, loss='mean_squared_error', metrics=['mean_squared_error', tf.keras.metrics.RootMeanSquaredError()])
    
    # Train the fine-tuning model with early stopping
    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=int(0.10 * num_epochs), min_delta=0.001, verbose=1, restore_best_weights=True)
    fine_tune_transfer_model_history = fine_tune_transfer_model.fit(X_train, Y_train, validation_data=(X_val, Y_val), batch_size=num_batch_size, epochs=num_epochs, verbose=2, callbacks=[early_stopping])
    
    return fine_tune_transfer_model_history, fine_tune_transfer_model
