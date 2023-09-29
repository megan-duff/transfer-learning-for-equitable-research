#!/bin/env python3

import tensorflow as tf
import pandas_plink as ppl
from pandas_plink import read_plink1_bin
from tensorflow import keras
from keras import backend as K
from tensorflow.keras import layers, regularizers
from tensorflow.keras.layers import experimental
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

##########################################################

p = psutil.Process(os.getpid())

print("Memory Usage (before import):", p.memory_info().rss/1024/1024, "MB")
print("Process is ... {}".format(p))
print("Start code")

# Read in PLINK files
(train_bim, train_fam, train_bed) = ppl.read_plink('whole_genome_GBR_1000GP_hm3_Phase3_120k_train.bed')
(val_bim, val_fam, val_bed) = ppl.read_plink('whole_genome_GBR_1000GP_hm3_Phase3_120k_val.bed')
(test_bim, test_fam, test_bed) = ppl.read_plink('whole_genome_GBR_1000GP_hm3_Phase3_120k_test.bed')

# Determine the size of the genotype matrix
print("Define variables")
num_train_snps, num_train_samples = train_bed.shape
num_val_snps, num_val_samples = val_bed.shape
num_test_snps, num_test_samples = test_bed.shape

############# To read in training genotypes

print("Create empty data frames and define number of chunks")
sample_wide_train_genotypes = []  # Initialize the list to store chunked data
num_chunks = 100
samples_per_chunk = num_train_samples // num_chunks

print("Start to read in data")

# Initialize empty lists to store the chunked data
for i in range(num_chunks):  # Iterate up to num_chunks
  print(f"Start chunk {i}")
  # Calculate the start and end indices for the current chunk
  start_index = i * samples_per_chunk
  end_index = (i + 1) * samples_per_chunk
  # Subset the data for the current chunk
  print("Starting to read genotypes...")
  chunk_train_compute = train_bed[:, start_index:end_index].compute()
  print(f"Got genotypes of chunk {i}")
  # Transpose the chunked genotypes
  chunk_train_genotypes = np.transpose(chunk_train_compute)
  print(f"Transposed chunk {i}")
  # Append the chunked data to the list
  sample_wide_train_genotypes.append(chunk_train_genotypes)
  print(f"Appended chunk {i}")
  # Remove variables from memory
  del chunk_train_genotypes
  del chunk_train_compute
  print("-------------------")

# empty dataframe to combine subset dataframes
train_genotypes=np.empty((num_train_samples, num_train_snps), dtype=np.int8)
print("start to read in training genotypes")

for i in range(num_chunks):
    print(f"Start {i}")
    start_person=i*samples_per_chunk
    print(f"{start_person}")
    end_person=(i+1)*samples_per_chunk
    print(f"{end_person}")
    print(f"{sample_wide_train_genotypes[i].shape}")
    print(f"{train_genotypes[start_person:end_person,:].shape}")
    train_genotypes[start_person:end_person,:]=sample_wide_train_genotypes[i]

# Print the shape of the genotype matrix
print(train_genotypes.shape)

# Read in the training sets allele frequencies
afreq_file_path = "train_WG_af.afreq"
# Read the .afreq file into a DataFrame
af_df = pd.read_csv(afreq_file_path, delimiter='\t')
train_freq_vector = 1-np.array(af_df.iloc[:, 4])
scale_factor = np.array(np.sqrt(train_freq_vector * (1 - train_freq_vector) * 2))

print("Memory Usage (before scaling):", p.memory_info().rss/1024/1024, "MB") 

scaled_train_genotypes=np.empty((num_train_samples, num_train_snps))

for i in range(num_chunks):
      print(f"start batch {i}")
      start_person=i*samples_per_chunk
      end_person=(i+1)*samples_per_chunk
      numerator = train_genotypes[start_person:end_person,:] - train_freq_vector
      scaled_train_chunk_genotypes = np.divide(numerator, scale_factor)
      scaled_train_genotypes[start_person:end_person,:] = scaled_train_chunk_genotypes

print("Mean of first SNP: " + str(np.mean(scaled_train_genotypes[:, 0])))
print("SD of first SNP: " + str(np.std(scaled_train_genotypes[:, 0])))

np.save("scaled_train_genotypes.npy", scaled_train_genotypes)
del train_genotypes
del scaled_train_genotypes
      
############ To read in validation genotypes 
print("Start validation genotypes...")
sample_wide_val_genotypes = []
num_chunks = 10
samples_per_chunk = num_val_samples // num_chunks

for i in range(num_chunks):  # Iterate up to num_chunks
  print(f"Start chunk {i}")
  # Calculate the start and end indices for the current chunk
  start_index = i * samples_per_chunk
  end_index = (i + 1) * samples_per_chunk
  # Subset the data for the current chunk
  print("Starting to read genotypes...")
  chunk_val_compute = val_bed[:, start_index:end_index].compute()
  print(f"Got genotypes of chunk {i}")
  # Transpose the chunked genotypes
  chunk_val_genotypes = np.transpose(chunk_val_compute)
  print(f"Transposed chunk {i}")
  # Append the chunked data to the list
  sample_wide_val_genotypes.append(chunk_val_genotypes)
  print(f"Appended chunk {i}")
  # Remove variables from memory
  del chunk_val_genotypes
  del chunk_val_compute
  print("-------------------")
  
val_genotypes=np.empty((num_val_samples, num_val_snps), dtype=np.int8)

for i in range(num_chunks):
  print(f"Start {i}")
  start_person=i*samples_per_chunk
  print(f"{start_person}")
  end_person=(i+1)*samples_per_chunk
  print(f"{end_person}")
  print(f"{sample_wide_val_genotypes[i].shape}")
  print(f"{val_genotypes[start_person:end_person,:].shape}")
  val_genotypes[start_person:end_person,:]=sample_wide_val_genotypes[i]
         
# Print the shape of the genotype matrix
print(val_genotypes.shape)

# Standardize genotype counts for val_genotypes
scaled_val_genotypes = np.divide(val_genotypes - 2 * train_freq_vector, scale_factor, out=np.zeros_like(val_genotypes), where=train_freq_vector != 0, casting='unsafe')
scaled_val_genotypes[np.isinf(scaled_val_genotypes)] = np.nan
print("VALIDATION: Mean of first SNP: " + str(np.mean(scaled_val_genotypes[:, 0])))
print("VALIDATION: SD of first SNP: " + str(np.std(scaled_val_genotypes[:, 0])))
np.save("scaled_val_genotypes.npy", scaled_val_genotypes)
    
############ To read in test genotypes 
print("Start test genotypes...")
sample_wide_test_genotypes = []
num_chunks = 10
samples_per_chunk = num_test_samples // num_chunks

for i in range(num_chunks):  # Iterate up to num_chunks
  print(f"Start chunk {i}")
  # Calculate the start and end indices for the current chunk
  start_index = i * samples_per_chunk
  end_index = (i + 1) * samples_per_chunk
  # Subset the data for the current chunk
  print("Starting to read genotypes...")
  chunk_test_compute = test_bed[:, start_index:end_index].compute()
  print(f"Got genotypes of chunk {i}")
  # Transpose the chunked genotypes
  chunk_test_genotypes = np.transpose(chunk_test_compute)
  print(f"Transposed chunk {i}")
  # Append the chunked data to the list
  sample_wide_test_genotypes.append(chunk_test_genotypes)
  print(f"Appended chunk {i}")
  # Remove variables from memory
  del chunk_test_genotypes
  del chunk_test_compute
  print("-------------------")
  
test_genotypes=np.empty((num_test_samples, num_test_snps), dtype=np.int8)
for i in range(num_chunks):
  print(f"Start {i}")
  start_person=i*samples_per_chunk
  print(f"{start_person}")
  end_person=(i+1)*samples_per_chunk
  print(f"{end_person}")
  print(f"{sample_wide_test_genotypes[i].shape}")
  print(f"{test_genotypes[start_person:end_person,:].shape}")
  test_genotypes[start_person:end_person,:]=sample_wide_test_genotypes[i]
         
# Print the shape of the genotype matrix
print(test_genotypes.shape)
# Standardize genotype counts for test_genotypes
scaled_test_genotypes = np.divide(test_genotypes - 2 * train_freq_vector, scale_factor, out=np.zeros_like(test_genotypes), where=train_freq_vector != 0, casting='unsafe')
scaled_test_genotypes[np.isinf(scaled_test_genotypes)] = np.nan
print("Mean of first SNP: " + str(np.mean(scaled_test_genotypes[:, 0])))
np.save("scaled_test_genotypes.npy", scaled_test_genotypes)
    
print("Finished scaling!")
print("---------------------------")

