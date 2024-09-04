#!/bin/env python

import argparse
from utils import *

def main():
    parser = argparse.ArgumentParser(description='Train a neural network on genetic data.')
    parser.add_argument('--train_prefix', type=str, required=True, help='Training set prefix')
    parser.add_argument('--val_prefix', type=str, required=True, help='Validation set prefix')
    parser.add_argument('--test_prefix', type=str, required=True, help='Test set prefix')
    parser.add_argument('--snp_file', type=str, required=True, help='SNP file')
    parser.add_argument('--p_value_threshold', type=float, default=1e-04, help='P-value threshold for SNP selection')
    parser.add_argument('--train_af_file', type=str, required=True, help='Training allele frequency file')
    parser.add_argument('--train_pheno_file', type=str, required=True, help='Training phenotype file')
    parser.add_argument('--val_pheno_file', type=str, required=True, help='Validation phenotype file')
    parser.add_argument('--test_pheno_file', type=str, required=True, help='Test phenotype file')
    parser.add_argument('--epochs', type=int, default=100, help='Number of epochs to run neural network training')
    parser.add_argument('--batch-size', type=int, default=64, help='Number of batch size to run neural network training')
    parser.add_argument('--output_dir', type=str, default=os.getcwd(), help='Directory to save outputs')
    parser.add_argument('--seed', type=int, default= None, help='Set seed - if not set specified a random seed is set')
    args = parser.parse_args()

    seed = set_seed(args.seed)

    timestamp = get_timestamp()
    start_time = datetime.datetime.now()
    
    # Read PLINK files
    train_data, val_data, test_data = read_plink_files(args.train_prefix, args.val_prefix, args.test_prefix)
    
    # Subset data based on GWAS results
    (
        filtered_train_bim, filtered_train_bed, filtered_train_fam,
        filtered_val_bim, filtered_val_bed, filtered_val_fam,
        filtered_test_bim, filtered_test_bed, filtered_test_fam,
        num_train_snps, num_train_samples, num_val_snps, num_val_samples,
        num_test_snps, num_test_samples, matching_rows
    ) = load_gwas_and_subset(args.snp_file, args.p_value_threshold, train_data, val_data, test_data)
    
    print("matching_rows length: " + str(len(matching_rows)), flush=True)
    
    # Read genotypes
    train_genotypes_modified = read_genotypes(filtered_train_bed, "train")
    val_genotypes_modified = read_genotypes(filtered_val_bed, "validation")
    test_genotypes_modified = read_genotypes(filtered_test_bed, "test")
    
    # Standardize genotypes
    scaled_train_genotypes, scaled_val_genotypes, scaled_test_genotypes = standardize_and_prepare_data(
        args.train_af_file, matching_rows, train_genotypes_modified, val_genotypes_modified, test_genotypes_modified
    )
    
    # Standardize phenotypes
    scaled_train_pheno, scaled_val_pheno, scaled_test_pheno = read_and_standardize_phenotypes(
        args.train_pheno_file, args.val_pheno_file, args.test_pheno_file
    )
    
    # Prepare tensors
    X_train_tensor, Y_train_tensor, X_val_tensor, Y_val_tensor, X_test_tensor, Y_test_tensor = prepare_tensors(
        scaled_train_genotypes, scaled_train_pheno, scaled_val_genotypes, scaled_val_pheno, scaled_test_genotypes, scaled_test_pheno
    )
    
    # Check shapes
    print("X_train shape:", X_train_tensor.shape)
    print("Y_train shape:", Y_train_tensor.shape)
    print("X_val shape:", X_val_tensor.shape)
    print("Y_val shape:", Y_val_tensor.shape)
    print("X_test shape:", X_test_tensor.shape)
    print("Y_test shape:", Y_test_tensor.shape)

    nn_directory_path = os.path.join(args.output_dir, 'nn_training')
    
    if not os.path.exists(nn_directory_path):
        os.makedirs(nn_directory_path)
        print(f"Directory '{nn_directory_path}' created successfully.")
    else:
        print(f"Directory '{nn_directory_path}' already exists.")

    nn_project_name="CNN"

    tuner = setup_tuner(args.epochs, nn_directory_path, nn_project_name, X_train_tensor)
    print("Start hyperparameter search!", flush=True)
    search_best_hyperparameters(tuner, X_train_tensor, Y_train_tensor, X_val_tensor, Y_val_tensor, args.batch_size, args.epochs)
    save_model_parameter_file = os.path.join(args.output_dir, "CNN_best_hp_model_params.obj")
    save_model_file = os.path.join(args.output_dir, "CNN_best_hp_model.obj")
    best_model_params = save_best_hyperparameters_and_model(tuner, save_model_parameter_file, save_model_file)
    
    print("################################### Neural Network Training ########################################")
    print("Fitting the best model...")
    
    early_stopping = setup_early_stopping(args.epochs)
    
    history = train_best_model(
        best_model_params,
        X_train_tensor,
        Y_train_tensor,
        X_val_tensor,
        Y_val_tensor,
        args.epochs,
        args.batch_size,
        early_stopping
    )
    
    print("Plotting results...")
    loss_plot_figure_file = os.path.join(args.output_dir, "CNN_loss_vs_epochs.png")
    plot_training_history(history, loss_plot_figure_file)
    
    print("Saving best model...")
    save_trained_best_model_file = os.path.join(args.output_dir, "CNN_trained_best_model.keras")
    save_best_model(best_model_params, save_trained_best_model_file)
    
    print("Viewing NN architecture...")
    display_model_summary(best_model_params)

    print("###################################### Neural Network Testing #####################################")
    
    Y_train_pred, Y_val_pred, Y_test_pred = get_predictions(best_model_params, X_train_tensor, X_val_tensor, X_test_tensor)

    train_correlation_NN, train_mse_NN = calculate_metrics(scaled_train_pheno, Y_train_pred)
    val_correlation_NN, val_mse_NN = calculate_metrics(scaled_val_pheno, Y_val_pred)
    test_correlation_NN, test_mse_NN = calculate_metrics(scaled_test_pheno, Y_test_pred)

    train_r_squared = train_correlation_NN**2
    val_r_squared = val_correlation_NN**2
    test_r_squared = test_correlation_NN**2

    test_evaluate = evaluate_model(best_model_params, X_test_tensor, Y_test_tensor)
    print("Test evaluate: ", flush=True)
    print(test_evaluate)

    csv_file = os.path.join(args.output_dir, "CNN_results.csv")
    save_results_to_csv(
        csv_file,
        train_correlation_NN,
        val_correlation_NN,
        test_correlation_NN,
        train_mse_NN,
        val_mse_NN,
        test_mse_NN,
        train_r_squared,
        val_r_squared,
        test_r_squared
    )

    print_results(
        train_correlation_NN,
        val_correlation_NN,
        test_correlation_NN,
        train_mse_NN,
        val_mse_NN,
        test_mse_NN,
        train_r_squared,
        val_r_squared,
        test_r_squared
    )

    end_time = datetime.datetime.now()
    time_elapsed = end_time - start_time
    print("Total time elapsed (seconds):", time_elapsed.total_seconds())


if __name__ == '__main__':
    main()
