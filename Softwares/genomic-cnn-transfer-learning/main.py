#!/bin/env python

import argparse
from utils import *

def main():
    parser = argparse.ArgumentParser(description='Train a neural network on genetic data.')
    parser.add_argument('--train_prefix', type=str, required=True, help='Training set prefix')
    parser.add_argument('--val_prefix', type=str, required=True, help='Validation set prefix')
    parser.add_argument('--test_prefix', type=str, required=True, help='Test set prefix')
    parser.add_argument('--source_snp_file', type=str, required=True, help='SNP file')
    parser.add_argument('--p_value_threshold', type=float, default=1e-04, help='P-value threshold for SNP selection')
    parser.add_argument('--train_af_file', type=str, required=True, help='Training allele frequency file')
    parser.add_argument('--train_pheno_file', type=str, required=True, help='Training phenotype file')
    parser.add_argument('--val_pheno_file', type=str, required=True, help='Validation phenotype file')
    parser.add_argument('--test_pheno_file', type=str, required=True, help='Test phenotype file')
    parser.add_argument('--source_model_best_params_file', type=str, required=True, help='Best parameters from HP optimization of source neural network model')
    parser.add_argument('--source_model_keras_file', type=str, required=True, help='Source neural network model (.keras format)')
    parser.add_argument('--epochs', type=int, default=100, help='Number of epochs to run neural network training')
    parser.add_argument('--batch_size', type=int, default=64, help='Number of batch size to run neural network training')
    parser.add_argument('--output_dir', type=str, default=os.getcwd(), help='Directory to save outputs')
    parser.add_argument('--seed', type=int, default= None, help='Set seed - if not set specified a random seed is set')
    args = parser.parse_args()

    seed = set_seed(args.seed)

    timestamp = get_timestamp()
    start_time = datetime.datetime.now()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Read PLINK files
    train_data, val_data, test_data = read_plink_files(args.train_prefix, args.val_prefix, args.test_prefix)
    
    # Subset data based on GWAS results
    (
        filtered_train_bim, filtered_train_bed, filtered_train_fam,
        filtered_val_bim, filtered_val_bed, filtered_val_fam,
        filtered_test_bim, filtered_test_bed, filtered_test_fam,
        num_train_snps, num_train_samples, num_val_snps, num_val_samples,
        num_test_snps, num_test_samples, matching_rows
    ) = load_gwas_and_subset(args.source_snp_file, args.p_value_threshold, train_data, val_data, test_data)
    
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

    pretrained_optimizer = load_best_hyperparameters(args.source_model_best_params_file)

    freeze_start_time = datetime.datetime.now()
    model_path = load_pretrained_model(args.source_model_keras_file)

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
        num_batch_size=args.batch_size,
        num_epochs=args.epochs
    )
    freeze_end_time = datetime.datetime.now()

    ft_start_time = datetime.datetime.now()

    layer_count = len(freeze_transfer_model.layers)
    best_correlation = 0
    best_layer = 0

    for layer_value in range(2, layer_count):
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
            num_batch_size=args.batch_size,
            num_epochs=args.epochs
        )
        print(f"Test correlation: {test_correlation}", flush=True)
        if test_correlation > best_correlation:
            best_correlation = test_correlation
            best_layer = layer_value

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
        num_batch_size=args.batch_size,
        num_epochs=args.epochs
    )
    ft_end_time = datetime.datetime.now()

    TL_freeze_loss_plot_figure_file=os.path.join(args.output_dir, "CNN_TL_freeze_loss_vs_epochs_plot.png")
    TL_fine_tune_loss_plot_figure_file=os.path.join(args.output_dir, "CNN_TL_fine_tune_loss_vs_epochs_plot.png")

    plot_training_validation_loss(freeze_transfer_model_history, 'TL Freeze Training and Validation Loss over Epochs', TL_freeze_loss_plot_figure_file)
    plot_training_validation_loss(best_fine_tune_model_history, 'TL Fine Tune Training and Validation Loss over Epochs', TL_fine_tune_loss_plot_figure_file)

    save_trained_best_model_TL_freeze_file=os.path.join(args.output_dir, "best_CNN_TL_freeze_model.keras")
    save_trained_best_model_TL_fine_tune_file=os.path.join(args.output_dir, "best_CNN_TL_fine_tune_model.keras")
    save_model(freeze_transfer_model, save_trained_best_model_TL_freeze_file)
    save_model(best_fine_tune_model, save_trained_best_model_TL_fine_tune_file)

    freeze_metrics = calculate_performance_metrics(
        freeze_transfer_model, X_train_tensor, X_val_tensor, X_test_tensor, Y_train_tensor, Y_val_tensor, Y_test_tensor
    )
    fine_tune_metrics = calculate_performance_metrics(
        best_fine_tune_model, X_train_tensor, X_val_tensor, X_test_tensor, Y_train_tensor, Y_val_tensor, Y_test_tensor
    )

    print_results("Freeze", freeze_metrics)
    print_results("Fine-Tune", fine_tune_metrics)

    freeze_time = freeze_end_time - freeze_start_time
    ft_time = ft_end_time - ft_start_time

    data = [
        ["Approach", "Train_Correlation", "Validation_Correlation", "Test_Correlation", "Train_MSE", "Validation_MSE", "Test_MSE", "Train_R2", "Validation_R2", "Test_R2", "p_value_threshold", "time"],
        ["Freeze", freeze_metrics['train'][0], freeze_metrics['val'][0], freeze_metrics['test'][0], freeze_metrics['train'][1], freeze_metrics['val'][1], freeze_metrics['test'][1], freeze_metrics['train'][2], freeze_metrics['val'][2], freeze_metrics['test'][2], args.p_value_threshold, freeze_time.total_seconds()],
        ["Fine_Tune", fine_tune_metrics['train'][0], fine_tune_metrics['val'][0], fine_tune_metrics['test'][0], fine_tune_metrics['train'][1], fine_tune_metrics['val'][1], fine_tune_metrics['test'][1], fine_tune_metrics['train'][2], fine_tune_metrics['val'][2], fine_tune_metrics['test'][2], args.p_value_threshold, ft_time.total_seconds()]
    ]

    csv_file=os.path.join(args.output_dir, "CNN_TL_results.csv")
    save_results_to_csv(data, csv_file)
    print(f"Results saved to {csv_file}")


if __name__ == '__main__':
    main()
