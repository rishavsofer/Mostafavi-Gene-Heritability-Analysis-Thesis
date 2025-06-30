#!/usr/bin/env python3
# extract_feature_batch.py - Extracts multiple features in one pass

import pandas as pd
import numpy as np
import sys
import os
import json
import traceback
import gc

def extract_features_batch(feature_file, feature_columns_str, gencode_file, output_dir):
    """Extract multiple features in one pass"""
    
    print(f"Reading feature file: {feature_file}")
    print(f"Processing feature columns: {feature_columns_str}")
    print(f"Reading gencode file: {gencode_file}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    def exit_with_error(message):
        print(f"ERROR: {message}")
        traceback.print_exc()
        sys.exit(1)
    
    # Parse feature columns
    feature_indices = [int(col) for col in feature_columns_str.split(',')]
    print(f"Processing {len(feature_indices)} features")
    
    # Read header to get column names
    try:
        header_df = pd.read_csv(feature_file, sep='\t', nrows=0)
        all_columns = header_df.columns.tolist()
        
        if 'ensg' not in all_columns:
            exit_with_error("Required 'ensg' column not found in the feature file.")
        
        # Validate column indices
        for idx in feature_indices:
            if idx >= len(all_columns):
                exit_with_error(f"Column index {idx} is out of range. Header has {len(all_columns)} columns.")
        
        # Get column positions
        ensg_col_idx = all_columns.index('ensg')
        usecols = [ensg_col_idx] + feature_indices
        
        del header_df
        gc.collect()
    except Exception as e:
        exit_with_error(f"Failed to read feature file header: {str(e)}")
    
    # Read data once for all features
    try:
        print(f"Reading columns: {usecols}")
        feature_data = pd.read_csv(feature_file, sep='\t', usecols=usecols)
        print(f"Feature data read successfully. Found {len(feature_data)} genes.")
    except Exception as e:
        exit_with_error(f"Failed to read feature data: {str(e)}")
    
    # Read gencode data once
    try:
        print("Reading gencode data...")
        gencode_data = pd.read_csv(
            gencode_file,
            sep='\t',
            header=None,
            names=['chrom', 'start', 'end', 'gene_id']
        )
        gencode_data['ensg'] = gencode_data['gene_id'].str.extract(r'(ENSG\d+)')
        print(f"Gencode data read successfully. Found {len(gencode_data)} entries.")
    except Exception as e:
        exit_with_error(f"Failed to read gencode file: {str(e)}")
    
    # Process each feature
    results = []
    
    for feat_idx in feature_indices:
        try:
            # Get feature name
            original_feature_name = all_columns[feat_idx]
            feature_name = original_feature_name.replace(':', '_').replace('/', '_').replace('\\', '_').replace(' ', '_')
            
            print(f"\nProcessing feature {feat_idx}: {original_feature_name} -> {feature_name}")
            
            # Get feature data for this column
            feature_subset = feature_data[['ensg', original_feature_name]].copy()
            feature_subset.rename(columns={original_feature_name: feature_name}, inplace=True)
            
            # Remove NaN values
            feature_subset = feature_subset.dropna(subset=[feature_name])
            print(f"After removing NaN values: {len(feature_subset)} genes.")
            
            # Analyze feature values
            unique_values = feature_subset[feature_name].unique()
            print(f"Found {len(unique_values)} unique values")
            
            is_binary = False
            
            if len(unique_values) <= 3:
                value_counts = feature_subset[feature_name].value_counts()
                print("Value counts:")
                print(value_counts)
                
                if set(value_counts.index).issubset({0, 1, 0.0, 1.0}):
                    print("Binary feature detected (0/1)")
                    is_binary = True
                elif len(unique_values) == 2 and set(value_counts.index).issubset({-1, 1, -1.0, 1.0}):
                    print("Binary feature detected (-1/1)")
                    is_binary = True
            
            # Select genes based on feature values
            if is_binary:
                max_value = feature_subset[feature_name].max()
                selected_genes = feature_subset[feature_subset[feature_name] == max_value]
            else:
                # Quantile-based selection for continuous features
                total_genes = len(feature_subset)
                min_required = max(2000, int(total_genes * 0.2))
                
                quantile_threshold = 0.9
                selected_genes = None
                
                while quantile_threshold > 0.1:
                    threshold_value = feature_subset[feature_name].quantile(quantile_threshold)
                    mask = feature_subset[feature_name] >= threshold_value
                    selected_count = mask.sum()
                    
                    if selected_count >= min_required:
                        selected_genes = feature_subset[mask].copy()
                        print(f"Selected {len(selected_genes)} genes using quantile {quantile_threshold}")
                        break
                    
                    quantile_threshold -= 0.05
                
                if selected_genes is None:
                    top_indices = feature_subset[feature_name].values.argsort()[-min_required:]
                    selected_genes = feature_subset.iloc[top_indices].copy()
            
            print(f"Final selection: {len(selected_genes)} genes")
            
            # Merge with gencode data
            merged_data = pd.merge(
                selected_genes,
                gencode_data,
                on='ensg',
                how='inner'
            )
            
            print(f"Found {len(merged_data)} genes matching between feature and gencode files.")
            
            if len(merged_data) == 0:
                print(f"WARNING: No matching genes found for feature {feature_name}")
                continue
            
            # Write outputs
            feature_output_dir = os.path.join(output_dir, feature_name)
            os.makedirs(feature_output_dir, exist_ok=True)
            
            gene_windows_file = os.path.join(feature_output_dir, f"{feature_name}_gene_windows.bed")
            merged_data[['chrom', 'start', 'end', 'gene_id']].to_csv(
                gene_windows_file, 
                sep='\t', 
                index=False,
                header=False
            )
            
            feature_file_output = os.path.join(feature_output_dir, f"{feature_name}_feature.tsv")
            merged_data[['ensg', feature_name]].to_csv(feature_file_output, sep='\t', index=False)
            
            results.append({
                'column': feat_idx,
                'feature_name': feature_name,
                'gene_windows_file': gene_windows_file,
                'num_genes': len(merged_data),
                'is_binary': is_binary
            })
            
            print(f"Feature {feature_name} processed successfully")
            
        except Exception as e:
            print(f"ERROR processing feature {feat_idx}: {str(e)}")
            traceback.print_exc()
            continue
    
    # Write manifest of processed features
    manifest_file = os.path.join(output_dir, 'batch_manifest.json')
    with open(manifest_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nBatch processing complete. Processed {len(results)} features.")
    print(f"Manifest written to: {manifest_file}")
    
    return results

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python extract_feature_batch.py <feature_file> <feature_columns> <gencode_file> <output_dir>")
        print("  feature_columns: comma-separated list of column indices")
        sys.exit(1)
    
    feature_file = sys.argv[1]
    feature_columns = sys.argv[2]
    gencode_file = sys.argv[3]
    output_dir = sys.argv[4]
    
    try:
        extract_features_batch(feature_file, feature_columns, gencode_file, output_dir)
    except Exception as e:
        print(f"Fatal error: {str(e)}")
        traceback.print_exc()
        sys.exit(1)

