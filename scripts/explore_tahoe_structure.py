#!/usr/bin/env python3
"""
Quick exploration of Tahoe-100M dataset structure.
This script explores the metadata tables to understand the data format
before running the full extraction.
"""

from datasets import load_dataset
import pandas as pd

print("="*60)
print("Exploring Tahoe-100M Dataset Structure")
print("="*60)

# 1. Load and explore cell line metadata
print("\n1. CELL LINE METADATA")
print("-"*40)
cell_line_ds = load_dataset(
    "tahoebio/Tahoe-100M", 
    "cell_line_metadata",
    split="train"
)
cell_line_df = cell_line_ds.to_pandas()
print(f"Shape: {cell_line_df.shape}")
print(f"Columns: {cell_line_df.columns.tolist()}")
print("\nAll cell lines:")
print(cell_line_df)

# Search for A549
print("\n--- Searching for A549 ---")
for col in cell_line_df.columns:
    if cell_line_df[col].dtype == object:
        matches = cell_line_df[cell_line_df[col].str.contains('A549', case=False, na=False)]
        if len(matches) > 0:
            print(f"\nFound 'A549' in column '{col}':")
            print(matches)

# 2. Load gene metadata
print("\n\n2. GENE METADATA")
print("-"*40)
gene_ds = load_dataset(
    "tahoebio/Tahoe-100M",
    "gene_metadata", 
    split="train"
)
gene_df = gene_ds.to_pandas()
print(f"Shape: {gene_df.shape}")
print(f"Columns: {gene_df.columns.tolist()}")
print("\nFirst 10 genes:")
print(gene_df.head(10))

# 3. Load sample metadata to understand sample structure
print("\n\n3. SAMPLE METADATA")
print("-"*40)
sample_ds = load_dataset(
    "tahoebio/Tahoe-100M",
    "sample_metadata",
    split="train"
)
sample_df = sample_ds.to_pandas()
print(f"Shape: {sample_df.shape}")
print(f"Columns: {sample_df.columns.tolist()}")
print("\nFirst 10 samples:")
print(sample_df.head(10))

# Check if sample metadata can tell us about A549
print("\n--- Searching for A549 in sample metadata ---")
for col in sample_df.columns:
    if sample_df[col].dtype == object:
        matches = sample_df[sample_df[col].str.contains('A549', case=False, na=False)]
        if len(matches) > 0:
            print(f"\nFound {len(matches)} rows with 'A549' in column '{col}'")
            print(matches.head(5))

# 4. Load obs_metadata
print("\n\n4. OBS METADATA")
print("-"*40)
try:
    obs_ds = load_dataset(
        "tahoebio/Tahoe-100M",
        "obs_metadata",
        split="train"
    )
    obs_df = obs_ds.to_pandas()
    print(f"Shape: {obs_df.shape}")
    print(f"Columns: {obs_df.columns.tolist()}")
    print("\nFirst 10 rows:")
    print(obs_df.head(10))
    
    # Search for A549
    print("\n--- Searching for A549 in obs metadata ---")
    for col in obs_df.columns:
        if obs_df[col].dtype == object:
            sample_vals = obs_df[col].head(1000)
            matches = sample_vals[sample_vals.str.contains('A549', case=False, na=False)]
            if len(matches) > 0:
                print(f"\nFound 'A549' in column '{col}' (first 1000 rows)")
                print(matches.head(5))
except Exception as e:
    print(f"Could not load obs_metadata: {e}")

# 5. Look at a few records from expression data
print("\n\n5. EXPRESSION DATA STRUCTURE (first 5 records)")
print("-"*40)
ds = load_dataset("tahoebio/Tahoe-100M", streaming=True, split="train")
for i, record in enumerate(ds):
    if i >= 5:
        break
    print(f"\nRecord {i}:")
    print(f"  Keys: {record.keys()}")
    for key, value in record.items():
        if isinstance(value, (list, tuple)):
            print(f"  {key}: list of length {len(value)} - first 5: {value[:5]}")
        else:
            print(f"  {key}: {value}")

# 6. Check if there's cell line info in expression records
print("\n\n6. CHECKING CELL LINE INFO IN EXPRESSION DATA")
print("-"*40)
print("Scanning first 1000 records to find A549...")
a549_count = 0
cell_lines_found = set()

for i, record in enumerate(ds):
    if i >= 1000:
        break
    
    # Check all string fields for cell line info
    for key, value in record.items():
        if isinstance(value, str):
            if 'A549' in value.upper():
                a549_count += 1
                if a549_count <= 3:
                    print(f"Found A549 in record {i}, key '{key}': {value}")
            if 'cell' in key.lower() or key == 'cell_line_id':
                cell_lines_found.add(value)

print(f"\nA549 records in first 1000: {a549_count}")
print(f"Unique cell line values found: {cell_lines_found}")

print("\n" + "="*60)
print("Exploration complete!")
print("="*60)
