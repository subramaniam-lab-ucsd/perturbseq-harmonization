#!/usr/bin/env python3
"""
Find A549 cell_line_id and extract the data.
"""

from datasets import load_dataset
import pandas as pd

print("="*60)
print("Finding A549 Cell Line ID")
print("="*60)

# Load cell line metadata
print("\n1. Loading cell_line_metadata...")
cell_line_ds = load_dataset(
    "tahoebio/Tahoe-100M", 
    "cell_line_metadata",
    split="train"
)
cell_line_df = cell_line_ds.to_pandas()
print(f"Shape: {cell_line_df.shape}")
print(f"Columns: {cell_line_df.columns.tolist()}")
print("\nAll cell lines:")
print(cell_line_df.to_string())

# Find A549
print("\n\n2. Searching for A549...")
for col in cell_line_df.columns:
    if cell_line_df[col].dtype == object:
        matches = cell_line_df[cell_line_df[col].str.contains('A549', case=False, na=False)]
        if len(matches) > 0:
            print(f"\nFound 'A549' in column '{col}':")
            print(matches.to_string())
            
# Also look at obs_metadata to find the cell_line_id for A549
print("\n\n3. Loading obs_metadata to find A549 cell_line mapping...")
obs_ds = load_dataset(
    "tahoebio/Tahoe-100M",
    "obs_metadata",
    split="train"
)
# Convert to pandas (this is a large table, might take a moment)
print("Converting to pandas (may take a moment)...")
obs_df = obs_ds.to_pandas()
print(f"OBS metadata shape: {obs_df.shape}")

# Find A549 cells and their cell_line info
a549_obs = obs_df[obs_df['cell_name'] == 'A549']
print(f"\nA549 cells in obs_metadata: {len(a549_obs)}")

if len(a549_obs) > 0:
    print("\nA549 sample (first 5 rows):")
    print(a549_obs.head().to_string())
    
    # Check if there's a cell_line_id column
    if 'cell_line' in a549_obs.columns:
        print(f"\nA549 cell_line values: {a549_obs['cell_line'].unique()}")
    
    # Get unique drugs for A549
    if 'drug' in a549_obs.columns:
        print(f"\nUnique drugs tested in A549: {a549_obs['drug'].nunique()}")
        print(f"First 20 drugs: {a549_obs['drug'].unique()[:20].tolist()}")

# Check the link to cell_line_metadata
print("\n\n4. Checking cell_line_metadata for mapping...")
print(cell_line_df.to_string())

print("\nDone!")
