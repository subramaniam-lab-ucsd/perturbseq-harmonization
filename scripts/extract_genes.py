#!/usr/bin/env python3
"""Extract all genes from Tahoe-100M and save as CSV."""

from datasets import load_dataset
import pandas as pd
from pathlib import Path

print('Loading gene metadata from Tahoe-100M...')
gene_ds = load_dataset('tahoebio/Tahoe-100M', 'gene_metadata', split='train')
gene_df = gene_ds.to_pandas()

print(f'Shape: {gene_df.shape}')
print(f'Columns: {gene_df.columns.tolist()}')
print(f'\nFirst 5 genes:')
print(gene_df.head())

# Save to CSV
output_dir = Path('/home/s5srinivasan/perturbseq-harmonization/data')
output_dir.mkdir(parents=True, exist_ok=True)
output_file = output_dir / 'tahoe_100m_genes.csv'
gene_df.to_csv(output_file, index=False)

print(f'\nSaved {len(gene_df)} genes to: {output_file}')
print(f'File size: {output_file.stat().st_size / 1024:.1f} KB')
