#!/usr/bin/env python3
"""
Add Transcription Factor information to Tahoe-100M genes.

Uses:
- CollecTRI - Literature-curated TF-target database
- DoRothEA - Curated TF regulons with confidence levels
"""

import pandas as pd
import numpy as np
from pathlib import Path
import decoupler as dc

print("="*60)
print("Adding TF Information to Tahoe-100M Genes")
print("="*60)

# Paths
DATA_DIR = Path('/home/s5srinivasan/perturbseq-harmonization/data')

# 1. Load gene list
print("\n1. Loading Tahoe-100M gene list...")
genes_df = pd.read_csv(DATA_DIR / 'tahoe_100m_genes.csv')
print(f"   Total genes: {len(genes_df)}")

# 2. Load TF databases
print("\n2. Loading TF databases...")

print("   Loading CollecTRI...")
collectri = dc.get_collectri(organism='human')
print(f"   CollecTRI: {collectri['source'].nunique()} TFs, {len(collectri)} interactions")

print("   Loading DoRothEA (levels A-D)...")
dorothea = dc.get_dorothea(organism='human', levels=['A', 'B', 'C', 'D'])
print(f"   DoRothEA: {dorothea['source'].nunique()} TFs, {len(dorothea)} interactions")

# 3. Get TF and target sets
collectri_tfs = set(collectri['source'].unique())
dorothea_tfs = set(dorothea['source'].unique())
all_tfs = collectri_tfs | dorothea_tfs

collectri_targets = set(collectri['target'].unique())
dorothea_targets = set(dorothea['target'].unique())
all_targets = collectri_targets | dorothea_targets

print(f"\n3. TF Statistics:")
print(f"   Combined unique TFs: {len(all_tfs)}")
print(f"   Combined unique targets: {len(all_targets)}")

# 4. Build regulon info
print("\n4. Building regulon information...")

def build_regulon_info(network_df, prefix):
    """Build regulon information for each TF."""
    regulon_info = []
    
    for tf in network_df['source'].unique():
        tf_network = network_df[network_df['source'] == tf]
        all_targets = tf_network['target'].tolist()
        
        if 'weight' in tf_network.columns:
            activated = tf_network[tf_network['weight'] > 0]['target'].tolist()
            repressed = tf_network[tf_network['weight'] < 0]['target'].tolist()
        else:
            activated = all_targets
            repressed = []
        
        regulon_info.append({
            'gene_symbol': tf,
            f'{prefix}_n_targets': len(all_targets),
            f'{prefix}_n_activated': len(activated),
            f'{prefix}_n_repressed': len(repressed),
            f'{prefix}_targets': ','.join(sorted(all_targets)),
        })
    
    return pd.DataFrame(regulon_info)

collectri_regulons = build_regulon_info(collectri, 'collectri')
dorothea_regulons = build_regulon_info(dorothea, 'dorothea')

print(f"   CollecTRI regulons: {len(collectri_regulons)}")
print(f"   DoRothEA regulons: {len(dorothea_regulons)}")

# 5. Count regulators per target
print("\n5. Counting regulators for target genes...")
collectri_regulator_counts = collectri.groupby('target')['source'].nunique().reset_index()
collectri_regulator_counts.columns = ['gene_symbol', 'collectri_n_regulators']

dorothea_regulator_counts = dorothea.groupby('target')['source'].nunique().reset_index()
dorothea_regulator_counts.columns = ['gene_symbol', 'dorothea_n_regulators']

# 6. Merge all information
print("\n6. Merging information...")
genes_enriched = genes_df.copy()

# Add TF flags
genes_enriched['is_tf_collectri'] = genes_enriched['gene_symbol'].isin(collectri_tfs)
genes_enriched['is_tf_dorothea'] = genes_enriched['gene_symbol'].isin(dorothea_tfs)
genes_enriched['is_tf'] = genes_enriched['gene_symbol'].isin(all_tfs)

# Add target flags
genes_enriched['is_target_collectri'] = genes_enriched['gene_symbol'].isin(collectri_targets)
genes_enriched['is_target_dorothea'] = genes_enriched['gene_symbol'].isin(dorothea_targets)
genes_enriched['is_target'] = genes_enriched['gene_symbol'].isin(all_targets)

# Merge regulon info
genes_enriched = genes_enriched.merge(collectri_regulons, on='gene_symbol', how='left')
genes_enriched = genes_enriched.merge(dorothea_regulons, on='gene_symbol', how='left')

# Merge regulator counts
genes_enriched = genes_enriched.merge(collectri_regulator_counts, on='gene_symbol', how='left')
genes_enriched = genes_enriched.merge(dorothea_regulator_counts, on='gene_symbol', how='left')

# Fill NaN values for numeric columns
for col in ['collectri_n_targets', 'collectri_n_activated', 'collectri_n_repressed',
            'dorothea_n_targets', 'dorothea_n_activated', 'dorothea_n_repressed',
            'collectri_n_regulators', 'dorothea_n_regulators']:
    if col in genes_enriched.columns:
        genes_enriched[col] = genes_enriched[col].fillna(0).astype(int)

# 7. Print summary
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"Total genes: {len(genes_enriched)}")
print(f"TFs (any database): {genes_enriched['is_tf'].sum()} ({genes_enriched['is_tf'].mean()*100:.2f}%)")
print(f"Targets (any database): {genes_enriched['is_target'].sum()} ({genes_enriched['is_target'].mean()*100:.2f}%)")

# Top TFs
print("\nTop 10 TFs by number of targets (CollecTRI):")
top_tfs = genes_enriched[genes_enriched['is_tf']].nlargest(10, 'collectri_n_targets')[
    ['gene_symbol', 'collectri_n_targets', 'collectri_n_activated', 'collectri_n_repressed']
]
print(top_tfs.to_string(index=False))

# 8. Save files
print("\n" + "="*60)
print("Saving files...")
print("="*60)

# Full version
full_output = DATA_DIR / 'tahoe_100m_genes_tf_enriched_full.csv'
genes_enriched.to_csv(full_output, index=False)
print(f"\nFull version: {full_output}")
print(f"  Size: {full_output.stat().st_size / (1024*1024):.2f} MB")

# Compact version (no target lists)
compact_cols = [
    'gene_symbol', 'ensembl_id', 'token_id',
    'is_tf', 'is_tf_collectri', 'is_tf_dorothea',
    'is_target', 'is_target_collectri', 'is_target_dorothea',
    'collectri_n_targets', 'collectri_n_activated', 'collectri_n_repressed',
    'dorothea_n_targets', 'dorothea_n_activated', 'dorothea_n_repressed',
    'collectri_n_regulators', 'dorothea_n_regulators'
]
genes_compact = genes_enriched[[c for c in compact_cols if c in genes_enriched.columns]]

compact_output = DATA_DIR / 'tahoe_100m_genes_tf_enriched.csv'
genes_compact.to_csv(compact_output, index=False)
print(f"\nCompact version: {compact_output}")
print(f"  Size: {compact_output.stat().st_size / 1024:.1f} KB")

# TFs only
tf_output = DATA_DIR / 'tahoe_100m_transcription_factors.csv'
tf_only = genes_enriched[genes_enriched['is_tf']]
tf_only.to_csv(tf_output, index=False)
print(f"\nTFs only: {tf_output}")
print(f"  TFs: {len(tf_only)}")
print(f"  Size: {tf_output.stat().st_size / (1024*1024):.2f} MB")

print("\n" + "="*60)
print("COMPLETE!")
print("="*60)
