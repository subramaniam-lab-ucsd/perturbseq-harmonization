<p align="center">
  <img src="https://img.shields.io/badge/Perturb--seq-Harmonization-blueviolet?style=for-the-badge&logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIyNCIgaGVpZ2h0PSIyNCIgdmlld0JveD0iMCAwIDI0IDI0IiBmaWxsPSJub25lIiBzdHJva2U9IndoaXRlIiBzdHJva2Utd2lkdGg9IjIiPjxjaXJjbGUgY3g9IjEyIiBjeT0iMTIiIHI9IjEwIi8+PHBhdGggZD0iTTggMTRzMS41IDIgNCAyczQtMiA0LTIiLz48bGluZSB4MT0iOSIgeTE9IjkiIHgyPSI5LjAxIiB5Mj0iOSIvPjxsaW5lIHgxPSIxNSIgeTE9IjkiIHgyPSIxNS4wMSIgeTI9IjkiLz48L3N2Zz4="/>
</p>

<h1 align="center">🧬 Perturb-seq Harmonization</h1>

<p align="center">
  <strong>A comprehensive toolkit for analyzing and harmonizing large-scale perturbation single-cell datasets</strong>
</p>

<p align="center">
  <a href="#-overview">Overview</a> •
  <a href="#-features">Features</a> •
  <a href="#-installation">Installation</a> •
  <a href="#-quick-start">Quick Start</a> •
  <a href="#-project-structure">Structure</a> •
  <a href="#-tutorials">Tutorials</a> •
  <a href="#-data">Data</a>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/python-3.8+-3776AB?style=flat-square&logo=python&logoColor=white"/>
  <img src="https://img.shields.io/badge/license-MIT-green?style=flat-square"/>
  <img src="https://img.shields.io/badge/Tahoe--100M-Dataset-orange?style=flat-square"/>
  <img src="https://img.shields.io/badge/scRNA--seq-Analysis-purple?style=flat-square"/>
</p>

---

## 🔬 Overview

**Perturb-seq Harmonization** is a research project designed to work with large-scale perturbation single-cell RNA sequencing datasets, with a particular focus on the **Tahoe-100M** dataset—one of the largest publicly available perturbation atlases.

This repository provides:

- 🎯 **Data extraction utilities** for working with massive perturbation datasets
- 🧠 **Transcription factor enrichment** using curated regulatory databases
- 📊 **Analysis pipelines** for perturbation response characterization
- 📚 **Comprehensive tutorials** covering state-of-the-art methods

---

## ✨ Features

### 🔹 Tahoe-100M Integration

Seamlessly work with the Tahoe-100M dataset containing **100+ million cells** across multiple cell lines and perturbations.

```python
from datasets import load_dataset

# Stream the massive dataset efficiently
ds = load_dataset("tahoebio/Tahoe-100M", streaming=True, split="train")
```

### 🔹 Transcription Factor Enrichment

Annotate gene lists with comprehensive TF information from two gold-standard databases:

| Database | Description | Coverage |
|----------|-------------|----------|
| **CollecTRI** | Literature-curated TF-target interactions | ~1,200 TFs |
| **DoRothEA** | Curated regulons with confidence levels (A-D) | ~1,400 TFs |

```python
import decoupler as dc

# Load TF databases
collectri = dc.get_collectri(organism='human')
dorothea = dc.get_dorothea(organism='human', levels=['A', 'B', 'C', 'D'])
```

### 🔹 Cell Line-Specific Analysis

Extract and analyze specific cell lines from the dataset:

- **A549** (lung adenocarcinoma)
- **HepG2** (hepatocellular carcinoma)
- **K562** (chronic myelogenous leukemia)
- And more...

---

## 📦 Installation

### Prerequisites

- Python 3.8+
- pip or conda

### Setup

```bash
# Clone the repository
git clone https://github.com/subramaniam-lab-ucsd/perturbseq-harmonization.git
cd perturbseq-harmonization

# Create virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Dependencies

```
datasets      # HuggingFace datasets library
anndata       # Annotated data matrices for scRNA-seq
scipy         # Scientific computing
pandas        # Data manipulation
pubchempy     # PubChem compound queries
decoupler     # TF activity inference (optional)
```

---

## 🚀 Quick Start

### 1. Extract Genes from Tahoe-100M

```python
from datasets import load_dataset

# Load gene metadata
gene_ds = load_dataset("tahoebio/Tahoe-100M", "gene_metadata", split="train")
gene_df = gene_ds.to_pandas()

print(f"Total genes: {len(gene_df)}")
# Output: Total genes: ~60,000
```

### 2. Enrich with TF Information

```bash
# Run the TF enrichment pipeline
python scripts/add_tf_info.py
```

This will generate:
- `tahoe_100m_genes_tf_enriched.csv` - Compact version with TF flags and counts
- `tahoe_100m_genes_tf_enriched_full.csv` - Full version with target gene lists
- `tahoe_100m_transcription_factors.csv` - TF-only subset

### 3. Explore the Data Structure

```bash
# Understand the Tahoe-100M dataset organization
python scripts/explore_tahoe_structure.py
```

---

## 📁 Project Structure

```
perturbseq-harmonization/
├── 📂 data/                          # Processed data files
│   ├── tahoe_100m_genes.csv          # Base gene list
│   ├── tahoe_100m_genes_tf_enriched.csv
│   ├── tahoe_100m_transcription_factors.csv
│   └── tahoe_a549/                   # A549-specific extracts
│       ├── a549_obs_metadata.csv
│       ├── cell_line_metadata.csv
│       └── gene_metadata.csv
│
├── 📂 notebooks/                     # Analysis notebooks
│   ├── add_tf_info_to_genes.ipynb    # TF enrichment workflow
│   └── tahoe_100m_gene_extraction.ipynb
│
├── 📂 scripts/                       # Python utilities
│   ├── add_tf_info.py                # TF annotation pipeline
│   ├── explore_tahoe_structure.py    # Dataset exploration
│   ├── extract_genes.py              # Gene extraction utilities
│   └── find_a549_id.py               # Cell line identification
│
├── 📂 tutorial/                      # External tutorials & references
│   ├── STATE_*.ipynb                 # STATE model tutorials
│   ├── Tahoe-*.ipynb                 # Tahoe dataset tutorials
│   ├── scPerturb_Tutorial_E-distance.ipynb
│   └── ...
│
├── 📂 docs/                          # Documentation
├── requirements.txt
└── README.md
```

---

## 📚 Tutorials

The `tutorial/` directory contains comprehensive notebooks for various analysis methods:

### 🏔️ Tahoe-100M Tutorials

| Notebook | Description |
|----------|-------------|
| `Tahoe-100M loading data.ipynb` | Load and convert Tahoe data to AnnData format |
| `tahoe_100M_download_tutorial.ipynb` | Download and manage the dataset |
| `Tahoe-x1_Clustering_Tutorial.ipynb` | Clustering analysis of perturbation responses |
| `Tahoe-x1_Training_Tutorial.ipynb` | Train models on Tahoe data |

### 🧪 STATE Model Tutorials

| Notebook | Description |
|----------|-------------|
| `STATE_Replogle_Nadig_Training.ipynb` | Train STATE on Replogle/Nadig data |
| `STATE_Run_Tahoe_Inference.ipynb` | Run STATE inference on Tahoe |
| `STATE_Use_State_Embedding.ipynb` | Work with STATE embeddings |
| `STATE_for_Virtual_Cell_Challenge.ipynb` | Virtual cell challenge applications |

### 🔬 Analysis Methods

| Notebook | Description |
|----------|-------------|
| `scPerturb_Tutorial_E-distance.ipynb` | E-distance for perturbation comparison |
| `scBaseCount_Tutorial.ipynb` | Base count normalization methods |
| `DrugReflector_Tutorial.ipynb` | Drug response prediction |

---

## 📊 Data

### Generated Datasets

After running the enrichment pipeline, you'll have access to:

#### `tahoe_100m_genes_tf_enriched.csv`

| Column | Description |
|--------|-------------|
| `gene_symbol` | HGNC gene symbol |
| `ensembl_id` | Ensembl gene ID |
| `token_id` | Tokenized ID for ML models |
| `is_tf` | Whether gene is a transcription factor |
| `is_tf_collectri` | TF in CollecTRI database |
| `is_tf_dorothea` | TF in DoRothEA database |
| `is_target` | Whether gene is a TF target |
| `collectri_n_targets` | Number of targets (CollecTRI) |
| `dorothea_n_targets` | Number of targets (DoRothEA) |
| `*_n_activated` | Number of activated targets |
| `*_n_repressed` | Number of repressed targets |
| `*_n_regulators` | Number of TFs regulating this gene |

### Example: Top Transcription Factors

```
gene_symbol  collectri_n_targets  collectri_n_activated  collectri_n_repressed
-----------  -------------------  ---------------------  ---------------------
ESR1                       2912                   1521                   1391
STAT3                      1823                   1089                    734
TP53                       1789                    892                    897
NFKB1                      1654                   1123                    531
MYC                        1612                   1098                    514
```

---

## 🛠️ Usage Examples

### Extract A549 Cell Line Data

```python
from datasets import load_dataset
import pandas as pd

# Load cell line metadata to find A549 identifier
cell_line_ds = load_dataset("tahoebio/Tahoe-100M", "cell_line_metadata", split="train")
cell_line_df = cell_line_ds.to_pandas()

# Find A549
a549_info = cell_line_df[cell_line_df['cell_line'].str.contains('A549', case=False)]
print(a549_info)
```

### Create AnnData from Tahoe Records

```python
from datasets import load_dataset
from scipy.sparse import csr_matrix
import anndata
import pandas as pd

def create_anndata_from_generator(generator, gene_vocab, sample_size=None):
    """Convert Tahoe records to AnnData format."""
    sorted_vocab_items = sorted(gene_vocab.items())
    token_ids, gene_names = zip(*sorted_vocab_items)
    
    data, indices, indptr = [], [], [0]
    obs_data = []
    
    for i, cell in enumerate(generator):
        if sample_size and i >= sample_size:
            break
            
        genes = cell['genes']
        expressions = cell['expressions']
        
        for gene, expr in zip(genes, expressions):
            data.append(expr)
            indices.append(gene)
        indptr.append(len(data))
        
        obs_entry = {k: v for k, v in cell.items() 
                     if k not in ['genes', 'expressions']}
        obs_data.append(obs_entry)
    
    expr_matrix = csr_matrix((data, indices, indptr), 
                              shape=(len(indptr) - 1, len(gene_names)))
    obs_df = pd.DataFrame(obs_data)
    
    adata = anndata.AnnData(X=expr_matrix, obs=obs_df)
    adata.var.index = pd.Index(gene_names, name='ensembl_id')
    
    return adata
```

---

## 📖 References

### Datasets

- **Tahoe-100M**: [HuggingFace Dataset](https://huggingface.co/datasets/tahoebio/Tahoe-100M)

### TF Databases

- **CollecTRI**: Müller-Dott, S., et al. (2023). *Nucleic Acids Research*
- **DoRothEA**: Garcia-Alonso, L., et al. (2019). *Genome Research*

### Related Tools

- [decoupler](https://github.com/saezlab/decoupler-py) - TF activity inference
- [scanpy](https://scanpy.readthedocs.io/) - Single-cell analysis
- [anndata](https://anndata.readthedocs.io/) - Annotated data matrices

---

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 👥 Acknowledgments

- **Tahoe Bio** for providing the Tahoe-100M dataset
- **Saez Lab** for the decoupler package and TF databases
- The broader single-cell and perturbation genomics community

---

<p align="center">
  <sub>Built with ❤️ by the Subramaniam Lab @ UCSD</sub>
</p>