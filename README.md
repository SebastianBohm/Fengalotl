# Fengalotl 🦎

[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)
[![Shiny](https://img.shields.io/badge/shiny-python-green.svg)](https://shiny.posit.co/py/)

An interactive Shiny application for exploring spatial transcriptomics data from the axolotl brain.

## 📋 Overview

Fengalotl provides an interactive web interface to explore spatial gene expression data from adult axolotl brain regions, including:

| Brain Region | Replicates |
|--------------|------------|
| 🧠 Metencephalon (hindbrain) | 2 |
| 👃 Olfactory bulb | 2 |
| 🔴 Pituitary | 1 |
| 🧩 Telencephalon (forebrain) | 2 |
| 🔷 Thalamencephalon (diencephalon) | 3 |

## ✨ Features

- **Spatial visualization**: View gene expression patterns in their spatial context
- **Gene expression plots**: Visualize expression levels across cell clusters
- **Differential expression**: Identify marker genes for each cluster
- **Annotated gene names**: Gene IDs mapped to Axolotl Tanaka annotations (~8,200 genes)
- **Fast loading**: Cached data for responsive interactions after initial load

---

## 🚀 Installation

### Requirements
- Python 3.12 or higher
- A package manager: [conda](https://docs.conda.io/) or [mamba](https://mamba.readthedocs.io/) (recommended)

### Setup Instructions

**a. On a remote server:** Connect via SSH  
**b. On a local machine (MacOS):** Open Terminal  
**c. On a local machine (Windows):** Press `Windows key` + `X`, select Windows Terminal

```bash
# Create a new conda environment
mamba create -n fengalotl python=3.12

# Activate the environment
mamba activate fengalotl

# Clone the repository
git clone --branch main https://github.com/quadbio/fengalotl.git

# Navigate to the directory
cd fengalotl

# Install the package
pip install -e .
```

---

## 📊 Data Setup

### Required Files

Place the `.h5ad` data files in the `data/` directory:

```
data/
├── Adult_metencephalon_rep1_2_DP8400015234BL_B1-2_region_ann.h5ad
├── Adult_metencephalon_rep3_DP8400015234BL_A3-1_region_ann.h5ad
├── Adult_olfactory_bulb_rep1_DP8400015234BL_A1-1_region_ann.h5ad
├── Adult_olfactory_bulb_rep2_DP8400015234BL_A2-2_region_ann.h5ad
├── Adult_pituitary_rep1_2_DP8400015234BL_B1-2_region_ann.h5ad
├── Adult_telencephalon_rep1_DP8400015234BL_A2-1_region_ann.h5ad
├── Adult_telencephalon_rep3_DP8400015234BL_A4-1_region_ann.h5ad
├── Adult_thalamencephalon_rep1_DP8400015234BL_A5-1_region_ann.h5ad
├── Adult_thalamencephalon_rep2_DP8400015234BL_A5-2_region_ann.h5ad
├── Adult_thalamencephalon_rep3_DP8400015234BL_A6-1_region_ann.h5ad
├── Adult_meta_DGE_markers.csv
├── genes.npy
└── samples.npy
```

### Gene Annotations

Gene annotations are automatically loaded from `data/Adult_meta_DGE_markers.csv`.

---

## 🖥️ Running the App

### On a Local Machine

```bash
# Activate the environment
mamba activate fengalotl

# Navigate to the project directory
cd fengalotl

# Run the Shiny app
shiny run src/fengalotl/app.py
```

Open your browser: **http://localhost:8000**

### On a Remote Server

1. Connect with port forwarding:
```bash
ssh -L 12345:localhost:8000 username@server
```

2. On the server, run:
```bash
mamba activate fengalotl
cd fengalotl
shiny run src/fengalotl/app.py --port 8000
```

3. Access locally at: **http://localhost:12345**

---

## 🎮 Usage Guide

1. **Select a dataset** from the dropdown menu
2. **Choose clustering** (Leiden clustering, Structure annotation, or Seurat clusters)
3. **Toggle cluster visualization** with the "Show clusters" switch
4. **Search for a gene** using the annotated gene names (e.g., "GLUL", "GAD1")
5. **Enable expression plotting** with the "Plot gene expression" switch
6. **Adjust visualization** using the dot size sliders
7. **Explore markers** in the differential expression accordion panel

---

## 📁 Project Structure

```
Fengalotl/
├── data/                       # H5AD data files
├── src/fengalotl/
│   ├── __init__.py
│   ├── _constants.py           # Configuration & gene annotations
│   ├── app.py                  # Main Shiny app entry point
│   ├── fct/
│   │   ├── expression.py       # Gene expression plotting
│   │   ├── load.py             # Data loading with caching
│   │   ├── spatial_widget.py   # Spatial plot functions
│   │   └── umap_widget.py      # PCA/UMAP plot functions
│   ├── js/
│   │   └── _format.py          # Dropdown formatting
│   └── mod/
│       ├── server.py           # Shiny server logic
│       └── ui.py               # Shiny UI definition
├── scripts/
│   └── create_tarball.sh       # Data packaging script
├── setup.py
├── pyproject.toml
└── README.md
```

---

## 🔧 Dependencies

| Package | Purpose |
|---------|---------|
| [Shiny for Python](https://shiny.posit.co/py/) | Web application framework |
| [Scanpy](https://scanpy.readthedocs.io/) | Single-cell analysis |
| [Plotly](https://plotly.com/python/) | Interactive visualizations |
| [Glasbey](https://github.com/lmcinnes/glasbey) | Color palette generation |
| [Pandas](https://pandas.pydata.org/) | Data manipulation |
| [NumPy](https://numpy.org/) | Numerical computing |

---

## 🙏 Acknowledgments

- **Adnan** for the template

---
