# 🧬 Ensemble Docking Score Workflow with Protein-Ligand Interaction Analysis

This repository contains a modular workflow for analyzing **ensemble docking scores** and evaluating protein-ligand interactions using the **AutoDock Vina scoring function**. The pipeline processes multiple protein conformations (frames), performs docking, scores ligand binding, and uses a custom PLI (Protein-Ligand Interaction) package to interpret atomic-level interactions.

---

## 🧠 Workflow Purpose

- **Goal**: Identify high-affinity ligands across an ensemble of protein structures.
- **Approach**: Use AutoDock Vina to dock ligands to multiple protein frames, then analyze interactions using a custom PLI package.
- **Applications**: Drug screening, flexible docking, structure-based design.

---

## 🔧 Workflow Steps

### `0.obabel.sh`
Converts ligand files to `.sdf` format using **Open Babel** for standardized input.

### `1.locol_add_sdf_score.py`
Performs **AutoDock Vina docking** across ensemble protein frames.  
- Calculates binding scores for each ligand-frame pair  
- Appends scores to ligand metadata

### `2.frame_pli.py`
Uses the **PLI package** to analyze atomic interactions for each docked pose.  
- Identifies bond types (hydrogen, ionic, hydrophobic, etc.)  
- Maps ligand atoms to protein groove geometry

### `3.h5_cvt_csv.py`
Converts intermediate `.h5` data to `.csv` format for inspection and downstream analysis.

### `4.cp_zeroless_csv.sh`
Filters out entries with zero or missing scores to clean the dataset.

### `5.cp_pli_csv.sh`
Copies and organizes final protein-ligand interaction `.csv` files for modeling or visualization.

---

## 📦 Requirements

- Python 3.8+
- Open Babel (`obabel`)
- AutoDock Vina
- NumPy, Pandas, h5py
- `PLI` package

---

## 🚀 Usage

```bash
bash 0.obabel.sh
python 1.locol_add_sdf_score.py
python 2.frame_pli.py
python 3.h5_cvt_csv.py
bash 4.cp_zeroless_csv.sh
bash 5.cp_pli_csv.sh
