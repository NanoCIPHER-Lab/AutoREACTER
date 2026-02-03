# Running Example Notebooks with VS Code + WSL

This guide explains how to run the example notebooks in this folder using **VS Code + WSL**.

---

## 1. Create and Activate the Conda Environment

```bash
conda create -n mupt313 -y -c conda-forge python=3.13 numpy pandas rdkit ipykernel
conda activate mupt313
```

## 2. Install the Project (Editable Mode)

```bash
cd /mnt/c/Users/<USER>/Documents/GitHub/reaction_lammps_mupt
python -m pip install -U pip
python -m pip install -e .
```

## 3. Register the Jupyter Kernel

```bash
python -m ipykernel install --user --name mupt313 --display-name "Python (mupt313)"
```

## 4. Open the Project in VS Code

```bash
cd /mnt/c/Users/<USER>/Documents/GitHub/reaction_lammps_mupt
code .
```