# Running Example Notebooks with VS Code + WSL

This guide explains how to run the example notebooks in this folder using **VS Code + WSL**.

---

## 1. Create and Activate the Conda Environment

```bash
conda create -n autoRX -y -c conda-forge python=3.13 numpy pandas rdkit ipykernel networkx
conda activate autoRX
```

## 2. Install the Project (Editable Mode)

```bash
cd ....../reaction_lammps_mupt
python -m pip install -U pip
python -m pip install -e .
```

## 3. Register the Jupyter Kernel

```bash
python -m ipykernel install --user --name autoRX --display-name "Python (autoRX)"
```

## 4. Open the Project in VS Code

```bash
code .
```