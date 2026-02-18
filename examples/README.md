# Running Example Notebooks with VS Code + WSL

This guide explains how to run the example notebooks in this folder using **VS Code + WSL**.

---
## 1. Installation
Clone the repository and install:

```bash
git clone https://github.com/NanoCIPHER-Lab/AutoREACTER
cd AutoREACTER
pip install -e .
```

## 2. Create and Activate the Conda Environment

```bash
conda create -n autoRX -y -c conda-forge python=3.13 numpy pandas rdkit ipykernel networkx
conda activate autoRX
```

## 3. Install the Project (Editable Mode)

```bash
cd AutoREACTER
python -m pip install -U pip
python -m pip install -e .
```

## 4. Register the Jupyter Kernel

```bash
python -m ipykernel install --user --name autoRX --display-name "Python (autoRX)"
```

## 5. Open the Project in VS Code

```bash
code .
```
