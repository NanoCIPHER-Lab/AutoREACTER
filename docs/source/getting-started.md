## Getting Started

This guide explains how to install AutoREACTER from PyPI and run your first workflow.

### Step 1: Create a Python Virtual Environment

We recommend using a virtual environment so AutoREACTER and its dependencies do not interfere with your system Python installation.

```bash
python -m venv autoRX
```

Activate the environment:

On Linux/macOS:

```bash
source autoRX/bin/activate
```

On Windows:

```bash
autoRX\Scripts\activate
```

### Step 2: Install AutoREACTER

```bash
python -m pip install -U pip
python -m pip install AutoREACTER
```

This will install AutoREACTER and its required Python dependencies.

### Step 3: Download and Prepare LUNAR

AutoREACTER requires LUNAR for atom typing and REACTER file preparation.

Download LUNAR from:

https://github.com/CMMRLab/LUNAR

Keep track of the full directory path where LUNAR is saved on your computer. During your first run, AutoREACTER will prompt you to enter this path.

### Step 4: Run AutoREACTER

You can either run by downloading [`AutoREACTER.py`](https://github.com/NanoCIPHER-Lab/AutoREACTER/blob/main/AutoREACTER.py) with [`example_1_inputs_count_mode.json`](https://github.com/NanoCIPHER-Lab/AutoREACTER/blob/main/examples/example_1_inputs_count_mode.json):

```bash
AutoREACTER.py -i example_1_inputs_count_mode.json
```
**Note:** Replace `examples/example_1_inputs_ratio_mode.json` with the actual relative path of a JSON file in your computer.

or you can run [`example_1.ipynb`](https://github.com/NanoCIPHER-Lab/AutoREACTER/blob/main/examples/example_1.ipynb) with the input file.
