## Getting Started - Pip Installation Guide

This guide explains how to install AutoREACTER from PyPI and run your first workflow.

### Step 1: Create a Python Virtual Environment

We recommend using a virtual environment so AutoREACTER and its dependencies do not interfere with your system Python installation.

```bash
python -m venv arx_env
```

Activate the environment:

On Linux/macOS:

```bash
source arx_env/bin/activate
```

On Windows:

```bash
arx_env\Scripts\activate
```

### Step 2: Install AutoREACTER

```bash
python -m pip install -U pip
python -m pip install AutoREACTER
```

This will install AutoREACTER and its required Python dependencies.

#### Step 2.1: Download and Prepare LUNAR (Prerequisite)

AutoREACTER requires the LUNAR package to handle atom typing. You must have this downloaded before running any examples.

Download LUNAR from: [https://github.com/CMMRLab/LUNAR](https://github.com/CMMRLab/LUNAR)

Note the Path: Keep track of the full directory path where LUNAR is saved on your computer (e.g., /home/user/software/LUNAR). During your first run, AutoREACTER will prompt you to enter this directory path. The path is then saved locally within the AutoREACTER.

Keep track of the full directory path where LUNAR is saved on your computer. During your first run, AutoREACTER will prompt you to enter this path.

### Step 3: Run AutoREACTER

You can either run by downloading [run_AutoREACTER.py](https://github.com/NanoCIPHER-Lab/AutoREACTER/blob/main/examples/run_AutoREACTER.py) with [example_1_inputs_count_mode.json](https://github.com/NanoCIPHER-Lab/AutoREACTER/blob/main/examples/example_1_inputs_count_mode.json):

```bash
python run_AutoREACTER.py -i example_1_inputs_count_mode.json
```
**Note:** Replace `example_1_inputs_count_mode.json` with the actual relative path of a JSON file in your computer.

or you can run [example_1.ipynb](https://github.com/NanoCIPHER-Lab/AutoREACTER/blob/main/examples/example_1.ipynb) with the input file.
