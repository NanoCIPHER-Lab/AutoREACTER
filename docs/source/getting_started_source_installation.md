
## Getting Started - Source Installation Guide

This guide will guide you through setting up your environment, installing the required dependencies, and running reaction modeling setup for REACTER.

AutoREACTER relies heavily on cheminformatics libraries like RDKit, and numeric computing libraries like Pandas, so we strongly recommend using **Conda** to manage your Python environment.

---

### Step 1: Clone the Repository

Download the AutoREACTER source code to your computer using Git. Open your terminal and run:

```bash
git clone https://github.com/NanoCIPHER-Lab/AutoREACTER.git
cd AutoREACTER
```

### Step 2: Set Up the Conda Environment

You need to create an environment containing Python 3.13(recommend) and the python libraries AutoREACTER needs to function.

*Run the following command to create a new Conda environment:*

```bash
conda create -n arx_env -y -c conda-forge python=3.13
```

*Once the installation is complete, activate your new environment:*

```bash
conda activate arx_env
```

**Note**: You must run conda activate arx_env every time you open a new terminal to use AutoREACTER.

*Install the required Python libraries using `requirements.txt`:*

```bash
python -m pip install -U pip
python -m pip install -r requirements.txt
```

#### Step 2.1: Download and Prepare LUNAR (Prerequisite)

AutoREACTER requires the LUNAR package to handle atom typing. You must have this downloaded before running any examples.

Download LUNAR from: [https://github.com/CMMRLab/LUNAR](https://github.com/CMMRLab/LUNAR)

Note the Path: Keep track of the full directory path where LUNAR is saved on your computer (e.g., /home/user/software/LUNAR). During your first run, AutoREACTER will prompt you to enter this directory path. The path is then saved locally within the AutoREACTER.

### Step 3: Run Your First Example

AutoREACTER provides two different ways to build your LAMMPS reaction files: 
- An interactive, visual Jupyter Notebook. 
- A fast, Command-Line Interface (CLI)

#### Option A: The Interactive Notebook.

If you want to see exactly how AutoREACTER detects functional groups, maps templates, and handles non-reactive monomers, the Jupyter Notebook is the best place to start.

*Install the project in editable mode*

```bash
python -m pip install -e .
```

*Register the Jupyter Kernel*

```bash
python -m ipykernel install --user --name arx_env --display-name "Python (arx_env)"
```

*Open the project in VS Code*

```bash
code .
```

Run the cells sequentially. The notebook will guide you step-by-step.

#### Option B: The Automated CLI

If want to generate the LAMMPS files quickly, you can run AutoREACTER directly from the terminal.

Ensure you are in the root AutoREACTER directory, then run the [`run_AutoREACTER.py`](https://github.com/NanoCIPHER-Lab/AutoREACTER/blob/main/examples/run_AutoREACTER.py) script with relative path to an your `input.json` file:

```bash
python examples/run_AutoREACTER.py -i examples/example_1_inputs_count_mode.json
```

**Note:** Replace `examples/example_1_inputs_count_mode.json` with the actual relative path of a JSON file in your computer.

AutoREACTER parses the JSON, processes the chemistry, and exports all LAMMPS scripts to a new directory named after your specific simulation.

### Step 4: Run AutoREACTER

You can either run by downloading [run_AutoREACTER.py](https://github.com/NanoCIPHER-Lab/AutoREACTER/blob/main/examples/run_AutoREACTER.py) with [`example_1_inputs_count_mode.json`](https://github.com/NanoCIPHER-Lab/AutoREACTER/blob/main/examples/example_1_inputs_count_mode.json):

```bash
python run_AutoREACTER.py -i example_1_inputs_count_mode.json
```

**Note:** Replace `example_1_inputs_count_mode.json` with the actual relative path of a JSON file in your computer.

or you can run [example_1.ipynb](https://github.com/NanoCIPHER-Lab/AutoREACTER/blob/main/examples/example_1.ipynb) with the input file.
