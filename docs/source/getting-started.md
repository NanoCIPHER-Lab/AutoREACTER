![AutoREACTER Logo](_static/logo.png)
# Getting Started 

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
conda create -n autoRX -y -c conda-forge python=3.13 numpy pandas rdkit ipykernel networkx
```

*Once the installation is complete, activate your new environment:*

```bash
conda activate autoRX
```
**Note**: You must run conda activate autoRX every time you open a new terminal to use AutoREACTER.

#### Step 2.1: Download and Prepare LUNAR (Prerequisite)
AutoREACTER requires the LUNAR package to handle atom typing. You must have this downloaded before running any examples. 

Download LUNAR: Clone or download the repository from [https://github.com/CMMRLab/LUNAR](https://github.com/CMMRLab/LUNAR).

Note the Path: Keep track of the full directory path where LUNAR is saved on your computer (e.g., /home/user/software/LUNAR). During your first run, AutoREACTER will prompt you to enter this directory path. The path is then saved locally within the AutoREACTER.

### Step 3: Run Your First Example
AutoREACTER provides two different ways to build your LAMMPS reaction files: \
&emsp; 1. An interactive, visual Jupyter Notebook. \
&emsp; 2. A fast, Command-Line Interface (CLI)
   
#### Option A: The Interactive Notebook (Recommended for Beginners)

If you want to see exactly how AutoREACTER detects functional groups, maps templates, and handles non-reactive monomers, the Jupyter Notebook is the best place to start.

*Navigate to the examples directory:*

```bash
cd examples
```

*Install the project in editable mode*

```bash
python -m pip install -U pip
python -m pip install -e ..
```
*Register the Jupyter Kernel*

```bash
python -m ipykernel install --user --name autoRX --display-name "Python (autoRX)"
```

*Open the project in VS Code*

```bash
code .
```

Run the cells sequentially. The notebook will guide you step-by-step.

#### Option B: The Automated CLI

If want to generate the LAMMPS files quickly, you can run AutoREACTER directly from the terminal.

Ensure you are in the root AutoREACTER directory, then run the `AutoREACTER.py` script with relative path to an your `input.json` file:

```bash
python AutoREACTER.py -i examples/example_1_inputs_ratio_mode.json
```

**Note:** Replace `examples/example_1_inputs_ratio_mode.json` with the actual relative path of a JSON file in your computer.

AutoREACTER will read the JSON, process the chemistry, output your LAMMPS setup scripts into a newly generated, timestamped run directory.

**Note:** A timestamped run directory means that AutoREACTER creates a folder at:

`root_directory/cache/<today>/N`

- `<today>` is the current date in `yyyy-mm-dd` format  
- `N` is an incrementing integer starting from 1 (e.g., 1, 2, 3, ...)

