![AutoREACTER Logo](_static/logo.png)
## Input Configuration.

The `input.json` file is the one and only input for your AutoREACTER workflow.  
It tells the system:

1. What molecules you are using  
2. Which force field to apply  
3. How many atoms or molecules should be in the system  
4. Density of the simulation  
5. Temperature of the simulation

AutoREACTER v0.2 introduces a highly flexible replica system. You can define multiple distinct simulation setups  within a single run. You must explicitly define parameters like `temperature` and `density` for **each replica**. This method ensures AutoREACTER generates only the exact LAMMPS files you need.

There are two primary ways to define your system composition: **Ratio Mode** and **Count Mode**.

---

### 1. Global Settings & Monomers

Regardless of which mode you use, every `input.json` needs global settings and a list of your chemical building blocks.

* **`simulation_name`**: A string used to name your output directories and files.
* **`force_field`**: *(Optional)* Specifies the force field LUNAR should use to assign atom types and charges (e.g., `"PCFF"`). If omitted, AutoREACTER defaults to `"PCFF"`.

---

### 2. Defining Replicas: Ratio Mode vs. Count Mode

When defining the `replicas` array, you have to decide how you want to calculate the number of molecules in the simulation box. 

#### Option A: Ratio Mode (Target Atom Count)
Use Ratio Mode when you know the total size of the simulation you want to run (e.g., ~10,000 atoms) and the ratio of your molecules, but you don't want to calculate the exact number of individual molecules by hand.

AutoREACTER will automatically calculate the correct number of molecules to hit your `total_atoms` target based on your `monomer_ratios`.

**Example `input.json` (Ratio Mode):**
```json
{
    "simulation_name": "Example_Ratio_Mode",
    "force_field": "PCFF",
    "replicas": [
        {
            "tag": "10k_base",
            "temperature": 300,
            "density": 0.8,
            "total_atoms": 10000,
            "monomer_ratios": {
                "tmc": 1.0,
                "mpd": 1.0,
                "ethanol": 0.5   
            }
        },
        {
            "tag": "100k_base",
            "temperature": 400,
            "density": 0.8,
            "total_atoms": 100000,
            "monomer_ratios": {
                "tmc": 1.0,
                "mpd": 1.0,
                "ethanol": 0.5
            }
        }
    ],
    "monomers": [
        { 
            "name": "tmc", 
            "smiles": "ClC(=O)c1cc(cc(c1)C(Cl)=O)C(Cl)=O" },
        {   
            "name": "mpd", 
            "smiles": "C1=CC(=CC(=C1)N)N"
        },
        { 
            "name": "ethanol", 
            "smiles": "CCO" 
        }
    ]
}
```

#### Option B: Count Mode (Exact Molecule Count)
Use Count Mode when you know number of molecules instead of number of atoms. Instead of providing a total atom target, you explicitly define the exact number of each molecule using `monomer_counts`.

**Example `input.json` (Count Mode):**
```json
{
    "simulation_name": "Example_Count_Mode",
    "force_field": "PCFF",
    "replicas": [
        {
            "tag": "10k",
            "temperature": 300,
            "density": 0.8,
            "monomer_counts": {
                "tmc": 220,
                "mpd": 220,
                "ethanol": 110
            }
        },
        {
            "tag": "100k",
            "temperature": 400,
            "density": 0.8,
            "monomer_counts": {
                "tmc": 2200,
                "mpd": 2200,
                "ethanol": 1100
            }
        }
    ],
    "monomers": [
        { 
            "name": "tmc", 
            "smiles": "ClC(=O)c1cc(cc(c1)C(Cl)=O)C(Cl)=O" },
        {   
            "name": "mpd", 
            "smiles": "C1=CC(=CC(=C1)N)N"
        },
        { 
            "name": "ethanol", 
            "smiles": "CCO" 
        }
    ]
}
```
---

### 3. Replica Parameters Breakdown

For each object inside the `replicas` list, you must define:

* **`tag`**: A unique label for the replica (e.g., `"10k_base"`). This will be used to name the output subdirectories.
* **`temperature`**: The system temperature in Kelvin. AutoREACTER uses this to configure the LAMMPS input scripts.
* **`density`**: The initial target density of the simulation box (in g/cm³).
* **Composition Setup:** Either provide `total_atoms` alongside `monomer_ratios` (Ratio Mode) **OR** provide `monomer_counts` (Count Mode). **Do not mix them in the same replica**.
---
### 4. Specifying the monomers or molecules.

In the monomers section, you define each molecule as a dictionary entry.
Each monomer (or molecule) must include:
1. A name field
2. A **valid SMILES string (smiles)**

```json
"monomers": [
        { 
            "name": "tmc", 
            "smiles": "ClC(=O)c1cc(cc(c1)C(Cl)=O)C(Cl)=O" },
        {   
            "name": "mpd", 
            "smiles": "C1=CC(=CC(=C1)N)N"
        },
        { 
            "name": "ethanol", 
            "smiles": "CCO" 
        }
    ]
```

**<u>IMPORTANT</u>**: The `monomers` section must remain consistent with the `monomer_counts` defined in each replica. All name tags must match exactly, and every monomer listed must have a corresponding count in each replica otherwise AutoREACTER will **raise an error** before proceeding with the chemistry. 

**Note:** You can use [SMILES Generator / Checker](https://www.cheminfo.org/flavor/malaria/Utilities/SMILES_generator___checker/index.html) to generate valid SMILES strings. If SMILES strings are incorrect AutoREACTER will **raise an error** before proceeding.

