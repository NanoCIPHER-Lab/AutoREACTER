# AutoREACTER

**AutoREACTER** is a Python-based toolkit for automating reaction-based molecular system generation for LAMMPS workflows. It is developed as part of the **Multiscale Polymer Toolkit (MuPT)**.

> **Status:** AutoREACTER is currently in **v0.2-beta** and under active development. APIs, configuration schemas, reaction libraries, and core functionality may change without notice.

## Documentation

Full documentation is available at:

**[https://nanocipher-lab.github.io/AutoREACTER/](https://nanocipher-lab.github.io/AutoREACTER/)**

The documentation includes installation instructions, input configuration, supported reactions, supported force fields, cleanup utilities, and developer API references.

## Current reaction support

AutoREACTER currently supports beta-stage step-growth **polycondensation** workflows for:

### Polyesterification

- Hydroxy-carboxylic acid polycondensation
- Hydroxy acid halide polycondensation
- Diol + di-acid halide polycondensation
- Diol + di-carboxylic acid polycondensation

### Polyamidation

- Amino acid polycondensation
- Diamine + di-carboxylic acid polycondensation
- Diamine + di-carboxylic acid halide polycondensation

For the detailed functional-group mapping and reaction rules, see the [supported reactions documentation](https://nanocipher-lab.github.io/AutoREACTER/supported-reactions.html).

## Installation

Clone the repository:

```bash
git clone https://github.com/NanoCIPHER-Lab/AutoREACTER.git
cd AutoREACTER
```

Create and activate the recommended Conda environment:

```bash
conda create -n autoRX -y -c conda-forge python=3.13 numpy pandas rdkit ipykernel networkx
conda activate autoRX
```

Install AutoREACTER in editable mode:

```bash
python -m pip install -U pip
python -m pip install -e .
```

AutoREACTER also requires **LUNAR** for atom typing. See the [getting started documentation](https://nanocipher-lab.github.io/AutoREACTER/getting-started.html) for setup details.

## Quick start

Run AutoREACTER with a JSON input file:

```bash
python AutoREACTER.py -i path/to/input.json
```

or:

```bash
python AutoREACTER.py --input path/to/input.json
```

Example:

```bash
python AutoREACTER.py -i examples/example_1_inputs_ratio_mode.json
```

View available commands and options:

```bash
python AutoREACTER.py --help
```

## Interactive notebook workflow

AutoREACTER can also be used through Jupyter notebooks for an interactive, visual, step-by-step workflow. This mode is recommended for inspecting monomers, functional groups, reaction templates, and generated LAMMPS setup files before running larger workflows.

See the examples directory for notebooks and usage notes:

**[examples/README.md](https://github.com/NanoCIPHER-Lab/AutoREACTER/blob/main/examples/README.md)**

## Cleanup utility

Delete cached runs older than a given number of days:

```bash
python AutoREACTER.py --cleanup 30
```

Delete all cached runs:

```bash
python AutoREACTER.py --cleanup all
```

Short flags are also available:

```bash
python AutoREACTER.py -c 30
python AutoREACTER.py -c all
```

## Help and support

If you find a bug, need a new reaction type, or want to request additional force-field support, please open an issue:

**[AutoREACTER Issues](https://github.com/NanoCIPHER-Lab/AutoREACTER/issues)**

## License

AutoREACTER is released under the **MIT License**. See [LICENSE](https://github.com/NanoCIPHER-Lab/AutoREACTER/blob/main/LICENSE.md) for details.