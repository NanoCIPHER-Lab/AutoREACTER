<p align="center">
  <img src="https://raw.githubusercontent.com/NanoCIPHER-Lab/AutoREACTER/main/docs/source/_static/logo.png" alt="AutoREACTER logo" width="220">
</p>

<h1 align="center">AutoREACTER</h1>

<p align="center"><b>Automated generation of LAMMPS/REACTER-ready reaction-template workflows.*</b></p>

> **Status:** AutoREACTER is currently in **v0.2-beta** and under active development. APIs, configuration schemas, reaction libraries, and core functionality may change without notice.
> Please refer to the [changelog](https://autoreacter.org/change_log.html) for the latest updates.

## Documentation

Full documentation is available at:

**[autoreacter.org](https://autoreacter.org/)**

The documentation includes installation instructions, input configuration, supported reactions, supported force fields, cleanup utilities, and developer API references.

For detailed functional-group mapping and reaction rules, see the [supported reactions documentation](https://autoreacter.org/supported-reactions.html).

## Installation

AutoREACTER can be installed directly from PyPI:

```bash
python -m pip install AutoREACTER
```

For users who want to modify the source code or run the latest development version, AutoREACTER can also be installed from source:

```bash
git clone https://github.com/NanoCIPHER-Lab/AutoREACTER.git
cd AutoREACTER
python -m pip install -e .
```

AutoREACTER also requires **LUNAR** for atom typing. See the [getting started documentation](https://autoreacter.org/getting-started.html) for the full setup guide.

## Quick start

Run AutoREACTER with a JSON input file:

```bash
python examples/run_AutoREACTER.py -i examples/example_1_inputs_count_mode.json
```

or:

```bash
python examples/run_AutoREACTER.py --input examples/example_1_inputs_count_mode.json
```

View available commands and options:

```bash
python examples/run_AutoREACTER.py --help
```

## Interactive notebook workflow

AutoREACTER can also be used through Jupyter notebooks for an interactive, visual, step-by-step workflow. This mode is recommended for inspecting monomers, functional groups, reaction templates, and generated LAMMPS setup files before running larger workflows.

See the examples directory for notebooks and usage notes:

**[examples/README.md](https://autoreacter.org/getting_started_source_installation.html)**

## Help and support

If you find a bug, need a new reaction type, or want to request additional force-field support, please open an issue:

**[AutoREACTER Issues](https://github.com/NanoCIPHER-Lab/AutoREACTER/issues)**

## License

AutoREACTER is released under the **MIT License**. See [LICENSE](https://github.com/NanoCIPHER-Lab/AutoREACTER/blob/main/LICENSE.md) for details.
