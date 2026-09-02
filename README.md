
<p align="center">
  <img src="https://raw.githubusercontent.com/NanoCIPHER-Lab/AutoREACTER/main/docs/source/_static/logo.png" alt="AutoREACTER logo" width="220">
</p>

<p align="center"><b>Automated generation of LAMMPS/REACTER-ready reaction-template workflows.</b></p>

> **Status:** AutoREACTER is currently in **v0.3** and under active development.
> APIs, configuration schemas, reaction libraries, and core functionality may change.

## Documentation

Full documentation is available at:

**[autoreacter.org](https://autoreacter.org/)**

The documentation covers installation, input configuration, supported reactions,
force fields, workflow options, and API usage.

## Installation

Install AutoREACTER from PyPI:

```bash
python -m pip install AutoREACTER
````

AutoREACTER also requires **LUNAR** for atom typing. See the
[Getting Started documentation](https://autoreacter.org/getting-started.html)
for setup instructions.

For development or source installation:

```bash
git clone https://github.com/NanoCIPHER-Lab/AutoREACTER.git
cd AutoREACTER
python -m pip install -e .
```

## Quick Start

After installing AutoREACTER, create a Python script such as `run_autoreacter.py`:

```python
import AutoREACTER as arx

arx.run("input.json")

arx.select_reactions()
arx.select_non_reactants()

arx.prepare_reactions()

session = arx.session()
print(f"Review generated images in: {session.images_dir}")

input("Press Enter to continue...")

arx.process()
```

Run it with:

```bash
python run_autoreacter.py
```

Example JSON input files and complete workflows are available in the
[`examples`](https://github.com/NanoCIPHER-Lab/AutoREACTER/tree/main/examples)
directory.

For users working directly from the source repository, the included example
runner can also be used:

```bash
python examples/example_1.py -i examples/polyamide_count_mode_basic.json
```

## Help and Support

For bugs, reaction requests, or force-field support requests, please open an issue:

**[AutoREACTER Issues](https://github.com/NanoCIPHER-Lab/AutoREACTER/issues)**

## License

AutoREACTER is released under the **MIT License**.
See [LICENSE](https://github.com/NanoCIPHER-Lab/AutoREACTER/blob/main/LICENSE.md).
