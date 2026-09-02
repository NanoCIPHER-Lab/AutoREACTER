# Advanced Options

AutoREACTER gives users control over advanced workflow behavior through the input JSON file. These options are intended for users who want to control reaction-product iteration, LAMMPS wildcard template generation, and output-file placement.

The following user-facing advanced options are currently available:

<table>
  <thead>
    <tr>
      <th style="text-align:left;">Option</th>
      <th style="text-align:left;">Default</th>
      <th style="text-align:left;">Purpose</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="white-space:nowrap;"><a href="#reaction_iteration_depth" style="color:black;"><strong>reaction_iteration_depth</strong></a></td>
      <td style="color:black;">5</td>
      <td>Controls the maximum number of reaction-product iterations used during reaction progression.</td>
    </tr>
    <tr>
      <td style="white-space:nowrap;"><a href="#wildcards" style="color:black;"><strong>wildcards</strong></a></td>
      <td style="color:black;">true</td>
      <td>Generates LAMMPS wildcard map sections to reduce redundant edge-template cases.</td>
    </tr>
    <tr>
      <td style="white-space:nowrap;"><a href="#output_dir" style="color:black;"><strong>output_dir</strong></a></td>
      <td style="color:black;">AutoREACTER_outputs/&lt;simulation_name&gt;</td>
      <td>Controls where AutoREACTER writes generated output files.</td>
    </tr>
  </tbody>
</table>

Default behavior:

```json
{
  "reaction_iteration_depth": 5,
  "wildcards": true
}
```

By default, `output_dir` does not need to be provided. If it is omitted, AutoREACTER writes outputs next to the input JSON file using:

```text
AutoREACTER_outputs/<simulation_name>
```

<p style="color:red;"><strong>Important: These options can change the number of generated reaction products, reaction templates, LAMMPS map files, and output locations.</strong></p>
---

(reaction_iteration_depth)=
## reaction_iteration_depth

Inside AutoREACTER, reaction generation can be performed iteratively. First, AutoREACTER detects all possible reactions from the initial monomers and generates the corresponding products without looping. After that, if reaction iteration is enabled, AutoREACTER takes the generated products and checks whether they can react with any of the other monomers or products.

This process continues until the maximum number of iterations is reached or until no new products are generated from any combination.

The default value is `5`, which means AutoREACTER will attempt to generate products through up to five reaction iterations.

This step is implemented because reaction templates are generated based on both monomers and products. In multistage reactions, products from the first stage can react with other monomers or products to generate new products and new reaction templates.

If this option is disabled, AutoREACTER will only generate products from the initial monomers and will not check whether those products can react further. If products are not checked for further reactions, some reaction templates may be missed. In the case of small molecules or monomers, the length of the template can be limited by the size of the monomer, so polymerization propagation may be limited to dimers.

You can set this variable to `false` so AutoREACTER will not enter the reaction-iteration loop. You can also use a non-negative integer to set the maximum number of iterations. This helps control reaction explosion. If five iterations are not enough to generate all possible products, you can increase the number of iterations.

<p style="color:red;"><strong>Important: Setting reaction_iteration_depth too high can generate a large number of products and reaction templates.</strong></p>

The examples below show this behavior for an epoxy-amine reaction and glycine polymerization.

For the epoxy-amine example, the input monomers are:

```json
{
  "monomers": [
    {
      "name": "bisphenol_A_diglycidyl_ether",
      "smiles": "CC(C)(c1ccc(OCC2CO2)cc1)c1ccc(OCC2CO2)cc1"
    },
    {
      "name": "1,5-diaminopentane",
      "smiles": "NCCCCCN"
    }
  ]
}
```

```{image} _static/loop_figures/loop_monomers.png
:alt: reaction_iteration_epoxy_amine_monomers
:width: 100%
:align: center
```

Without `reaction_iteration_depth` enabled, the reaction templates will be generated as shown below.

```{image} _static/loop_figures/wo_loop_templates_template.png
:alt: reaction_iteration_without_loop_templates
:width: 100%
:align: center
```

With `reaction_iteration_depth` enabled, the reaction templates will be generated as shown below.

```{image} _static/loop_figures/loop_templates_template.png
:alt: reaction_iteration_with_loop_templates
:width: 100%
:align: center
```

For glycine polymerization, the input monomer is:

```json
{
  "monomers": [
    {
      "name": "glycine",
      "smiles": "NCC(=O)O"
    }
  ]
}
```

```{image} _static/loop_figures/glycine_templates_template.png
:alt: glycine_polymerization_reaction_templates
:width: 100%
:align: center
```

---

(wildcards)=
## wildcards

<p style="color:red;"><strong>Important: wildcards requires LAMMPS 22 July 2025 or later. The wildcard feature is not available in earlier versions of LAMMPS.</strong></p>

```{image} _static/wild_cards_figures/REACTER_wildcards_Diels-Alder-example_ed.avif
:alt: REACTER_wildcards_Diels-Alder-example
:width: 100%
:align: center
```

More information about REACTER wildcards can be found here:

- [Type Label Framework for Bonded Force Fields in LAMMPS](https://doi.org/10.1021/acs.jpcb.3c08419)
- [REACTER website](https://www.reacter.org/)

In LAMMPS bond/react, wildcards allow LAMMPS to infer the fourth atom used for dihedral typing from the first three atoms. This is useful for reactions where the fourth atom or edge atom is allowed to vary.

When `wildcards` is enabled, AutoREACTER generates reaction templates that include the fourth atoms or edge atoms as wildcards in the map file. This reduces the number of templates generated for reactions where the edge environment varies.

We highly recommend using wildcards for most reactions because it can reduce the number of generated templates and make the simulation more efficient.

However, if you want to control specific reactions during the reactive simulation, you may need to use templates without wildcards. In that case, keeping explicit edge atoms can help you control specific reaction sites more directly.

For normal vinyl reactions, AutoREACTER may need three templates without wildcards. With wildcard support enabled, the same simulation can often be performed using only two templates.

See the styrene vinyl reaction example below.

The following templates are generated without wildcards.

```{image} _static/wild_cards_figures/prior_templates_template.png
:alt: styrene_templates_without_wildcards
:width: 100%
:align: center
```

The following templates are generated with wildcards enabled.

```{image} _static/wild_cards_figures/wild_cards-templates_template.png
:alt: styrene_templates_with_wildcards
:width: 100%
:align: center
```

---

(output_dir)=
## output_dir

`output_dir` is an optional top-level input setting that controls where AutoREACTER writes generated output files.

If `output_dir` is not provided, AutoREACTER writes outputs next to the input JSON file using the default folder structure:

```text
AutoREACTER_outputs/<simulation_name>
```

Example:

```json
{
  "simulation_name": "Polyamide_Count_Mode_Basic",
  "force_field": "PCFF",
  "output_dir": "my_outputs",
  "simulations": [
    {
      "tag": "10k_300K",
      "temperature": 300,
      "density": 0.8,
      "monomer_counts": {
        "trimesoyl_chloride": 220,
        "m_phenylenediamine": 330
      }
    }
  ],
  "monomers": [
    {
      "name": "trimesoyl_chloride",
      "smiles": "ClC(=O)c1cc(cc(c1)C(Cl)=O)C(Cl)=O"
    },
    {
      "name": "m_phenylenediamine",
      "smiles": "C1=CC(=CC(=C1)N)N"
    }
  ]
}
```

For this example, AutoREACTER writes output to:

```text
my_outputs
```

Relative paths are resolved relative to the input JSON file location. Absolute paths are used directly.

Examples:

```json
{
  "output_dir": "AutoREACTER_outputs/polyamide_test"
}
```

```json
{
  "output_dir": "/home/user/AutoREACTER_runs/polyamide_test"
}
```

On Windows or WSL, Windows-style paths can also be used:

```json
{
  "output_dir": "C:/Users/Janitha/Documents/AutoREACTER_runs/polyamide_test"
}
```

<p style="color:red;"><strong>Important: If the output directory already exists, AutoREACTER may overwrite or clear generated workflow files from the previous run.</strong></p>