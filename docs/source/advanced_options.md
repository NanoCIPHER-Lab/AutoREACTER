# Advanced Options

AutoREACTER gives users control over advanced reaction-generation behavior through the input JSON file. These options are intended for users who want to control reaction-product iteration and LAMMPS wildcard template generation.

The following advanced options are available for users:

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
      <td>Controls how many reaction-product iterations are attempted.</td>
    </tr>
    <tr>
      <td style="white-space:nowrap;"><a href="#wildcards" style="color:black;"><strong>wildcards</strong></a></td>
      <td style="color:black;">true</td>
      <td>Uses LAMMPS wildcard maps to reduce edge-template duplication.</td>
    </tr>
  </tbody>
</table>

Default JSON block:

```json
{
  "reaction_iteration_depth": 5,
  "wildcards": true
}
```

<p style="color:red;"><strong>Important: These options can change how many reaction products and LAMMPS reaction templates are generated.</strong></p>

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