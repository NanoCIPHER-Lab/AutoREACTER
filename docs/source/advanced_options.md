# Advanced Options

AutoREACTER gives users control over reaction generation through the input JSON file. These options allow users to tune how reaction products are iteratively searched, how reaction templates are compared, how duplicate templates are handled, how LAMMPS wildcard templates are generated, and whether one or two LAMMPS reaction stages are written.

```json
{
  "simulation_name": "Polyamide_Deep_Search_Demo",
  "force_field": "PCFF",
  "deep_search": true,
  "reaction_iteration_depth": false,
  "wildcards": false,
  "deduplicate_reaction_templates": true,
  "write_second_reaction_stage": true,

  "simulations": [
    {
      "tag": "100k_300K",
      "temperature": 300,
      "density": 0.8,
      "monomer_counts": {
        "4-methylheptanedioyl_dichloride": 2380,
        "m-phenylenediamine": 3570,
        "heptanedioyl_dichloride": 100
      }
    }
  ],

  "monomers": [
    {
      "name": "4-methylheptanedioyl_dichloride",
      "smiles": "CC(CCC(=O)Cl)CCC(=O)Cl"
    },
    {
      "name": "heptanedioyl_dichloride",
      "smiles": "O=C(Cl)CCCCCC(=O)Cl"
    },
    {
      "name": "m-phenylenediamine",
      "smiles": "C1=CC(=CC(=C1)N)N"
    }
  ]
}
```

The following advanced options are currently available:

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
      <td style="white-space:nowrap;"><a href="#deep_search" style="color:black;"><strong>deep_search</strong></a></td>
      <td style="color:black;">true</td>
      <td>Controls stricter RDKit/NetworkX template deduplication before atom typing.</td>
    </tr>
    <tr>
      <td style="white-space:nowrap;"><a href="#deduplicate_reaction_templates" style="color:black;"><strong>deduplicate_reaction_templates</strong></a></td>
      <td style="color:black;">true</td>
      <td>Removes duplicate reaction templates before writing LAMMPS files.</td>
    </tr>
    <tr>
      <td style="white-space:nowrap;"><a href="#wildcards" style="color:black;"><strong>wildcards</strong></a></td>
      <td style="color:black;">true</td>
      <td>Uses LAMMPS wildcard maps to reduce edge-template duplication.</td>
    </tr>
    <tr>
      <td style="white-space:nowrap;"><a href="#write_second_reaction_stage" style="color:black;"><strong>write_second_reaction_stage</strong></a></td>
      <td style="color:black;">false</td>
      <td>Controls whether the second LAMMPS reaction stage is written.</td>
    </tr>
  </tbody>
</table>

Default JSON block:

```text
{
  "reaction_iteration_depth": 5,
  "deep_search": true,
  "deduplicate_reaction_templates": true,
  "wildcards": true,
  "write_second_reaction_stage": false
}
```

<p style="color:red;"><strong>Important: These options affect template generation and LAMMPS input generation, so changing them can change the number of generated reaction templates and the final LAMMPS workflow.</strong></p>

---

(reaction_iteration_depth)=
## reaction_iteration_depth

Inside AutoREACTER, reaction generation can be performed iteratively. First, AutoREACTER detects all possible reactions from the initial monomers and generates the corresponding products without looping. After that, if reaction iteration is enabled, AutoREACTER takes the generated products and checks whether they can react with any of the other monomers or products.

This process continues until the maximum number of iterations is reached or until no new products are generated from any combination.

The default value is `5`, which means AutoREACTER will attempt to generate products through up to five reaction iterations.

This step is implemented because reaction templates are generated based on both monomers and products. In multistage reactions, products from the first stage can react with other monomers or products to generate new products and new reaction templates.

If this option is disabled, AutoREACTER will only generate products from the initial monomers and will not check whether those products can react further. If products are not checked for further reactions, some reaction templates may be missed. In the case of small molecules or monomers, the length of the template can be limited by the size of the monomer, so polymerization propagation may be limited to dimers.

You can set this variable to `false` so AutoREACTER will not enter the reaction-iteration loop. You can also use a non-negative integer to set the maximum number of iterations. This helps control reaction explosion. If five iterations are not enough to generate all possible products, you can increase the number of iterations.

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

(deep_search)=
## deep_search

`deep_search` controls how strictly AutoREACTER compares reaction templates during RDKit/NetworkX-based deduplication.

During early reaction-template generation, AutoREACTER has not yet assigned force-field atom types. At this stage, RDKit only knows the chemical graph, not the final force-field atom types. Because of this, two templates can look identical to RDKit even though they may later receive different atom types after force-field assignment.

When `deep_search` is enabled, AutoREACTER extends the graph comparison by one additional neighbor beyond the normal template edge atoms. In practical terms, this means the comparison also checks the “5th atom” outside the normal reaction-template distance. This extra check helps prevent two templates from being incorrectly treated as duplicates when their immediate reaction core is the same but their outer chemical environment is different.

This option is useful because atom typing is not performed during RDKit reaction detection. The extra graph environment acts as a substitute for missing atom-type information during deduplication.

For example, two reaction templates may have the same reacting atoms and the same four-atom dihedral-distance template core, but differ at the edge, as shown in the example below.

Using the following polyamidation reaction as an example:

```json
{
  "monomers": [
    {
      "name": "4-methylheptanedioyl_dichloride",
      "smiles": "CC(CCC(=O)Cl)CCC(=O)Cl"
    },
    {
      "name": "heptanedioyl_dichloride",
      "smiles": "O=C(Cl)CCCCCC(=O)Cl"
    },
    {
      "name": "m-phenylenediamine",
      "smiles": "C1=CC(=CC(=C1)N)N"
    }
  ]
}
```

The input monomers are shown below.

```{image} _static/deep_search_figures/input_monomers.png
:alt: deep_search_demo_monomers
:width: 100%
:align: center
```

Without deep search, the reaction between `4-methylheptanedioyl_dichloride` and `m-phenylenediamine` can be treated as a duplicate of the reaction between `heptanedioyl_dichloride` and `m-phenylenediamine`, because the reaction core up to the normal dihedral-distance template region is identical.

```{image} _static/deep_search_figures/w_o_deep_search.png
:alt: without_deep_search_template_deduplication
:width: 100%
:align: center
```

With deep search enabled, the additional graph-environment check detects that the two templates are chemically different and does not treat them as duplicates.

```{image} _static/deep_search_figures/w_deep_search.png
:alt: with_deep_search_template_deduplication
:width: 100%
:align: center
```

However, there can still be edge cases, as shown in the example below.

```json
{
  "monomers": [
    {
      "name": "branched_alkenyl_diacid_chloride",
      "smiles": "C=CCC(CCCCCCC(CC=C)C(=O)Cl)C(=O)Cl"
    },
    {
      "name": "m-phenylenediamine",
      "smiles": "C1=CC(=CC(=C1)N)N"
    },
    {
      "name": "extended_branched_alkenyl_diacid_chloride",
      "smiles": "C/C=C/CC(CCCCCCC(C/C=C/C)C(=O)Cl)C(=O)Cl"
    }
  ]
}
```

```{image} _static/deep_search_figures/edge_monomers.png
:alt: deep_search_edge_case_monomers
:width: 100%
:align: center
```

In this case, the template contains double-bonded carbons. In one molecule, the double bond is part of a vinyl group, and the carbon bonded to the alkene chain will be typed as `c=1`, which is part of the template. In the other molecule, the corresponding carbon is bonded to another alkene chain. In this case, both carbons will be typed as `c=2`, including the carbon in the template.

Even with deep search, AutoREACTER will only generate the template shown below.

```{image} _static/deep_search_figures/edge_templates_template.png
:alt: deep_search_edge_case_generated_template
:width: 100%
:align: center
```

The following template will be missed.

```{image} _static/deep_search_figures/missed_templates_template.png
:alt: deep_search_edge_case_missed_template
:width: 100%
:align: center
```

In this situation, deep search is not enough to detect the difference between the two templates, so the templates can still be treated as duplicates.

With wildcard search enabled, this issue can be eliminated. See the `wildcards` section below.

---

(deduplicate_reaction_templates)=
## deduplicate_reaction_templates

AutoREACTER first generates all possible reactions and then generates LAMMPS reaction templates for those reactions. Some of these templates can be duplicates. The `deduplicate_reaction_templates` option controls whether AutoREACTER removes those duplicate templates.

This deduplication step is one of the final steps in the AutoREACTER workflow. If `deduplicate_reaction_templates` is set to `false`, AutoREACTER will keep all generated templates, including duplicates.

We recommend keeping this option set to `true` so AutoREACTER only writes unique templates.

<p style="color:red;"><strong>Important: If deduplicate_reaction_templates is set to false, AutoREACTER may write many redundant LAMMPS reaction templates.</strong></p>

Example output:

```text
Preparing templates and map file for reaction ID: 1
Preparing templates and map file for reaction ID: 2
Preparing templates and map file for reaction ID: 4
Preparing templates and map file for reaction ID: 8
Preparing templates and map file for reaction ID: 9
Duplicate template disabled: RXN_8
Duplicate template disabled: RXN_9
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

(write_second_reaction_stage)=
## write_second_reaction_stage

`write_second_reaction_stage` controls whether AutoREACTER writes the second LAMMPS reaction-stage input file.

When this option is set to `true`, AutoREACTER writes both reaction stages:

```text
3_reaction_first_stage
4_reaction_second_stage
```

The first reaction stage usually uses a shorter reaction cutoff and is intended to capture close-contact reactions. The second reaction stage uses a larger reaction cutoff and allows additional reactions to occur after the first-stage structure has already reacted.

When this option is set to `false`, AutoREACTER skips the second reaction stage and writes only the first reaction stage. In that case, the post-equilibration stage will use the output from the first reaction stage as its input.

Use `write_second_reaction_stage: true` when you want a two-stage reaction workflow. Use `write_second_reaction_stage: false` when you want a shorter workflow or when the first reaction stage is sufficient for the system.

Example:

```json
{
  "write_second_reaction_stage": true
}
```