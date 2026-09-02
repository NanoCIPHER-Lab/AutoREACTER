# Template Deduplication Options

This page describes AutoREACTER options related to reaction-template comparison and duplicate-template removal. These options are mainly useful for developers or advanced users who need to understand how AutoREACTER decides whether two generated reaction templates are equivalent.

The following template-deduplication options are available:

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
      <td style="white-space:nowrap;"><a href="#deep_search" style="color:black;"><strong>deep_search</strong></a></td>
      <td style="color:black;">true</td>
      <td>Controls stricter RDKit/NetworkX template deduplication before atom typing.</td>
    </tr>
    <tr>
      <td style="white-space:nowrap;"><a href="#deduplicate_reaction_templates" style="color:black;"><strong>deduplicate_reaction_templates</strong></a></td>
      <td style="color:black;">true</td>
      <td>Removes duplicate reaction templates before writing LAMMPS files.</td>
    </tr>
  </tbody>
</table>

Default JSON block:

```json
{
  "deep_search": true,
  "deduplicate_reaction_templates": true
}
```

<p style="color:red;"><strong>Important: These options affect reaction-template comparison and duplicate-template removal before LAMMPS files are written.</strong></p>

---

(deep_search)=
## deep_search

`deep_search` controls how strictly AutoREACTER compares reaction templates during RDKit/NetworkX-based deduplication.

During early reaction-template generation, AutoREACTER has not yet assigned force-field atom types. At this stage, RDKit only knows the chemical graph, not the final force-field atom types. Because of this, two templates can look identical to RDKit even though they may later receive different atom types after force-field assignment.

When `deep_search` is enabled, AutoREACTER extends the graph comparison by one additional neighbor beyond the normal template edge atoms. In practical terms, this means the comparison also checks the “5th atom” outside the normal reaction-template distance. This extra check helps prevent two templates from being incorrectly treated as duplicates when their immediate reaction core is the same but their outer chemical environment is different.

This option is useful because atom typing is not performed during RDKit reaction detection. The extra graph environment acts as a substitute for missing atom-type information during deduplication.

<p style="color:red;"><strong>Important: deep_search helps distinguish templates before atom typing, but it is not a full replacement for force-field atom typing.</strong></p>

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

With wildcard search enabled, this issue can be eliminated. See the {doc}`advanced_options` page.

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
