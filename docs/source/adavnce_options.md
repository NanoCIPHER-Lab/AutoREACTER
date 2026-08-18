## Advanced Options

AutoREACTER gives users control over reaction generation through the input JSON file. These options allow users to tune how aggressively reactions are searched, how many reaction-growth iterations are attempted, how reaction templates are simplified, and whether one or two reaction stages are written for LAMMPS.

### Deep search

`deep_search` controls whether AutoREACTER performs iterative reaction discovery beyond the first set of directly detected reactions.

When `deep_search` is disabled, AutoREACTER detects reactions only from the functional groups present in the initial input monomers. This is useful for simple systems where only the first reaction event is needed.

When `deep_search` is enabled, AutoREACTER uses newly generated reaction products as possible reactants for later reaction searches. This allows the workflow to discover chain-growth, crosslinking, and multi-generation polymerization reactions that are not visible from the starting monomers alone.

For example, in a vinyl polymerization workflow, the first reaction may generate an active chain-end radical. That new radical product can then react with another vinyl monomer in the next search iteration. Without deep search, AutoREACTER may only detect the initial reaction. With deep search enabled, AutoREACTER can continue discovering propagation-style reactions from newly generated products.

```json
{
  "deep_search": true,
  "reaction_iteration_depth": 5
}