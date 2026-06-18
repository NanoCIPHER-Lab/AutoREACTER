# API Reference

This documentation provides users with the public APIs that they can use with AutoREACTER. AutoREACTER has **10 public user-facing API functions**. From them, users only have to run **5 required functions**, while the remaining **5 functions are optional**.

Users must run the required functions in order. Some optional functions also have to be run only after a certain required function. The optional functions are there for users to visualize their molecules, functional groups, reactions, reaction templates, and non-reactants inside their system.

AutoREACTER accepts an input file path, initializes the session, and exposes methods to step through — or run end-to-end — the reaction-detection and simulation-preparation workflow. Intermediate images such as molecules, functional groups, reactions, reaction templates, and non-reactants are saved automatically to the session image directory.

---

### Import AutoREACTER

```python
import AutoREACTER as arx
```

Users can import AutoREACTER as `arx` so the code is shorter and the long module name is eliminated.

---

## Required Functions in Order

<span style="color:red">**Important**: Users must run these functions in the following order.</span> 

### 1. Run AutoREACTER

```python
arx.run("input.json")
```

This function does the initial processing of AutoREACTER. It clears out the necessary directories, initializes the session, calculates what functional groups are available, and detects what reactions are available based on the functional groups that are available.

This step also initializes the AutoREACTER environment and prepares the working directory for the current run.

After this function, users can optionally visualize molecules, functional groups, and detected reactions.

Optional functions available after `arx.run(...)`:

```python
arx.show_molecules()
arx.show_functional_groups()
arx.show_reactions()
```

---

### 2. Select Reactions

```python
arx.select_reactions()
```

This function lets users select which reactions they want to process with AutoREACTER.

**Users must run this function even if only one reaction is detected.** If only one reaction is available, AutoREACTER will automatically proceed with that reaction after this function is called. If more than one reaction is found, the user is prompted to select which reactions to proceed with.

This function marks the `select_reactions` waterfall stage as complete.

After this function, users can optionally visualize the detected non-reactants.

Optional function available after `arx.select_reactions()`:

```python
arx.show_non_reactants()
```
---

### 3. Select Non-Reactants

```python
arx.select_non_reactants()
```

This function lets users select which non-reactant molecules they want to include in the AutoREACTER workflow.

**Users must run this function even if no non-reactants are detected.** If non-reactants are found, the user is prompted to select which non-reactants they want to proceed with. If no non-reactants are found, this step is skipped automatically after the function is called.

A molecule can be treated as a non-reactant if:

* The molecule does not qualify as a monomer.
* The molecule qualifies as a monomer, but the user chooses not to include it in any selected reaction.



### 4. Prepare Reactions

```python
arx.prepare_reactions()
```

This function prepares reaction templates from the selected reactions for downstream processing.

This function marks the reaction-template preparation stage as complete.

After this function, users can optionally visualize reaction templates.

Optional function available after `arx.prepare_reactions()`:

```python
arx.show_reaction_templates()
```

---

### 5. Process

```python
arx.process()
```

This function executes the back-half of the AutoREACTER pipeline in one shot.

According to the docstring, this method runs reaction template preparation, 3D geometry setup, force-field generation through the LUNAR API, REACTER file building, and LAMMPS simulation writing.

This function should be run only after the required stages are completed:

```python
arx.select_reactions()
arx.select_non_reactants()
arx.prepare_reactions()
```

If the required stages are not completed, AutoREACTER raises a `RuntimeError`.

---

## Full Required Workflow

```python
import AutoREACTER as arx

arx.run("input.json")

arx.select_reactions()

arx.select_non_reactants()

arx.prepare_reactions()

arx.process()
```

---

## Optional Visualization Functions

The following functions are optional. They are provided so users can visualize molecules, functional groups, reactions, non-reactants, and reaction templates.

---

### Show Molecules

Available after:

```python
arx.run("input.json")
```

Usage:

```python
arx.show_molecules()
```

This function returns a PIL/RDKit image grid of the initial molecules or monomers.

Returns:

```text
Image
```

---

### Show Functional Groups

Available after:

```python
arx.run("input.json")
```

Usage:

```python
arx.show_functional_groups()
```

This function returns an image grid with functional groups highlighted on each molecule.

If functional-group detection has not run yet, this function triggers functional-group detection automatically.

Returns:

```text
Image
```

---

### Show Reactions

Available after:

```python
arx.run("input.json")
```

Usage:

```python
arx.show_reactions()
```

This function returns an image grid showing the detected reactions.

If reaction detection has not run yet, this function triggers reaction detection automatically. Since reaction detection depends on functional-group detection, functional-group detection is also triggered if needed.

Returns:

```text
Image
```

---

### Show Non-Reactants

Available after:

```python
arx.select_reactions()
```

Usage:

```python
arx.show_non_reactants()
```

This function returns an image visualizing detected non-reactant species.

If non-reactant detection has not run yet, this function triggers the detection automatically. It also marks the non-reactant detection waterfall stage as complete.

Returns:

```text
Image
```

---

### Show Reaction Templates

Available after:

```python
arx.prepare_reactions()
```

Usage:

```python
arx.show_reaction_templates()
```

This function returns an image grid visualizing the reaction templates.

The `highlight_type` parameter can be used to visualize different parts of the reaction templates. Default type is `"template"`

```python
arx.show_reaction_templates(highlight_type="template")
arx.show_reaction_templates(highlight_type="edge")
arx.show_reaction_templates(highlight_type="delete")
arx.show_reaction_templates(highlight_type="initiators")
```

#### Reaction Template Visualization Options

Visualize reaction templates with different highlighting options by setting the `highlight_type` parameter to one of the following values:

* `template`: Highlights all structural changes in the reaction templates.
* `edge`: Highlights edge atoms of the templates.
* `delete`: Highlights removed components, if applicable.
* `initiators`: Highlights reaction initiator atoms.

The default value is:

```python
highlight_type="template"
```

Returns:

```text
Image
```


## Summary of Public APIs

| API                             | Required or Optional | When to Run                        |
| ------------------------------- | -------------------: | ---------------------------------- |
| `arx.run("input.json")`         |             Required | First                              |
| `arx.select_reactions()`        |             Required | After `arx.run(...)`               |
| `arx.select_non_reactants()`    |             Required | After `arx.select_reactions()`     |
| `arx.prepare_reactions()`       |             Required | After `arx.select_non_reactants()` |
| `arx.process()`                 |             Required | Final required step                |
| `arx.show_molecules()`          |             Optional | After `arx.run(...)`               |
| `arx.show_functional_groups()`  |             Optional | After `arx.run(...)`               |
| `arx.show_reactions()`          |             Optional | After `arx.run(...)`               |
| `arx.show_non_reactants()`      |             Optional | After `arx.select_reactions()`     |
| `arx.show_reaction_templates()` |             Optional | After `arx.prepare_reactions()`    |
