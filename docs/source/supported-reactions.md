![AutoREACTER Logo](_static/logo.png)
## Supported Reactions

AutoREACTER is currently in **v0.2-beta**. At this stage of development, the reaction library is limited, specifically step-growth **polycondensation** reactions. 

The core `Detector` module automatically identifies the following functional groups and maps them to their respective reaction pathways. 

**Important:** If your `input.json` contains monomers with functional groups outside of this list, AutoREACTER will classify them as *non-reactive molecules* (which you can choose to retain as solvent/additives or discard).

---

### 1. Polyesterification
These reactions form ester linkages (`-COO-`) and typically release water (`H₂O`) or hydrogen halides (e.g., `HCl`) as byproducts. AutoREACTER tracks and handles these byproducts natively.

* **Hydroxy–Carboxylic Acid Polycondensation**
  * *Reactants:* `-OH` + `-COOH`
* **Hydroxy Acid Halide Polycondensation**
  * *Reactants:* `-OH` + `-COX` (where X is typically Cl)
* **Diol + Di-Acid Halide Polycondensation**
  * *Reactants:* Two `-OH` groups + Two `-COX` groups
* **Diol + Di-Carboxylic Acid Polycondensation**
  * *Reactants:* Two `-OH` groups + Two `-COOH` groups

---

### 2. Polyamidation
These reactions form amide linkages (`-CONH-`), which are critical for synthesizing nylons, Kevlar, and structural proteins. 

* **Amino Acid Polycondensation**
  * *Reactants:* `-NH₂` + `-COOH`
* **Diamine + Di-Carboxylic Acid Polycondensation**
  * *Reactants:* Two `-NH₂` groups + Two `-COOH` groups
* **Diamine + Di-Carboxylic Acid Halide Polycondensation**
  * *Reactants:* Two `-NH₂` groups + Two `-COX` groups

---

## Future Expansions

The NanoCIPHER team is actively working to expand this library.  
Future versions of AutoREACTER will include support for additional reactions.

If you would like support for a specific reaction, please open an issue on  
[AutoREACTER GitHub Repository](https://github.com/NanoCIPHER-Lab/AutoREACTER).
