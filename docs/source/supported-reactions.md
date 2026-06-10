## Supported Reactions

AutoREACTER is currently in **v0.2.2-beta**. At this stage of development, the reaction library is limited to selected step-growth polymerization reactions, including **polycondensation**, **transesterification**, and **polyaddition** reactions.

The core `Detector` module automatically identifies the following functional groups and maps them to their respective reaction pathways.

**Important:** If your `input.json` contains monomers with functional groups outside of this list, AutoREACTER will classify them as *non-reactive molecules* (which you can choose to retain as solvents/additives or discard).

**NOTE:** Certain force fields do not support all atom types; for example, iodine ``(I)`` is sometimes unsupported.

---

### 1. Polyesterification

These reactions form ester linkages (`-COO-`) and typically release water (`H₂O`), alcohols (`R-OH`), or hydrogen halides (e.g., `HCl`) as byproducts.

* **Hydroxy–Carboxylic Acid Polycondensation**

  * *Reactants:* `-OH` + `-COOH`

* **Hydroxy Acid Halide Polycondensation**

  * *Reactants:* `-OH` + `-COX` where `X = Cl, Br, I`

* **Diol + Di-Carboxylic Acid Polycondensation**

  * *Reactants:* Two `-OH` groups + Two `-COOH` groups

* **Diol + Di-Acid Halide Polycondensation**

  * *Reactants:* Two `-OH` groups + Two `-COX` groups where `X = Cl, Br, I`

* **Diol + Di-Carboxylic Ester Transesterification**

  * *Reactants:* Two `-OH` groups + Two ester groups (`-COOR`)

---

### 2. Polyamidation

These reactions form amide linkages (`-CONH-`) and typically release water (`H₂O`) or hydrogen halides (e.g., `HCl`) as byproducts.

* **Amino Acid Polycondensation**

  * *Reactants:* `-NH₂` / `-NH-` + `-COOH`

* **Amino Acid + Amino Acid Polycondensation**

  * *Reactants:* `-NH₂` / `-NH-` + `-COOH`

* **Diamine + Di-Carboxylic Acid Polycondensation**

  * *Reactants:* Two amine groups (`-NH₂` / `-NH-`) + Two `-COOH` groups

* **Diamine + Di-Carboxylic Acid Halide Polycondensation**

  * *Reactants:* Two amine groups (`-NH₂` / `-NH-`) + Two `-COX` groups where `X = Cl, Br, I`

---

### 3. Polyanhydride Formation

These reactions form anhydride linkages (`-CO-O-CO-`) and typically release hydrogen halides (e.g., `HCl`) as byproducts.

* **Carboxylic Acid + Acid Halide Polycondensation**

  * *Reactants:* `-COOH` + `-COX` where `X = Cl, Br, I`

---

### 4. Polythioesterification

These reactions form thioester linkages (`-COS-`) and typically release water (`H₂O`) or hydrogen halides (e.g., `HCl`) as byproducts.

* **Dithiol + Di-Carboxylic Acid Polycondensation**

  * *Reactants:* Two `-SH` groups + Two `-COOH` groups

* **Dithiol + Di-Carboxylic Acid Halide Polycondensation**

  * *Reactants:* Two `-SH` groups + Two `-COX` groups where `X = Cl, Br, I`

---

### 5. Mixed Polyester/Polythioester Formation

These reactions are supported for hydroxy–thiol monomers reacting with acid halides. Depending on the reacting group, either an ester or thioester linkage can be formed.

* **Hydroxy–Thiol + Di-Carboxylic Acid Halide through Hydroxy Group**

  * *Reactants:* `-OH` + `-COX` where `X = Cl, Br, I`

* **Hydroxy–Thiol + Di-Carboxylic Acid Halide through Thiol Group**

  * *Reactants:* `-SH` + `-COX` where `X = Cl, Br, I`

---

### 6. Polyurethane Formation

These reactions form urethane linkages (`-O-CO-NH-`). Unlike most polycondensation reactions, this reaction is a **polyaddition** reaction.

* **Diol + Di-Isocyanate Polyaddition**

  * *Reactants:* Two `-OH` groups + Two isocyanate groups (`-NCO`)

---

NOTE: If you would like support for a specific reaction, please open an issue on
 [AutoREACTER GitHub Repository](https://github.com/NanoCIPHER-Lab/AutoREACTER).
