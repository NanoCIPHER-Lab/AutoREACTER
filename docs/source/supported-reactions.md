## Supported Reactions

AutoREACTER is currently in **v{{ autoreacter_version }}**. At this stage of development, the reaction library supports a broad range of step-growth and chain-growth polymerization reactions, including **polycondensation**, **transesterification**, **polyaddition**, **hydrolysis initiation**, and **addition polymerization**.

The core `Detector` module automatically identifies the following functional groups and maps them to their respective reaction pathways.

**Important:** If your `input.json` contains monomers with functional groups outside of this list, AutoREACTER will classify them as *non-reactive molecules* (which you can choose to retain as solvents/additives or discard).

**NOTE:** Certain force fields do not support all atom types; for example, iodine `(I)` is sometimes unsupported.

---

### 1. Polyesterification

These reactions form ester linkages (`-COO-`) and typically release water (`H₂O`), alcohols (`R-OH`), or hydrogen halides (e.g., `HCl`) as byproducts.

* **Hydroxy Carboxylic Acid Polycondensation**
* *Reactants:* `-OH` + `-COOH`


* **Hydroxy Carboxylic Acid and Hydroxy Carboxylic Acid Polycondensation**
* *Reactants:* `-OH` + `-COOH` (Intermolecular)


* **Hydroxy Acid Halides Polycondensation**
* *Reactants:* `-OH` + `-COX` where `X = Cl, Br, I`


* **Hydroxy Acid Halides Hydroxy Acid Halides Polycondensation**
* *Reactants:* `-OH` + `-COX` where `X = Cl, Br, I` (Intermolecular)


* **Diol and Di-Carboxylic Acid Polycondensation**
* *Reactants:* Two `-OH` groups + Two `-COOH` groups


* **Diol and Di-Acid Halide Polycondensation**
* *Reactants:* Two `-OH` groups + Two `-COX` groups where `X = Cl, Br, I`


* **Diol and Di-Carboxylic Ester Polycondensation (Transesterification)**
* *Reactants:* Two `-OH` groups + Two ester groups (`-COOR`)



---

### 2. Polyamidation

These reactions form amide linkages (`-CONH-`) and typically release water (`H₂O`) or hydrogen halides (e.g., `HCl`) as byproducts.

* **Amino Acid Polycondensation**
* *Reactants:* `-NH₂` / `-NH-` + `-COOH`


* **Amino Acid and Amino Acid Polycondensation**
* *Reactants:* `-NH₂` / `-NH-` + `-COOH` (Intermolecular)


* **Di-Amine and Di-Carboxylic Acid Polycondensation**
* *Reactants:* Two amine groups (`-NH₂` / `-NH-`) + Two `-COOH` groups


* **Di-Amine and Di-Carboxylic Acid Halide Polycondensation**
* *Reactants:* Two amine groups (`-NH₂` / `-NH-`) + Two `-COX` groups where `X = Cl, Br, I`


* **Hydrolytic Initiation of Caprolactam**
* *Reactants:* Water (`H₂O`) + Lactam ring opening



---

### 3. Polyanhydride Formation

These reactions form anhydride linkages (`-CO-O-CO-`) and typically release hydrogen halides (e.g., `HCl`) as byproducts.

* **Carboxylic Acid and Acid Halide Polycondensation**
* *Reactants:* `-COOH` + `-COX` where `X = Cl, Br, I`


* **Carboxylic Acid and Acid Halide Copolycondensation**
* *Reactants:* Mixed `-COOH` + `-COX` copolymerization systems



---

### 4. Polythioesterification

These reactions form thioester linkages (`-COS-`) and typically release water (`H₂O`) or hydrogen halides (e.g., `HCl`) as byproducts.

* **Dithiol and Di-Carboxylic Acid Halide Polycondensation**
* *Reactants:* Two `-SH` groups + Two `-COX` groups where `X = Cl, Br, I`


* **Dithiol and Di-Carboxylic Acid Polycondensation**
* *Reactants:* Two `-SH` groups + Two `-COOH` groups



---

### 5. Mixed Polyester/Polythioester Formation

These reactions are supported for hydroxy–thiol monomers reacting with acid halides. Depending on the reacting group, either an ester or thioester linkage can be formed.

* **Hydroxy-Thiol and Di-Carboxylic Acid Halide Polycondensation through Hydroxy Group**
* *Reactants:* `-OH` + `-COX` where `X = Cl, Br, I`


* **Hydroxy-Thiol and Di-Carboxylic Acid Halide Polycondensation through Thiol Group**
* *Reactants:* `-SH` + `-COX` where `X = Cl, Br, I`



---

### 6. Polyurethane, Polythiourethane, and Polyurea Formation

These reactions form urethane, thiourethane, or urea linkages via **polyaddition** pathways.

* **Diol and Di-Isocyanate Polyaddition (Polyurethane Formation)**
* *Reactants:* Two `-OH` groups + Two isocyanate groups (`-NCO`)


* **Dithiol and Di-Isocyanate Polyaddition (Polythiourethane Formation)**
* *Reactants:* Two `-SH` groups + Two isocyanate groups (`-NCO`)


* **Di-Amine and Di-Isocyanate Polyaddition (Polyurea Formation)**
* *Reactants:* Two amine groups (`-NH₂` / `-NH-`) + Two isocyanate groups (`-NCO`)



---

### 7. Epoxy-Amine Addition and Crosslinking

These reactions model step-growth/network formation between amine curing agents and epoxy rings.

* **Primary Amine and Epoxide Polyaddition (First Addition)**
* *Reactants:* Primary amine (`-NH₂`) + Epoxide ring


* **Secondary Amine and Epoxide Polyaddition (Second Addition / Crosslink)**
* *Reactants:* Secondary amine (`-NH-`) + Epoxide ring



---

### 8. Vinyl and Fluoropolymer Addition Polymerization

These chain-growth pathways model radical initiation, propagation, and copolymerization of vinyl and fluorinated monomers.

* **Vinyl Addition Polymerization Initiation**
* *Reactants:* Vinyl double bonds (`-CH=C-`)


* **Vinyl Addition Polymerization Propagation**
* *Reactants:* Vinyl monomer + Chain-end radical


* **Vinyl Copolymerization**
* *Reactants:* Mixed vinyl monomer systems


* **Tetrafluoroethylene Addition Polymerization Initiation**
* *Reactants:* Tetrafluoroethylene (`TFE`) self-initiation


* **Tetrafluoroethylene Addition Polymerization Propagation**
* *Reactants:* Tetrafluoroethylene monomer + TFE radical chain-end



---

### 9. Polycarbonate Formation

These reactions build carbonate linkages via condensation or transcarbonation.

* **Diol and Phosgene Polycondensation (Polycarbonate Formation)**
* *Reactants:* Two `-OH` groups + Phosgene (`COCl₂`)


* **Diol and Diphenyl Carbonate Polycondensation (Transcarbonation)**
* *Reactants:* Two `-OH` groups + Diphenyl carbonate



---

### 10. Polysiloxane Formation

These pathways handle hydrolysis of chlorosilanes and condensation of silanols into silicone chains.

* **Dichlorosilane Hydrolysis to Silanol**
* *Reactants:* Dichlorosilane (`-Si-Cl`) + Water (`H₂O`)


* **Silanediol Polycondensation (Polysiloxane Formation)**
* *Reactants:* Silanediols (`-Si-OH`)


* **Silanediol and Silanediol Copolycondensation (Polysiloxane Formation)**
* *Reactants:* Mixed silanediol systems



---

### 11. Thiol-Ene Click Polymerization

* **Dithiol and Diene Thiol-Ene Click Polymerization**
* *Reactants:* Dithiol (`-SH`) + Diene (`-C=C-`)



---

NOTE: If you would like support for a specific reaction, please open an issue on [AutoREACTER GitHub Repository](https://github.com/NanoCIPHER-Lab/AutoREACTER).