# **LEDAW-GUI: LED Analysis Wizard with GUI**  

## **A Python-Based Program Package with GUI for Automating Local Energy Decomposition Analysis Using ORCA Outputs**  

### **Author**  
**Prof. Dr. Ahmet Altun**  
Max-Planck-Institut für Kohlenforschung  
Department of Molecular Theory and Spectroscopy  

---

## **Features**

- **LEDAW** automates all types of LED interaction energy analyses, including:
  - Interactions between arbitrary numbers of fragments (e.g., water cluster formation).
  - Interactions of single- or multi-fragment systems with other single- or multi-fragment systems (e.g., lattice energy calculations, duplex DNA formation with multiple fragments per strand).
- Features a **user-friendly, self-explanatory GUI** with built-in info buttons and help messages, providing guidance at every step.
- Provides **example Python input files** for code-oriented users who prefer script-based workflows.
- Calculates **N-body, two-body, and cooperativity LED interaction energy matrices** for both standard and fragment pairwise (fp)-LED schemes from ORCA output files — independent of the number of fragments in the supersystem and its subsystems, within seconds.
- Performs **Complete PNO Space (CPS)** and **Complete Basis Set (CBS) extrapolations** based on unextrapolated LED terms from ORCA outputs and generates corresponding energy matrices.
- **Automatically standardizes fragment labels** to match those in the supersystem file if they differ between supersystem and subsystem ORCA output files.
- Allows **relabeling of fragments** if the user wishes to adjust the fragment labeling in the supersystem ORCA output.
- Supports specifying an **alternative file** if the primary ORCA output file lacks certain required energy terms.
- Collects LED terms **method-specifically** (for DLPNO-CCSD(T), DLPNO-CCSD, and HFLD), including terms like London dispersion.
- Detects the use of **implicit solvation schemes** (CPCM, SMD, etc.) and distributes dielectric contributions across pairwise terms.
- Detects automatically if **BSSE-correction** is requested and proceeds subsystem files accordingly.
- Writes **standard and fp-LED interaction energy matrices** into separate Excel files, with each matrix on a separate sheet.
- Provides **heatmaps** of all interaction energy matrices for convenient data interpretation and presentation.

---

## **How to Run**

### **Downloading**
- Download the `ledaw_package` directory along with:
  - `main.py`
  - `LEDAW.spec`
  - Example input Python files:  
    - `water-dimer.py`  
    - `crystal.py`  
    - `boat.py`  
  - Example ORCA output files directory: `ORCA-OUT`
- Place these into your working directory.

---

### **Installing LEDAW-GUI (Executable Generation)**

- **LEDAW-GUI** is pre-configured to generate an executable using **PyInstaller**.
- Python must be installed on your system.
- To create the executable:

```bash
cd /path/to/the/downloaded/directory
pyinstaller LEDAW.spec
```

- Then, run the following command:

``bash
pyinstaller LEDAW.spec
```

- The generated executable will work without requiring the downloaded directory or Python installation.

### **Running LEDAW-GUI Directly with Python**

- As an alternative to generating the executable, you can run LEDAW-GUI directly with Python:

```bash
cd /path/to/the/downloaded/directory
python main.py
```

### **Running LEDAW Without GUI (Script-Based Workflow)**

For more code-oriented users, three example Python input scripts are provided:

- **`water-dimer.py`**:  
  Demonstrates BSSE-corrected and BSSE-uncorrected LED analyses for a water dimer (a two-fragment system).

- **`crystal.py`**:  
  Designed for performing **N-body**, **two-body**, and **cooperativity** HFLD/LED analysis of the interaction between a central monomer and its environment in a crystal.

- **`boat.py`**:  
  A comprehensive example that runs **all modules** of LEDAW, including CPS and CBS extrapolations on the DLPNO-CCSD(T)/LED terms for the boat conformer of a water hexamer.

---

#### **Personalizing the Example Scripts**

- To adapt the provided scripts for your own system:
  - Specify the **path and filenames** of your ORCA output files.
  - Specify the **output path** where LEDAW should write the results.
  - If certain parts of the analysis (e.g., CBS extrapolation, CPS extrapolation, cooperativity analysis) are not needed, simply **comment out** or **remove** the relevant sections.

---

#### **Recommended Workflow**

- It is recommended to first explore **`crystal.py`** to understand the basic LEDAW logic and workflow.
- Use **`boat.py`** if you want to explore **all available modules**, including CPS and CBS extrapolations.  
  This script includes multiple file path specifications and several calls to the engine functions, demonstrating the full flexibility of LEDAW.

> **Tip:**  
> If your ORCA outputs do not contain certain energy terms, the corresponding parts of the example scripts can be commented out or removed.

---

## **References**

If you use any part of this code or its results in your research, in addition to the original LED, CPS, and CBS studies, please cite:

- **LEDAW GitHub Repository:**  
  [https://github.com/ahmetaltunfatih/LEDAW](https://github.com/ahmetaltunfatih/LEDAW)  
  DOI: [https://zenodo.org/doi/10.5281/zenodo.13756704](https://zenodo.org/doi/10.5281/zenodo.13756704)

- **Preprint on ChemRxiv:**  
  [https://chemrxiv.org/engage/chemrxiv/article-details/6698104e01103d79c547414c](https://chemrxiv.org/engage/chemrxiv/article-details/6698104e01103d79c547414c)

---

## **License**

### **Academic Use License (with Commercial Restriction)**

Copyright (c) 2025  
**Ahmet Altun**

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction for **academic and non-commercial use**, including without limitation the rights to use, copy, modify, merge, publish, and distribute copies of the Software, subject to the following conditions:

> **Commercial use, including redistribution or incorporation into commercial products, is not permitted without prior written permission from the author.  
> Please contact the author for commercial licensing.**

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

> **Disclaimer:**  
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHOR OR COPYRIGHT HOLDER BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

---
