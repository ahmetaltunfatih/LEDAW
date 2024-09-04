#  LEDAW - LED Analysis Wizard   
## A Python-Based Program Package for Automating Local Energy Decomposition Analysis Using Data in ORCA Outputs 
## Author

- Prof. Dr. Ahmet Altun, Max-Planck-Institut für Kohlenforshung, Department of Molecular Theory and Spectroscopy, https://github.com/ahmetaltunfatih/


## Features

- LEDAW automizes all kinds of LED interaction energy analyses such as the interaction of arbitrary number of fragments as in water cluster formation, and interaction of a single or multifragment system with other arbitrary number of single or multifragment systems as in lattice energy computations and duplex DNA formation with multiple fragments on each strand.
- LEDAW calculates N-body, two-body, and cooperativity LED interaction energy matrices for both standard and fragment pairwise (fp)-LED schemes from ORCA output files irrespective of the number of fragments in the adduct (supersystem) and its subsystems in seconds.
- LEDAW performs complete PNO space (CPS) and complete basis set (CBS) extrapolations based on the unextrapolated LED terms in ORCA output files and provides corresponding matrices.
- LEDAW standardizes fragment labels to those in the supersystem file automatically if they are different in supersystem and subsystem ORCA output files.
- LEDAW allows relabeling of fragments if you are not satisfied with the fragment labeling in supersystem ORCA output file.
- LEDAW allows specifying an alternative file in the case that an ORCA output file does not contain all necessary energy terms. 
- Method-specific (DLPNO-CCSD(T), DLPNO-CCSD, and HFLD) collection of LED terms, such as London dispersion, are done automatically.
- Standard and fp-LED interaction energy maps are written for separate excel files.
- LEDAW plots finally heat maps for all the LED interaction energy matrices.

## How to Run
- Download ledaw_package directory, example input python files crystal.py and boat.py, and example ORCA output files directory ORCA-OUT into your working directory.
- Run crystal.py to perform N-body, two-body, and cooperativity HFLD/LED analysis of the interaction of a central monomer of a crsytal with its environment.
- Run boat.py to perform additionaly CPS and CBS extrapolations
- The generic crystal.py and boat.py can be easily personalzed for the user's use case. In addition to the workflow, they contain detailed explanations of all necessary parameters. If you do not want to perform some parts of the analyses, the related parts can be commented or removed. 
## Reference
- If you use any part of this code, in addition to            original LED, CPS, and CBS studies, please cite: 
    - https://github.com/ahmetaltunfatih/LEDAW
    - https://chemrxiv.org/engage/chemrxiv/article-details/6698104e01103d79c547414c
