#  LEDAW - LED Analysis Wizard   
## A Python-Based Program Package for Automating Local Energy Decomposition Analysis Using Data in ORCA Outputs 
## Author

- Prof. Dr. Ahmet Altun, Max-Planck-Institut für Kohlenforshung, Department of Molecular Theory and Spectroscopy

## Features

- LEDAW automizes all kinds of LED interaction energy analyses such as the interaction of arbitrary number of fragments as in water cluster formation, and interaction of a single or multifragment system with other arbitrary number of single or multifragment systems as in lattice energy computations and duplex DNA formation with multiple fragments on each strand.
- LEDAW calculates N-body, two-body, and cooperativity LED interaction energy matrices for both standard and fragment pairwise (fp)-LED schemes from ORCA output files irrespective of the number of fragments in the adduct (supersystem) and its subsystems in seconds.
- LEDAW performs complete PNO space (CPS) and complete basis set (CBS) extrapolations based on the unextrapolated LED terms in ORCA output files and provides corresponding matrices.
- LEDAW standardizes fragment labels to those in the supersystem file automatically if they are different in supersystem and subsystem ORCA output files.
- LEDAW allows relabeling of fragments if the user is not satisfied with the fragment labeling in supersystem ORCA output file.
- LEDAW allows specifying an alternative file in the case that an ORCA output file does not contain all necessary energy terms. 
- Method-specific (DLPNO-CCSD(T), DLPNO-CCSD, and HFLD) collection of LED terms, such as London dispersion, are done automatically.
- Standard and fp-LED interaction energy maps are written on separate excel files (each on a seperate sheet).
- LEDAW plots finally heat maps for all the LED interaction energy matrices.

## How to Run
- Download the ledaw_package directory, the example input Python files (crystal.py and boat.py), and the example ORCA output files directory (ORCA-OUT) into your working directory.
- The generic crystal.py and boat.py scripts can be easily personalized for your specific use case. In addition to the workflow, they contain detailed explanations of all necessary parameters. If you do not wish to perform some parts of the analysis, the related sections can be commented or removed.
- In pesonalizing the generic crystal.py and boat.py scripts, you typically just need to specify the path/name of ORCA output files and of where to write LEDAW output files. The parts that you do not have necessary ORCA output files can be commented or deleted.
- To get used to the logic of LEDAW, it is best to start with "crystal.py". It is for performing N-body, two-body, and cooperativity HFLD/LED analysis of the interaction of a central monomer in a crystal with its environment.
- "boat.py" is for running all modules of LEDAW, and thus it is a bit more crowded, due to several file path specifications and multiple calls of engine functions. It performs CPS and CBS extrapolations on the DLPNO-CCSD(T)/LED terms for the interaction energy of boat conformer of water hexamer in addition to N-body, two-body, and cooperativity analyses from standard settings in ORCA output files.

## MIT License

Copyright (c) 2024 Ahmet Altun

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## References
- If you use any part of this code, in addition to original LED, CPS, and CBS studies, please cite: 
    - https://github.com/ahmetaltunfatih/LEDAW (DOI: 10.13140/RG.2.2.35160.30729)
    - https://chemrxiv.org/engage/chemrxiv/article-details/6698104e01103d79c547414c
