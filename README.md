#  LEDAW - LED Analysis Wizard   
## A Python-Based Program Package for Automating Local Energy Decomposition Analysis Using ORCA Outputs 
## Author

- Prof. Dr. Ahmet Altun, Max-Planck-Institut für Kohlenforshung, Department of Molecular Theory and Spectroscopy

## Features

- LEDAW automizes all kinds of LED interaction energy analyses such as the interaction of arbitrary number of fragments as in water cluster formation, and interaction of a single or multifragment system with other arbitrary number of single or multifragment systems as in lattice energy computations and duplex DNA formation with multiple fragments on each strand.
- LEDAW calculates N-body, two-body, and cooperativity LED interaction energy matrices for both standard and fragment pairwise (fp)-LED schemes from ORCA output files irrespective of the number of fragments in the adduct (supersystem) and its subsystems in seconds.
- LEDAW performs complete PNO space (CPS) and complete basis set (CBS) extrapolations based on the unextrapolated LED terms in ORCA output files and provides corresponding matrices.
- LEDAW standardizes fragment labels to those in the supersystem file automatically if they are different in supersystem and subsystem ORCA output files.
- LEDAW allows relabeling of fragments if the user is not satisfied with the fragment labeling in supersystem ORCA output file.
- LEDAW allows specifying an alternative file in the case that an ORCA output file does not contain all necessary energy terms. 
- Method-specific (DLPNO-CCSD(T), DLPNO-CCSD, and HFLD) collection of LED terms, such as London dispersion, are done automatically. If an implicit solvation scheme (CPCM, SMD, etc.) is used, the code recognizes this automatically and distribute dielectric contribution to pairwise terms. 
- Standard and fp-LED interaction energy maps are written on separate excel files (each on a seperate sheet).
- LEDAW plots finally heat maps for all the LED interaction energy matrices.

## How to Run
- Download the ledaw_package directory, the example input Python files (water-dimer.py, crystal.py, and boat.py), and the example ORCA output files directory (ORCA-OUT) into your working directory.
- As most LED applications in literature are for two-fragment systems, as an example, a simple script for BSSE-corrected and BSSE-uncorrected LED analyses is provided for water dimer (water-dimer.py).
- The generic crystal.py and boat.py scripts can be easily personalized for your specific use case. In addition to the workflow, they contain detailed explanations of all necessary parameters. If you do not wish to perform some parts of the analysis, the related sections can be commented or removed.  
- In pesonalizing the generic crystal.py and boat.py scripts, you typically just need to specify the path/name of ORCA output files and of where to write LEDAW output files. The parts that you do not have necessary ORCA output files can be commented or deleted.
- To get used to the logic of LEDAW, it is best to start with "crystal.py". It is for performing N-body, two-body, and cooperativity HFLD/LED analysis of the interaction of a central monomer in a crystal with its environment.
- "boat.py" is for running all modules of LEDAW, and thus it is a bit more crowded, due to several file path specifications and multiple calls of engine functions. It performs CPS and CBS extrapolations on the DLPNO-CCSD(T)/LED terms for the interaction energy of boat conformer of water hexamer in addition to N-body, two-body, and cooperativity analyses from standard settings in ORCA output files.

## References
- If you use any part of this code, in addition to original LED, CPS, and CBS studies, please cite: 
    - https://github.com/ahmetaltunfatih/LEDAW (DOI: https://zenodo.org/doi/10.5281/zenodo.13756704)
    - https://chemrxiv.org/engage/chemrxiv/article-details/6698104e01103d79c547414c
