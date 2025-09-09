########################################################################################
#                                                                                      #
#                           LEDAW - LED Analysis Wizard                                #
#                                    Ledip Module                                      #
#                                                                                      #
#                 Automate Local Energy Decomposition Analysis                         #
#                           Using Data in ORCA Outputs                                 #
#                                                                                      #
#                            written by Ahmet Altun                                    #
#                    Max-Planck-Institut für Kohlenforshung                            #
#                Department of Molecular Theory and Spectroscopy                       #
#                                                                                      #
#                                                                                      #
#                                 Citation                                             #
#               If you use any part of this code, in addition to                       #
#               original LED, CPS, and CBS studies, please cite:                       #
#                                                                                      #
#                 1) https://pubs.acs.org/doi/full/10.1021/acs.jcim.5c01561            #
#                    (J. Chem. Inf. Model. 65/17, 2025, 8917–8923)                     #
#                 2) https://github.com/ahmetaltunfatih/LEDAW                          #
#                 3) https://doi.org/10.1002/anie.202421922                            #
#                    (Angew. Chemie. Int. Ed. 64/12, 2025, e202421922)                 #
#                                                                                      #
#                                    License                                           #
#                             Free for academic use.                                   #
#           For commercial use or redistribution, contact the author.                  #
#             The author provides this code as-is, without warranty.                   #
#                                                                                      #
########################################################################################


import sys
sys.path.insert(0,"../..")
from ledip_package import *


xyzfile = "wt.xyz"

cutoff = 1.62

cut_bonds = [
    [0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [10, 11], 
    [12, 13], [14, 15], [16, 17], [18, 19], [20, 21], 
    [22, 23], [24, 25], [26, 27]
]

nonstandard_coordination_numbers = {
#    "Bi": 4,
#    "Fe": 0,
#    "7-9,11,13": 3
}


try:
    fragmentation_engine(xyzfile, cutoff=cutoff, cut_bonds=cut_bonds,
                         nonstandard_coordination_numbers=nonstandard_coordination_numbers)
    print("Fragmentation engine completed successfully!\n")
except Exception as e:
    print(f"Error: {e}")
    print("Fragmentation engine failed to complete.\n")


if sys.platform.startswith('win'):
    os.system("pause")
else:
    input("Press Enter to exit...") # For Linux/macOS
