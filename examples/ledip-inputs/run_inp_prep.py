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


xyzfile = "wt_labeled.xyz"

output_dir = "LED-INPS/DNA"

subsystems = [[1,9],[2,3,4,5,6,7,8],[10,11,12,13,14,15,16]]

orca_inp_heading = """
! DLPNO-CCSD(T) RIJCOSX def2-tzvp(-f) def2/J def2-tzvp/C VeryTightSCF CPCM(water) DefGrid3 NormalPNO LED

%pal nprocs 16 end
%maxcore 6000
%scf maxiter 999 end

%mdci TCutPairs 1e-5 printlevel 3 end

*xyz 0 1
"""


try:
    led_input_prep_engine(xyzfile, subsystems, orca_inp_heading, output_dir)
    print("LED input preparation engine completed successfully!\n")
except Exception as e:
    print(f"Error: {e}")
    print("LED input preparation engine failed to complete.")


if sys.platform.startswith('win'):
    os.system("pause")
else:
    input("Press Enter to exit...") # For Linux/macOS
