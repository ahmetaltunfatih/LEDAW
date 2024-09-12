########################################################################################
#                                                                                      #
#                           LEDAW - LED Analysis Wizard                                #
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
#                 1) https://github.com/ahmetaltunfatih/LEDAW                          #
#   2) https://chemrxiv.org/engage/chemrxiv/article-details/6698104e01103d79c547414c   #
#                                                                                      #
########################################################################################


from ledaw_package import *

### Choose the method
# Note 1: Place the proper method your ORCA outputs belong to. Choices: DLPNO-CCSD(T), DLPNO-CCSD, HFLD
# Note 2: Method name is not case sensitive.
# Note 3: The code does not check for the method name in the output files. It is necessary for properly
# combining decomposed LED terms.

method = "DLPNO-CCSD(T)"

# In HFLD calculations, for one-body subsystems, RHF energies rather than reference E(0) energies can be 
# used in conjunction with LED in calculating electronic preparation energies with negligibly small variations.
# If you want to proceed this way, set use_ref_as_rhf_in_hfld = True (default: False or None). Then, if LEDAW
# could not find E(0) energies, it will use RHF energies of one-body subsystems in HFLD/LED calculations.
# CAUTION: This approach cannot use used in open-shell calculations as E(0) corresponds to QRO energy.

use_ref_as_rhf_in_hfld = False

### Convert Hartree to kcal/mol
conversion_factor = 627.5095  # for kj/mol, use: 2625.5

# If you are not happy with the fragment labeling in ORCA supersystem (adduct) LED file
# you can reorder the labels with relabel_mapping variable (default: None, equivalent to [1,2,3,..,the_number_of_fragments])
# As an example, relabel_mapping=[5,3,2,1,4] swaps label 1 with label 5; 2 with 3; 3 with 2; 4 with 1; and 5 with 4.
# Note: This is different than inconsistent labelling the user introduced in supersystem and subsystem
# ORCA jobs. The code automatically standardize subsystem fragment labelings to those in supersystem file.
# Since the present example has two fragments, this is irrelevant.  

relabel_mapping = None


###############################################################################################################################################################
### SINCE THE SUPERSYSTEM (ADDUCT) IS COMPOSED OF TWO FRAGMENTS, ANY OF THE engine_LED_N_body AND engine_LED_two_body FUNCTIONS CAN BE USED IN THIS EXAMPLE:  #
# To demonstrate how to use both functions, in the following, as an example, engine_LED_N_body function is called for BSSE-corrected interaction energy       #
# and engine_LED_two_body function is called for BSSE-uncorrected interaction energy.                                                                         #
############################################################################################################################################################### 


nbody_title='''
##############################################################################################################
#                         BSSE-CORRECTED INTERACTION ENERGY USING N-BODY ENGINE                              #             
##############################################################################################################
'''

### Specify ORCA Output FilesNeeded For N-Body LED Analysis
# Below specify main and alternative filenames with their paths as lists. 
# Note 1: Always specify the output file for the entire system (super system) first in main_filenames and
# alternative_filenames lists. Where necessary, the code differentiate the super system and its subsystems 
# from these lists.
# Note 2: All subsystem output files (even if they do not have LED data as they are one-fragment subsystem)
# must be listed. 
# Note 3: If a part of needed LED data is not present in your main file but present in an another file 
# (such as due to "a crash", or "DoLEDHF False", or absence of "triples"), an alternative file with
# needed data can be specified. Hence, for example, if your method choice is DLPNO-CCSD(T), 
# DLPNO-CCSD output can  still be provided in the main or alternative file name.
# Note 4: The corresponding main and alternative ORCA output files must be given at the same index in the file lists.
# One letter or one number file names may sometimes cause problems. Hence, awoid using such short names, e.g., "A.mpi4.out".
# If you are on windows and using back slash, it is safer to specify all file paths as raw string, e.g., r'.\ORCA-OUT\ADDUCT.mpi4.out'
# Note 5: You do not have to specify any file name for alternative_filenames. But an empty list with the length of main_filenames must be initiated. 
# Note 6: If you specified the same fragment with different labels in supersystem and subsystem ORCA output files, 
# the code automatically labels subsystem fragment labels as in the supersystem output file.
# Note 7: The code only accepts seperate ORCA output files for each system (supersystem and subsystems) and computational setting 
# If you have compound job output, you need to split this file to seperate files for each job.  
# Note 8: Optionally, you can perform Complete PNO Space (CPS) and Complete Basis Set (CBS) extrapolations.
# In this case, you need to specify the path of each ORCA output file with each computational setting.
# For CPS extrapolation, you need two calculations with a looser (LPNO) and a tighter (TPNO) TCutPNO values, 
# such as 1e-6 and 1e-7, respectively.
# For CBS extrapolation, you need two calculations with a smaller (SB) and a larger (LB) basis sets, 
# such as aug-cc-pVTZ and aug-cc-pVQZ, respectively.
# In the list variables, for example, SB_LPNO corresponds to the calculations with smaller basis set and looser TCutPNO setting.
# The code does not check basis set and TCutPNO values from output files but relies on the user input in the lists below.

print(nbody_title)

main_filenames = [r'./ORCA-OUT/WATER-DIMER/DIMER/DIMER.out', 
                          r'./ORCA-OUT/WATER-DIMER/MONO-CP/BSSE1.out',
                          r'./ORCA-OUT/WATER-DIMER/MONO-CP/BSSE2.out',
                 ]

alternative_filenames = ['', '', '']


### Specify the Directories where LEDAW will write N-Body LED Matrices, which will be
# read from ORCA output files
LEDAW_output_path = r'./LEDAW-OUT/WATER-DIMER/withBSSE'

### Run N-Body LED engine.
# Standard and fp-LED N-body matrices will be written to excel files in specified LEDAW output directory
engine_LED_N_body(main_filenames=main_filenames, 
                  alternative_filenames=alternative_filenames, 
                  conversion_factor=conversion_factor, 
                  method=method,
                  relabel_mapping=relabel_mapping,
                  use_ref_as_rhf_in_hfld=use_ref_as_rhf_in_hfld,
                  LEDAW_output_path=LEDAW_output_path)


twobody_title='''
#####################################################################################################
#                      BSSE-UNCORRECTED INTERACTION ENERGY USING TWO-BODY ENGINE                    #
#####################################################################################################
'''
print(twobody_title)

### Speciy the one-body ORCA output files
# The order of fragments in the list matters. In the follwing, we list consistent with the fragment labeling sequence in the adduct file.
# However, there  are other ways. To get a feeling, see other examples.

one_body_orcaout_filenames = [r'./ORCA-OUT/WATER-DIMER/MONO/MONO1.out',
                                r'./ORCA-OUT/WATER-DIMER/MONO/MONO2.out']


### Specify the two-body ORCA output file directory
# Note: This directory must contain only the necessary two-body ORCA output file.
two_body_orcaout_directory = r'./ORCA-OUT/WATER-DIMER/DIMER'


### Specify the Directory where LEDAW will write Two-Body LED Matrices
LEDAW_output_path_twobody = r'./LEDAW-OUT/WATER-DIMER/withoutBSSE'


### Run two-body LED engine
# Standard and fp-LED two-body matrices will be written to excel files in specified LEDAW two-body output directory
engine_LED_two_body(one_body_orcaout_filenames=one_body_orcaout_filenames,
                    two_body_orcaout_directory=two_body_orcaout_directory,
                    conversion_factor=conversion_factor, 
                    method=method,
                    LEDAW_output_path_two_body=LEDAW_output_path_twobody,
                    use_ref_as_rhf_in_hfld=use_ref_as_rhf_in_hfld,
                    relabel_mapping=relabel_mapping)