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
#    1) https://github.com/ahmetaltunfatih/LEDAW (DOI: 10.13140/RG.2.2.35160.30729)    #
#   2) https://chemrxiv.org/engage/chemrxiv/article-details/6698104e01103d79c547414c   # 
#                                                                                      #
########################################################################################

### The following is an automated LEDAW procedure for obtaining LED interactoin energy matrices and heat maps 
### on the interaction energy of the boat conformer of water hexamer from ORCA output files.
### For more details, see https://chemrxiv.org/engage/chemrxiv/article-details/6698104e01103d79c547414c
### N-body, two-body, and cpooperativty LED maps are generated within standard and fp-LED schemes for each 
### computational setting. At the same time, LED terms are extrapolated to CPS and CBS limits. 
### Finally, standard and fp-LED LED heat maps are generated for all cases.
 
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

use_ref_as_rhf_in_hfld = None

### Convert Hartree to kcal/mol
conversion_factor = 627.5095  # for kj/mol, use: 2625.5

# If you are not happy with the fragment labeling in ORCA supersystem (adduct) LED file
# you can reorder the labels with relabel_mapping variable (default: None, equivalent to [1,2,3,..,the_number_of_fragments])
# As an example, relabel_mapping=[5,3,2,1,4] swaps label 1 with label 5; 2 with 3; 3 with 2; 4 with 1; and 5 with 4.
# Note: This is different than inconsistent labelling the user introduced in supersystem and subsystem
# ORCA jobs. The code automatically standardize subsystem fragment labelings to those in supersystem file.  

relabel_mapping = [1,3,2,4,6,5]

### Extrapolation coefficient F for E_extr = E_X + F*(E_Y - E_X)
# where E_X corresponds to the energy with smaller basis set (in CBS) or looser TCutPNO (in CPS)and
# E_Y corresponds to the energy with larger basis set (in CBS) or tighter TCutPNO (in CPS).
# For more details, see the notes in: https://github.com/ORCAQuantumChemistry/CompoundScripts/blob/main/EnergyExtrapolation/extrapolate_CPS_CBS.cmp 

# F for CBS(3/4) using (aug)-cc-pV(T/Q)Z: 
F_ref_cbs = 1.30130392358216 
F_corr_cbs = 1.71188948355961

# F for CPS
# Reference part is not affected from the CPS extrapolation (E_Y - E_X = 0). As for CPS and CBS the same function is used, a redundant F for reference needs to be specified. 
# Below the same F will be specified in CPS. Anyway, the results do not depend on the choice.
F = 1.5  

nbody_title='''
##############################################################################################################
#                                               N-BODY LED                                                   #             
##############################################################################################################
'''
print(nbody_title)


### Specify ORCA Output Files Needed For N-Body LED Analysis
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
# If you are on windows and using back slash, it is safer to specify all file paths as raw string, e.g., r'.\ORCA-OUT\ADDUCT.out'
# Note 5: You do not have to specify any file name for alternative_filenames. But an empty list with the length of main_filenames must be initiated. 
# Note 6: If you specified the same fragment with different labels in supersystem and subsystem ORCA output files, 
# the code automatically labels subsystem fragment labels as in the supersystem output file.
# Note 7: The code only accepts separate ORCA output files for each system (supersystem and subsystems) and computational setting 
# If you have compound job output, you need to split this file to separate files for each job.  
# Note 8: Optionally, you can perform Complete PNO Space (CPS) and Complete Basis Set (CBS) extrapolations.
# In this case, you need to specify the path of each ORCA output file with each computational setting.
# For CPS extrapolation, you need two calculations with a looser (LPNO) and a tighter (TPNO) TCutPNO values, 
# such as 1e-6 and 1e-7, respectively.
# For CBS extrapolation, you need two calculations with a smaller (SB) and a larger (LB) basis sets, 
# such as aug-cc-pVTZ and aug-cc-pVQZ, respectively.
# In the list variables, for example, SB_LPNO corresponds to the calculations with smaller basis set and looser TCutPNO setting.
# The code does not check basis set and TCutPNO values from output files but relies on the user input in the lists below.

main_filenames_SB_LPNO = [r'./ORCA-OUT/BOAT/aTZ/PNO6/HEXAMER/ADDUCT.mpi16.out', 
                          r'./ORCA-OUT/BOAT/aTZ/PNO6/ONEBODY/frag1.mpi4.out',
                          r'./ORCA-OUT/BOAT/aTZ/PNO6/ONEBODY/frag2.mpi4.out',
                          r'./ORCA-OUT/BOAT/aTZ/PNO6/ONEBODY/frag3.mpi4.out',
                          r'./ORCA-OUT/BOAT/aTZ/PNO6/ONEBODY/frag4.mpi4.out',
                          r'./ORCA-OUT/BOAT/aTZ/PNO6/ONEBODY/frag5.mpi4.out',
                          r'./ORCA-OUT/BOAT/aTZ/PNO6/ONEBODY/frag6.mpi4.out',
                 ]

alternative_filenames_SB_LPNO = ['', '', '', '', '', '' ,'']


main_filenames_SB_TPNO = [r'./ORCA-OUT/BOAT/aTZ/PNO7/HEXAMER/ADDUCT.mpi16.out', 
                          r'./ORCA-OUT/BOAT/aTZ/PNO7/ONEBODY/frag1.mpi4.out',
                          r'./ORCA-OUT/BOAT/aTZ/PNO7/ONEBODY/frag2.mpi4.out',
                          r'./ORCA-OUT/BOAT/aTZ/PNO7/ONEBODY/frag3.mpi4.out',
                          r'./ORCA-OUT/BOAT/aTZ/PNO7/ONEBODY/frag4.mpi4.out',
                          r'./ORCA-OUT/BOAT/aTZ/PNO7/ONEBODY/frag5.mpi4.out',
                          r'./ORCA-OUT/BOAT/aTZ/PNO7/ONEBODY/frag6.mpi4.out',
                 ]

alternative_filenames_SB_TPNO = ['', '', '', '', '', '' ,'']

main_filenames_LB_LPNO = [r'./ORCA-OUT/BOAT/aQZ/PNO6/HEXAMER/ADDUCT.mpi16.out', 
                          r'./ORCA-OUT/BOAT/aQZ/PNO6/ONEBODY/frag1.mpi4.out',
                          r'./ORCA-OUT/BOAT/aQZ/PNO6/ONEBODY/frag2.mpi4.out',
                          r'./ORCA-OUT/BOAT/aQZ/PNO6/ONEBODY/frag3.mpi4.out',
                          r'./ORCA-OUT/BOAT/aQZ/PNO6/ONEBODY/frag4.mpi4.out',
                          r'./ORCA-OUT/BOAT/aQZ/PNO6/ONEBODY/frag5.mpi4.out',
                          r'./ORCA-OUT/BOAT/aQZ/PNO6/ONEBODY/frag6.mpi4.out',
                 ]

alternative_filenames_LB_LPNO = ['', '', '', '', '', '', '']


main_filenames_LB_TPNO = [r'./ORCA-OUT/BOAT/aQZ/PNO7/HEXAMER/ADDUCT.mpi16.out', 
                          r'./ORCA-OUT/BOAT/aQZ/PNO7/ONEBODY/frag1.mpi4.out',
                          r'./ORCA-OUT/BOAT/aQZ/PNO7/ONEBODY/frag2.mpi4.out',
                          r'./ORCA-OUT/BOAT/aQZ/PNO7/ONEBODY/frag3.mpi4.out',
                          r'./ORCA-OUT/BOAT/aQZ/PNO7/ONEBODY/frag4.mpi4.out',
                          r'./ORCA-OUT/BOAT/aQZ/PNO7/ONEBODY/frag5.mpi4.out',
                          r'./ORCA-OUT/BOAT/aQZ/PNO7/ONEBODY/frag6.mpi4.out',
                 ]

alternative_filenames_LB_TPNO = ['', '', '', '', '', '', '']


### Specify the Directories where LEDAW will write N-Body LED Matrices, which will be
# read from ORCA output files
LEDAW_output_path_SB_LPNO = r'./LEDAW-OUT/BOAT/SB-LPNO/NBODY'
LEDAW_output_path_SB_TPNO = r'./LEDAW-OUT/BOAT/SB-TPNO/NBODY'
LEDAW_output_path_LB_LPNO = r'./LEDAW-OUT/BOAT/LB-LPNO/NBODY'
LEDAW_output_path_LB_TPNO = r'./LEDAW-OUT/BOAT/LB-TPNO/NBODY'
# extrapolated to CPS and CBS
LEDAW_output_path_SB_CPS = r'./LEDAW-OUT/BOAT/SB-CPS/NBODY'
LEDAW_output_path_LB_CPS = r'./LEDAW-OUT/BOAT/LB-CPS/NBODY'
LEDAW_output_path_CBS_LPNO = r'./LEDAW-OUT/BOAT/CBS-LPNO/NBODY'
LEDAW_output_path_CBS_TPNO = r'./LEDAW-OUT/BOAT/CBS-TPNO/NBODY'
LEDAW_output_path_CBS_CPS = r'./LEDAW-OUT/BOAT/CBS-CPS/NBODY'


## Run N-Body LED engine for all computational settings.
#Standard and fp-LED N-body matrices will be written to excel files in specified LEDAW output directory

# for smaller basis set and looser TCutPNO setting
engine_LED_N_body(main_filenames=main_filenames_SB_LPNO, 
                  alternative_filenames=alternative_filenames_SB_LPNO, 
                  conversion_factor=conversion_factor, 
                  method=method,
                  relabel_mapping=relabel_mapping,
                  use_ref_as_rhf_in_hfld=use_ref_as_rhf_in_hfld,
                  LEDAW_output_path=LEDAW_output_path_SB_LPNO)

# for smaller basis set and tighter TCutPNO setting
engine_LED_N_body(main_filenames=main_filenames_SB_TPNO, 
                  alternative_filenames=alternative_filenames_SB_TPNO, 
                  conversion_factor=conversion_factor, 
                  method=method, 
                  relabel_mapping=relabel_mapping,
                  use_ref_as_rhf_in_hfld=use_ref_as_rhf_in_hfld,
                  LEDAW_output_path=LEDAW_output_path_SB_TPNO)

# for larger basis set and looser TCutPNO setting
engine_LED_N_body(main_filenames=main_filenames_LB_LPNO, 
                  alternative_filenames=alternative_filenames_LB_LPNO, 
                  conversion_factor=conversion_factor, 
                  method=method, 
                  relabel_mapping=relabel_mapping,
                  use_ref_as_rhf_in_hfld=use_ref_as_rhf_in_hfld,
                  LEDAW_output_path=LEDAW_output_path_LB_LPNO)

# for larger basis set and tighter TCutPNO setting
engine_LED_N_body(main_filenames=main_filenames_LB_TPNO, 
                  alternative_filenames=alternative_filenames_LB_TPNO, 
                  conversion_factor=conversion_factor, 
                  method=method, 
                  relabel_mapping=relabel_mapping,
                  use_ref_as_rhf_in_hfld=use_ref_as_rhf_in_hfld,
                  LEDAW_output_path=LEDAW_output_path_LB_TPNO)


### Extrapolate LPNO and TPNO N-body results to CPS limit.
# CPS for smaller basis set (SB) between looser and tighter TCutPNO settings (LPNO, TPNO)
extrapolate_engine(standard_LED_summary_file_X = LEDAW_output_path_SB_LPNO + r'/Summary_Standard_LED_matrices.xlsx', 
                   standard_LED_summary_file_Y = LEDAW_output_path_SB_TPNO + r'/Summary_Standard_LED_matrices.xlsx', 
                   fp_LED_summary_file_X = LEDAW_output_path_SB_LPNO + r'/Summary_fp-LED_matrices.xlsx', 
                   fp_LED_summary_file_Y = LEDAW_output_path_SB_TPNO + r'/Summary_fp-LED_matrices.xlsx', 
                   LEDAW_output_path = LEDAW_output_path_SB_CPS, 
                   F_ref = F, F_corr = F, method = method)

# CPS for larger basis set (LB) between looser and tighter TCutPNO settings (LPNO, TPNO)
extrapolate_engine(standard_LED_summary_file_X = LEDAW_output_path_LB_LPNO + r'/Summary_Standard_LED_matrices.xlsx', 
                   standard_LED_summary_file_Y = LEDAW_output_path_LB_TPNO + r'/Summary_Standard_LED_matrices.xlsx', 
                   fp_LED_summary_file_X = LEDAW_output_path_LB_LPNO + r'/Summary_fp-LED_matrices.xlsx', 
                   fp_LED_summary_file_Y = LEDAW_output_path_LB_TPNO + r'/Summary_fp-LED_matrices.xlsx', 
                   LEDAW_output_path = LEDAW_output_path_LB_CPS, 
                   F_ref = F, F_corr = F, method = method)

### Extrapolate (aug-)cc-pVTZ and (aug-)cc-pVQZ N-body energies to CBS limit.
# CBS for LPNO energies
extrapolate_engine(standard_LED_summary_file_X = LEDAW_output_path_SB_LPNO + r'/Summary_Standard_LED_matrices.xlsx', 
                   standard_LED_summary_file_Y = LEDAW_output_path_LB_TPNO + r'/Summary_Standard_LED_matrices.xlsx', 
                   fp_LED_summary_file_X = LEDAW_output_path_SB_LPNO + r'/Summary_fp-LED_matrices.xlsx', 
                   fp_LED_summary_file_Y = LEDAW_output_path_LB_LPNO + r'/Summary_fp-LED_matrices.xlsx', 
                   LEDAW_output_path = LEDAW_output_path_CBS_LPNO, 
                   F_ref = F_ref_cbs, F_corr = F_corr_cbs, method = method)

# CBS for TPNO energies
extrapolate_engine(standard_LED_summary_file_X = LEDAW_output_path_SB_TPNO + r'/Summary_Standard_LED_matrices.xlsx', 
                   standard_LED_summary_file_Y = LEDAW_output_path_LB_TPNO + r'/Summary_Standard_LED_matrices.xlsx', 
                   fp_LED_summary_file_X = LEDAW_output_path_SB_TPNO + r'/Summary_fp-LED_matrices.xlsx', 
                   fp_LED_summary_file_Y = LEDAW_output_path_LB_TPNO + r'/Summary_fp-LED_matrices.xlsx', 
                   LEDAW_output_path = LEDAW_output_path_CBS_TPNO, 
                   F_ref = F_ref_cbs, F_corr = F_corr_cbs, method = method)

# CBS for CPS-extrpolated energies. These are the final N-body energies to be considered.
extrapolate_engine(standard_LED_summary_file_X = LEDAW_output_path_SB_CPS + r'/Summary_Standard_LED_matrices.xlsx', 
                   standard_LED_summary_file_Y = LEDAW_output_path_LB_CPS + r'/Summary_Standard_LED_matrices.xlsx', 
                   fp_LED_summary_file_X = LEDAW_output_path_SB_CPS + r'/Summary_fp-LED_matrices.xlsx', 
                   fp_LED_summary_file_Y = LEDAW_output_path_LB_CPS + r'/Summary_fp-LED_matrices.xlsx', 
                   LEDAW_output_path = LEDAW_output_path_CBS_CPS, 
                   F_ref = F_ref_cbs, F_corr = F_corr_cbs, method = method)


twobody_title='''
#####################################################################################################
#                                          TWO-BODY LED                                             #
#####################################################################################################
'''
print(twobody_title)


### Specifying the one-body ORCA output files and their labels

# If onebody ORA output directory is specified without any file name, the code compares the labels of the supersystem file with 
# onebody files and automatically standardize the labeling of monomers to the original supersystem labeling.
# In this example, as relabel_mapping is initiated at the beginning of this file ([1,3,2,4,6,5]), labels will then be reordered.
# Note: onebody_out_directory directory must contain only the necessary one-body ORCA output files.
# In the following it is assumed that fragment labeling in all supersystem files are the same.
# Otherwise, you need to specify all the supersystem files that has differing labeling schemes.

supersystem_file = r'./ORCA-OUT/BOAT/aTZ/PNO6/HEXAMER/ADDUCT.mpi16.out' 
onebody_out_directory_SB_LPNO = r'./ORCA-OUT/BOAT/aTZ/PNO6/ONEBODY'
onebody_out_directory_SB_TPNO = r'./ORCA-OUT/BOAT/aTZ/PNO7/ONEBODY'
onebody_out_directory_LB_LPNO = r'./ORCA-OUT/BOAT/aQZ/PNO6/ONEBODY'
onebody_out_directory_LB_TPNO = r'./ORCA-OUT/BOAT/aQZ/PNO7/ONEBODY'

one_body_orcaout_filenames_SB_LPNO = extract_one_body_orcaout_filenames(supersystem_file, onebody_out_directory=onebody_out_directory_SB_LPNO)
one_body_orcaout_filenames_SB_TPNO = extract_one_body_orcaout_filenames(supersystem_file, onebody_out_directory=onebody_out_directory_SB_TPNO)
one_body_orcaout_filenames_LB_LPNO = extract_one_body_orcaout_filenames(supersystem_file, onebody_out_directory=onebody_out_directory_LB_LPNO)
one_body_orcaout_filenames_LB_TPNO = extract_one_body_orcaout_filenames(supersystem_file, onebody_out_directory=onebody_out_directory_LB_TPNO)


# Alternatively, you can enter the list of ORCA output files manually. But, in this approach, the order you specify the file names
# matters: (a) if it is as in the supersystem output file, i.e.,:
# one_body_orcaout_filenames_SB_LPNO = [r'./ORCA-OUT/BOAT/aTZ/PNO6/ONEBODY/frag1.mpi4.out',
#                                       r'./ORCA-OUT/BOAT/aTZ/PNO6/ONEBODY/frag3.mpi4.out',
#                                       r'./ORCA-OUT/BOAT/aTZ/PNO6/ONEBODY/frag2.mpi4.out',
#                                       r'./ORCA-OUT/BOAT/aTZ/PNO6/ONEBODY/frag4.mpi4.out',
#                                       r'./ORCA-OUT/BOAT/aTZ/PNO6/ONEBODY/frag6.mpi4.out',
#                                       r'./ORCA-OUT/BOAT/aTZ/PNO6/ONEBODY/frag5.mpi4.out',
#                                      ]
# ... (for the other computational settings, onebody orcaout file names must be given as separate lists.)   
# As relabel_mapping is initiated at the beginning of this file, these will then be reordered.
# (b) If you provide the order of the files consistent with the outcome of relabel_mapping intiated at the beginning of this file,
# then you need to reset relabel_mapping to None or False (avoid double reordering) to get exact same result as those obtained with
# the above two ways.


### Specify the two-body ORCA output file directories for each computational setting.
# Note: These directories must contain only the necessary two-body ORCA output files.
# The code will automatically read the files in this directory and label the fragments
# consistent first with the labeling in "one_body_orcaout_filenames" variables above and
# then with that in supersystem files. 
two_body_orcaout_directory_SB_LPNO = r'./ORCA-OUT/BOAT/aTZ/PNO6/TWOBODY'
two_body_orcaout_directory_SB_TPNO = r'./ORCA-OUT/BOAT/aTZ/PNO7/TWOBODY'
two_body_orcaout_directory_LB_LPNO = r'./ORCA-OUT/BOAT/aQZ/PNO6/TWOBODY'
two_body_orcaout_directory_LB_TPNO = r'./ORCA-OUT/BOAT/aQZ/PNO7/TWOBODY'


### Specify the directories where LEDAW will write Two-Body LED Matrices, which will be
# read from ORCA output files
LEDAW_output_path_twobody_SB_LPNO = r'./LEDAW-OUT/BOAT/SB-LPNO/TWOBODY'
LEDAW_output_path_twobody_SB_TPNO = r'./LEDAW-OUT/BOAT/SB-TPNO/TWOBODY'
LEDAW_output_path_twobody_LB_LPNO = r'./LEDAW-OUT/BOAT/LB-LPNO/TWOBODY'
LEDAW_output_path_twobody_LB_TPNO = r'./LEDAW-OUT/BOAT/LB-TPNO/TWOBODY'
# extrapolated to CPS and CBS
LEDAW_output_path_twobody_SB_CPS = r'./LEDAW-OUT/BOAT/SB-CPS/TWOBODY'
LEDAW_output_path_twobody_LB_CPS = r'./LEDAW-OUT/BOAT/LB-CPS/TWOBODY'
LEDAW_output_path_twobody_CBS_LPNO = r'./LEDAW-OUT/BOAT/CBS-LPNO/TWOBODY'
LEDAW_output_path_twobody_CBS_TPNO = r'./LEDAW-OUT/BOAT/CBS-TPNO/TWOBODY'
LEDAW_output_path_twobody_CBS_CPS = r'./LEDAW-OUT/BOAT/CBS-CPS/TWOBODY'


### Run two-body LED engine for all computational settings.
# Standard and fp-LED two-body matrices will be written to excel files in specified LEDAW two-body output directory

# for smaller basis set and looser TCutPNO setting
engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_SB_LPNO,
                    two_body_orcaout_directory = two_body_orcaout_directory_SB_LPNO,
                    conversion_factor = conversion_factor, 
                    method = method,
                    relabel_mapping=relabel_mapping,
                    use_ref_as_rhf_in_hfld=use_ref_as_rhf_in_hfld,
                    LEDAW_output_path_two_body = LEDAW_output_path_twobody_SB_LPNO)

# for smaller basis set and tighter TCutPNO setting
engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_SB_TPNO,
                    two_body_orcaout_directory = two_body_orcaout_directory_SB_TPNO,
                    conversion_factor = conversion_factor, 
                    method = method,
                    relabel_mapping=relabel_mapping,
                    use_ref_as_rhf_in_hfld=use_ref_as_rhf_in_hfld,
                    LEDAW_output_path_two_body = LEDAW_output_path_twobody_SB_TPNO)

# for larger basis set and looser TCutPNO setting
engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_LB_LPNO,
                    two_body_orcaout_directory = two_body_orcaout_directory_LB_LPNO,
                    conversion_factor = conversion_factor, 
                    method = method,
                    relabel_mapping=relabel_mapping,
                    use_ref_as_rhf_in_hfld=use_ref_as_rhf_in_hfld,
                    LEDAW_output_path_two_body = LEDAW_output_path_twobody_LB_LPNO)

# for larger basis set and tighter TCutPNO setting
engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_LB_TPNO,
                    two_body_orcaout_directory = two_body_orcaout_directory_LB_TPNO,
                    conversion_factor = conversion_factor, 
                    method = method,
                    relabel_mapping=relabel_mapping,
                    use_ref_as_rhf_in_hfld=use_ref_as_rhf_in_hfld,
                    LEDAW_output_path_two_body = LEDAW_output_path_twobody_LB_TPNO)


### Extrapolate LPNO and TPNO two-body results to CPS limit.
# CPS for smaller basis set (SB) between looser and tighter TCutPNO settings (LPNO, TPNO)
extrapolate_engine(standard_LED_summary_file_X = LEDAW_output_path_twobody_SB_LPNO + r'/Summary_Standard_LED_matrices.xlsx', 
                   standard_LED_summary_file_Y = LEDAW_output_path_twobody_SB_TPNO + r'/Summary_Standard_LED_matrices.xlsx', 
                   fp_LED_summary_file_X = LEDAW_output_path_twobody_SB_LPNO + r'/Summary_fp-LED_matrices.xlsx', 
                   fp_LED_summary_file_Y = LEDAW_output_path_twobody_SB_TPNO + r'/Summary_fp-LED_matrices.xlsx', 
                   LEDAW_output_path = LEDAW_output_path_twobody_SB_CPS, 
                   F_ref = F, F_corr = F, method = method)

# CPS for larger basis set (LB) between looser and tighter TCutPNO settings (LPNO, TPNO)
extrapolate_engine(standard_LED_summary_file_X = LEDAW_output_path_twobody_LB_LPNO + r'/Summary_Standard_LED_matrices.xlsx', 
                   standard_LED_summary_file_Y = LEDAW_output_path_twobody_LB_TPNO + r'/Summary_Standard_LED_matrices.xlsx', 
                   fp_LED_summary_file_X = LEDAW_output_path_twobody_LB_LPNO + r'/Summary_fp-LED_matrices.xlsx', 
                   fp_LED_summary_file_Y = LEDAW_output_path_twobody_LB_TPNO + r'/Summary_fp-LED_matrices.xlsx', 
                   LEDAW_output_path = LEDAW_output_path_twobody_LB_CPS, 
                   F_ref = F, F_corr = F, method = method)

### Extrapolate (aug-)cc-pVTZ and (aug-)cc-pVQZ two-body energies to CBS limit.
# CBS for LPNO energies
extrapolate_engine(standard_LED_summary_file_X = LEDAW_output_path_twobody_SB_LPNO + r'/Summary_Standard_LED_matrices.xlsx', 
                   standard_LED_summary_file_Y = LEDAW_output_path_twobody_LB_TPNO + r'/Summary_Standard_LED_matrices.xlsx', 
                   fp_LED_summary_file_X = LEDAW_output_path_twobody_SB_LPNO + r'/Summary_fp-LED_matrices.xlsx', 
                   fp_LED_summary_file_Y = LEDAW_output_path_twobody_LB_LPNO + r'/Summary_fp-LED_matrices.xlsx', 
                   LEDAW_output_path = LEDAW_output_path_twobody_CBS_LPNO, 
                   F_ref = F_ref_cbs, F_corr = F_corr_cbs, method = method)

# CBS for TPNO energies
extrapolate_engine(standard_LED_summary_file_X = LEDAW_output_path_twobody_SB_TPNO + r'/Summary_Standard_LED_matrices.xlsx', 
                   standard_LED_summary_file_Y = LEDAW_output_path_twobody_LB_TPNO + r'/Summary_Standard_LED_matrices.xlsx', 
                   fp_LED_summary_file_X = LEDAW_output_path_twobody_SB_TPNO + r'/Summary_fp-LED_matrices.xlsx', 
                   fp_LED_summary_file_Y = LEDAW_output_path_twobody_LB_TPNO + r'/Summary_fp-LED_matrices.xlsx', 
                   LEDAW_output_path = LEDAW_output_path_twobody_CBS_TPNO, 
                   F_ref = F_ref_cbs, F_corr = F_corr_cbs, method = method)

# CBS for CPS-extrpolated energies. These are the final N-body energies to be considered.
extrapolate_engine(standard_LED_summary_file_X = LEDAW_output_path_twobody_SB_CPS + r'/Summary_Standard_LED_matrices.xlsx', 
                   standard_LED_summary_file_Y = LEDAW_output_path_twobody_LB_CPS + r'/Summary_Standard_LED_matrices.xlsx', 
                   fp_LED_summary_file_X = LEDAW_output_path_twobody_SB_CPS + r'/Summary_fp-LED_matrices.xlsx', 
                   fp_LED_summary_file_Y = LEDAW_output_path_twobody_LB_CPS + r'/Summary_fp-LED_matrices.xlsx', 
                   LEDAW_output_path = LEDAW_output_path_twobody_CBS_CPS, 
                   F_ref = F_ref_cbs, F_corr = F_corr_cbs, method = method)


coop_title='''
##################################################################################################
#                                      COOPERATIVITY                                             #
##################################################################################################
'''
print(coop_title)


# Calculate cooperativity LED matrices from the excel files under the base_path/*/NBODY and base_path/*/TWOBODY directories
# and save them to the base_path/*/COOPERATIVITY directory generated (directory_level=2)
# Default directory_level = 1
cooperativity_engine(base_path=r'./LEDAW-OUT/BOAT', nbody_dir_name='NBODY', twobody_dir_name='TWOBODY', directory_level=2)


heat_map_title='''
##################################################################################################
#                                       HEAT MAP PLOT                                            #
##################################################################################################
'''
print(heat_map_title)


# Define plot parameters for standard LED heat maps
plot_params_for_std_led_matrices = {
    "figsize": (6, 6),
    "vmin": -100,
    "vmax": 100,
    "fig_format": "tif",
    "set_dpi": 400,
    "cutoff_annot": None, # The values smaller than absolute value of cutoff_annot are not shown on the heat maps. Default: None
    "submatrix_coords_to_be_highlighted": None, # to enclose the border of specified submatrix with black box. See docstring
    "display_heatmap": False # If true, the plot is also displayed in the environment (python  IDE, VS, jupyter notebook, etc) where the code is executed.
}


# Define plot parameters for fp-LED heat maps
plot_params_for_fp_led_matrices = {
    "figsize": (6, 6),
    "vmin": -15,
    "vmax": 15,
    "fig_format": "tif",
    "set_dpi": 400,
    "cutoff_annot": None, # The values smaller than absolute value of cutoff_annot are not shown on the heat maps. Default: None
    "submatrix_coords_to_be_highlighted": None, # to enclose the border of a specified submatrix with black box. See docstring
    "display_heatmap": False # If true, the plot is displayed in the environment (python  IDE, VS, jupyter notebook, etc) where the code is executed.
}

# Generate heat maps from Summary_Standard_LED_matrices.xlsx and Summary_fp_LED_matrices.xlsx.
# Save them to HEAT-MAP-STD-LED and HEAT-MAP-fp-LED directories generated in base_path (if directory_level=0),
# or its first or second level subdirectories (directory_level= 1 or 2): Default = 0
heatmap_plot_engine(base_path=r'./LEDAW-OUT/BOAT', 
                    plot_params_for_std_led_matrices=plot_params_for_std_led_matrices,
                    plot_params_for_fp_led_matrices=plot_params_for_fp_led_matrices,
                    show_diag_cells_for_fp_led=False,
                    delete_existing_heatmap_directories_first=True,
                    directory_level=2)
