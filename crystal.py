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

### The following is an automated LEDAW procedure for obtaining LED interactoin energy matrices and heat maps 
### on the interaction of a central monomer in a crystal with its environment from ORCA output files.
### For more details, see https://chemrxiv.org/engage/chemrxiv/article-details/6698104e01103d79c547414c
### N-body, two-body, and cooperativty LED maps are generated within standard and fp-LED schemes,
### followed by the plot of the corresponding heat maps.

from ledaw_package import *

### Choose the method
# Note 1: Place the proper method your ORCA outputs belong to. Choices: DLPNO-CCSD(T), DLPNO-CCSD, HFLD
# Note 2: Method name is not case sensitive.
# Note 3: The code does not check for the method name in the output files. It is necessary for properly
# combining decomposed LED terms.

method = "HFLD"

# In HFLD calculations, for one-body subsystems, RHF energies rather than reference E(0) energies can be 
# used in conjunction with LED in calculating electronic preparation energies with negligibly small variations.
# If you want to proceed this way, set use_ref_as_rhf_in_hfld = True (default: False or None). In this case, 
# if LEDAW could not find E(0) energies, it will use RHF energies of one-body subsystems in HFLD/LED calculations.
# CAUTION: This approach cannot use used in open-shell calculations as E(0) corresponds to QRO energy.

use_ref_as_rhf_in_hfld = True

### Convert Hartree to kcal/mol
conversion_factor = 627.5095  # for kj/mol, use: 2625.5

# If you are not happy with the fragment labeling in ORCA supersystem (adduct) LED file
# you can reorder the labels with relabel_mapping variable (default: None, equivalent to [1,2,3,..,the_number_of_fragments])
# As an example, relabel_mapping=[5,3,2,1,4] swaps label 1 with label 5; 2 with 3; 3 with 2; 4 with 1; and 5 with 4.
# Note: This is different than inconsistent fragment labelling the user introduced in supersystem and subsystem
# ORCA jobs. The code automatically standardize subsystem fragment labelings to those in supersystem file.  

relabel_mapping = [1,6,4,5,3,2,7,8,9,10,11] 
# Applied for consistency with https://chemrxiv.org/engage/chemrxiv/article-details/6698104e01103d79c547414c


nbody_title='''
##############################################################################################################
#                                               N-BODY LED                                                   #             
##############################################################################################################
'''

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

main_filenames = [r'./ORCA-OUT/CRYSTAL/MULTIFRAG/dimer.mpi4.out', 
                          r'./ORCA-OUT/CRYSTAL/MULTIFRAG/env.mpi2.out',
                          r'./ORCA-OUT/CRYSTAL/ONEBODY/mono1.mpi16.out',
                 ]

alternative_filenames = ['', '', '']


### Specify the Directories where LEDAW will write N-Body LED Matrices, which will be
# read from ORCA output files
LEDAW_output_path = r'./LEDAW-OUT/CRYSTAL/NBODY'

### Run N-Body LED engine for all computational settings.
# Standard and fp-LED N-body matrices will be written to excel files in specified LEDAW output directory

# for smaller basis set and looser TCutPNO setting
engine_LED_N_body(main_filenames=main_filenames, 
                  alternative_filenames=alternative_filenames, 
                  conversion_factor=conversion_factor, 
                  method=method,
                  relabel_mapping=relabel_mapping,
                  use_ref_as_rhf_in_hfld=use_ref_as_rhf_in_hfld,
                  LEDAW_output_path=LEDAW_output_path)

#-----------------------------------------------------------------------------------------------------------------
#                                 TYPICAL LED ANALYSIS JOB ENDED HERE. 
# If you further want assess cooperative effect of the interactions, proceed with the corresponding ORCA output 
# files as in the following "TWO-BOY LED" and "COOPREATIVITY". blocks.
# Otherwise commnet/delete them until "HEAT MAP PLOT "
#-----------------------------------------------------------------------------------------------------------------


twobody_title='''
#####################################################################################################
#                                          TWO-BODY LED                                             #
#####################################################################################################
'''
print(twobody_title)

### Get the one-body ORCA output files and their labels

## First Way ##
# Compare the labels of supersystem file with onebody files and automatically 
# standardize the labelling of monomers to the original supersystem labelling.
# In this example, as relabel_mapping is initiated at the beginning of this file ([1,6,4,5,3,2,7,8,9,10,11]),
# labels will be reordered.
# Note: onebody_out_directory directory must contain only the necessary one-body ORCA output files.

supersystem_file = r'./ORCA-OUT/CRYSTAL/MULTIFRAG/dimer.mpi4.out'
onebody_out_directory = r'./ORCA-OUT/CRYSTAL/ONEBODY'
one_body_orcaout_filenames = extract_one_body_orcaout_filenames(supersystem_file, onebody_out_directory)

## Second Way ##
# You can manually specify one-boy output files with the order consistent with the original supersystem labelling.
# relabel_mapping from N-body ([1,6,4,5,3,2,7,8,9,10,11]) is still active. Thus final LED maps will have reordered labels. 
#
# one_body_orcaout_filenames = [r'./ORCA-OUT/CRYSTAL/ONEBODY/mono1.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono2.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono3.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono4.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono5.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono6.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono7.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono8.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono9.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono10.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono11.mpi16.out',
#                              ]

## Third Way ##
# Specify one-body fragments in the order consistent with the reordered N-body labels.
# Then you need to set relabel_mapping=None in order not to reorder twice.

# relabel_mapping=None
# one_body_orcaout_filenames = [r'./ORCA-OUT/CRYSTAL/ONEBODY/mono1.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono6.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono5.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono3.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono4.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono2.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono7.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono8.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono9.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono10.mpi16.out',
#                               r'./ORCA-OUT/CRYSTAL/ONEBODY/mono11.mpi16.out',
#                              ]


### Specify the two-boy ORCA output file directory
# Note: This directory must contain only the necessary two-body ORCA output files.
# The code will automatically read the files in this directory and label the fragments consistent with the order you specified fragments in "one_body_orcaout_filenames" variables above. 
two_body_orcaout_directory = r'./ORCA-OUT/CRYSTAL/TWOBODY'


### Specify the Directory where LEDAW will write Two-Body LED Matrices
LEDAW_output_path_twobody = r'./LEDAW-OUT/CRYSTAL/TWOBODY'


### Run two-body LED engine
# Standard and fp-LED two-body matrices will be written to excel files in specified LEDAW two-body output directory
engine_LED_two_body(one_body_orcaout_filenames=one_body_orcaout_filenames,
                    two_body_orcaout_directory=two_body_orcaout_directory,
                    conversion_factor=conversion_factor, 
                    method=method,
                    LEDAW_output_path_two_body=LEDAW_output_path_twobody,
                    use_ref_as_rhf_in_hfld=use_ref_as_rhf_in_hfld,
                    relabel_mapping=relabel_mapping)


coop_title='''
# ##################################################################################################
# #                                      COOPERATIVITY                                             #
# ##################################################################################################
# '''

print(coop_title)

# Calculate cooperativity LED matrices from the excel files under the base_path/NBODY and base_path/TWOBODY directories
# and save them to the base_path/COOPERATIVITY directory generated (if directory_level=1),
# or its second or third level subdirectories (directory_level= 2 or 3): Default = 1
cooperativity_engine(base_path=r'./LEDAW-OUT/CRYSTAL', nbody_dir_name='NBODY', twobody_dir_name='TWOBODY', directory_level=1)



heat_map_title='''
##################################################################################################
#                                       HEAT MAP PLOT                                            #
##################################################################################################
'''

print(heat_map_title)

# Define plot parameters for standard LED heat maps
plot_params_for_std_led_matrices = {
    "figsize": (12, 12),
    "vmin": -60,
    "vmax": 60,
    "fig_format": "tif",
    "set_dpi": 400,
    "cutoff_annot": None, # The values smaller than absolute value of cutoff_annot are not shown on the heat maps. Default: None
    "submatrix_coords_to_be_highlighted": ((0,1),(0,10)), # to enclose the border of a specified submatrix with black box. See docstring
    "display_heatmap": False # If true, the plot is also displayed in the environment (python  IDE, VS, jupyter notebook, etc) where the code is executed.
}


# Define plot parameters for fp-LED heat maps
plot_params_for_fp_led_matrices = {
    "figsize": (12, 12),
    "vmin": -30,
    "vmax": 30,
    "fig_format": "tif",
    "set_dpi": 400,
    "cutoff_annot": None, # The values smaller than absolute value of cutoff_annot are not shown on the heat maps. Default: None
    "submatrix_coords_to_be_highlighted": ((0,1),(0,10)), # to enclose the border of specified submatrix with black box. See docstring
    "display_heatmap": False # If True, the plot is displayed in the environment (python  IDE, VS, jupyter notebook, etc) where the code is executed.
}

# Generate heat maps from Summary_Standard_LED_matrices.xlsx and Summary_fp_LED_matrices.xlsx.
# Save them to HEAT-MAP-STD-LED and HEAT-MAP-fp-LED directories generated in base_path (if directory_level=0),
# or its first or second level subdirectories (directory_level= 1 or 2): Default = 0
heatmap_plot_engine(base_path=r'./LEDAW-OUT/CRYSTAL', 
                    plot_params_for_std_led_matrices=plot_params_for_std_led_matrices,
                    plot_params_for_fp_led_matrices=plot_params_for_fp_led_matrices,
                    show_diag_cells_for_fp_led=False,
                    delete_existing_heatmap_directories_first=True,
                    directory_level=1)
