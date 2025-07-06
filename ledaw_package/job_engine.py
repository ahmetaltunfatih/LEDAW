import os
import shutil
from .nbody_engine import engine_LED_N_body
from .nbody_to_twobody_engine import extract_one_body_orcaout_filenames, get_reduced_relabel_mapping
from .twobody_engine import engine_LED_two_body
from .cooperativity_engine import cooperativity_engine
from .extrapolate_engine import extrapolate_engine
from .plot_engine import heatmap_plot_engine
from PySide6.QtCore import QThread, Signal


def normalize_twobody_path(path):
    """Normalize the path by converting backslashes to forward slashes and removing trailing slashes."""
    # Replace backslashes with forward slashes for consistency
    normalized_path = path.replace("\\", "/")
    return normalized_path.rstrip("/")  # remove any trailing slashes


def run_cps_cbs(app_instance):
    # Manage  N-body and/or two-body jobs with CPS and CBS extrapolations with optional plot

    ##########################################  N-BODY LED  #########################################

    if app_instance.perform_nbody:

        if app_instance.worker.is_cancelled:
            return
        
        # N-body LEDAW output directories
        LEDAW_output_path_SB_LPNO = os.path.join(app_instance.ledaw_out_root, "SB-LPNO", "NBODY")
        LEDAW_output_path_SB_TPNO = os.path.join(app_instance.ledaw_out_root, "SB-TPNO", "NBODY")
        LEDAW_output_path_LB_LPNO = os.path.join(app_instance.ledaw_out_root, "LB-LPNO", "NBODY")
        LEDAW_output_path_LB_TPNO = os.path.join(app_instance.ledaw_out_root, "LB-TPNO", "NBODY")
        LEDAW_output_path_SB_CPS = os.path.join(app_instance.ledaw_out_root, "SB-CPS", "NBODY")
        LEDAW_output_path_LB_CPS = os.path.join(app_instance.ledaw_out_root, "LB-CPS", "NBODY")
        LEDAW_output_path_CBS_LPNO = os.path.join(app_instance.ledaw_out_root, "CBS-LPNO", "NBODY")
        LEDAW_output_path_CBS_TPNO = os.path.join(app_instance.ledaw_out_root, "CBS-TPNO", "NBODY")
        LEDAW_output_path_CBS_CPS = os.path.join(app_instance.ledaw_out_root, "CBS-CPS", "NBODY")

        if app_instance.worker.is_cancelled:
            return

        ## Run N-Body LED engine for all computational settings.

        # for smaller basis set and looser TCutPNO setting
        engine_LED_N_body(main_filenames=app_instance.main_filenames_SB_LPNO, 
                          alternative_filenames=app_instance.alternative_filenames_SB_LPNO, 
                          conversion_factor=app_instance.conversion_factor, 
                          method=app_instance.method.currentText(),
                          relabel_mapping=app_instance.relabel_mappings["SB_LPNO"],
                          use_ref_as_rhf_in_hfld=app_instance.use_ref_as_rhf_in_hfld,
                          LEDAW_output_path=LEDAW_output_path_SB_LPNO)

        if app_instance.worker.is_cancelled:
            return
        
        # for smaller basis set and tighter TCutPNO setting
        engine_LED_N_body(main_filenames=app_instance.main_filenames_SB_TPNO, 
                          alternative_filenames=app_instance.alternative_filenames_SB_TPNO, 
                          conversion_factor=app_instance.conversion_factor, 
                          method=app_instance.method.currentText(),
                          relabel_mapping=app_instance.relabel_mappings["SB_TPNO"],
                          use_ref_as_rhf_in_hfld=app_instance.use_ref_as_rhf_in_hfld,
                          LEDAW_output_path=LEDAW_output_path_SB_TPNO)

        if app_instance.worker.is_cancelled:
            return

        # for larger basis set and looser TCutPNO setting
        engine_LED_N_body(main_filenames=app_instance.main_filenames_LB_LPNO, 
                          alternative_filenames=app_instance.alternative_filenames_LB_LPNO, 
                          conversion_factor=app_instance.conversion_factor, 
                          method=app_instance.method.currentText(),
                          relabel_mapping=app_instance.relabel_mappings["LB_LPNO"],
                          use_ref_as_rhf_in_hfld=app_instance.use_ref_as_rhf_in_hfld,
                          LEDAW_output_path=LEDAW_output_path_LB_LPNO)

        if app_instance.worker.is_cancelled:
            return

        # for larger basis set and tighter TCutPNO setting
        engine_LED_N_body(main_filenames=app_instance.main_filenames_LB_TPNO, 
                          alternative_filenames=app_instance.alternative_filenames_LB_TPNO, 
                          conversion_factor=app_instance.conversion_factor, 
                          method=app_instance.method.currentText(),
                          relabel_mapping=app_instance.relabel_mappings["LB_TPNO"],
                          use_ref_as_rhf_in_hfld=app_instance.use_ref_as_rhf_in_hfld,
                          LEDAW_output_path=LEDAW_output_path_LB_TPNO)

        if app_instance.worker.is_cancelled:
            return


        ### Extrapolate LPNO and TPNO N-body results to CPS limit.
        # CPS for smaller basis set (SB) between looser and tighter TCutPNO settings (LPNO, TPNO)
        extrapolate_engine(standard_LED_summary_file_X = os.path.join(LEDAW_output_path_SB_LPNO, 'Summary_Standard_LED_matrices.xlsx'), 
                           standard_LED_summary_file_Y = os.path.join(LEDAW_output_path_SB_TPNO, 'Summary_Standard_LED_matrices.xlsx'), 
                           fp_LED_summary_file_X = os.path.join(LEDAW_output_path_SB_LPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           fp_LED_summary_file_Y = os.path.join(LEDAW_output_path_SB_TPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           LEDAW_output_path = LEDAW_output_path_SB_CPS, 
                           F_ref = app_instance.F, F_corr = app_instance.F, method = app_instance.method.currentText())

        if app_instance.worker.is_cancelled:
            return

        # CPS for larger basis set (LB) between looser and tighter TCutPNO settings (LPNO, TPNO)
        extrapolate_engine(standard_LED_summary_file_X = os.path.join(LEDAW_output_path_LB_LPNO, 'Summary_Standard_LED_matrices.xlsx'), 
                           standard_LED_summary_file_Y = os.path.join(LEDAW_output_path_LB_TPNO, 'Summary_Standard_LED_matrices.xlsx'), 
                           fp_LED_summary_file_X = os.path.join(LEDAW_output_path_LB_LPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           fp_LED_summary_file_Y = os.path.join(LEDAW_output_path_LB_TPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           LEDAW_output_path = LEDAW_output_path_LB_CPS, 
                           F_ref = app_instance.F, F_corr = app_instance.F, method = app_instance.method.currentText())

        if app_instance.worker.is_cancelled:
            return

        ### Extrapolate (aug-)cc-pVTZ and (aug-)cc-pVQZ N-body energies to CBS limit.
        # CBS for LPNO energies
        extrapolate_engine(standard_LED_summary_file_X = os.path.join(LEDAW_output_path_SB_LPNO, 'Summary_Standard_LED_matrices.xlsx'),
                           standard_LED_summary_file_Y = os.path.join(LEDAW_output_path_LB_LPNO, 'Summary_Standard_LED_matrices.xlsx'),
                           fp_LED_summary_file_X = os.path.join(LEDAW_output_path_SB_LPNO, 'Summary_fp-LED_matrices.xlsx'),
                           fp_LED_summary_file_Y = os.path.join(LEDAW_output_path_LB_LPNO, 'Summary_fp-LED_matrices.xlsx'),
                           LEDAW_output_path = LEDAW_output_path_CBS_LPNO, 
                           F_ref = app_instance.f_ref, F_corr = app_instance.f_corr, method = app_instance.method.currentText())

        if app_instance.worker.is_cancelled:
            return

        # CBS for TPNO energies
        extrapolate_engine(standard_LED_summary_file_X = os.path.join(LEDAW_output_path_SB_TPNO, 'Summary_Standard_LED_matrices.xlsx'),
                           standard_LED_summary_file_Y = os.path.join(LEDAW_output_path_LB_TPNO, 'Summary_Standard_LED_matrices.xlsx'),
                           fp_LED_summary_file_X = os.path.join(LEDAW_output_path_SB_LPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           fp_LED_summary_file_Y = os.path.join(LEDAW_output_path_LB_LPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           LEDAW_output_path = LEDAW_output_path_CBS_TPNO, 
                           F_ref = app_instance.f_ref, F_corr = app_instance.f_corr, method = app_instance.method.currentText())

        if app_instance.worker.is_cancelled:
            return

        # CBS for CPS-extrpolated energies. These are the final N-body energies to be considered.
        extrapolate_engine(standard_LED_summary_file_X = os.path.join(LEDAW_output_path_SB_CPS, 'Summary_Standard_LED_matrices.xlsx'), 
                           standard_LED_summary_file_Y = os.path.join(LEDAW_output_path_LB_CPS, 'Summary_Standard_LED_matrices.xlsx'), 
                           fp_LED_summary_file_X = os.path.join(LEDAW_output_path_SB_CPS, 'Summary_fp-LED_matrices.xlsx'), 
                           fp_LED_summary_file_Y = os.path.join(LEDAW_output_path_LB_CPS, 'Summary_fp-LED_matrices.xlsx'), 
                           LEDAW_output_path = LEDAW_output_path_CBS_CPS, 
                           F_ref = app_instance.f_ref, F_corr = app_instance.f_corr, method = app_instance.method.currentText())    

        if app_instance.worker.is_cancelled:
            return
    
    ##########################################  TWO-BODY LED  #########################################
    if app_instance.perform_twobody:
 
        if app_instance.worker.is_cancelled:
            return

        # Two-body LEDAW output directories
        LEDAW_output_path_twobody_SB_LPNO = os.path.join(app_instance.ledaw_out_root, "SB-LPNO", "TWOBODY")
        LEDAW_output_path_twobody_SB_TPNO = os.path.join(app_instance.ledaw_out_root, "SB-TPNO", "TWOBODY")
        LEDAW_output_path_twobody_LB_LPNO = os.path.join(app_instance.ledaw_out_root, "LB-LPNO", "TWOBODY")
        LEDAW_output_path_twobody_LB_TPNO = os.path.join(app_instance.ledaw_out_root, "LB-TPNO", "TWOBODY")
        LEDAW_output_path_twobody_SB_CPS = os.path.join(app_instance.ledaw_out_root, "SB-CPS", "TWOBODY")
        LEDAW_output_path_twobody_LB_CPS = os.path.join(app_instance.ledaw_out_root, "LB-CPS", "TWOBODY")
        LEDAW_output_path_twobody_CBS_LPNO = os.path.join(app_instance.ledaw_out_root, "CBS-LPNO", "TWOBODY")
        LEDAW_output_path_twobody_CBS_TPNO = os.path.join(app_instance.ledaw_out_root, "CBS-TPNO", "TWOBODY")
        LEDAW_output_path_twobody_CBS_CPS = os.path.join(app_instance.ledaw_out_root, "CBS-CPS", "TWOBODY")

        if app_instance.worker.is_cancelled:
            return

        if app_instance.perform_nbody:
            
            if app_instance.worker.is_cancelled:
                return

            # One-body ORCA output files whose fragment labelings to be syncronized with the corresponding N-body main supersystem files.    
            one_body_orcaout_filenames_SB_LPNO = extract_one_body_orcaout_filenames(app_instance.main_filenames_SB_LPNO[0], app_instance.onebody_orcaout_directory_SB_LPNO)
            one_body_orcaout_filenames_SB_TPNO = extract_one_body_orcaout_filenames(app_instance.main_filenames_SB_TPNO[0], app_instance.onebody_orcaout_directory_SB_TPNO)
            one_body_orcaout_filenames_LB_LPNO = extract_one_body_orcaout_filenames(app_instance.main_filenames_LB_LPNO[0], app_instance.onebody_orcaout_directory_LB_LPNO)
            one_body_orcaout_filenames_LB_TPNO = extract_one_body_orcaout_filenames(app_instance.main_filenames_LB_TPNO[0], app_instance.onebody_orcaout_directory_LB_TPNO)

            if app_instance.worker.is_cancelled:
                return

            reduced_relabel_mapping_SB_LPNO = get_reduced_relabel_mapping(supersystem_file=app_instance.main_filenames_SB_LPNO[0], relabel_mapping=app_instance.relabel_mappings["SB_LPNO"])
            reduced_relabel_mapping_SB_TPNO = get_reduced_relabel_mapping(supersystem_file=app_instance.main_filenames_SB_TPNO[0], relabel_mapping=app_instance.relabel_mappings["SB_TPNO"])
            reduced_relabel_mapping_LB_LPNO = get_reduced_relabel_mapping(supersystem_file=app_instance.main_filenames_LB_LPNO[0], relabel_mapping=app_instance.relabel_mappings["LB_LPNO"])
            reduced_relabel_mapping_LB_TPNO = get_reduced_relabel_mapping(supersystem_file=app_instance.main_filenames_LB_TPNO[0], relabel_mapping=app_instance.relabel_mappings["LB_TPNO"])

            ### Run two-body LED engine for all computational settings.

            # for smaller basis set and looser TCutPNO setting
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_SB_LPNO,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory_SB_LPNO,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=reduced_relabel_mapping_SB_LPNO,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody_SB_LPNO,
                                LEDAW_output_path_nbody=LEDAW_output_path_SB_LPNO)


            if app_instance.worker.is_cancelled:
                return

            # for smaller basis set and tighter TCutPNO setting
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_SB_TPNO,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory_SB_TPNO,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=reduced_relabel_mapping_SB_TPNO,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody_SB_TPNO,
                                LEDAW_output_path_nbody=LEDAW_output_path_SB_TPNO)

            if app_instance.worker.is_cancelled:
                return

            # for larger basis set and looser TCutPNO setting
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_LB_LPNO,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory_LB_LPNO,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=reduced_relabel_mapping_LB_LPNO,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody_LB_LPNO,
                                LEDAW_output_path_nbody=LEDAW_output_path_LB_LPNO)

            if app_instance.worker.is_cancelled:
                return

            # for larger basis set and tighter TCutPNO setting
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_LB_TPNO,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory_LB_TPNO,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=reduced_relabel_mapping_LB_TPNO,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody_LB_TPNO,
                                LEDAW_output_path_nbody=LEDAW_output_path_LB_TPNO)

            if app_instance.worker.is_cancelled:
                return

        if not app_instance.perform_nbody:
            one_body_orcaout_filenames_SB_LPNO = app_instance.onebody_orcaout_files_SB_LPNO  
            one_body_orcaout_filenames_SB_TPNO = app_instance.onebody_orcaout_files_SB_TPNO 
            one_body_orcaout_filenames_LB_LPNO = app_instance.onebody_orcaout_files_LB_LPNO  
            one_body_orcaout_filenames_LB_TPNO = app_instance.onebody_orcaout_files_LB_TPNO 

            if app_instance.worker.is_cancelled:
                return

            ### Run two-body LED engine for all computational settings.

            # for smaller basis set and looser TCutPNO setting
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_SB_LPNO,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory_SB_LPNO,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=None,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody_SB_LPNO)

            if app_instance.worker.is_cancelled:
                return

            # for smaller basis set and tighter TCutPNO setting
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_SB_TPNO,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory_SB_TPNO,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=None,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody_SB_TPNO)

            if app_instance.worker.is_cancelled:
                return

            # for larger basis set and looser TCutPNO setting
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_LB_LPNO,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory_LB_LPNO,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=None,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody_LB_LPNO)

            if app_instance.worker.is_cancelled:
                return

            # for larger basis set and tighter TCutPNO setting
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_LB_TPNO,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory_LB_TPNO,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=None,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody_LB_TPNO)

            if app_instance.worker.is_cancelled:
                return

        ### Extrapolate LPNO and TPNO two-body results to CPS limit.
        # CPS for smaller basis set (SB) between looser and tighter TCutPNO settings (LPNO, TPNO)
        extrapolate_engine(standard_LED_summary_file_X = os.path.join(LEDAW_output_path_twobody_SB_LPNO, 'Summary_Standard_LED_matrices.xlsx'), 
                           standard_LED_summary_file_Y = os.path.join(LEDAW_output_path_twobody_SB_TPNO, 'Summary_Standard_LED_matrices.xlsx'), 
                           fp_LED_summary_file_X = os.path.join(LEDAW_output_path_twobody_SB_LPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           fp_LED_summary_file_Y = os.path.join(LEDAW_output_path_twobody_SB_TPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           LEDAW_output_path = LEDAW_output_path_twobody_SB_CPS, 
                           F_ref = app_instance.F, F_corr = app_instance.F, method = app_instance.method.currentText())

        if app_instance.worker.is_cancelled:
            return

        # CPS for larger basis set (LB) between looser and tighter TCutPNO settings (LPNO, TPNO)
        extrapolate_engine(standard_LED_summary_file_X = os.path.join(LEDAW_output_path_twobody_LB_LPNO, 'Summary_Standard_LED_matrices.xlsx'),  
                           standard_LED_summary_file_Y = os.path.join(LEDAW_output_path_twobody_LB_TPNO, 'Summary_Standard_LED_matrices.xlsx'),  
                           fp_LED_summary_file_X = os.path.join(LEDAW_output_path_twobody_SB_LPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           fp_LED_summary_file_Y = os.path.join(LEDAW_output_path_twobody_SB_TPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           LEDAW_output_path = LEDAW_output_path_twobody_LB_CPS, 
                           F_ref = app_instance.F, F_corr = app_instance.F, method = app_instance.method.currentText())
        if app_instance.worker.is_cancelled:
            return

        ### Extrapolate SB and LB two-body energies to CBS limit.
        # CBS for LPNO energies
        extrapolate_engine(standard_LED_summary_file_X = os.path.join(LEDAW_output_path_twobody_SB_LPNO, 'Summary_Standard_LED_matrices.xlsx'), 
                           standard_LED_summary_file_Y = os.path.join(LEDAW_output_path_twobody_LB_LPNO, 'Summary_Standard_LED_matrices.xlsx'), 
                           fp_LED_summary_file_X = os.path.join(LEDAW_output_path_twobody_SB_LPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           fp_LED_summary_file_Y = os.path.join(LEDAW_output_path_twobody_LB_LPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           LEDAW_output_path = LEDAW_output_path_twobody_CBS_LPNO, 
                           F_ref = app_instance.f_ref, F_corr = app_instance.f_corr, method = app_instance.method.currentText())

        if app_instance.worker.is_cancelled:
            return

        # CBS for TPNO energies
        extrapolate_engine(standard_LED_summary_file_X = os.path.join(LEDAW_output_path_twobody_SB_TPNO, 'Summary_Standard_LED_matrices.xlsx'), 
                           standard_LED_summary_file_Y = os.path.join(LEDAW_output_path_twobody_LB_TPNO, 'Summary_Standard_LED_matrices.xlsx'), 
                           fp_LED_summary_file_X = os.path.join(LEDAW_output_path_twobody_SB_TPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           fp_LED_summary_file_Y = os.path.join(LEDAW_output_path_twobody_LB_TPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           LEDAW_output_path = LEDAW_output_path_twobody_CBS_TPNO, 
                           F_ref = app_instance.f_ref, F_corr = app_instance.f_corr, method = app_instance.method.currentText())
        if app_instance.worker.is_cancelled:
            return

        # CBS for CPS-extrpolated energies. These are the final two-body energies to be considered.
        extrapolate_engine(standard_LED_summary_file_X = os.path.join(LEDAW_output_path_twobody_SB_CPS, 'Summary_Standard_LED_matrices.xlsx'), 
                           standard_LED_summary_file_Y = os.path.join(LEDAW_output_path_twobody_LB_CPS, 'Summary_Standard_LED_matrices.xlsx'), 
                           fp_LED_summary_file_X = os.path.join(LEDAW_output_path_twobody_SB_CPS, 'Summary_fp-LED_matrices.xlsx'), 
                           fp_LED_summary_file_Y = os.path.join(LEDAW_output_path_twobody_LB_CPS, 'Summary_fp-LED_matrices.xlsx'), 
                           LEDAW_output_path = LEDAW_output_path_twobody_CBS_CPS, 
                           F_ref = app_instance.f_ref, F_corr = app_instance.f_corr, method = app_instance.method.currentText())

        if app_instance.worker.is_cancelled:
            return

    ##########################################  COOPERATIVITY  #########################################

    if app_instance.perform_nbody  and app_instance.perform_twobody:
        
        if app_instance.worker.is_cancelled:
            return
        
        cooperativity_engine(base_path=app_instance.ledaw_out_root, nbody_dir_name='NBODY', twobody_dir_name='TWOBODY', directory_level=2)

        if app_instance.worker.is_cancelled:
            return

     #############################################  PLOT  #############################################       


    if app_instance.perform_plot:

        if app_instance.worker.is_cancelled:
            return           
        
        heatmap_plot_engine(app_instance, base_path=app_instance.ledaw_out_root,
                            plot_params_for_std_led_matrices=app_instance.plot_params_for_std_led_matrices,
                            plot_params_for_fp_led_matrices=app_instance.plot_params_for_fp_led_matrices,
                            show_diag_cells_for_fp_led=app_instance.show_diag_cells_for_fp_led,
                            delete_existing_heatmap_directories_first=app_instance.delete_old_plots,
                            directory_level=2)

        if app_instance.worker.is_cancelled:
            return           
            

def run_cps(app_instance):
      # Manage combined N-body or two-body job with CPS extrapolation with optional plot
    
    ##########################################  N-BODY LED  #########################################

    if app_instance.perform_nbody:
        
        if app_instance.worker.is_cancelled:
            return
    
        # N-body LEDAW output directories
        LEDAW_output_path_LPNO = os.path.join(app_instance.ledaw_out_root, "LPNO", "NBODY")
        LEDAW_output_path_TPNO = os.path.join(app_instance.ledaw_out_root, "TPNO", "NBODY")
        LEDAW_output_path_CPS = os.path.join(app_instance.ledaw_out_root, "CPS", "NBODY")

        if app_instance.worker.is_cancelled:
            return

        # for looser TCutPNO setting
        engine_LED_N_body(main_filenames=app_instance.main_filenames_LPNO, 
                          alternative_filenames=app_instance.alternative_filenames_LPNO, 
                          conversion_factor=app_instance.conversion_factor, 
                          method=app_instance.method.currentText(),
                          relabel_mapping=app_instance.relabel_mappings["LPNO"],
                          use_ref_as_rhf_in_hfld=app_instance.use_ref_as_rhf_in_hfld,
                          LEDAW_output_path=LEDAW_output_path_LPNO)

        if app_instance.worker.is_cancelled:
            return
        
        # for tighter TCutPNO setting
        engine_LED_N_body(main_filenames=app_instance.main_filenames_TPNO, 
                          alternative_filenames=app_instance.alternative_filenames_TPNO, 
                          conversion_factor=app_instance.conversion_factor, 
                          method=app_instance.method.currentText(),
                          relabel_mapping=app_instance.relabel_mappings["TPNO"],
                          use_ref_as_rhf_in_hfld=app_instance.use_ref_as_rhf_in_hfld,
                          LEDAW_output_path=LEDAW_output_path_TPNO)
        if app_instance.worker.is_cancelled:
            return

        ### Extrapolate LPNO and TPNO N-body results to CPS limit.
        extrapolate_engine(standard_LED_summary_file_X = os.path.join(LEDAW_output_path_LPNO, 'Summary_Standard_LED_matrices.xlsx'), 
                           standard_LED_summary_file_Y = os.path.join(LEDAW_output_path_TPNO, 'Summary_Standard_LED_matrices.xlsx'), 
                           fp_LED_summary_file_X = os.path.join(LEDAW_output_path_LPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           fp_LED_summary_file_Y = os.path.join(LEDAW_output_path_TPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           LEDAW_output_path = LEDAW_output_path_CPS, 
                           F_ref = app_instance.F, F_corr = app_instance.F, method = app_instance.method.currentText())

        if app_instance.worker.is_cancelled:
            return

    ##########################################  TWO-BODY LED  #########################################

    if app_instance.perform_twobody:
        if app_instance.worker.is_cancelled:
            return
        
        # Two-body LEDAW output directories
        LEDAW_output_path_twobody_LPNO = os.path.join(app_instance.ledaw_out_root, "LPNO", "TWOBODY")
        LEDAW_output_path_twobody_TPNO = os.path.join(app_instance.ledaw_out_root, "TPNO", "TWOBODY")
        LEDAW_output_path_twobody_CPS = os.path.join(app_instance.ledaw_out_root, "CPS", "TWOBODY")

        if app_instance.perform_nbody:

            if app_instance.worker.is_cancelled:
                return

            # One-body ORCA output files whose fragment labelings to be syncronized with the corresponding N-body main supersystem files.    
            one_body_orcaout_filenames_LPNO = extract_one_body_orcaout_filenames(app_instance.main_filenames_LPNO[0], app_instance.onebody_orcaout_directory_LPNO)
            one_body_orcaout_filenames_TPNO = extract_one_body_orcaout_filenames(app_instance.main_filenames_TPNO[0], app_instance.onebody_orcaout_directory_TPNO)

            reduced_relabel_mapping_LPNO = get_reduced_relabel_mapping(supersystem_file=app_instance.main_filenames_LPNO[0], relabel_mapping=app_instance.relabel_mappings["LPNO"])
            reduced_relabel_mapping_TPNO = get_reduced_relabel_mapping(supersystem_file=app_instance.main_filenames_TPNO[0], relabel_mapping=app_instance.relabel_mappings["TPNO"])

            ### Run two-body LED engine for all computational settings.

            if app_instance.worker.is_cancelled:
                return

            # for looser TCutPNO setting
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_LPNO,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory_LPNO,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=reduced_relabel_mapping_LPNO,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody_LPNO,
                                LEDAW_output_path_nbody=LEDAW_output_path_LPNO)

            if app_instance.worker.is_cancelled:
                return

            # for tighter TCutPNO setting
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_TPNO,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory_TPNO,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=reduced_relabel_mapping_TPNO,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody_TPNO,
                                LEDAW_output_path_nbody=LEDAW_output_path_TPNO)

            if app_instance.worker.is_cancelled:
                return

        if not app_instance.perform_nbody:
            
            if app_instance.worker.is_cancelled:
                return

            one_body_orcaout_filenames_LPNO = app_instance.onebody_orcaout_files_LPNO  
            one_body_orcaout_filenames_TPNO = app_instance.onebody_orcaout_files_TPNO 

            ### Run two-body LED engine for all computational settings.

            if app_instance.worker.is_cancelled:
                return

            # for looser TCutPNO setting
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_LPNO,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory_LPNO,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=None,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody_LPNO)

            if app_instance.worker.is_cancelled:
                return

            # for tighter TCutPNO setting
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_TPNO,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory_TPNO,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=None,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody_TPNO)

            if app_instance.worker.is_cancelled:
                return

        ### Extrapolate LPNO and TPNO two-body results to CPS limit
        extrapolate_engine(standard_LED_summary_file_X = os.path.join(LEDAW_output_path_twobody_LPNO, 'Summary_Standard_LED_matrices.xlsx'), 
                           standard_LED_summary_file_Y = os.path.join(LEDAW_output_path_twobody_TPNO, 'Summary_Standard_LED_matrices.xlsx'), 
                           fp_LED_summary_file_X = os.path.join(LEDAW_output_path_twobody_LPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           fp_LED_summary_file_Y = os.path.join(LEDAW_output_path_twobody_TPNO, 'Summary_fp-LED_matrices.xlsx'), 
                           LEDAW_output_path = LEDAW_output_path_twobody_CPS, 
                           F_ref = app_instance.F, F_corr = app_instance.F, method = app_instance.method.currentText())

        if app_instance.worker.is_cancelled:
            return

    ##########################################  COOPERATIVITY  #########################################

    if app_instance.perform_nbody  and app_instance.perform_twobody:

        if app_instance.worker.is_cancelled:
            return

        cooperativity_engine(base_path=app_instance.ledaw_out_root, nbody_dir_name='NBODY', twobody_dir_name='TWOBODY', directory_level=2)

        if app_instance.worker.is_cancelled:
            return

     #############################################  PLOT  #############################################       


    if app_instance.perform_plot:

        if app_instance.worker.is_cancelled:
            return           
        
        heatmap_plot_engine(app_instance, base_path=app_instance.ledaw_out_root,
                            plot_params_for_std_led_matrices=app_instance.plot_params_for_std_led_matrices,
                            plot_params_for_fp_led_matrices=app_instance.plot_params_for_fp_led_matrices,
                            show_diag_cells_for_fp_led=app_instance.show_diag_cells_for_fp_led,
                            delete_existing_heatmap_directories_first=app_instance.delete_old_plots,
                            directory_level=2)

        if app_instance.worker.is_cancelled:
            return           

        
def run_cbs(app_instance):
       # Manage combined N-body or two-body job with CBS extrapolation with optional plot
    
    ##########################################  N-BODY LED  #########################################

    if app_instance.perform_nbody:
        
        if app_instance.worker.is_cancelled:
            return
    
        # N-body LEDAW output directories
        LEDAW_output_path_SB = os.path.join(app_instance.ledaw_out_root, "SB", "NBODY")
        LEDAW_output_path_LB = os.path.join(app_instance.ledaw_out_root, "LB", "NBODY")
        LEDAW_output_path_CBS = os.path.join(app_instance.ledaw_out_root, "CBS", "NBODY")
        
        if app_instance.worker.is_cancelled:
            return

        # for smaller basis set
        engine_LED_N_body(main_filenames=app_instance.main_filenames_SB, 
                          alternative_filenames=app_instance.alternative_filenames_SB, 
                          conversion_factor=app_instance.conversion_factor, 
                          method=app_instance.method.currentText(),
                          relabel_mapping=app_instance.relabel_mappings["SB"],
                          use_ref_as_rhf_in_hfld=app_instance.use_ref_as_rhf_in_hfld,
                          LEDAW_output_path=LEDAW_output_path_SB)
        
        if app_instance.worker.is_cancelled:
            return

        # for larger basis set
        engine_LED_N_body(main_filenames=app_instance.main_filenames_LB, 
                          alternative_filenames=app_instance.alternative_filenames_LB, 
                          conversion_factor=app_instance.conversion_factor, 
                          method=app_instance.method.currentText(),
                          relabel_mapping=app_instance.relabel_mappings["LB"],
                          use_ref_as_rhf_in_hfld=app_instance.use_ref_as_rhf_in_hfld,
                          LEDAW_output_path=LEDAW_output_path_LB)
        
        if app_instance.worker.is_cancelled:
            return

        ### Extrapolate SB and LB N-body results to CBS limit.
        extrapolate_engine(standard_LED_summary_file_X = os.path.join(LEDAW_output_path_SB, 'Summary_Standard_LED_matrices.xlsx'), 
                           standard_LED_summary_file_Y = os.path.join(LEDAW_output_path_LB, 'Summary_Standard_LED_matrices.xlsx'), 
                           fp_LED_summary_file_X = os.path.join(LEDAW_output_path_SB, 'Summary_fp-LED_matrices.xlsx'), 
                           fp_LED_summary_file_Y = os.path.join(LEDAW_output_path_LB, 'Summary_fp-LED_matrices.xlsx'), 
                           LEDAW_output_path = LEDAW_output_path_CBS, 
                           F_ref = app_instance.f_ref, F_corr = app_instance.f_corr, method = app_instance.method.currentText())
        
        if app_instance.worker.is_cancelled:
            return
    
    ##########################################  TWO-BODY LED  #########################################

    if app_instance.perform_twobody:
        
        if app_instance.worker.is_cancelled:
            return
        
        # Two-body LEDAW output directories
        LEDAW_output_path_twobody_SB = os.path.join(app_instance.ledaw_out_root, "SB", "TWOBODY")
        LEDAW_output_path_twobody_LB = os.path.join(app_instance.ledaw_out_root, "LB", "TWOBODY")
        LEDAW_output_path_twobody_CBS = os.path.join(app_instance.ledaw_out_root, "CBS", "TWOBODY")
        
        if app_instance.worker.is_cancelled:
            return

        if app_instance.perform_nbody:
                    
            if app_instance.worker.is_cancelled:
                return

            # One-body ORCA output files whose fragment labelings to be syncronized with the corresponding N-body main supersystem files.    
            one_body_orcaout_filenames_SB = extract_one_body_orcaout_filenames(app_instance.main_filenames_SB[0], app_instance.onebody_orcaout_directory_SB)
            one_body_orcaout_filenames_LB = extract_one_body_orcaout_filenames(app_instance.main_filenames_LB[0], app_instance.onebody_orcaout_directory_LB)

            reduced_relabel_mapping_SB = get_reduced_relabel_mapping(supersystem_file=app_instance.main_filenames_SB[0], relabel_mapping=app_instance.relabel_mappings["SB"])
            reduced_relabel_mapping_LB = get_reduced_relabel_mapping(supersystem_file=app_instance.main_filenames_LB[0], relabel_mapping=app_instance.relabel_mappings["LB"])


            ### Run two-body LED engine for all computational settings.
            if app_instance.worker.is_cancelled:
                return

            # for smaller basis set
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_SB,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory_SB,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=reduced_relabel_mapping_SB,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody_SB,
                                LEDAW_output_path_nbody=LEDAW_output_path_SB)

            if app_instance.worker.is_cancelled:
                return

            # for larger basis set
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_LB,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory_LB,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=reduced_relabel_mapping_LB,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody_LB,
                                LEDAW_output_path_nbody=LEDAW_output_path_LB)

            if app_instance.worker.is_cancelled:
                return

        if not app_instance.perform_nbody:

            if app_instance.worker.is_cancelled:
                return

            one_body_orcaout_filenames_SB = app_instance.onebody_orcaout_files_SB
            one_body_orcaout_filenames_LB = app_instance.onebody_orcaout_files_LB
            
            ### Run two-body LED engine for all computational settings.

            if app_instance.worker.is_cancelled:
                return

            # for smaller basis set
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_SB,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory_SB,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=None,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody_SB)

            if app_instance.worker.is_cancelled:
                return

            # for larger basis set
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames_LB,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory_LB,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=None,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody_LB)

            if app_instance.worker.is_cancelled:
                return

        ### Extrapolate SB and LB two-body results to CBS limit
        extrapolate_engine(standard_LED_summary_file_X = os.path.join(LEDAW_output_path_twobody_SB, 'Summary_Standard_LED_matrices.xlsx'), 
                           standard_LED_summary_file_Y = os.path.join(LEDAW_output_path_twobody_LB, 'Summary_Standard_LED_matrices.xlsx'), 
                           fp_LED_summary_file_X = os.path.join(LEDAW_output_path_twobody_SB, 'Summary_fp-LED_matrices.xlsx'), 
                           fp_LED_summary_file_Y = os.path.join(LEDAW_output_path_twobody_LB, 'Summary_fp-LED_matrices.xlsx'), 
                           LEDAW_output_path = LEDAW_output_path_twobody_CBS, 
                           F_ref = app_instance.f_ref, F_corr = app_instance.f_corr, method = app_instance.method.currentText())

        if app_instance.worker.is_cancelled:
            return

    ##########################################  COOPERATIVITY  #########################################

    if app_instance.perform_nbody  and app_instance.perform_twobody:
        
        if app_instance.worker.is_cancelled:
            return

        cooperativity_engine(base_path=app_instance.ledaw_out_root, nbody_dir_name='NBODY', twobody_dir_name='TWOBODY', directory_level=2)

        if app_instance.worker.is_cancelled:
            return

     #############################################  PLOT  #############################################       


    if app_instance.perform_plot:

        if app_instance.worker.is_cancelled:
            return           
        
        heatmap_plot_engine(app_instance, base_path=app_instance.ledaw_out_root,
                            plot_params_for_std_led_matrices=app_instance.plot_params_for_std_led_matrices,
                            plot_params_for_fp_led_matrices=app_instance.plot_params_for_fp_led_matrices,
                            show_diag_cells_for_fp_led=app_instance.show_diag_cells_for_fp_led,
                            delete_existing_heatmap_directories_first=app_instance.delete_old_plots,
                            directory_level=2)

        if app_instance.worker.is_cancelled:
            return           


def run_without_cps_cbs(app_instance):
       # Manage N-body and/or two-body job without CPS and CBS extrapolations with optional plot
    
    ##########################################  N-BODY LED  #########################################

    if app_instance.perform_nbody:

        if app_instance.worker.is_cancelled:
            return

        # N-body LEDAW output directories
        LEDAW_output_path = os.path.join(app_instance.ledaw_out_root, "NBODY")

        # Run N-body LED engine        
        engine_LED_N_body(main_filenames=app_instance.main_filenames, 
                          alternative_filenames=app_instance.alternative_filenames, 
                          conversion_factor=app_instance.conversion_factor, 
                          method=app_instance.method.currentText(),
                          relabel_mapping=app_instance.relabel_mappings["Standard"],
                          use_ref_as_rhf_in_hfld=app_instance.use_ref_as_rhf_in_hfld,
                          LEDAW_output_path=LEDAW_output_path)

        if app_instance.worker.is_cancelled:
            return

    ##########################################  TWO-BODY LED  #########################################

    if app_instance.perform_twobody:
        
        if app_instance.worker.is_cancelled:
            return
        
        # Two-body LEDAW output directories
        LEDAW_output_path_twobody = os.path.join(app_instance.ledaw_out_root, "TWOBODY")

        if app_instance.perform_nbody:
            if app_instance.worker.is_cancelled:
                return

            # One-body ORCA output files whose fragment labelings to be syncronized with the corresponding N-body main supersystem files.    
            one_body_orcaout_filenames = extract_one_body_orcaout_filenames(app_instance.main_filenames[0], app_instance.onebody_orcaout_directory)

            if app_instance.worker.is_cancelled:
                return

            reduced_relabel_mapping = get_reduced_relabel_mapping(supersystem_file=app_instance.main_filenames[0], relabel_mapping=app_instance.relabel_mappings["Standard"])

            # Run two-body LED engine
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory,
                                conversion_factor=app_instance.conversion_factor, 
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=reduced_relabel_mapping,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody,
                                LEDAW_output_path_nbody=LEDAW_output_path)

            if app_instance.worker.is_cancelled:
                return

        if not app_instance.perform_nbody:
            one_body_orcaout_filenames = app_instance.onebody_orcaout_files
            
            if app_instance.worker.is_cancelled:
                return

           # Run two-body LED engine
            engine_LED_two_body(one_body_orcaout_filenames = one_body_orcaout_filenames,
                                two_body_orcaout_directory = app_instance.twobody_orcaout_directory,
                                conversion_factor=app_instance.conversion_factor,
                                method=app_instance.method.currentText(),
                                reduced_relabel_mapping=None,
                                LEDAW_output_path_two_body=LEDAW_output_path_twobody)
            if app_instance.worker.is_cancelled:
                return


    ##########################################  COOPERATIVITY  #########################################

    if app_instance.perform_nbody  and app_instance.perform_twobody:

        if app_instance.worker.is_cancelled:
            return

        cooperativity_engine(base_path=app_instance.ledaw_out_root, nbody_dir_name='NBODY', twobody_dir_name='TWOBODY', directory_level=1)

        if app_instance.worker.is_cancelled:
            return

     #############################################  PLOT  #############################################       


    if app_instance.perform_plot:

        if app_instance.worker.is_cancelled:
            return           
        
        heatmap_plot_engine(app_instance, base_path=app_instance.ledaw_out_root,
                            plot_params_for_std_led_matrices=app_instance.plot_params_for_std_led_matrices,
                            plot_params_for_fp_led_matrices=app_instance.plot_params_for_fp_led_matrices,
                            show_diag_cells_for_fp_led=app_instance.show_diag_cells_for_fp_led,
                            delete_existing_heatmap_directories_first=app_instance.delete_old_plots,
                            directory_level=1)

        if app_instance.worker.is_cancelled:
            return


############################################  LED JOB MANAGER  ###########################################

def manage_led_jobs(app_instance):
    try:
        if app_instance.F is not None and app_instance.f_ref is not None and app_instance.f_corr is not None:
            run_cps_cbs(app_instance)
            if app_instance.worker.is_cancelled:
                return
        elif app_instance.F is not None:
            run_cps(app_instance)
            if app_instance.worker.is_cancelled:
                return
        elif app_instance.f_ref is not None and app_instance.f_corr is not None:
            run_cbs(app_instance)
            if app_instance.worker.is_cancelled:
                return
        elif app_instance.F is None and app_instance.f_ref is None and app_instance.f_corr is None:
            run_without_cps_cbs(app_instance)
            if app_instance.worker.is_cancelled:
                return
        else:
            raise ValueError("Check the supplied parameters and files/directories and try again.")
    except Exception as e:
        print(f"Error in manage_led_jobs: {e}")
        raise


def manage_plot_jobs(app_instance):
        
    if app_instance.F is not None or app_instance.f_ref is not None or app_instance.f_corr is not None:
        directory_level=2
    else:
        directory_level=1
        
    try:
        if app_instance.worker.is_cancelled:
            return

        if app_instance.plot_locked:
            if app_instance.worker.is_cancelled:
                return

            heatmap_plot_engine(app_instance, base_path=app_instance.ledaw_out_root,
                                plot_params_for_std_led_matrices=app_instance.plot_params_for_std_led_matrices,
                                plot_params_for_fp_led_matrices=app_instance.plot_params_for_fp_led_matrices,
                                show_diag_cells_for_fp_led=app_instance.show_diag_cells_for_fp_led,
                                delete_existing_heatmap_directories_first=app_instance.delete_old_plots,
                                directory_level=directory_level)

            if app_instance.worker.is_cancelled:
                return

    except Exception as e:
        print(f"Error in manage_jobs: {e}")
        raise

class JobWorker(QThread):
    """Worker class for managing the LEDAW job in a separate thread"""

    # Define signal for status updates
    status_signal = Signal(str)

    def __init__(self, app_instance, parent=None):
        super().__init__(parent)
        self.app_instance = app_instance
        self.is_cancelled = False

    def run(self):
        try:
            # Emit "In progress" status to update UI
            self.status_signal.emit("In progress...")

            if self.is_cancelled:
                self.status_signal.emit("Cancelled\n\nIf needed, press Run LEDAW button to resubmit the job.") 
                return

            # Perform the LEDAW job, but allow for periodic cancellation checks
            if not self.is_cancelled:
                manage_led_jobs(self.app_instance)

                if self.is_cancelled:
                    self.status_signal.emit("Cancelled\n\nIf needed, press Run LEDAW button to resubmit the job.") 
                # Emit "Completed" status if the job finishes successfully and delete two-body tmp directory
                elif not self.is_cancelled and not self.app_instance.perform_plot:
                    self.status_signal.emit("Completed\n\nLED matrices were written to the LEDAW output directory.")
                    self.app_instance.delete_tmp_files()
                elif not self.is_cancelled and self.app_instance.perform_plot:
                    self.status_signal.emit('Completed\n\nLED matrices and heat maps were written to the LEDAW output directory.\n\nPlease check heat maps plotted. If necessary, adjust and then confirm plot parameters in the Plot tab to rerun the plot job.')
                    self.app_instance.delete_tmp_files()

        except Exception as e:
            # Emit "Failed" status
            self.status_signal.emit("Failed.\nPossible Reason: User specified files do not correspond to actual ORCA output files, or there is an issue with file permissions.")
            # Ensure cleanup happens even if an error occurs
            self.app_instance.delete_tmp_files()

    def cancel(self):
        """Call this method to cancel the job."""
        self.is_cancelled = True


class PlotJobWorker(QThread):
    """Worker class for managing the plot job in a separate thread"""

    # Define signal for status updates
    status_signal = Signal(str)

    def __init__(self, app_instance, parent=None):
        super().__init__(parent)
        self.app_instance = app_instance
        self.is_cancelled = False

    def run(self):
        try:
            # Emit "In progress" status to update UI
            self.status_signal.emit("In progress...")

            if self.is_cancelled:
                self.status_signal.emit("Cancelled\n\nIf needed, adjust and then confirm plot parameters in the Plot tab to resubmit the plot job.") 
                return

            # Perform the Plot job, but allow for periodic cancellation checks
            if not self.is_cancelled:
                manage_plot_jobs(self.app_instance)

                # Emit "Cancelled" status if the job finishes successfully
                if self.is_cancelled:
                    self.status_signal.emit("Cancelled\n\nIf needed, adjust and then confirm plot parameters in the Plot tab to resubmit the plot job.") 
                # Emit "Completed" status if the job finishes successfully
                elif not self.is_cancelled:
                    self.status_signal.emit("Completed\n\nIf needed, adjust and then confirm plot parameters in the Plot tab to rerun the plot job.") 

        except Exception as e:
            # Emit "Failed" status
            self.status_signal.emit("Failed.\nPossible Reason: Necessary Excel files that include interaction energy matrices could not be found, or there is an issue with file permissions.")

    def cancel(self):
        """Call this method to cancel the job."""
        self.is_cancelled = True
