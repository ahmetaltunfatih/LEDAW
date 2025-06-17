import os
import pandas as pd
from .nbody_engine import normalize_path


def extrapolate_matrices(X, Y, F):
    """Extrapolate matrix A and B using factor F."""
    return X + F * (Y - X)


def categorize_sheets(sheet_name):
    """Categorize the sheets based on their names."""
    sheet_name_upper = sheet_name.upper()
    reference = ['REF', 'ELECTROSTAT', 'EXCHANGE', 'REF-EL-PREP', 'SOLV']
    correlation_prefixes = ['C-', 'DISP', 'INTER-NONDISP']
    total = ['TOTAL', 'CCSD-EL-PREP', 'CCSD(T)-EL-PREP']

    if any(sheet_name_upper.startswith(ref) for ref in reference):
        return 'Reference'
    elif any(sheet_name_upper.startswith(cor) for cor in correlation_prefixes):
        return 'Correlation'
    elif any(sheet_name_upper.startswith(tot) for tot in total):
        return 'Total'
    else:
        return None


def calculate_total(summary_sheets, method):
    """Calculate the TOTAL matrix based on the specified method."""
    ref = summary_sheets.get('REF', pd.DataFrame())
    solv = summary_sheets.get('SOLV', pd.DataFrame())

    if method.lower() == 'dlpno-ccsd(t)':
        c_ccsd_t = summary_sheets.get('C-CCSD(T)', pd.DataFrame())
        if solv.empty:
            total = ref + c_ccsd_t
        else:
            total = ref + c_ccsd_t + solv
    elif method.lower() == 'dlpno-ccsd':
        c_ccsd = summary_sheets.get('C-CCSD', pd.DataFrame())
        if solv.empty:
            total = ref + c_ccsd
        else:
            total = ref + c_ccsd + solv
    elif method.lower() == 'hfld':
        disp_hfld = summary_sheets.get('Disp HFLD', pd.DataFrame())
        if solv.empty:
            total = ref + disp_hfld
        else:
            total = ref + disp_hfld + solv
    else:
        total = pd.DataFrame()  # Empty DataFrame if method is unrecognized
    
    return total


def calculate_ccsd_elprep(summary_sheets, method):
    """Calculate CCSD-EL-PREP or CCSD(T)-EL-PREP based on the method."""
    ref_elprep = summary_sheets.get('REF-EL-PREP', pd.DataFrame())

    if method.lower() == 'dlpno-ccsd(t)':
        c_ccsd_t_elprep = summary_sheets.get('C-CCSD(T)-EL-PREP', pd.DataFrame())
        return ref_elprep + c_ccsd_t_elprep
    elif method.lower() == 'dlpno-ccsd':
        c_ccsd_elprep = summary_sheets.get('C-CCSD-EL-PREP', pd.DataFrame())
        return ref_elprep + c_ccsd_elprep
    else:
        return pd.DataFrame()  # Empty DataFrame if method is unrecognized


def extrapolate_engine(standard_LED_summary_file_X, standard_LED_summary_file_Y,
                       fp_LED_summary_file_X, fp_LED_summary_file_Y,
                       LEDAW_output_path, F_ref, F_corr, method):
    """
    Processes and extrapolates matrices from two sets of summary and solvent files,
    and saves the results to a new directory. Solvent file paths are automatically derived.
    """

    # Normalize the output path
    normalized_LEDAW_output_path = normalize_path(LEDAW_output_path)

    # Create the extrapolation output directory if it doesn't exist
    if not os.path.exists(normalized_LEDAW_output_path):
        os.makedirs(normalized_LEDAW_output_path)
    
    # Files to be written
    extrapolated_standard_LED_summary_file = os.path.join(normalized_LEDAW_output_path, 'Summary_Standard_LED_matrices.xlsx')
    extrapolated_fp_LED_summary_file = os.path.join(normalized_LEDAW_output_path, 'Summary_fp-LED_matrices.xlsx')
    
    # New files to be written for solvent data
    extrapolated_standard_LED_solvent_file_output = os.path.join(normalized_LEDAW_output_path, 'SOLV-STD.xlsx')
    extrapolated_fp_LED_solvent_file_output = os.path.join(normalized_LEDAW_output_path, 'SOLV-fp.xlsx')


    # --- Derive solvent file paths based on summary file paths ---
    standard_LED_solvent_file_X = None
    standard_LED_solvent_file_Y = None
    if standard_LED_summary_file_X and os.path.exists(standard_LED_summary_file_X):
        standard_LED_solvent_file_X = os.path.join(os.path.dirname(standard_LED_summary_file_X), 'SOLV-STD.xlsx')
    if standard_LED_summary_file_Y and os.path.exists(standard_LED_summary_file_Y):
        standard_LED_solvent_file_Y = os.path.join(os.path.dirname(standard_LED_summary_file_Y), 'SOLV-STD.xlsx')

    fp_LED_solvent_file_X = None
    fp_LED_solvent_file_Y = None
    if fp_LED_summary_file_X and os.path.exists(fp_LED_summary_file_X):
        fp_LED_solvent_file_X = os.path.join(os.path.dirname(fp_LED_summary_file_X), 'SOLV-fp.xlsx')
    if fp_LED_summary_file_Y and os.path.exists(fp_LED_summary_file_Y):
        fp_LED_solvent_file_Y = os.path.join(os.path.dirname(fp_LED_summary_file_Y), 'SOLV-fp.xlsx')


    # Process the standard LED summary files
    with pd.ExcelWriter(extrapolated_standard_LED_summary_file, engine='openpyxl') as writer_standard:
        # Check if both summary files exist before proceeding
        if os.path.exists(standard_LED_summary_file_X) and os.path.exists(standard_LED_summary_file_Y):
            xl1 = pd.ExcelFile(standard_LED_summary_file_X)
            xl2 = pd.ExcelFile(standard_LED_summary_file_Y)

            summary_sheets = {}
            for sheet_name in xl1.sheet_names:
                try:
                    X = pd.read_excel(xl1, sheet_name=sheet_name, index_col=0)
                    # Attempt to read from Y, expecting the sheet to exist for extrapolation
                    Y = pd.read_excel(xl2, sheet_name=sheet_name, index_col=0) 
                except ValueError:
                    print(f"Warning: Sheet '{sheet_name}' not found in '{standard_LED_summary_file_Y}'. Skipping extrapolation for this sheet in standard summary files.")
                    continue
                
                category = categorize_sheets(sheet_name)
                
                if category == 'Reference':
                    F = F_ref
                elif category == 'Correlation':
                    F = F_corr
                else:
                    continue  # Skip uncategorized sheets
                
                extrapolated_matrix = extrapolate_matrices(X, Y, F)
                summary_sheets[sheet_name] = extrapolated_matrix

            # Calculate and store the TOTAL sheet
            total_matrix = calculate_total(summary_sheets, method)
            if not total_matrix.empty:
                summary_sheets['TOTAL'] = total_matrix

            # Calculate and store the CCSD-EL-PREP or CCSD(T)-EL-PREP sheet
            ccsd_elprep_matrix = calculate_ccsd_elprep(summary_sheets, method)
            ccsd_elprep_name = 'CCSD(T)-EL-PREP' if method.lower() == 'dlpno-ccsd(t)' else 'CCSD-EL-PREP'
            if not ccsd_elprep_matrix.empty:
                summary_sheets[ccsd_elprep_name] = ccsd_elprep_matrix

            # Write all sheets to Excel
            sheet_order = [
                'TOTAL', 'SOLV', 'REF', 'Electrostat', 'Exchange',
                'C-CCSD' if method.lower() == 'dlpno-ccsd' else 'C-CCSD(T)',
                'Disp CCSD' if method.lower() == 'dlpno-ccsd' else 'Disp CCSD(T)' if method.lower() == 'dlpno-ccsd(t)' else 'Disp HFLD',
                'Inter-NonDisp-C-CCSD' if method.lower() == 'dlpno-ccsd' else 'Inter-NonDisp-C-CCSD(T)'
            ]
            
            # Ensure EL-PREP sheets are added if they exist and are not empty
            if not summary_sheets.get('REF-EL-PREP', pd.DataFrame()).empty:
                if 'REF-EL-PREP' not in sheet_order: sheet_order.append('REF-EL-PREP')
            if not summary_sheets.get(ccsd_elprep_name, pd.DataFrame()).empty:
                if ccsd_elprep_name not in sheet_order: sheet_order.append(ccsd_elprep_name)

            for sheet_name in sheet_order:
                if sheet_name in summary_sheets:
                    summary_sheets[sheet_name].to_excel(writer_standard, sheet_name=sheet_name)
        else:
            print(f"Standard LED summary files (X: {standard_LED_summary_file_X}, Y: {standard_LED_summary_file_Y}) not found. Skipping standard summary extrapolation.")

    # Process the fp-LED summary files
    with pd.ExcelWriter(extrapolated_fp_LED_summary_file, engine='openpyxl') as writer_fp:
        # Check if both summary files exist before proceeding
        if os.path.exists(fp_LED_summary_file_X) and os.path.exists(fp_LED_summary_file_Y):
            xl1 = pd.ExcelFile(fp_LED_summary_file_X)
            xl2 = pd.ExcelFile(fp_LED_summary_file_Y)

            summary_sheets = {}
            for sheet_name in xl1.sheet_names:
                try:
                    X = pd.read_excel(xl1, sheet_name=sheet_name, index_col=0)
                    # Attempt to read from Y, expecting the sheet to exist for extrapolation
                    Y = pd.read_excel(xl2, sheet_name=sheet_name, index_col=0) 
                except ValueError:
                    print(f"Warning: Sheet '{sheet_name}' not found in '{fp_LED_summary_file_Y}'. Skipping extrapolation for this sheet in fp summary files.")
                    continue
                
                category = categorize_sheets(sheet_name)
                
                if category == 'Reference':
                    F = F_ref
                elif category == 'Correlation':
                    F = F_corr
                else:
                    continue  # Skip uncategorized sheets
                
                extrapolated_matrix = extrapolate_matrices(X, Y, F)
                summary_sheets[sheet_name] = extrapolated_matrix

            # Calculate and store the TOTAL sheet
            total_matrix = calculate_total(summary_sheets, method)
            if not total_matrix.empty:
                summary_sheets['TOTAL'] = total_matrix

            # Calculate and store the CCSD-EL-PREP or CCSD(T)-EL-PREP sheet
            ccsd_elprep_matrix = calculate_ccsd_elprep(summary_sheets, method)
            ccsd_elprep_name = 'CCSD(T)-EL-PREP' if method.lower() == 'dlpno-ccsd(t)' else 'CCSD-EL-PREP'
            if not ccsd_elprep_matrix.empty:
                summary_sheets[ccsd_elprep_name] = ccsd_elprep_matrix

            # Write all sheets to Excel
            sheet_order = [
                'TOTAL', 'SOLV', 'REF', 'Electrostat', 'Exchange', 'REF-EL-PREP',
                'C-CCSD' if method.lower() == 'dlpno-ccsd' else 'C-CCSD(T)',
                'Disp CCSD' if method.lower() == 'dlpno-ccsd' else 'Disp CCSD(T)' if method.lower() == 'dlpno-ccsd(t)' else 'Disp HFLD',
                'Inter-NonDisp-C-CCSD' if method.lower() == 'dlpno-ccsd' else 'Inter-NonDisp-C-CCSD(T)',
                'C-CCSD-EL-PREP' if method.lower() == 'dlpno-ccsd' else 'C-CCSD(T)-EL-PREP',
                'CCSD-EL-PREP' if method.lower() == 'dlpno-ccsd' else 'CCSD(T)-EL-PREP'
            ]
            for sheet_name in sheet_order:
                if sheet_name in summary_sheets:
                    summary_sheets[sheet_name].to_excel(writer_fp, sheet_name=sheet_name)
        else:
            print(f"fp-LED summary files (X: {fp_LED_summary_file_X}, Y: {fp_LED_summary_file_Y}) not found. Skipping fp summary extrapolation.")


    # Process the standard LED solvent files
    print("\nAttempting to extrapolate standard solvent files...")
    if standard_LED_solvent_file_X and os.path.exists(standard_LED_solvent_file_X) and \
       standard_LED_solvent_file_Y and os.path.exists(standard_LED_solvent_file_Y):
        try:
            xl_solv_X = pd.ExcelFile(standard_LED_solvent_file_X)
            xl_solv_Y = pd.ExcelFile(standard_LED_solvent_file_Y)

            with pd.ExcelWriter(extrapolated_standard_LED_solvent_file_output, engine='openpyxl') as writer_solv_std:
                for sheet_name in xl_solv_X.sheet_names:
                    try:
                        X = pd.read_excel(xl_solv_X, sheet_name=sheet_name, index_col=0)
                        Y = pd.read_excel(xl_solv_Y, sheet_name=sheet_name, index_col=0)
                        
                        extrapolated_matrix = extrapolate_matrices(X, Y, F_ref) # Use F_ref for solvent sheets
                        extrapolated_matrix.to_excel(writer_solv_std, sheet_name=sheet_name)
                        print(f"Extrapolated standard LED solvent '{sheet_name}' term written to '{normalize_path(extrapolated_standard_LED_solvent_file_output)}'.")
                    except ValueError:
                        print(f"Warning: '{sheet_name}' not found in '{standard_LED_solvent_file_Y}'. Skipping extrapolation for this standard LED solvent term.")
                    except Exception as e:
                        print(f"Error processing sheet '{sheet_name}' in standard solvent files: {e}. Skipping.")
            print(f"  Extrapolation for standard LED solvent files completed and saved to '{normalize_path(extrapolated_standard_LED_solvent_file_output)}'.")
        except Exception as e:
            print(f"Error opening/processing standard LED solvent files: {e}. Skipping standard LED solvent extrapolation.")
    else:
        print(f"Standard solvent files not found or accessed. Skipping corresponding extrapolation.")

    # Process the fp-LED solvent files
    print("\nAttempting to extrapolate fp-LED solvent files...")
    if fp_LED_solvent_file_X and os.path.exists(fp_LED_solvent_file_X) and \
       fp_LED_solvent_file_Y and os.path.exists(fp_LED_solvent_file_Y):
        try:
            xl_solv_fp_X = pd.ExcelFile(fp_LED_solvent_file_X)
            xl_solv_fp_Y = pd.ExcelFile(fp_LED_solvent_file_Y)

            with pd.ExcelWriter(extrapolated_fp_LED_solvent_file_output, engine='openpyxl') as writer_solv_fp:
                for sheet_name in xl_solv_fp_X.sheet_names:
                    try:
                        X = pd.read_excel(xl_solv_fp_X, sheet_name=sheet_name, index_col=0)
                        Y = pd.read_excel(xl_solv_fp_Y, sheet_name=sheet_name, index_col=0)
                        
                        extrapolated_matrix = extrapolate_matrices(X, Y, F_ref) # Use F_ref for solvent sheets
                        extrapolated_matrix.to_excel(writer_solv_fp, sheet_name=sheet_name)
                        print(f"Extrapolated fp-LED solvent contribution '{sheet_name}' written to '{normalize_path(extrapolated_fp_LED_solvent_file_output)}'.")
                    except ValueError:
                        print(f"Warning: '{sheet_name}' not found in '{fp_LED_solvent_file_Y}'. Skipping corresponding extrapolation.")
                    except Exception as e:
                        print(f"Error processing '{sheet_name}' fp-LED solvent term: {e}. Skipping.")
            print(f"  Extrapolation for fp-LED solvent terms completed and saved to '{normalize_path(extrapolated_fp_LED_solvent_file_output)}'.")
        except Exception as e:
            print(f"Error opening/processing fp-LED solvent files: {e}. Skipping fp-LED solvent extrapolation.")
    else:
        print("fp-LED solvent files not found or accessed. Skipping corresponding extrapolation.")

    print(f"\n  Extrapolation job was terminated NORMALLY. Extrapolated standard and fp-LED matrices are at{normalized_LEDAW_output_path}")
  