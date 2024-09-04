import os
import pandas as pd


def extrapolate_matrices(X, Y, F):
    """Extrapolate matrix A and B using factor F."""
    return X + F * (Y - X)


def categorize_sheets(sheet_name):
    """Categorize the sheets based on their names."""
    sheet_name_upper = sheet_name.upper()
    reference = ['REF', 'ELECTROSTAT', 'EXCHANGE', 'REF-EL-PREP', 'DIEL']
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
    diel = summary_sheets.get('DIEL', pd.DataFrame())

    if method.lower() == 'dlpno-ccsd(t)':
        c_ccsd_t = summary_sheets.get('C-CCSD(T)', pd.DataFrame())
        if diel.empty:
            total = ref + c_ccsd_t
        else:
            total = ref + c_ccsd_t + diel
    elif method.lower() == 'dlpno-ccsd':
        c_ccsd = summary_sheets.get('C-CCSD', pd.DataFrame())
        if diel.empty:
            total = ref + c_ccsd
        else:
            total = ref + c_ccsd + diel
    elif method.lower() == 'hfld':
        disp_hfld = summary_sheets.get('Disp HFLD', pd.DataFrame())
        if diel.empty:
            total = ref + disp_hfld
        else:
            total = ref + disp_hfld + diel
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


def extrapolate_engine(standard_LED_summary_file_X, standard_LED_summary_file_Y, fp_LED_summary_file_X, fp_LED_summary_file_Y, LEDAW_output_path, F_ref, F_corr, method):
    """Process and extrapolate matrices from two directories and save to a new directory."""

    # Create the extrapolation output directory if it doesn't exist
    if not os.path.exists(LEDAW_output_path):
        os.makedirs(LEDAW_output_path)
    
    # Files to be written
    extrapolated_standard_LED_summary_file = os.path.join(LEDAW_output_path, 'Summary_Standard_LED_matrices.xlsx')
    extrapolated_fp_LED_summary_file = os.path.join(LEDAW_output_path, 'Summary_fp-LED_matrices.xlsx')
    
    # Process the standard LED files
    with pd.ExcelWriter(extrapolated_standard_LED_summary_file, engine='openpyxl') as writer_standard:
        for file1, file2 in [(standard_LED_summary_file_X, standard_LED_summary_file_Y)]:
            xl1 = pd.ExcelFile(file1)
            xl2 = pd.ExcelFile(file2)

            summary_sheets = {}
            for sheet_name in xl1.sheet_names:
                X = pd.read_excel(xl1, sheet_name=sheet_name, index_col=0)
                Y = pd.read_excel(xl2, sheet_name=sheet_name, index_col=0)
                
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
                'TOTAL', 'REF', 'Electrostat', 'Exchange',
                'C-CCSD' if method.lower() == 'dlpno-ccsd' else 'C-CCSD(T)',
                'Disp CCSD' if method.lower() == 'dlpno-ccsd' else 'Disp CCSD(T)' if method.lower() == 'dlpno-ccsd(t)' else 'Disp HFLD',
                'Inter-NonDisp-C-CCSD' if method.lower() == 'dlpno-ccsd' else 'Inter-NonDisp-C-CCSD(T)'
            ]
            for sheet_name in sheet_order:
                if sheet_name in summary_sheets:
                    summary_sheets[sheet_name].to_excel(writer_standard, sheet_name=sheet_name)

    # Process the fp-LED files
    with pd.ExcelWriter(extrapolated_fp_LED_summary_file, engine='openpyxl') as writer_fp:
        for file1, file2 in [(fp_LED_summary_file_X, fp_LED_summary_file_Y)]:
            xl1 = pd.ExcelFile(file1)
            xl2 = pd.ExcelFile(file2)

            summary_sheets = {}
            for sheet_name in xl1.sheet_names:
                X = pd.read_excel(xl1, sheet_name=sheet_name, index_col=0)
                Y = pd.read_excel(xl2, sheet_name=sheet_name, index_col=0)
                
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
                'TOTAL', 'DIEL', 'REF', 'Electrostat', 'Exchange', 'REF-EL-PREP',
                'C-CCSD' if method.lower() == 'dlpno-ccsd' else 'C-CCSD(T)',
                'Disp CCSD' if method.lower() == 'dlpno-ccsd' else 'Disp CCSD(T)' if method.lower() == 'dlpno-ccsd(t)' else 'Disp HFLD',
                'Inter-NonDisp-C-CCSD' if method.lower() == 'dlpno-ccsd' else 'Inter-NonDisp-C-CCSD(T)',
                'C-CCSD-EL-PREP' if method.lower() == 'dlpno-ccsd' else 'C-CCSD(T)-EL-PREP',
                'CCSD-EL-PREP' if method.lower() == 'dlpno-ccsd' else 'CCSD(T)-EL-PREP'
            ]
            for sheet_name in sheet_order:
                if sheet_name in summary_sheets:
                    summary_sheets[sheet_name].to_excel(writer_fp, sheet_name=sheet_name)

    print(f"  Extrapolation job was terminated NORMALLY. Extrapolated standard and fp-LED matrices are at {LEDAW_output_path}")
  