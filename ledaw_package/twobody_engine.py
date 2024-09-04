import os
import re
import numpy as np
import pandas as pd
from .classes import *


def get_two_body_filenames(two_body_orcaout_directory):
    """Reads all files from the specified directory and saves their full paths to a list."""
    # Normalize the directory path to remove any trailing slashes
    two_body_orcaout_directory = os.path.normpath(two_body_orcaout_directory)

    # Initialize an empty list to store file paths
    two_body_orcaout_filenames = []

    # Walk through the directory and collect file paths
    for root, _, files in os.walk(two_body_orcaout_directory):
        for file in files:
            # Combine the root path with the file name to get the full file path
            full_path = os.path.normpath(os.path.join(root, file))
            two_body_orcaout_filenames.append(full_path)

    return two_body_orcaout_filenames


def extract_labels_from_one_body_files(one_body_orcaout_filenames, relabel_mapping=None, tolerance=1e-3):
    """Extract labels from one-body files based on the first line after *xyz pattern and store numerical values."""
    
    # Check if relabel_mapping was explicitly passed as an argument
    if relabel_mapping is None:
        # If relabel_mapping is None, generate a default mapping [1, 2, 3, ..., n]
        relabel_mapping = list(range(1, len(one_body_orcaout_filenames) + 1))
    
    label_mapping = {}
    
    for idx, filename in enumerate(one_body_orcaout_filenames, start=1):
        with open(filename, 'r') as file:
            content = file.read()
            match = re.search(r'\*\s*xyz[^\n]*\n(.*)', content)
            if match:
                next_line = match.group(1)
                numbers = re.findall(r'-?\d+\.\d+', next_line)
                numbers = [float(num) for num in numbers]
                original_label = idx  # idx corresponds to the original file label (1-based)
                remapped_label = relabel_mapping[original_label - 1]  # Adjusted for 0-based index
                label_mapping[remapped_label] = numbers
                
    return label_mapping


def extract_labels_from_two_body_files(two_body_orcaout_directory, label_mapping, tolerance=1e-3):
    """Extract labels from two-body files using coordinate comparison with a specified tolerance."""
    
    two_body_labels = {}

    # Iterate over all files in the directory
    for file_name in os.listdir(two_body_orcaout_directory):
        file_path = os.path.join(two_body_orcaout_directory, file_name)
        if not os.path.isfile(file_path):
            continue  # Skip if not a file
        
        with open(file_path, 'r') as file:
            lines = file.readlines()
        
        labels = []
        fragment_count = 0
        
        for line in lines:
            line = line.strip()
            
            if line.startswith("FRAGMENT"):
                fragment_count += 1
            
            if fragment_count > 0 and fragment_count <= 2:
                parts = line.split()
                if len(parts) == 4:  # Check for atomic symbol and three coordinates
                    try:
                        coordinates = list(map(float, parts[1:4]))
                    except ValueError:
                        continue  # Skip lines that don't have valid floats for coordinates
                    
                    for key, value in label_mapping.items():
                        # Compare coordinates with tolerance
                        if all(abs(c1 - c2) <= tolerance for c1, c2 in zip(coordinates, value)):
                            labels.append(key)
                            break
            
            if fragment_count == 2 and len(labels) == 2:
                break  # Stop after finding two labels
        
        if len(labels) == 2:
            two_body_labels[file_name] = labels
    
    return two_body_labels


def extract_onebody_dielectric_values(one_body_orcaout_filenames, dielectric_pattern, relabel_mapping=None):
    """Extract dielectric values from one-body files based on the relabel_mapping."""
    
    # Check if relabel_mapping was explicitly passed as an argument
    if relabel_mapping is None:
        # If relabel_mapping is None, generate a default mapping [1, 2, 3, ..., n]
        relabel_mapping = list(range(1, len(one_body_orcaout_filenames) + 1))
    
    dielectric_values = {}
    
    # Assume relabel_mapping is a list of labels corresponding to one_body_orcaout_filenames
    for filename, label in zip(one_body_orcaout_filenames, relabel_mapping):
        with open(filename, 'r') as file:
            content = file.read()
            match = re.search(dielectric_pattern, content)
            if match:
                dielectric_value = float(match.group(1))
                dielectric_values[label] = dielectric_value
    
    return dielectric_values


def extract_twobody_dielectric_values(two_body_orcaout_directory, two_body_labels, dielectric_pattern):
    """Extract dielectric values from two-body files."""
    dielectric_values = {}
    
    for filename, (label1, label2) in two_body_labels.items():
        file_path = os.path.join(two_body_orcaout_directory, filename)
        if not os.path.isfile(file_path):
            continue  # Skip if not a file
        
        with open(file_path, 'r') as file:
            content = file.read()
            match = re.search(dielectric_pattern, content)
            if match:
                dielectric_value = float(match.group(1))
                dielectric_values[filename] = (label1, label2, dielectric_value)
    
    return dielectric_values


def populate_twobody_dielectric_matrices(one_body_orcaout_filenames, two_body_orcaout_directory, conversion_factor, two_body_labels, LEDAW_output_path_two_body, relabel_mapping=None):
    """Populate the two-body dielectric interaction matrix and save it to an Excel file."""
    
    # Define the dielectric extraction pattern
    dielectric_pattern = r'CPCM Dielectric\s*:\s*(-?\d+\.\d+)\s*Eh'

    # Extract dielectric values from one-body files using relabel_mapping
    one_body_dielectric_values = extract_onebody_dielectric_values(one_body_orcaout_filenames, dielectric_pattern, relabel_mapping)
    
    # Extract dielectric values from two-body files
    two_body_dielectric_values = extract_twobody_dielectric_values(two_body_orcaout_directory, two_body_labels, dielectric_pattern)
    
    # Initialize the matrix
    n = len(one_body_orcaout_filenames)
    matrix = np.zeros((n, n))

    # Ensure relabel_mapping is not None
    if relabel_mapping is None:
        relabel_mapping = list(range(1, len(one_body_orcaout_filenames) + 1))

    # Check if any values are missing
    if len(one_body_dielectric_values) < len(relabel_mapping) or len(two_body_dielectric_values) < len(two_body_labels):
        print("At least one of the output files does not contain dielectric contribution.\nDielectric contribution to the interaction energy is taken as 0")
        matrix[:] = 0  # Set entire matrix to zero
    else:
        # Populate the matrix with dielectric interaction values
        for filename, (label1, label2, two_body_value) in two_body_dielectric_values.items():
            one_body_value_i = one_body_dielectric_values.get(label1, 0)
            one_body_value_j = one_body_dielectric_values.get(label2, 0)
            
            # Calculate the adjusted dielectric interaction value
            adjusted_value = (two_body_value - one_body_value_i - one_body_value_j) * conversion_factor

            # Only fill non-diagonal elements, set diagonal to NaN
            if label1 != label2:
                matrix[label1 - 1, label2 - 1] = adjusted_value
                matrix[label2 - 1, label1 - 1] = adjusted_value
            else:
                matrix[label1 - 1, label2 - 1] = np.nan

        # Set values below the diagonal to zero
        for i in range(1, n):
            for j in range(i):
                matrix[i, j] = 0
    
    # Convert to DataFrame
    df = pd.DataFrame(matrix, index=range(1, n + 1), columns=range(1, n + 1))

    # Write the matrix to an Excel file in the specified output directory
    if LEDAW_output_path_two_body:
        output_file_path = os.path.join(LEDAW_output_path_two_body, 'DIEL.xlsx')
        df.to_excel(output_file_path, sheet_name='DIEL')
        print(f"Two-body dielectric LED interaction energy matrix was written to {output_file_path}")


def populate_twobody_inter_matrices(one_body_orcaout_filenames, two_body_orcaout_directory, conversion_factor, relabel_mapping=None, two_body_labels=None, LEDAW_output_path_two_body=None):
    results = {}
    patterns = Patterns()
    patterns.PATTERNS['corr_inter_full'] = (
        r"Interaction correlation for Fragments\s+(\d+)\s+and\s+(\d+):\s+[-]+\s+"
        r"Inter strong pairs\s+([-]?\d*\.\d+)\s+\(.*?\)\s+"
        r"(Inter triples\s+([-]?\d*\.\d+)\s+\(.*?\)\s+)?"
        r"Inter weak pairs\s+([-]?\d*\.\d+)\s+\(.*?\)\s+[-]+\s+Total interaction"
    )

    property_mapping = {
        'Electrostat': {'pattern': 'ref_inter', 'group': 3},
        'Exchange': {'pattern': 'ref_inter', 'group': 4},
        'Inter SP': {'pattern': 'corr_inter_full', 'group': 3},
        'Inter T': {'pattern': 'corr_inter_full', 'group': 5},
        'Inter WP': {'pattern': 'corr_inter_full', 'group': 6},
        'Disp SP': {'pattern': 'dispersion_strong_pairs', 'group': 3},
    }

    # Check if relabel_mapping was explicitly passed as an argument
    if relabel_mapping is None:
        # If relabel_mapping is None, generate a default mapping [1, 2, 3, ..., n]
        relabel_mapping = list(range(1, len(one_body_orcaout_filenames) + 1))

    # Initialize matrices and populate them based on extracted data
    for sheet_name, props in property_mapping.items():
        n = len(one_body_orcaout_filenames)
        matrix = np.zeros((n, n))
        pattern = patterns.PATTERNS[props['pattern']]
        group_num = props['group']

        for filename, (label1, label2) in two_body_labels.items():
            file_path = os.path.join(two_body_orcaout_directory, filename)
            with open(file_path, 'r') as file:
                content = file.read()
                match = re.search(pattern, content)
                if match and len(match.groups()) >= group_num:
                    value = match.group(group_num)
                    if value:
                        value = value.strip()
                        numeric_value = float(value) * conversion_factor
                        i, j = label1 - 1, label2 - 1
                        if i < j:
                            matrix[i, j] = numeric_value
                        elif j < i:
                            matrix[j, i] = numeric_value

        df = pd.DataFrame(matrix, index=range(1, n+1), columns=range(1, n+1))
        results[sheet_name] = df

    # Save the results to an Excel file
    output_file_path = os.path.join(LEDAW_output_path_two_body, 'INTER.xlsx')
    with pd.ExcelWriter(output_file_path) as writer:
        for sheet_name, df in results.items():
            df.to_excel(writer, sheet_name=sheet_name)

    print(f"Two-body interfragment LED interaction energy matrices were written to {output_file_path}")


import re

def extract_onebody_values(one_body_filenames, patterns, method, use_ref_as_rhf_in_hfld=None, relabel_mapping=None):
    """Extract various values from one-body files using provided patterns and optionally a relabel_mapping."""
    one_body_values = {'ref': {}, 'strong_corr': {}, 'weak_corr': {}, 'triples_corr': {}}
    
    # If relabel_mapping is not provided, use enumerate to generate default labels
    if relabel_mapping is None:
        relabel_mapping = list(range(1, len(one_body_filenames) + 1))
    
    for filename, label in zip(one_body_filenames, relabel_mapping):
        if filename is None:
            continue  # Skip None filenames
        
        try:
            with open(filename, 'r') as file:
                content = file.read()

                # Attempt to extract the E(0) value for Intra REF
                ref_match = re.search(patterns.pattern_e0, content)
                if ref_match:
                    one_body_values['ref'][label] = float(ref_match.group(1))
                else:
                    # Fallback for HFLD method if E(0) is not found and use_ref_as_rhf_in_hfld is True
                    if method.lower() == 'hfld' and use_ref_as_rhf_in_hfld:
                        total_energy_match = re.search(r"Total Energy\s+:\s+([-]?\d+\.\d+)", content)
                        dielectric_match = re.search(r"CPCM Dielectric\s+:\s+([-]?\d+\.\d+)", content)
                        
                        if total_energy_match:
                            total_energy = float(total_energy_match.group(1))
                            dielectric_value = float(dielectric_match.group(1)) if dielectric_match else 0.0
                            e0_value = total_energy - dielectric_value
                            one_body_values['ref'][label] = e0_value
                        else:
                            raise ValueError(f"E(0) and fallback patterns not found in file: {filename}")

                # Extract the strong pairs correction
                strong_corr_match = re.search(patterns.pattern_strong_corr, content)
                if strong_corr_match:
                    one_body_values['strong_corr'][label] = float(strong_corr_match.group(1))

                # Extract the weak pairs correction
                weak_corr_match = re.search(patterns.pattern_weak_corr, content)
                if weak_corr_match:
                    one_body_values['weak_corr'][label] = float(weak_corr_match.group(1))

                # Extract the triples correction
                triples_corr_match = re.search(patterns.pattern_triples_corr, content)
                if triples_corr_match:
                    one_body_values['triples_corr'][label] = float(triples_corr_match.group(1))

        except FileNotFoundError:
            raise FileNotFoundError(f"File {filename} does not exist. Please check the path.")
        except Exception as e:
            raise RuntimeError(f"An error occurred while processing file {filename}: {e}")
    
    return one_body_values


def extract_twobody_values(two_body_directory, two_body_labels, pattern_key, patterns):
    """Process two-body files and extract values based on the provided pattern key."""
    extracted_values = {}

    for filename in os.listdir(two_body_directory):
        file_path = os.path.join(two_body_directory, filename)
        with open(file_path, 'r') as file:
            content = file.read()

            # Extract the section using the appropriate pattern
            pattern = patterns.PATTERNS[pattern_key]
            matches = re.findall(pattern, content)

            if matches and filename in two_body_labels:
                label1, label2 = two_body_labels[filename]

                # Default values for all potential corrections
                value1 = float(matches[0][0].strip())
                value2 = float(matches[0][1].strip()) if len(matches[0]) > 1 else 0.0  # Handle absence of "Inter T"
                value3 = float(matches[0][2].strip()) if len(matches[0]) > 2 else 0.0  # Handle absence of "Inter WP"

                extracted_values[filename] = (label1, label2, value1, value2, value3)

    return extracted_values


def populate_twobody_elprep_matrices(one_body_orcaout_filenames, two_body_orcaout_directory, conversion_factor, method, two_body_labels, LEDAW_output_path_two_body, use_ref_as_rhf_in_hfld=None, relabel_mapping=None):
    """Populate matrices for different intra-fragment properties and save them to an Excel file."""
    
    # Patterns class instance
    patterns = Patterns()
    
    # Extract values from one-body files using the potentially provided relabel_mapping
    one_body_values = extract_onebody_values(one_body_orcaout_filenames, patterns, method, use_ref_as_rhf_in_hfld, relabel_mapping)
    
    # Define the properties to extract and process
    properties = {
        'REF': 'intra_ref',
        'SP': 'intra_strong_pairs',
        'WP': 'intra_weak_pairs',
        'T': 'intra_triples',
        'Singles': 'singles_contribution'
    }
    
    # Store results in a dictionary
    results = {}

    for prop_name, pattern_key in properties.items():
        # Extract values from two-body files for the given property
        two_body_values = extract_twobody_values(two_body_orcaout_directory, two_body_labels, pattern_key, patterns)
        
        # Initialize the matrix
        n = len(one_body_orcaout_filenames)
        matrix = np.zeros((n, n))
        
        # Determine which one-body values to use
        if prop_name == 'REF':
            one_body_key = 'ref'
        elif prop_name == 'SP':
            one_body_key = 'strong_corr'
        elif prop_name == 'WP':
            one_body_key = 'weak_corr'
        elif prop_name == 'T':
            one_body_key = 'triples_corr'
        else:
            one_body_key = None  # For Singles, we consider the one-body value as zero

        # Process extracted two-body values
        for filename, (label1, label2, value1, value2, value3) in two_body_values.items():
            if one_body_key:
                one_body_value1 = one_body_values[one_body_key].get(label1, 0.0)
                one_body_value2 = one_body_values[one_body_key].get(label2, 0.0)
            else:
                one_body_value1 = one_body_value2 = 0.0  # Singles Contribution

            # Calculate adjusted values
            adjusted_value1 = (value1 - one_body_value1) * conversion_factor
            adjusted_value2 = (value2 - one_body_value2) * conversion_factor

            # Populate the matrix
            matrix[label1 - 1, label1 - 1] += adjusted_value1
            matrix[label2 - 1, label2 - 1] += adjusted_value2
            total_value = adjusted_value1 + adjusted_value2
            if label1 < label2:
                matrix[label1 - 1, label2 - 1] = total_value
            else:
                matrix[label2 - 1, label1 - 1] = total_value
        
        # Store the resulting matrix as a DataFrame
        df = pd.DataFrame(matrix, index=range(1, n + 1), columns=range(1, n + 1))
        results[prop_name] = df

    # Write results to Excel, each property to a separate sheet
    output_file_path = os.path.join(LEDAW_output_path_two_body, 'ELPREP.xlsx')
    with pd.ExcelWriter(output_file_path) as writer:
        for sheet_name, df in results.items():
            df.to_excel(writer, sheet_name=sheet_name)

    print(f"Two-body LED electronic preparation matrices were written to {output_file_path}")
    print("Standard (diagonal) and fp-LED (nondiagonal) electronic preparation terms are shown on the same matrix")


def extract_elprep_diagonals(elprep_file_path):
    """Extract diagonals from each sheet in the ELPREP.xlsx file."""
    elprep_sheets = pd.read_excel(elprep_file_path, sheet_name=None, index_col=0)
    diagonals = {}
    for sheet_name, df in elprep_sheets.items():
        diagonal = np.diag(df)
        diagonals[sheet_name] = diagonal
    return diagonals


def combine_diagonals(diagonals, method):
    """Combine diagonals according to the specified method."""
    
    # Calculate REF diagonals
    ref_diagonal = diagonals.get('REF', np.zeros_like(next(iter(diagonals.values()))))

    combined_diagonals = {}

    if method.lower() == 'dlpno-ccsd':
        c_ccsd_diagonal = (
            diagonals.get('SP', np.zeros_like(ref_diagonal)) +
            diagonals.get('WP', np.zeros_like(ref_diagonal)) +
            diagonals.get('Singles', np.zeros_like(ref_diagonal))
        )

        total_diagonal = c_ccsd_diagonal + ref_diagonal

        combined_diagonals['REF'] = ref_diagonal
        combined_diagonals['C-CCSD'] = c_ccsd_diagonal
        combined_diagonals['TOTAL'] = total_diagonal

    elif method.lower() == 'dlpno-ccsd(t)':
        c_ccsd_t_diagonal = (
            diagonals.get('SP', np.zeros_like(ref_diagonal)) +
            diagonals.get('WP', np.zeros_like(ref_diagonal)) +
            diagonals.get('T', np.zeros_like(ref_diagonal)) +
            diagonals.get('Singles', np.zeros_like(ref_diagonal))
        )

        total_ccsd_t_diagonal = c_ccsd_t_diagonal + ref_diagonal

        combined_diagonals['REF'] = ref_diagonal
        combined_diagonals['C-CCSD(T)'] = c_ccsd_t_diagonal
        combined_diagonals['TOTAL'] = total_ccsd_t_diagonal

    elif method.lower() == 'hfld':
        combined_diagonals['REF'] = ref_diagonal
        combined_diagonals['TOTAL'] = ref_diagonal

    else:
        raise ValueError("Unknown method provided.")

    return combined_diagonals


def set_diag_belowdiag_nan(df):
    """Set all diagonal and below-diagonal elements to NaN for Electrostat and Exchange matrices."""
    df = df.astype(float)  # Convert the DataFrame to float to allow NaN values
    df = df.copy()  # Make a copy of the DataFrame to avoid modifying it in place
    for i in range(df.shape[0]):
        for j in range(i+1):  # Ensure all elements on the diagonal and below it are set to NaN
            df.iat[i, j] = np.nan
    return df


def set_belowdiag_nan(df):
    """Set the values below the diagonal in the matrix to NaN."""
    df = df.astype(float)  # Convert the DataFrame to float to allow NaN values
    mask = np.tril(np.ones(df.shape), k=-1).astype(bool)
    df.values[mask] = np.nan
    return df


def set_diag_to_nan_if_zero(df):
    """Set diagonal elements to NaN if they are all zero."""
    df = df.astype(float)  # Convert the DataFrame to float to allow NaN values
    diagonal = np.diag(df)
    
    # Check if all diagonal elements are zero
    if np.all(diagonal == 0):
        np.fill_diagonal(df.values, np.nan)
    
    return df


def finalize_els_exch_matrices_for_writing(summary_sheets):
    """Finalize matrices by setting diagonals and below-diagonal elements to NaN where necessary."""
    for key in summary_sheets.keys():
        if key in ['Electrostat', 'Exchange']:
            summary_sheets[key] = set_diag_belowdiag_nan(summary_sheets[key])
        else:
            summary_sheets[key] = set_belowdiag_nan(summary_sheets[key])
            summary_sheets[key] = set_diag_to_nan_if_zero(summary_sheets[key])
    return summary_sheets


def calculate_twobody_standard_LED_summary_matrices(LEDAW_output_path_two_body, method, use_ref_as_rhf_in_hfld):
    # File paths
    elprep_file = os.path.join(LEDAW_output_path_two_body, 'ELPREP.xlsx')
    inter_file = os.path.join(LEDAW_output_path_two_body, 'INTER.xlsx')
    summary_file = os.path.join(LEDAW_output_path_two_body, 'Summary_Standard_LED_matrices.xlsx')
    
    # Load sheets from INTER.xlsx
    inter_sheets = pd.read_excel(inter_file, sheet_name=None, index_col=0)
    
    # Extract diagonals from ELPREP.xlsx
    diagonals = extract_elprep_diagonals(elprep_file)

    # Combine the diagonals
    combined_diagonals = combine_diagonals(diagonals, method)
    
    # Create a dictionary to hold summary results
    summary_sheets = {}

    # Transfer Electrostat and Exchange sheets
    summary_sheets['Electrostat'] = inter_sheets['Electrostat']
    summary_sheets['Exchange'] = inter_sheets['Exchange']

    # Calculate REF sheet (Electrostat + Exchange) and then add the diagonal of Intra REF
    ref_matrix = inter_sheets['Electrostat'] + inter_sheets['Exchange']
    
    # Add diagonal corrections (only once)
    np.fill_diagonal(ref_matrix.values, np.diag(ref_matrix) + combined_diagonals['REF'])
    summary_sheets['REF'] = ref_matrix

    # Calculate Dispersion sheets
    if method.lower() == 'dlpno-ccsd':
        disp_ccsd = inter_sheets['Disp SP'] + inter_sheets['Inter WP']
        summary_sheets['Disp CCSD'] = set_diag_belowdiag_nan(disp_ccsd)  # Set diag and lower diag to NaN

    elif method.lower() == 'hfld':
        disp_hfld = inter_sheets['Disp SP'] + inter_sheets['Inter WP']
        summary_sheets['Disp HFLD'] = set_diag_belowdiag_nan(disp_hfld)

    if method.lower() == 'dlpno-ccsd(t)':
        with np.errstate(divide='ignore', invalid='ignore'):
            disp_t = np.divide(inter_sheets['Disp SP'] * inter_sheets['Inter T'], inter_sheets['Inter SP'])
            disp_t[np.isnan(disp_t)] = 0
        
        disp_t = disp_t.where(np.triu(np.ones(disp_t.shape), k=0).astype(bool))
        inter_sheets['Disp T'] = disp_t  # Prepare to append Disp T to INTER.xlsx

        disp_ccsd_t = inter_sheets['Disp SP'] + inter_sheets['Inter WP'] + disp_t
        summary_sheets['Disp CCSD(T)'] = set_diag_belowdiag_nan(disp_ccsd_t)  # Set diag and lower diag to NaN

    # Calculate Inter-NonDisp-C matrices
    if method.lower() == 'dlpno-ccsd':
        inter_nondisp_ccsd = inter_sheets['Inter SP'] - inter_sheets['Disp SP']
        summary_sheets['Inter-NonDisp-C-CCSD'] = set_diag_belowdiag_nan(inter_nondisp_ccsd)  # Set diag and lower diag to NaN

    elif method.lower() == 'dlpno-ccsd(t)':
        inter_nondisp_ccsd_t = (inter_sheets['Inter SP'] + inter_sheets['Inter WP'] + 
                                inter_sheets['Inter T'] - summary_sheets['Disp CCSD(T)'])
        summary_sheets['Inter-NonDisp-C-CCSD(T)'] = set_diag_belowdiag_nan(inter_nondisp_ccsd_t)  # Set diag and lower diag to NaN

    # Correct diagonal elements for C-CCSD and related matrices
    if method.lower() == 'dlpno-ccsd':
        c_ccsd_matrix = inter_sheets['Inter SP'] + inter_sheets['Inter WP']
        np.fill_diagonal(c_ccsd_matrix.values, np.diag(c_ccsd_matrix) + combined_diagonals['C-CCSD'])
        summary_sheets['C-CCSD'] = c_ccsd_matrix

        total_ccsd = summary_sheets['REF'] + summary_sheets['C-CCSD']
        summary_sheets['TOTAL'] = total_ccsd

    elif method.lower() == 'dlpno-ccsd(t)':
        c_ccsd_t_matrix = inter_sheets['Inter SP'] + inter_sheets['Inter WP'] + inter_sheets['Inter T']
        np.fill_diagonal(c_ccsd_t_matrix.values, np.diag(c_ccsd_t_matrix) + combined_diagonals['C-CCSD(T)'])
        summary_sheets['C-CCSD(T)'] = c_ccsd_t_matrix

        total_ccsd_t = summary_sheets['REF'] + summary_sheets['C-CCSD(T)']
        summary_sheets['TOTAL'] = total_ccsd_t

    elif method.lower() == 'hfld':
        summary_sheets['TOTAL'] = summary_sheets['REF'] + summary_sheets['Disp HFLD']

    # Adjust labels of Intra sheets by adding 1
    adjusted_intra_sheets = {}
    for name, diagonal in combined_diagonals.items():
        adjusted_df = pd.DataFrame(np.diag(diagonal), index=range(1, len(diagonal) + 1), columns=range(1, len(diagonal) + 1))
        adjusted_intra_sheets[f"Intra {name}"] = adjusted_df

    # Sum the Intra sheets with their corresponding non-Intra sheets
    for name in adjusted_intra_sheets.keys():
        base_name = name.replace('Intra ', '')
        if base_name in summary_sheets:
            summary_sheets[base_name] = summary_sheets[base_name].add(adjusted_intra_sheets[name], fill_value=0)

    # Finalize matrices before writing to Excel
    summary_sheets = finalize_els_exch_matrices_for_writing(summary_sheets)

    # Save or overwrite Disp T to INTER.xlsx
    if 'Disp T' in inter_sheets:
        with pd.ExcelWriter(inter_file, mode='a', engine='openpyxl') as writer:
            workbook = writer.book
            if 'Disp T' in workbook.sheetnames:
                del workbook['Disp T']  # Remove existing 'Disp T' sheet
            inter_sheets['Disp T'].to_excel(writer, sheet_name='Disp T')

    # Write to Summary Excel file in the specified order
    sheet_order = [
        'TOTAL', 'REF', 'Electrostat', 'Exchange',
        'C-CCSD' if method.lower() == 'dlpno-ccsd' else 'C-CCSD(T)',
        'Disp CCSD' if method.lower() == 'dlpno-ccsd' else 'Disp CCSD(T)' if method.lower() == 'dlpno-ccsd(t)' else 'Disp HFLD',
        'Inter-NonDisp-C-CCSD' if method.lower() == 'dlpno-ccsd' else 'Inter-NonDisp-C-CCSD(T)'
    ]

    with pd.ExcelWriter(summary_file, engine='openpyxl') as writer:
        for sheet_name in sheet_order:
            if sheet_name in summary_sheets:
                summary_sheets[sheet_name].to_excel(writer, sheet_name=sheet_name)

        print(f"Standard LED two-body summary interaction energy matrices were written to '{summary_file}'")


def extract_upper_diagonal(ref_sheet):
    """Extract the upper-diagonal elements"""
    upper_diag = np.triu(ref_sheet, k=1)  # Extract upper diagonal elements (excluding diagonal)
    return pd.DataFrame(upper_diag, index=ref_sheet.index, columns=ref_sheet.columns)


def sum_corr_elprep_upperdiagonals(elprep_file_path, include_t=False):
    """Sum the upper-diagonal elements of correlation components of ELPREP."""
    
    # Load the specified sheets from ELPREP.xlsx
    sheet_names = ['SP', 'WP', 'Singles']
    if include_t:
        sheet_names.append('T')
        
    elprep_sheets = pd.read_excel(elprep_file_path, sheet_name=sheet_names, index_col=0)
    
    # Initialize the sum with SP sheet's upper diagonal elements
    summed_matrix = elprep_sheets['SP'].where(np.triu(np.ones(elprep_sheets['SP'].shape), k=1).astype(bool), 0)

    # Add WP sheet's upper diagonal elements
    if 'WP' in elprep_sheets:
        wp_upper_diag = elprep_sheets['WP'].where(np.triu(np.ones(elprep_sheets['WP'].shape), k=1).astype(bool), 0)
        summed_matrix += wp_upper_diag

    # Add Singles sheet's upper diagonal elements
    if 'Singles' in elprep_sheets:
        singles_upper_diag = elprep_sheets['Singles'].where(np.triu(np.ones(elprep_sheets['Singles'].shape), k=1).astype(bool), 0)
        summed_matrix += singles_upper_diag

    # Add T sheet's upper diagonal elements if required
    if include_t and 'T' in elprep_sheets:
        t_upper_diag = elprep_sheets['T'].where(np.triu(np.ones(elprep_sheets['T'].shape), k=1).astype(bool), 0)
        summed_matrix += t_upper_diag

    return summed_matrix


def calculate_twobody_fpLED_matrices(LEDAW_output_path_two_body, method):
    # Ensure the output directory exists
    if not os.path.exists(LEDAW_output_path_two_body):
        os.makedirs(LEDAW_output_path_two_body)

    # File paths
    elprep_file = os.path.join(LEDAW_output_path_two_body, 'ELPREP.xlsx')
    summary_standard_file = os.path.join(LEDAW_output_path_two_body, 'Summary_Standard_LED_matrices.xlsx')
    diel_file = os.path.join(LEDAW_output_path_two_body, 'DIEL.xlsx')
    summary_file = os.path.join(LEDAW_output_path_two_body, 'Summary_fp-LED_matrices.xlsx')
    
    # Load standard sheets from Summary_Standard_LED_matrices.xlsx
    standard_sheets = pd.read_excel(summary_standard_file, sheet_name=None, index_col=0)
    
    # Initialize the summary sheets dictionary
    summary_sheets = {}

    # Transfer Electrostat and Exchange sheets from the standard sheets
    summary_sheets['Electrostat'] = standard_sheets['Electrostat']
    summary_sheets['Exchange'] = standard_sheets['Exchange']

    # Load the REF sheet from ELPREP.xlsx
    elprep_sheets = pd.read_excel(elprep_file, sheet_name=['REF'], index_col=0)
    ref_sheet = elprep_sheets['REF']
    
    # Calculate REF-EL-PREP as the upper diagonal of the REF sheet
    summary_sheets['REF-EL-PREP'] = extract_upper_diagonal(ref_sheet)
    
    # Calculate REF by summing Electrostat, Exchange, and REF-EL-PREP
    summary_sheets['REF'] = summary_sheets['Electrostat'] + summary_sheets['Exchange'] + summary_sheets['REF-EL-PREP']
    
    # Calculate C-CCSD-EL-PREP or C-CCSD(T)-EL-PREP depending on the method
    include_t = 't' in method.lower()
    c_ccsd_elprep = sum_corr_elprep_upperdiagonals(elprep_file, include_t=include_t)
    elprep_key = 'C-CCSD-EL-PREP' if not include_t else 'C-CCSD(T)-EL-PREP'
    
    # Assign the summed matrix to the appropriate key
    summary_sheets[elprep_key] = c_ccsd_elprep
    
    # Transfer Disp and Inter-NonDisp sheets from the standard sheets
    if method.lower() == 'dlpno-ccsd':
        summary_sheets['Disp CCSD'] = standard_sheets['Disp CCSD']
        summary_sheets['Inter-NonDisp-C-CCSD'] = standard_sheets['Inter-NonDisp-C-CCSD']
        # Calculate C-CCSD
        summary_sheets['C-CCSD'] = summary_sheets['Disp CCSD'] + summary_sheets['Inter-NonDisp-C-CCSD'] + summary_sheets['C-CCSD-EL-PREP']
    elif method.lower() == 'dlpno-ccsd(t)':
        summary_sheets['Disp CCSD(T)'] = standard_sheets['Disp CCSD(T)']
        summary_sheets['Inter-NonDisp-C-CCSD(T)'] = standard_sheets['Inter-NonDisp-C-CCSD(T)']
        # Calculate C-CCSD(T)
        summary_sheets['C-CCSD(T)'] = summary_sheets['Disp CCSD(T)'] + summary_sheets['Inter-NonDisp-C-CCSD(T)'] + summary_sheets['C-CCSD(T)-EL-PREP']
    elif method.lower() == 'hfld':
        summary_sheets['Disp HFLD'] = standard_sheets['Disp HFLD']

    # Calculate CCSD-EL-PREP or CCSD(T)-EL-PREP as the sum of REF-EL-PREP and C-CCSD-EL-PREP (or C-CCSD(T)-EL-PREP)
    if method.lower() == 'dlpno-ccsd':
        summary_sheets['CCSD-EL-PREP'] = summary_sheets['REF-EL-PREP'] + summary_sheets['C-CCSD-EL-PREP']
    elif method.lower() == 'dlpno-ccsd(t)':
        summary_sheets['CCSD(T)-EL-PREP'] = summary_sheets['REF-EL-PREP'] + summary_sheets['C-CCSD(T)-EL-PREP']

    # Transfer DIEL sheet if it is not zero
    if os.path.exists(diel_file):
        diel_sheet = pd.read_excel(diel_file, sheet_name=None, index_col=0)
        diel_df = diel_sheet['DIEL']

        # Check if the DIEL sheet is not entirely zero
        if not np.all(diel_df.values == 0):
            # Set diagonal and lower-diagonal elements to NaN
            mask = np.tril(np.ones(diel_df.shape), k=0).astype(bool)
            diel_df = diel_df.mask(mask)
            summary_sheets['DIEL'] = diel_df

    # Calculate TOTAL
    if method.lower() == 'hfld':
        summary_sheets['TOTAL'] = summary_sheets['REF'] + summary_sheets['Disp HFLD']
        if 'DIEL' in summary_sheets:
            summary_sheets['TOTAL'] += summary_sheets['DIEL']
    elif method.lower() == 'dlpno-ccsd':
        summary_sheets['TOTAL'] = summary_sheets['REF'] + summary_sheets['C-CCSD']
        if 'DIEL' in summary_sheets:
            summary_sheets['TOTAL'] += summary_sheets['DIEL']
    elif method.lower() == 'dlpno-ccsd(t)':
        summary_sheets['TOTAL'] = summary_sheets['REF'] + summary_sheets['C-CCSD(T)']
        if 'DIEL' in summary_sheets:
            summary_sheets['TOTAL'] += summary_sheets['DIEL']

    # Set diagonal and below-diagonal elements to NaN for all matrices
    for sheet_name in summary_sheets:
        df = summary_sheets[sheet_name]
        mask = np.tril(np.ones(df.shape), k=0).astype(bool)
        df.values[mask] = np.nan
        summary_sheets[sheet_name] = df

    # Ensure that the file is created and write the results to the Excel file
    with pd.ExcelWriter(summary_file, engine='openpyxl') as writer:
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
                summary_sheets[sheet_name].to_excel(writer, sheet_name=sheet_name)

    print(f"fp-LED two-body summary interaction energy matrices were written to '{summary_file}'")


def engine_LED_two_body(one_body_orcaout_filenames, two_body_orcaout_directory, conversion_factor, method, LEDAW_output_path_two_body, use_ref_as_rhf_in_hfld=None, relabel_mapping=None):
    ''' Collect and process the files needed for two-body LED and provide standard and fp-LED matrices'''

    # Ensure the output directory exists
    if not os.path.exists(LEDAW_output_path_two_body):
        os.makedirs(LEDAW_output_path_two_body)

    # Get two-body filenames
    two_body_orcaout_filenames = get_two_body_filenames(two_body_orcaout_directory)
    print(f"one_body_orcaout_filenames: {one_body_orcaout_filenames}")

    # Extract labels from one-body files
    label_mapping = extract_labels_from_one_body_files(one_body_orcaout_filenames, relabel_mapping=relabel_mapping)
    
    # Extract labels from two-body files based on one-body labels
    two_body_labels = extract_labels_from_two_body_files(two_body_orcaout_directory, label_mapping)
    
    # Populate dielectric matrices
    populate_twobody_dielectric_matrices(one_body_orcaout_filenames, two_body_orcaout_directory,
                                         conversion_factor, two_body_labels, 
                                         LEDAW_output_path_two_body=LEDAW_output_path_two_body, 
                                         relabel_mapping=relabel_mapping)
    
    # Populate inter-fragment interaction matrices
    populate_twobody_inter_matrices(one_body_orcaout_filenames, two_body_orcaout_directory,
                                    conversion_factor, relabel_mapping=relabel_mapping, 
                                    two_body_labels=two_body_labels, 
                                    LEDAW_output_path_two_body=LEDAW_output_path_two_body)
    
    # Populate energy preparation matrices
    populate_twobody_elprep_matrices(one_body_orcaout_filenames, two_body_orcaout_directory,
                                     conversion_factor, two_body_labels=two_body_labels, 
                                     LEDAW_output_path_two_body=LEDAW_output_path_two_body, 
                                     method=method, relabel_mapping=relabel_mapping,
									 use_ref_as_rhf_in_hfld=use_ref_as_rhf_in_hfld
									 )
    
    # Calculate and store standard LED matrices
    calculate_twobody_standard_LED_summary_matrices(LEDAW_output_path_two_body=LEDAW_output_path_two_body,
                                                    method=method, use_ref_as_rhf_in_hfld=use_ref_as_rhf_in_hfld)

    
    # Calculate and store fp-LED matrices
    calculate_twobody_fpLED_matrices(LEDAW_output_path_two_body=LEDAW_output_path_two_body, method=method)
    
    # Print completion message
    print('\n')
    print('*' * 125)
    print(f"  Two-body LED analyses were terminated NORMALLY. Standard and fp-LED two-body summary matrices are at {LEDAW_output_path_two_body}")
    print('*' * 125)

	
