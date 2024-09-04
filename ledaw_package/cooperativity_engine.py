import os
import pandas as pd


def find_nbody_twobody_subdirectories(base_path, nbody_dir_name, twobody_dir_name, directory_level=1):
    """Find all directories under the given base path containing both NBODY and TWOBODY subdirectories."""
    directories = []
    if directory_level == 1:
        # When directory_level is 1, look directly under base_path
        nbody_path = os.path.join(base_path, nbody_dir_name)
        twobody_path = os.path.join(base_path, twobody_dir_name)
        if os.path.exists(nbody_path) and os.path.exists(twobody_path):
            directories.append(base_path)
    else:
        # When directory_level is > 1, look under subdirectories
        for dir_name in os.listdir(base_path):
            full_dir_path = os.path.join(base_path, dir_name)
            if os.path.isdir(full_dir_path):
                nbody_path = os.path.join(full_dir_path, nbody_dir_name)
                twobody_path = os.path.join(full_dir_path, twobody_dir_name)
                if os.path.exists(nbody_path) and os.path.exists(twobody_path):
                    directories.append(full_dir_path)

    print(f"Found {len(directories)} directories with {nbody_dir_name} and {twobody_dir_name} subdirectories.")
    return directories


def calculate_cooperativity_matrices(nbody_file, twobody_file, output_file):
    """Calculate NBODY minus TWOBODY matrices and write them to the output Excel file."""
    print(f"Processing files: {nbody_file} and {twobody_file}")
    nbody_excel = pd.ExcelFile(nbody_file)
    twobody_excel = pd.ExcelFile(twobody_file)
    
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for sheet_name in nbody_excel.sheet_names:
            if sheet_name in twobody_excel.sheet_names:
                nbody_df = pd.read_excel(nbody_excel, sheet_name=sheet_name, index_col=0)
                twobody_df = pd.read_excel(twobody_excel, sheet_name=sheet_name, index_col=0)
                
                # Calculate NBODY minus TWOBODY matrix
                coop_df = nbody_df - twobody_df
                
                # Write the result to the new file
                coop_df.to_excel(writer, sheet_name=sheet_name)


def cooperativity_engine(base_path, nbody_dir_name='NBODY', twobody_dir_name='TWOBODY', directory_level=1):
    """Create cooperativity Excel files for all directories"""
    # Pass nbody_dir_name, twobody_dir_name, and directory_level to the find_nbody_twobody_subdirectories function
    directories = find_nbody_twobody_subdirectories(base_path, nbody_dir_name, twobody_dir_name, directory_level)
    
    for directory in directories:
        if directory_level == 1:
            nbody_dir = os.path.join(directory, nbody_dir_name)
            twobody_dir = os.path.join(directory, twobody_dir_name)
        else:
            nbody_dir = os.path.join(directory, nbody_dir_name)
            twobody_dir = os.path.join(directory, twobody_dir_name)
        
        for nbody_file in os.listdir(nbody_dir):
            if nbody_file.endswith('.xlsx'):
                nbody_filepath = os.path.join(nbody_dir, nbody_file)
                twobody_filepath = os.path.join(twobody_dir, nbody_file)
                
                if os.path.exists(twobody_filepath):
                    # Create corresponding output directory in COOPERATIVITY
                    output_dir = os.path.join(directory, 'COOPERATIVITY')
                    os.makedirs(output_dir, exist_ok=True)
                    
                    # Create output file path
                    output_filepath = os.path.join(output_dir, nbody_file)
                    
                    # Process the files
                    calculate_cooperativity_matrices(nbody_filepath, twobody_filepath, output_filepath)
                    print(f"Processed and saved to: {output_dir}")
                else:
                    print(f"TWOBODY file not found for: {nbody_file}")
    
    print('\n')
    print('*'*100)
    print("Cooperativity job was terminated NORMALLY")
    print('*'*100)
