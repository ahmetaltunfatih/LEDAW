import os
import re


def extract_fragment_coords(supersystem_file, tolerance=1e-3):
    """Extracts fragment coordinates from the supersystem file."""
    label_coord_dict = {}
    with open(supersystem_file, 'r') as file:
        content = file.read()
        
        match = re.search(r'CARTESIAN COORDINATES OF FRAGMENTS \(ANGSTROEM\)(.*?)INTERNAL COORDINATES \(ANGSTROEM\)', content, re.S)
        if match:
            fragment_section = match.group(1)
            
            fragment_matches = re.findall(r'FRAGMENT (\d+)(.*?)\n\s*([A-Za-z]+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)', fragment_section, re.S)
            
            for frag_num, _, _, x, y, z in fragment_matches:
                frag_num = int(frag_num)
                coords = (float(x), float(y), float(z))
                label_coord_dict[frag_num] = coords
    
    return label_coord_dict


def extract_coordinates_from_onebody_file(file_path):
    """Extracts coordinates from one-body files, skips lines with colons for BSSE calculations, and stops at 'END OF INPUT'."""
    
    coords = []
    
    with open(file_path, 'r') as file:
        content = file.readlines()

    pattern_found = False

    for i, line in enumerate(content):
        # Check for the *xyz pattern
        if re.search(r'\*\s*xyz', line):
            pattern_found = True

            # Process lines after the *xyz pattern to find the first one without a colon
            for subsequent_line in content[i+1:]:
                # Stop if END OF INPUT is found
                if 'END OF INPUT' in subsequent_line:
                    break
                
                # Skip lines with a colon
                if ':' in subsequent_line:
                    continue  # to skip this line
                
                # Extract coordinates from lines without a colon
                match = re.findall(r'-?\d+\.\d+', subsequent_line)
                if match:
                    # Convert to float and format with 6 digits after the decimal point
                    coord_set = tuple(float(f"{float(num):.6f}") for num in match[:3])  # Extract the first 3 values
                    coords.append(coord_set)
            
            # Break the outer loop after processing the first *xyz block
            break

    return coords


def match_fragments_with_onebody_files(label_coord_dict, onebody_out_directory, tolerance=1e-3):
    """Matches fragment coordinates with the coordinates in the one-body files."""
    one_body_orcaout_filenames = [None] * len(label_coord_dict)
    
    for frag_num, frag_coords in label_coord_dict.items():
        for filename in os.listdir(onebody_out_directory):
            file_path = os.path.join(onebody_out_directory, filename)
            if not os.path.isfile(file_path):
                continue

            mono_coords_list = extract_coordinates_from_onebody_file(file_path)
            
            for coord_mono in mono_coords_list:
                if all(abs(c1 - c2) <= tolerance for c1, c2 in zip(coord_mono, frag_coords)):
                    one_body_orcaout_filenames[frag_num - 1] = file_path
                    break
            if one_body_orcaout_filenames[frag_num - 1] is not None:
                break
    
    if None in one_body_orcaout_filenames:
        print("Some fragments could not be matched. Here's the list of matched files:")
        print(one_body_orcaout_filenames)
    
    return one_body_orcaout_filenames


def extract_one_body_orcaout_filenames(supersystem_file, onebody_out_directory):
    """Wrapper function to extract and match fragment coordinates with one-body files."""
    label_coord_dict = extract_fragment_coords(supersystem_file, tolerance=1e-3)
    one_body_orcaout_filenames = match_fragments_with_onebody_files(label_coord_dict, onebody_out_directory, tolerance=1e-3)
    return one_body_orcaout_filenames
