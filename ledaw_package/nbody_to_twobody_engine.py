import os
import re
from .nbody_engine import normalize_path, extract_coords_from_line


def extract_labels_coords_from_supersys(supersystem_file, tolerance=1e-3):
    """Extract fragment labels and coordinates from the FRAGMENT X section of the supersystem file."""

    label_coord_dict = {}
    current_label = None
    current_coords = []
    recording = False

    with open(supersystem_file, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if re.search(r"CARTESIAN COORDINATES OF FRAGMENTS\s*\(ANGSTROEM\)", line):
            recording = True
            continue
        if re.search(r"INTERNAL COORDINATES\s*\(ANGSTROEM\)", line):
            recording = False
            continue

        frag_match = re.match(r"\s*FRAGMENT\s+(\d+)", line)
        if frag_match:
            if current_label is not None and current_coords:
                label_coord_dict[current_label] = current_coords
            current_label = int(frag_match.group(1))
            current_coords = []
            continue

        if recording and current_label is not None:
            coords = extract_coords_from_line(line)
            if coords:
                current_coords.append(tuple(round(x, 6) for x in coords))

    # Add last fragment
    if current_label is not None and current_coords:
        label_coord_dict[current_label] = current_coords

    return label_coord_dict


def get_reduced_relabel_mapping(supersystem_file, relabel_mapping):
    """Extracts fragment labels from a supersystem file, filters out labels that do not correspond to actual fragments,
    and reduces the relabel_mapping accordingly. The reduced mapping is then rewritten in terms of consecutive labels 
    ranging from 1 to the number of actual fragments in the supersystem."""

    label_coord_dict = extract_labels_coords_from_supersys(supersystem_file)
    original_labels = sorted(label_coord_dict.keys())

    if relabel_mapping is None:
        max_label = max(original_labels)
        relabel_mapping = list(range(1, max_label + 1))

    # Step 1: Get reduced list of remapped labels
    reduced_relabel_values = [relabel_mapping[orig - 1] for orig in original_labels]

    # Step 2: Sort and map to 1, 2, 3...
    sorted_values = sorted(reduced_relabel_values)

    value_to_new_label = {val: idx + 1 for idx, val in enumerate(sorted_values)}

    # Step 3: Remap reduced values to contiguous new labels
    reduced_relabel_mapping = [value_to_new_label[val] for val in reduced_relabel_values]

    return reduced_relabel_mapping


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
    """Matches fragment coordinates with the coordinates in the one-body files.
    Keeps the first matched file per fragment in the sequential order of supersystem file, then appends unmatched files at the end."""

    normalized_onebody_out_directory = normalize_path(onebody_out_directory)
    all_files = [os.path.join(normalized_onebody_out_directory, f)
                 for f in os.listdir(normalized_onebody_out_directory)
                 if os.path.isfile(os.path.join(normalized_onebody_out_directory, f))]

    matched_files_dict = {}
    used_files = []

    for frag_num, frag_coords in sorted(label_coord_dict.items()):
        for file_path in all_files:
            if file_path in used_files:
                continue

            mono_coords_list = extract_coordinates_from_onebody_file(file_path)
            for coord_mono in mono_coords_list:
                for frag_coord in frag_coords:
                    if all(abs(c1 - c2) <= tolerance for c1, c2 in zip(coord_mono, frag_coord)):
                        matched_files_dict[frag_num] = file_path
                        used_files.append(file_path)
                        break
                if frag_num in matched_files_dict:
                    break
            if frag_num in matched_files_dict:
                break

    unmatched_files = [f for f in all_files if f not in used_files]
    one_body_orcaout_filenames = [matched_files_dict[k] for k in sorted(matched_files_dict.keys())] + unmatched_files

    if len(matched_files_dict) < len(label_coord_dict):
        print("Some fragments could not be matched. Here's the partial match list:")
        print([matched_files_dict.get(k) for k in sorted(label_coord_dict.keys())])

    return one_body_orcaout_filenames


def extract_one_body_orcaout_filenames(supersystem_file, onebody_out_directory):
    """Wrapper function to extract and match fragment coordinates with one-body files."""
    label_coord_dict = extract_labels_coords_from_supersys(supersystem_file, tolerance=1e-3)
    one_body_orcaout_filenames = match_fragments_with_onebody_files(label_coord_dict, onebody_out_directory, tolerance=1e-3)
    return one_body_orcaout_filenames
