import os
import re
import gc
import copy
import glob
import time
import shutil
import psutil
import platform
import openpyxl
import numpy as np
import pandas as pd
from openpyxl import load_workbook, Workbook
from .classes import *


patterns = Patterns()


def normalize_path(path):
    return path.replace("\\", "/").rstrip("/")


def label_systems(filenames):
    """Assign SUPERSYS to the first file and SUBSYS1, SUBSYS2, ... to the rest."""
    labeled_filenames = {}
    for i, filename in enumerate(filenames):
        if filename:  # Ensure filename is not None or empty
            normalized_filename = normalize_path(filename)  # Normalize the path
            if i == 0:
                labeled_filenames[normalized_filename] = "SUPERSYS"
            else:
                labeled_filenames[normalized_filename] = f"SUBSYS{i}"
    return labeled_filenames


def extract_coords_from_line(line):
    parts = line.strip().split()
    try:
        return [float(parts[-3]), float(parts[-2]), float(parts[-1])]
    except (ValueError, IndexError):
        return None


def extract_real_coords_xyz_block(file_path, tag, debug_lines):
    coords = []
    in_xyz_block = False
    with open(file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if re.search(r"\*\s*xyz\b", line, re.IGNORECASE):
            in_xyz_block = True
            debug_lines.append(f"[{tag}] Found *xyz block at line {i}: {line.strip()}")
            continue
        if in_xyz_block:
            if re.match(r"\s*\*\s*$", line) or "END OF INPUT" in line:
                debug_lines.append(f"[{tag}] End of *xyz block at line {i}: {line.strip()}")
                break
            if ":" in line:
                debug_lines.append(f"[{tag}] Skipping ghost line {i}: {line.strip()}")
                continue
            coord = extract_coords_from_line(line)
            if coord:
                coords.append(coord)
                debug_lines.append(f"[{tag}] Extracted coord at line {i}: {coord}")
            else:
                debug_lines.append(f"[{tag}] Failed to extract coord at line {i}: {line.strip()}")
    return coords


def fragments_equal(frag1, frag2, tol=1e-3):
    """Compare two lists of coordinates to determine if they match within a tolerance. Order of atoms in different files does not matter."""

    # Check for the length of fragments
    if len(frag1) != len(frag2):
        return False

    # Sort both lists of atom coordinates to make the comparison order-independent.
    sorted_frag1 = sorted(frag1)
    sorted_frag2 = sorted(frag2)

    # Compare the sorted lists atom by atom.
    for atom1, atom2 in zip(sorted_frag1, sorted_frag2):
        # 'any' allows stopping on the first mismatch.
        if any(abs(coord_a - coord_b) > tol for coord_a, coord_b in zip(atom1, atom2)):
            return False # Mismatch found, fragments are not equal.

    # If the loop completes, all atoms match.
    return True


def extract_fragments(file_path, supersystem_coords=None, tol=1e-3):
    fragments = {}
    ghost_flags = {}
    bsse_found = False
    current_fragment = None
    current_atoms = []

    with open(file_path, 'r') as f:
        lines = f.readlines()

    fragment_mode = False
    fragment_found = False

    # Loop to detect and parse FRAGMENT blocks
    for line in lines:
        if re.search(r"CARTESIAN COORDINATES OF FRAGMENTS\s*\(ANGSTROEM\)", line):
            fragment_mode = True
            continue
        if re.search(r"INTERNAL COORDINATES\s*\(ANGSTROEM\)", line):
            fragment_mode = False
            continue

        frag_match = re.match(r"\s*FRAGMENT\s+(\d+)", line)
        if frag_match:
            if current_fragment is not None and current_atoms:
                fragments[current_fragment] = current_atoms
                ghost_flags[current_fragment] = False
            current_fragment = int(frag_match.group(1))
            current_atoms = []
            fragment_found = True
            continue

        if fragment_mode and current_fragment is not None:
            coord = extract_coords_from_line(line)
            if coord:
                current_atoms.append(tuple(round(x, 6) for x in coord))

    # If no FRAGMENT blocks found, treat the system as one fragment and extract coordinates from *xyz block 
    if not fragment_found:
        current_fragment = 1
        current_atoms = []
        for line in lines:
            if re.search(r"\*\s*xyz\b", line, re.IGNORECASE):
                for j in range(lines.index(line) + 1, len(lines)):
                    if re.match(r"\*", lines[j]) or "END OF INPUT" in lines[j]:
                        break
                    if ":" in lines[j]:
                        continue
                    coord = extract_coords_from_line(lines[j])
                    if coord:
                        current_atoms.append(tuple(round(x, 6) for x in coord)) 
                break
        if current_atoms:
            fragments[current_fragment] = current_atoms
            ghost_flags[current_fragment] = False

    if current_fragment is not None and current_atoms:
        fragments[current_fragment] = current_atoms
        ghost_flags[current_fragment] = False

    # If FRAGMENT and *xyz blocks do not exist, treat the system as one fragment and extract coordinates from "CARTESIAN COORDINATES (ANGSTROEM)" block.
    if not fragments:
        current_fragment = 1
        current_atoms = []
        in_cart_block = False
        for i, line in enumerate(lines):
            if re.search(r"CARTESIAN COORDINATES\s*\(ANGSTROEM\)", line, re.IGNORECASE):
                in_cart_block = True
                continue
            if in_cart_block:
                if re.match(r"\s*-{5,}", line):  # skip divider lines like "-----"
                    continue
                if line.strip() == "":
                    break  # empty line signals end of block
                coord = extract_coords_from_line(line)
                if coord:
                    current_atoms.append(tuple(round(x, 6) for x in coord))
        if current_atoms:
            fragments[current_fragment] = current_atoms
            ghost_flags[current_fragment] = False

    # Check *xyz block for ghost atoms
    for i, line in enumerate(lines):
        if re.search(r"\*\s*xyz\b", line, re.IGNORECASE):
            for j in range(i + 1, len(lines)):
                if re.match(r"\*", lines[j]) or "END OF INPUT" in lines[j]:
                    break
                if ':' not in lines[j]:
                    continue
                coord = extract_coords_from_line(lines[j])
                if coord is None:
                    continue
                coord = tuple(round(x, 6) for x in coord)

                bsse_found = True
                matched = False
                for frag_id, atom_list in fragments.items():
                    for atom in atom_list:
                        if np.allclose(coord, atom, atol=tol):
                            ghost_flags[frag_id] = True
                            matched = True
                            break
                    if matched:
                        break
            break

    return fragments, ghost_flags, bsse_found


def construct_label_mappings(main_filenames, alternative_filenames, LEDAW_output_path):
    if len(main_filenames) != len(alternative_filenames):
        raise ValueError("main_filenames and alternative_filenames must have the same length")

    supersystem_file = normalize_path(main_filenames[0])
    subsystem_files = [normalize_path(f) for f in main_filenames[1:]]
    alternative_files = [normalize_path(f) if f else "" for f in alternative_filenames[1:]]

    # === STEP 0: Check for BSSE one-fragment-per-subsystem case ===
    bsse_found_temp = False
    is_all_single_real = True

    for subsystem_file in subsystem_files:
        _, ghost_flags, bsse = extract_fragments(subsystem_file)
        bsse_found_temp |= bsse
        real_frags = [k for k, v in ghost_flags.items() if not v]
        if len(real_frags) != 1:
            is_all_single_real = False
            break

    if bsse_found_temp and is_all_single_real:
        main_label_mappings_ghost_free, alt_label_mappings_ghost_free = construct_label_mappings_singlefrag_bsse(
            main_filenames, alternative_filenames)

        print("======== Fragment Mappings ========")
        print("main_label_mappings_ghost_free:", main_label_mappings_ghost_free)
        print("alt_label_mappings_ghost_free:", alt_label_mappings_ghost_free)
        print("bsse_found:", True)

        return main_label_mappings_ghost_free, alt_label_mappings_ghost_free, True

    # === ELSE: Continue with multi-fragment logic ===
    supersystem_frags, super_ghost_flags, _ = extract_fragments(supersystem_file)

    main_label_mappings = {"SUPERSYS": {k: k for k in supersystem_frags}}
    alt_label_mappings = {"SUPERSYS": {}}
    main_label_mappings_ghost_free = {"SUPERSYS": {k: k for k, v in super_ghost_flags.items() if not v}}
    alt_label_mappings_ghost_free = {"SUPERSYS": {}}
    bsse_found = False

    # === MAIN SUBSYSTEMS ===
    for i, subsystem_file in enumerate(subsystem_files):
        subsys_label = f"SUBSYS{i+1}"
        subsystem_frags, sub_ghost_flags, bsse = extract_fragments(subsystem_file)
        bsse_found |= bsse

        forward_mapping = {}
        forward_mapping_ghost_free = {}

        real_fragments = [k for k, v in sub_ghost_flags.items() if not v]
        if len(real_fragments) == 1:
            real_label = real_fragments[0]
            real_coords = subsystem_frags[real_label]

            matched = False
            for super_label, super_atoms in supersystem_frags.items():
                if fragments_equal(real_coords, super_atoms):
                    forward_mapping[real_label] = super_label
                    forward_mapping_ghost_free[real_label] = super_label
                    matched = True
                    break
            if not matched:
                forward_mapping[real_label] = None
                forward_mapping_ghost_free[real_label] = None

            for frag_label in subsystem_frags:
                if frag_label != real_label:
                    forward_mapping[frag_label] = None
        else:
            for sub_label, sub_atoms in subsystem_frags.items():
                matched = False
                for super_label, super_atoms in supersystem_frags.items():
                    if fragments_equal(sub_atoms, super_atoms):
                        forward_mapping[sub_label] = super_label
                        if not sub_ghost_flags.get(sub_label, False):
                            forward_mapping_ghost_free[sub_label] = super_label
                        matched = True
                        break
                if not matched:
                    forward_mapping[sub_label] = None
                    if not sub_ghost_flags.get(sub_label, False):
                        forward_mapping_ghost_free[sub_label] = None

        main_label_mappings[subsys_label] = forward_mapping
        main_label_mappings_ghost_free[subsys_label] = forward_mapping_ghost_free

    # === ALT SUPERSYSTEM ===
    alt_supersystem_file = normalize_path(alternative_filenames[0])
    if alt_supersystem_file and os.path.isfile(alt_supersystem_file):
        alt_frags, alt_ghost_flags, _ = extract_fragments(alt_supersystem_file)

        for alt_label, alt_atoms in alt_frags.items():
            matched = False
            for super_label, super_atoms in supersystem_frags.items():
                if fragments_equal(alt_atoms, super_atoms):
                    alt_label_mappings["SUPERSYS"][alt_label] = super_label
                    if not alt_ghost_flags.get(alt_label, False):
                        alt_label_mappings_ghost_free["SUPERSYS"][alt_label] = super_label
                    matched = True
                    break
            if not matched:
                alt_label_mappings["SUPERSYS"][alt_label] = None
                if not alt_ghost_flags.get(alt_label, False):
                    alt_label_mappings_ghost_free["SUPERSYS"][alt_label] = None
    else:
        alt_label_mappings["SUPERSYS"] = main_label_mappings["SUPERSYS"].copy()
        alt_label_mappings_ghost_free["SUPERSYS"] = main_label_mappings_ghost_free["SUPERSYS"].copy()

    # === ALT SUBSYSTEMS ===
    for i, alt_file in enumerate(alternative_files):
        subsys_label = f"SUBSYS{i+1}"
        if not alt_file or not os.path.isfile(alt_file):
            alt_label_mappings[subsys_label] = main_label_mappings[subsys_label].copy()
            alt_label_mappings_ghost_free[subsys_label] = main_label_mappings_ghost_free[subsys_label].copy()
            continue

        alt_frags, alt_ghost_flags, _ = extract_fragments(alt_file)
        forward_mapping = {}
        forward_mapping_ghost_free = {}

        for alt_label, alt_atoms in alt_frags.items():
            matched = False
            for super_label, super_atoms in supersystem_frags.items():
                if fragments_equal(alt_atoms, super_atoms):
                    forward_mapping[alt_label] = super_label
                    if not alt_ghost_flags.get(alt_label, False):
                        forward_mapping_ghost_free[alt_label] = super_label
                    matched = True
                    break
            if not matched:
                forward_mapping[alt_label] = None
                if not alt_ghost_flags.get(alt_label, False):
                    forward_mapping_ghost_free[alt_label] = None

        alt_label_mappings[subsys_label] = forward_mapping
        alt_label_mappings_ghost_free[subsys_label] = forward_mapping_ghost_free

    print("======== Fragment Mappings ========")
    print("main_label_mappings_ghost_free:", main_label_mappings_ghost_free)
    print("alt_label_mappings_ghost_free:", alt_label_mappings_ghost_free)
    print("bsse_found:", bsse_found)

    return main_label_mappings_ghost_free, alt_label_mappings_ghost_free, bsse_found


def construct_label_mappings_singlefrag_bsse(main_filenames, alternative_filenames):
    """Construct label mappings for BSSE: one real fragment per subsystem."""

    supersystem_file = main_filenames[0]
    subsystem_files = main_filenames[1:]
    alt_subsystem_files = alternative_filenames[1:] if len(alternative_filenames) > 1 else []

    # SUPERSYS: take only non-ghost labels for the ghost-free identity
    supersys_frags, super_ghost_flags, _ = extract_fragments(supersystem_file)
    supersys_identity_ghost_free = {k: k for k, is_ghost in super_ghost_flags.items() if not is_ghost}

    def build_mapping(files):
        out = {}
        for idx, f in enumerate(files):
            label = f"SUBSYS{idx+1}"
            if not f or not os.path.isfile(f):
                out[label] = {}
                continue

            subsys_frags, ghost_flags, _ = extract_fragments(f)
            real_frags = [fid for fid, is_ghost in ghost_flags.items() if not is_ghost]
            mapping = {}

            if len(real_frags) == 1:
                real_id = real_frags[0]
                real_coords = subsys_frags[real_id]

                # Find the matched SUPERSYS fragment
                matched_super = None
                for super_id, super_coords in supersys_frags.items():
                    if fragments_equal(real_coords, super_coords):
                        matched_super = super_id
                        break

                # IMPORTANT: map fabricated index 1 (the 1x1 sheet) to the matched SUPERSYS label
                mapping[1] = matched_super  # may be None if no match
                # ghosts are intentionally ignored for the ghost-free mapping
            else:
                # Unexpected in BSSE single-frag path; leave empty
                mapping = {}

            out[label] = mapping
        return out

    main_label_mappings_ghost_free = {"SUPERSYS": supersys_identity_ghost_free}
    main_label_mappings_ghost_free.update(build_mapping(subsystem_files))

    alt_label_mappings_ghost_free = {"SUPERSYS": supersys_identity_ghost_free.copy()}
    alt_label_mappings_ghost_free.update(build_mapping(alt_subsystem_files))

    return main_label_mappings_ghost_free, alt_label_mappings_ghost_free


def subsystem_label_lists_alligned_to_main_supersystem(system_labels, main_label_mappings_ghost_free):
    """Returns lists of fragment labels of each subsystem alligned to the main SUPERSY labels."""
    alligned_subsystem_labels = []
    for label in system_labels[1:]:  # Skip SUPERSYS
        mapping = main_label_mappings_ghost_free.get(label, {})
        matched_labels = [super_label for super_label in mapping.values() if super_label is not None]
        matched_labels = sorted(set(matched_labels))
        alligned_subsystem_labels.append(matched_labels)
    return alligned_subsystem_labels


def compute_ref_diel_int_en(main_filenames, alternative_filenames, conversion_factor):
    """Extract REF DIEL interaction energy from CPCM Dielectric values."""
    ref_diel_values = {}
    subsystem_diel_values = {}

    pattern = r"CPCM Dielectric\s*:\s*([-+]?\d*\.\d+|\d+)"

    def extract_diel(filename):
        filename = normalize_path(filename)
        if not filename:  # Avoid error message for empty string or None in ALT file list
            return None
        try:
            with open(filename, 'r') as f:
                for line in f:
                    match = re.search(pattern, line)
                    if match:
                        return float(match.group(1))
        except FileNotFoundError:
            print(f"File not found: {filename}")
        return None

    ref_diel_values["SUPERSYS"] = extract_diel(main_filenames[0]) or extract_diel(alternative_filenames[0])

    for i in range(1, len(main_filenames)):
        label = f"SUBSYS{i}"
        diel = extract_diel(main_filenames[i]) or extract_diel(alternative_filenames[i])
        ref_diel_values[label] = diel
        if diel is not None:
            subsystem_diel_values[label] = diel

    if ref_diel_values["SUPERSYS"] is None or len(subsystem_diel_values) != len(main_filenames) - 1:
        print("REF dielectric values incomplete; REF DIEL int energy set to 0")
        return ref_diel_values, 0.0

    ref_diel_int_energy = (ref_diel_values["SUPERSYS"] - sum(subsystem_diel_values.values())) * conversion_factor
    return ref_diel_values, ref_diel_int_energy


def compute_ref_cds_int_en(main_filenames, alternative_filenames, conversion_factor):
    """Extract REF CDS interaction energy from SMD CDS values."""
    ref_cds_values = {}
    subsystem_cds_values = {}

    # Modified pattern to target SMD CDS (Gcds)
    pattern = r"SMD CDS \(Gcds\)\s*:\s*([-+]?\d*\.\d+)"

    def extract_cds(filename):
        filename = normalize_path(filename)
        if not filename:  # Avoid error message for empty string or None in ALT file list
            return None
        try:
            with open(filename, 'r') as f:
                for line in f:
                    match = re.search(pattern, line)
                    if match:
                        return float(match.group(1))
        except FileNotFoundError:
            print(f"File not found: {filename}")
        except Exception as e:
            print(f"Error reading {filename} for CDS: {e}")
        return None

    # Extract CDS for the SUPERSYS, with fallback
    ref_cds_values["SUPERSYS"] = extract_cds(main_filenames[0])
    if ref_cds_values["SUPERSYS"] is None and len(alternative_filenames) > 0 and alternative_filenames[0]:
        ref_cds_values["SUPERSYS"] = extract_cds(alternative_filenames[0])

    # Extract CDS for each SUBSYS, with fallback
    for i in range(1, len(main_filenames)):
        label = f"SUBSYS{i}"
        cds = extract_cds(main_filenames[i])
        if cds is None and len(alternative_filenames) > i and alternative_filenames[i]:
            cds = extract_cds(alternative_filenames[i])
            
        ref_cds_values[label] = cds
        if cds is not None:
            subsystem_cds_values[label] = cds

    # Check for completeness of extracted CDS values
    if ref_cds_values["SUPERSYS"] is None or len(subsystem_cds_values) != len(main_filenames) - 1:
        print("REF CDS values incomplete; REF CDS int energy set to 0")
        # Ensure all expected subsystem labels are in ref_cds_values, even if None
        for i in range(1, len(main_filenames)):
            label = f"SUBSYS{i}"
            if label not in ref_cds_values:
                ref_cds_values[label] = None
        return ref_cds_values, 0.0

    # Calculate REF CDS interaction energy
    ref_cds_int_energy = (ref_cds_values["SUPERSYS"] - sum(subsystem_cds_values.values())) * conversion_factor
    return ref_cds_values, ref_cds_int_energy


def compute_corr_diel_int_en(main_filenames, alternative_filenames, conversion_factor, method):
    """Extract CORR DIEL interaction energy from C-PCM corr. term."""
    corr_diel_values = {}
    subsystem_diel_values = {}

    pattern = r"C-PCM corr\. term \(included in E\(CORR\)\).*?([-+]?\d*\.\d+)"

    def extract_corr_diel(filename):
        filename = normalize_path(filename)
        if not filename:  # Avoid error message for empty string or None in ALT file list
            return None
        try:
            with open(filename, 'r') as f:
                for line in f:
                    match = re.search(pattern, line)
                    if match:
                        return float(match.group(1))
        except FileNotFoundError:
            print(f"File not found: {filename}")
        return None

    corr_diel_values["SUPERSYS"] = extract_corr_diel(main_filenames[0]) or extract_corr_diel(alternative_filenames[0])

    for i in range(1, len(main_filenames)):
        label = f"SUBSYS{i}"
        if method.lower() == "hfld":
            diel = 0.0
        else:
            diel = extract_corr_diel(main_filenames[i]) or extract_corr_diel(alternative_filenames[i]) or 0.0
        corr_diel_values[label] = diel
        subsystem_diel_values[label] = diel

    if corr_diel_values["SUPERSYS"] is None:
        print("CORR dielectric value for SUPERSYS missing; CORR DIEL int energy set to 0")
        return corr_diel_values, 0.0

    corr_diel_int_energy = (corr_diel_values["SUPERSYS"] - sum(subsystem_diel_values.values())) * conversion_factor
    return corr_diel_values, corr_diel_int_energy


def compute_bulk_solvation_contribution(main_filenames, alternative_filenames, conversion_factor, method):
    """Compute bulk solvent interaction energy as REF DIEL + CORR DIEL + SMD CDS.
	Besides this cumulative term, QM energy incorporates polarization contribution from DIEL."""
    ref_diel_values, ref_diel_int_energy = compute_ref_diel_int_en(main_filenames, alternative_filenames, conversion_factor)
    ref_cds_values, ref_cds_int_energy = compute_ref_cds_int_en(main_filenames, alternative_filenames, conversion_factor)
    corr_diel_values, corr_diel_int_energy = compute_corr_diel_int_en(main_filenames, alternative_filenames, conversion_factor, method)

    solv_values = {}
    for key in set(ref_diel_values) | set(corr_diel_values) | set(ref_cds_values):
        solv_values[key] = (ref_diel_values.get(key) or 0.0) + (corr_diel_values.get(key) or 0.0) + (ref_cds_values.get(key) or 0.0)

    solv_int_energy = ref_diel_int_energy + corr_diel_int_energy + ref_cds_int_energy
    return solv_values, solv_int_energy, ref_diel_int_energy, ref_cds_int_energy, corr_diel_int_energy


def check_local_energy_decomposition(filename, patterns):
    """Check if the file contains the 'LOCAL ENERGY DECOMPOSITION' section."""
    normalized_filename = normalize_path(filename)
    
    try:
        with open(normalized_filename, 'r') as file:
            content = file.read()
        return bool(re.search(patterns.local_energy_decomp_pattern, content))
    except FileNotFoundError:
        print(f"File {normalized_filename} not found. Skipping.")
        return False


def extract_generic_numbers_from_regex(regex_str, content):
    """Extracts a list of all numbers based on a direct regex string."""
    numbers = []
    if not regex_str:
        return numbers

    match = re.search(regex_str, content, re.DOTALL | re.UNICODE) 
    if match:
        number_pattern = r"([-]?\d*\.\d+)"
        numbers = [float(number) for number in re.findall(number_pattern, match.group(0))]
            
    return numbers


def extract_intra_ref_alt_numbers(content):
    """Extracts numbers specifically for the 'intra_ref_alt' pattern."""
    numbers = []
    pattern_regex = patterns.PATTERNS.get("intra_ref_alt")
    
    if not pattern_regex:
        return numbers

    for match in re.finditer(pattern_regex, content, re.DOTALL | re.UNICODE): 
        try:
            total_energy = float(match.group(2)) 
            numbers.append(total_energy)
        except (ValueError, IndexError):
            pass
            
    return numbers


def extract_numbers(regex_str_input, content):
    """Extracts numbers based on a provided regex string. Handles special fallback logic for the 'intra_ref' pattern."""
    if not regex_str_input:
        return []

    # Get the specific regex strings for intra_ref and intra_ref_alt from the Patterns instance.
    intra_ref_regex = patterns.PATTERNS.get("intra_ref")
    intra_ref_alt_regex = patterns.PATTERNS.get("intra_ref_alt") # Needed for fallback

    # Check if the input regex string is the specific 'intra_ref' pattern.
    if regex_str_input == intra_ref_regex:
        # First, try to extract using the intra_ref pattern itself.
        results = extract_generic_numbers_from_regex(regex_str_input, content)
        
        # If no results are found from intra_ref, then try intra_ref_alt as a fallback.
        if not results:
            results = extract_intra_ref_alt_numbers(content)
        return results
    
    # Use the generic regex extraction for any other regex string (e.g., intra_corr, intra_strong_pairs, etc.)
    else:
        return extract_generic_numbers_from_regex(regex_str_input, content)


def extract_first_match_from_file(filename, patterns, method, use_ref_as_rhf_in_hfld=None):
    """Extract the first matches for E(0), strong pairs, weak pairs, and triples correction from the file."""
    normalized_filename = normalize_path(filename)

    try:
        with open(normalized_filename, 'r') as file:
            content = file.read()
        
        # Search for the energy patterns
        e_ref_match = re.search(patterns.pattern_e0, content)
        e_wp_match = re.search(patterns.pattern_weak_corr, content)
        e_t_match = re.search(patterns.pattern_triples_corr, content)
        e_sp_match = re.search(patterns.pattern_strong_corr, content)

        e_ref = float(e_ref_match.group(1)) if e_ref_match else 0.0
        e_wp = float(e_wp_match.group(1)) if e_wp_match else 0.0
        e_t = float(e_t_match.group(1)) if e_t_match else 0.0
        if e_sp_match:
            e_sp = float(e_sp_match.group(1))
            cpcm_corr_match = re.search(r"C-PCM corr\. term \(included in E\(CORR\)\).*?([-+]?\d*\.\d+)", content)
            if cpcm_corr_match:
                e_sp -= float(cpcm_corr_match.group(1))
        else:
            e_sp = 0.0

        # Subtract dielectric and cds either from REF or Total energy
        diel_match = re.search(r"CPCM Dielectric\s*:\s*([-+]?\d*\.\d+|\d+)", content)
        cds_match = re.search(r"SMD CDS \(Gcds\)\s*:\s*([-+]?\d*\.\d+)", content)

        if e_ref == 0.0 and method.lower() == 'hfld' and use_ref_as_rhf_in_hfld:
            total_match = re.search(r"Total Energy\s+:\s+([-]?\d+\.\d+)", content)
            if total_match:
                e_ref = float(total_match.group(1))
                if diel_match:
                    e_ref -= float(diel_match.group(1))
                if cds_match:
                    e_ref -= float(cds_match.group(1))
        elif diel_match or cds_match:
            if diel_match:
                e_ref -= float(diel_match.group(1))
            if cds_match:
                e_ref -= float(cds_match.group(1))

        return e_ref, e_sp, e_wp, e_t

    except FileNotFoundError:
        print(f"File {normalized_filename} not found. Skipping.")
        return 0.0, 0.0, 0.0, 0.0


def determine_matrix_size(filename, patterns):
    """Determine the size of the matrix based on the intra_ref pattern in the file."""
    normalized_filename = normalize_path(filename)

    try:
        with open(normalized_filename, 'r') as file:
            content = file.read()
        intra_ref = extract_numbers(patterns.PATTERNS["intra_ref"], content)
        if not intra_ref:
            intra_ref = extract_numbers(patterns.PATTERNS["intra_ref_alt"], content)
        return len(intra_ref)
    except FileNotFoundError:
        raise FileNotFoundError(f"File {normalized_filename} not found. Unable to determine matrix size.")


def extract_els_exch_matrices(content, primary_pattern, alternative_pattern):
    """Extracts electrostatics and exchange matrices based on the primary and alternative patterns."""
    matches = re.findall(primary_pattern, content)
    if not matches:
        matches = re.findall(alternative_pattern, content)
    if not matches:
        return None, None

    fragment_ids = set()
    for match in matches:
        fragment_ids.add(int(match[0]))
        fragment_ids.add(int(match[1]))
    matrix_size = max(fragment_ids)

    electrostatics_matrix = np.zeros((matrix_size, matrix_size))
    exchange_matrix = np.zeros((matrix_size, matrix_size))

    for match in matches:
        i = int(match[0]) - 1 
        j = int(match[1]) - 1 
        els_value = float(match[2])
        exch_value = float(match[3])
        electrostatics_matrix[i, j] = electrostatics_matrix[j, i] = els_value
        exchange_matrix[i, j] = exchange_matrix[j, i] = exch_value

    return electrostatics_matrix, exchange_matrix


def process_led_file(filename, patterns, intra_ref_list, intra_corr_list, intra_strong_pairs_list, intra_triples_list, intra_weak_pairs_list,
    singles_contribution_list, inter_ref_matrices, inter_corr_matrices, electrostat_matrices, exchange_matrices, inter_strong_pairs_matrices,
    inter_triples_matrices, inter_weak_pairs_matrices, dispersion_strong_pairs_matrices, system_label,prefix_suffix=""):

    if not filename:
        return

    normalized_filename = normalize_path(filename)

    try:
        with open(normalized_filename, 'r') as file:
            content = file.read()
    except FileNotFoundError:
        print(f"File {normalized_filename} not found. Skipping.")
        return

    # Extract intra-component values
    intra_ref = extract_numbers(patterns.PATTERNS["intra_ref"], content) # automatically check intra_ref_alt when intra_ref is not present
    intra_corr = extract_numbers(patterns.PATTERNS["intra_corr"], content)
    intra_strong_pairs = extract_numbers(patterns.PATTERNS["intra_strong_pairs"], content)
    intra_triples = extract_numbers(patterns.PATTERNS["intra_triples"], content)
    intra_weak_pairs = extract_numbers(patterns.PATTERNS["intra_weak_pairs"], content)
    singles_contribution = extract_numbers(patterns.PATTERNS["singles_contribution"], content)

    # Match patterns
    matches_combined = re.findall(patterns.PATTERNS["ref_corr_inter"], content)
    matches_inter_correlation_full = re.findall(patterns.PATTERNS["corr_inter_full"], content)
    matches_inter_correlation_partial = re.findall(patterns.PATTERNS["corr_inter_partial"], content)
    matches_dispersion = re.findall(patterns.PATTERNS["dispersion_strong_pairs"], content)

    # Determine matrix size dynamically from fragment IDs
    all_matches = matches_combined + matches_inter_correlation_full + matches_inter_correlation_partial + matches_dispersion
    fragment_ids = set()
    for match in all_matches:
        try:
            fragment_ids.add(int(match[0]))
            fragment_ids.add(int(match[1]))
        except (IndexError, ValueError):
            continue
    matrix_size = max(fragment_ids) if fragment_ids else 1

    # Initialize matrices
    inter_ref_matrix = np.zeros((matrix_size, matrix_size))
    inter_corr_matrix = np.zeros((matrix_size, matrix_size))
    electrostatics_matrix = np.zeros((matrix_size, matrix_size))
    exchange_matrix = np.zeros((matrix_size, matrix_size))
    inter_strong_pairs_matrix = np.zeros((matrix_size, matrix_size))
    inter_triples_matrix = np.zeros((matrix_size, matrix_size))
    inter_weak_pairs_matrix = np.zeros((matrix_size, matrix_size))
    dispersion_strong_pairs_matrix = np.zeros((matrix_size, matrix_size))

    # Electrostatics & Exchange
    electrostat_pattern = patterns.PATTERNS["ref_inter"]
    electrostat_alt_pattern = patterns.PATTERNS["ref_corr_inter"]
    electrostat_matrix, exch_matrix = extract_els_exch_matrices(content, electrostat_pattern, electrostat_alt_pattern)
    if electrostat_matrix is not None:
        electrostatics_matrix = electrostat_matrix
        exchange_matrix = exch_matrix

    # Matrix population helper
    def populate_matrix(matches, matrix, idx):
        for match in matches:
            try:
                i = int(match[0]) - 1
                j = int(match[1]) - 1
                if i < matrix_size and j < matrix_size and idx < len(match):
                    value = float(match[idx])
                    matrix[i, j] = matrix[j, i] = value
            except (IndexError, ValueError):
                continue

    populate_matrix(matches_combined, inter_ref_matrix, 2)
    populate_matrix(matches_combined, inter_corr_matrix, 3)
    if matches_inter_correlation_full:
        populate_matrix(matches_inter_correlation_full, inter_strong_pairs_matrix, 2)
        populate_matrix(matches_inter_correlation_full, inter_triples_matrix, 3)
        populate_matrix(matches_inter_correlation_full, inter_weak_pairs_matrix, 4)
    else:
        populate_matrix(matches_inter_correlation_partial, inter_strong_pairs_matrix, 2)
        populate_matrix(matches_inter_correlation_partial, inter_weak_pairs_matrix, 3)

    for match in matches_dispersion:
        try:
            i = int(match[0]) - 1
            j = int(match[1]) - 1
            value = float(match[2])
            if i < matrix_size and j < matrix_size:
                dispersion_strong_pairs_matrix[i, j] = dispersion_strong_pairs_matrix[j, i] = value
        except (IndexError, ValueError):
            continue

    file_prefix = system_label + prefix_suffix

    # Append results
    intra_ref_list.append((intra_ref, file_prefix))
    intra_corr_list.append((intra_corr, file_prefix))
    intra_strong_pairs_list.append((intra_strong_pairs, file_prefix))
    intra_triples_list.append((intra_triples, file_prefix))
    intra_weak_pairs_list.append((intra_weak_pairs, file_prefix))
    singles_contribution_list.append((singles_contribution, file_prefix))
    inter_ref_matrices.append((inter_ref_matrix, file_prefix))
    inter_corr_matrices.append((inter_corr_matrix, file_prefix))
    electrostat_matrices.append((electrostatics_matrix, file_prefix))
    exchange_matrices.append((exchange_matrix, file_prefix))
    inter_strong_pairs_matrices.append((inter_strong_pairs_matrix, file_prefix))
    inter_triples_matrices.append((inter_triples_matrix, file_prefix))
    inter_weak_pairs_matrices.append((inter_weak_pairs_matrix, file_prefix))
    dispersion_strong_pairs_matrices.append((dispersion_strong_pairs_matrix, file_prefix))


def write_matrices_to_excel(filename, intra_ref_list, intra_corr_list, intra_strong_pairs_list, intra_triples_list, intra_weak_pairs_list,
    singles_contribution_list, inter_ref_matrices, inter_corr_matrices, electrostat_matrices, exchange_matrices, inter_strong_pairs_matrices,
    inter_triples_matrices, inter_weak_pairs_matrices, dispersion_strong_pairs_matrices):
    """Writes collected data to an Excel file using each matrix’s actual shape."""

    normalized_filename = normalize_path(filename)

    with pd.ExcelWriter(normalized_filename, engine='xlsxwriter') as writer:

        def write_diag_matrix(data_list, name_prefix):
            for values, label in data_list:
                size = len(values)
                mat = np.zeros((size, size))
                for i in range(size):
                    mat[i, i] = values[i]
                df = pd.DataFrame(mat, index=range(1, size + 1), columns=range(1, size + 1))
                df.to_excel(writer, sheet_name=f'{name_prefix} {label}')

        def write_square_matrix(matrix_list, name_prefix):
            for matrix, label in matrix_list:
                size = matrix.shape[0]
                df = pd.DataFrame(matrix, index=range(1, size + 1), columns=range(1, size + 1))
                df.to_excel(writer, sheet_name=f'{name_prefix} {label}')

        write_diag_matrix(intra_ref_list, 'Intra REF')
        write_diag_matrix(intra_corr_list, 'Intra CORR')
        write_diag_matrix(intra_strong_pairs_list, 'Intra SP')
        write_diag_matrix(intra_triples_list, 'Intra T')
        write_diag_matrix(intra_weak_pairs_list, 'Intra WP')
        write_diag_matrix(singles_contribution_list, 'Singles')

        write_square_matrix(inter_ref_matrices, 'Inter REF')
        write_square_matrix(inter_corr_matrices, 'Inter CORR')
        write_square_matrix(electrostat_matrices, 'Electrostat')
        write_square_matrix(exchange_matrices, 'Exchange')
        write_square_matrix(inter_strong_pairs_matrices, 'Inter SP')
        write_square_matrix(inter_triples_matrices, 'Inter T')
        write_square_matrix(inter_weak_pairs_matrices, 'Inter WP')
        write_square_matrix(dispersion_strong_pairs_matrices, 'Disp SP')


def multifrag_system_processing(main_filenames, alternative_filenames, LEDAW_output_path, system_labels):
    # Initialize lists to store results for multiple files
    intra_ref_list = []
    intra_corr_list = []
    intra_strong_pairs_list = []
    intra_triples_list = []
    intra_weak_pairs_list = []
    singles_contribution_list = []
    inter_ref_matrices = []
    inter_corr_matrices = []
    electrostat_matrices = []
    exchange_matrices = []
    inter_strong_pairs_matrices = []
    inter_triples_matrices = []
    inter_weak_pairs_matrices = []
    dispersion_strong_pairs_matrices = []

    # Process each main and alternative file using unified function
    for main_file, alt_file, system_label in zip(main_filenames, alternative_filenames, system_labels):
        process_led_file(
            filename=main_file,
            patterns=patterns,
            intra_ref_list=intra_ref_list,
            intra_corr_list=intra_corr_list,
            intra_strong_pairs_list=intra_strong_pairs_list,
            intra_triples_list=intra_triples_list,
            intra_weak_pairs_list=intra_weak_pairs_list,
            singles_contribution_list=singles_contribution_list,
            inter_ref_matrices=inter_ref_matrices,
            inter_corr_matrices=inter_corr_matrices,
            electrostat_matrices=electrostat_matrices,
            exchange_matrices=exchange_matrices,
            inter_strong_pairs_matrices=inter_strong_pairs_matrices,
            inter_triples_matrices=inter_triples_matrices,
            inter_weak_pairs_matrices=inter_weak_pairs_matrices,
            dispersion_strong_pairs_matrices=dispersion_strong_pairs_matrices,
            system_label=system_label,
            prefix_suffix=""
        )

        if alt_file:
            process_led_file(
                filename=alt_file,
                patterns=patterns,
                intra_ref_list=intra_ref_list,
                intra_corr_list=intra_corr_list,
                intra_strong_pairs_list=intra_strong_pairs_list,
                intra_triples_list=intra_triples_list,
                intra_weak_pairs_list=intra_weak_pairs_list,
                singles_contribution_list=singles_contribution_list,
                inter_ref_matrices=inter_ref_matrices,
                inter_corr_matrices=inter_corr_matrices,
                electrostat_matrices=electrostat_matrices,
                exchange_matrices=exchange_matrices,
                inter_strong_pairs_matrices=inter_strong_pairs_matrices,
                inter_triples_matrices=inter_triples_matrices,
                inter_weak_pairs_matrices=inter_weak_pairs_matrices,
                dispersion_strong_pairs_matrices=dispersion_strong_pairs_matrices,
                system_label=system_label,
                prefix_suffix=" ALT"
            )

    # Ensure output path exists
    normalized_LEDAW_output_path = normalize_path(LEDAW_output_path)
    os.makedirs(normalized_LEDAW_output_path, exist_ok=True)

    output_file = os.path.join(normalized_LEDAW_output_path, 'tmp1.xlsx')

    # Write matrices to Excel (label mappings not needed for sizing anymore)
    write_matrices_to_excel(
        filename=output_file,
        intra_ref_list=intra_ref_list,
        intra_corr_list=intra_corr_list,
        intra_strong_pairs_list=intra_strong_pairs_list,
        intra_triples_list=intra_triples_list,
        intra_weak_pairs_list=intra_weak_pairs_list,
        singles_contribution_list=singles_contribution_list,
        inter_ref_matrices=inter_ref_matrices,
        inter_corr_matrices=inter_corr_matrices,
        electrostat_matrices=electrostat_matrices,
        exchange_matrices=exchange_matrices,
        inter_strong_pairs_matrices=inter_strong_pairs_matrices,
        inter_triples_matrices=inter_triples_matrices,
        inter_weak_pairs_matrices=inter_weak_pairs_matrices,
        dispersion_strong_pairs_matrices=dispersion_strong_pairs_matrices,
    )


def singlefrag_system_processing(labeled_main_filenames, labeled_alt_filenames, LEDAW_output_path, method, use_ref_as_rhf_in_hfld=None):

    matrices = {}

    for main_file, system_label in labeled_main_filenames.items():
        alt_file = None
        for af, af_label in labeled_alt_filenames.items():
            if af_label == system_label:
                alt_file = af
                break

        # Process main file
        if main_file:
            led_exists = check_local_energy_decomposition(main_file, patterns)
            if not led_exists:
                e_ref, e_sp, e_wp, e_t = extract_first_match_from_file(main_file, patterns, method, use_ref_as_rhf_in_hfld)
                size = 1 if e_ref != 0 else 0
                ref_matrix = np.zeros((size, size))
                sp_matrix = np.zeros((size, size))
                wp_matrix = np.zeros((size, size))
                t_matrix = np.zeros((size, size))

                if size == 1:
                    ref_matrix[0, 0] = e_ref
                    sp_matrix[0, 0] = e_sp
                    wp_matrix[0, 0] = e_wp
                    t_matrix[0, 0] = e_t

                matrices[f'Intra REF {system_label}'] = ref_matrix
                matrices[f'Intra SP {system_label}'] = sp_matrix
                matrices[f'Intra WP {system_label}'] = wp_matrix
                matrices[f'Intra T {system_label}'] = t_matrix

        # Process ALT file
        if alt_file:
            led_exists = check_local_energy_decomposition(alt_file, patterns)
            if not led_exists:
                e_ref, e_sp, e_wp, e_t = extract_first_match_from_file(alt_file, patterns, method, use_ref_as_rhf_in_hfld)
                size = 1 if e_ref != 0 else 0
                ref_matrix = np.zeros((size, size))
                sp_matrix = np.zeros((size, size))
                wp_matrix = np.zeros((size, size))
                t_matrix = np.zeros((size, size))

                if size == 1:
                    ref_matrix[0, 0] = e_ref
                    sp_matrix[0, 0] = e_sp
                    wp_matrix[0, 0] = e_wp
                    t_matrix[0, 0] = e_t

                file_prefix = system_label + ' ALT'
                matrices[f'Intra REF {file_prefix}'] = ref_matrix
                matrices[f'Intra SP {file_prefix}'] = sp_matrix
                matrices[f'Intra WP {file_prefix}'] = wp_matrix
                matrices[f'Intra T {file_prefix}'] = t_matrix

    # Ensure the output directory exists
    normalized_LEDAW_output_path = normalize_path(LEDAW_output_path)
    os.makedirs(normalized_LEDAW_output_path, exist_ok=True)

    output_file = os.path.join(normalized_LEDAW_output_path, 'tmp2.xlsx')

    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        for name, matrix in matrices.items():
            size = matrix.shape[0]
            df = pd.DataFrame(matrix, index=range(1, size + 1), columns=range(1, size + 1))
            df.to_excel(writer, sheet_name=name)

    return matrices


def combine_unprocessed_LED_data_fies(LEDAW_output_path):
    """Combines two temporary Excel files with multiple sheets into one file and deletes the source Excel files."""
    
    # Normalize the LEDAW_output_path
    normalized_LEDAW_output_path = normalize_path(LEDAW_output_path)
    
    file1 = os.path.join(normalized_LEDAW_output_path, "tmp1.xlsx")
    file2 = os.path.join(normalized_LEDAW_output_path, "tmp2.xlsx")
    write_to_excel_filename = os.path.join(normalized_LEDAW_output_path, "Unprocessed_LED_matrices0.xlsx")

    
    # Create a new Excel writer object
    with pd.ExcelWriter(write_to_excel_filename) as writer:
        # Load the first file and write its sheets to the output file
        with pd.ExcelFile(file1) as xls1:
            for sheet_name in xls1.sheet_names:
                df = pd.read_excel(xls1, sheet_name=sheet_name)
                # If the first column has the 'Unnamed' pattern, remove it, but only if the DataFrame is not empty
                if not df.empty and df.columns[0] == "Unnamed: 0":
                    df = df.rename(columns={"Unnamed: 0": ""})
                df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        # Load the second file and write its sheets to the output file
        with pd.ExcelFile(file2) as xls2:
            for sheet_name in xls2.sheet_names:
                df = pd.read_excel(xls2, sheet_name=sheet_name)
                # If the first column has the 'Unnamed' pattern, remove it, but only if the DataFrame is not empty
                if not df.empty and df.columns[0] == "Unnamed: 0":
                    df = df.rename(columns={"Unnamed: 0": ""})
                df.to_excel(writer, sheet_name=sheet_name, index=False)
    
    # Delete the source Excel files
    os.remove(file1)
    os.remove(file2)
    
    # Normalize the path before printing
    normalized_output_file = normalize_path(write_to_excel_filename)
    print(f"LED energy components from each file were written as matrices without any processing to {normalized_output_file}")


def unify_labels(system_labels, alternative_labels, main_label_mappings_ghost_free, alternative_label_mappings_ghost_free, LEDAW_output_path):
    """Relabels and resizes all matrix sheets in Excel based on mappings to main SUPERSYS fragment labels."""

    normalized_LEDAW_output_path = normalize_path(LEDAW_output_path)
    input_excel_file = os.path.join(normalized_LEDAW_output_path, 'Unprocessed_LED_matrices0.xlsx')
    output_excel_file = os.path.join(normalized_LEDAW_output_path, 'Unprocessed_LED_matrices1_withUnifiedFragmentLabeling.xlsx')

    xl = pd.ExcelFile(input_excel_file)

    # Determine full label range from SUPERSYS mapping
    supersys_labels = [
        value for value in main_label_mappings_ghost_free.get("SUPERSYS", {}).values()
        if value is not None
    ]
    if supersys_labels:
        max_label = max(supersys_labels)
        full_label_range = list(range(1, max_label + 1))  # Include all labels up to max
    else:
        full_label_range = []

    with pd.ExcelWriter(output_excel_file, engine='openpyxl') as writer:
        for sheet_name in xl.sheet_names:
            df = xl.parse(sheet_name, index_col=0)
            if df.empty:
                continue

            is_alt = 'ALT' in sheet_name
            # Extract the system label
            words = sheet_name.strip().split()
            if not words:
                matched_label = None
            else:
                if words[-1] == "ALT" and len(words) >= 2:
                    candidate_label = words[-2]
                else:
                    candidate_label = words[-1]
                matched_label = candidate_label if candidate_label in (system_labels + alternative_labels) else None

            if not matched_label:
                print(f"Skipping sheet '{sheet_name}': no system label matched.")
                df.to_excel(writer, sheet_name=sheet_name)
                continue

            mappings = alternative_label_mappings_ghost_free if is_alt else main_label_mappings_ghost_free
            label_mapping = mappings.get(matched_label, {})
            if not label_mapping:
                if df.shape[0] == 1 and matched_label.startswith("SUBSYS"):
                    # Handle single-fragment subsystem: map its single intra value
                    extended_df = pd.DataFrame(0.0, index=full_label_range, columns=full_label_range)

                    # Find the mapped SUPERSYS label (take the first available one)
                    super_labels = list(main_label_mappings_ghost_free.get("SUPERSYS", {}).values())
                    if super_labels:
                        new_label = super_labels[0]
                        extended_df.at[new_label, new_label] = df.iat[0, 0]

                    extended_df.to_excel(writer, sheet_name=sheet_name)
                    continue
                else:
                    print(f"Skipping sheet '{sheet_name}': no label mapping found for {matched_label}.")
                    continue

            try:
                # Ensure both index and columns are numeric
                df.index = pd.to_numeric(df.index, errors='coerce')
                df.columns = pd.to_numeric(df.columns, errors='coerce')

                # Prepare the extended matrix
                extended_df = pd.DataFrame(0.0, index=full_label_range, columns=full_label_range)

                # Safely remap and copy values
                for orig_row in df.index:
                    new_row = label_mapping.get(orig_row)
                    if new_row is None or new_row not in full_label_range:
                        continue
                    for orig_col in df.columns:
                        new_col = label_mapping.get(orig_col)
                        if new_col is None or new_col not in full_label_range:
                            continue
                        extended_df.at[new_row, new_col] = df.at[orig_row, orig_col]

                extended_df.to_excel(writer, sheet_name=sheet_name)

            except Exception as e:
                print(f"Error processing sheet '{sheet_name}': {e}")
                df.to_excel(writer, sheet_name=sheet_name)

    print(f"Matrices with unified labels written to '{normalize_path(output_excel_file)}'")


def combine_main_ALT(LEDAW_output_path):
    """Process the Excel file to compare and handle ALT sheets, then write results to a new file."""
    
    # Normalize the LEDAW_output_path
    normalized_LEDAW_output_path = normalize_path(LEDAW_output_path)
    
    # Static file names with normalized path
    input_excel_file = os.path.join(normalized_LEDAW_output_path, 'Unprocessed_LED_matrices1_withUnifiedFragmentLabeling.xlsx')
    output_excel_file = os.path.join(normalized_LEDAW_output_path, 'Unprocessed_LED_matrices2_withoutALTlabel.xlsx')
    
    difference_detected = False  # Flag to check if any differences were found
    differences = []  # List to store sheet names with differences
    sheets_to_remove = []
    sheets_to_rename = {}

    # Dictionary to temporarily store sheet data
    sheet_data = {}

    # Open the input Excel file and read all sheets once
    with pd.ExcelFile(input_excel_file) as xl:
        for sheet_name in xl.sheet_names:
            df = xl.parse(sheet_name, index_col=0)
            sheet_data[sheet_name] = df

            # Identify corresponding ALT sheets
            if 'ALT' in sheet_name:
                original_sheet_name = sheet_name.replace(' ALT', '')
                if original_sheet_name in xl.sheet_names:
                    df_original = sheet_data[original_sheet_name]
                    
                    if df_original.values.sum() == 0:
                        sheets_to_remove.append(original_sheet_name)
                        sheets_to_rename[sheet_name] = original_sheet_name
                    elif df.values.sum() == 0:
                        sheets_to_remove.append(sheet_name)
                    else:
                        diff = abs(df - df_original)
                        if (diff.values > 1e-4).any():
                            differences.append(original_sheet_name)
                        else:
                            sheets_to_remove.append(sheet_name)
                else:
                    sheets_to_rename[sheet_name] = original_sheet_name

    # Write the processed sheets to a new Excel file
    with pd.ExcelWriter(output_excel_file, engine='openpyxl') as writer:
        for sheet_name, df in sheet_data.items():
            if sheet_name not in sheets_to_remove:
                # Rename sheet if needed
                if sheet_name in sheets_to_rename:
                    sheet_name = sheets_to_rename[sheet_name]
                df.to_excel(writer, sheet_name=sheet_name)
    
    # Print message before raising the exception
    normalized_output_file = normalize_path(output_excel_file)
    print(f"Redundant matrices from main and alternative files were eliminated. Cleaned data were written to '{normalized_output_file}'")

    # Handle differences after writing to the new Excel file
    if differences:
        for sheet_name in differences:
            print(f"Matrices in '{sheet_name}' and its alternative are different")
        # Raise the exception after printing the message
        raise MatrixDifferenceError("Significant differences were found in the above listed matrices. Execution was continued with the data in the main file.")


def combine_intra_inter_matrices(LEDAW_output_path):
    """Combines Intra and Inter matrices from an Excel file and writes the results to a new file."""
    
    # Normalize the LEDAW_output_path
    normalized_LEDAW_output_path = normalize_path(LEDAW_output_path)
    
    # Static file names with normalized paths
    input_excel_file = os.path.join(normalized_LEDAW_output_path, 'Unprocessed_LED_matrices2_withoutALTlabel.xlsx')
    output_excel_file = os.path.join(normalized_LEDAW_output_path, 'Unprocessed_LED_matrices3_combined_INTRA-INTER.xlsx')
    
    combined_sheets = {}
    sheets_to_transfer = []

    # Open the input Excel file and process it
    with pd.ExcelFile(input_excel_file) as xl:
        # Group sheets by their base names without "Intra" or "Inter"
        sheet_groups = {}
        for sheet_name in xl.sheet_names:
            if 'Intra' in sheet_name:
                base_name = sheet_name.replace('Intra ', '')
                if base_name not in sheet_groups:
                    sheet_groups[base_name] = {'Intra': None, 'Inter': None}
                sheet_groups[base_name]['Intra'] = sheet_name
            elif 'Inter' in sheet_name:
                base_name = sheet_name.replace('Inter ', '')
                if base_name not in sheet_groups:
                    sheet_groups[base_name] = {'Intra': None, 'Inter': None}
                sheet_groups[base_name]['Inter'] = sheet_name
            else:
                # Non-Intra/Inter sheets are directly added for transfer
                sheets_to_transfer.append(sheet_name)

        # Combine Intra and Inter sheets
        for base_name, sheets in sheet_groups.items():
            intra_df = xl.parse(sheets['Intra'], index_col=0) if sheets['Intra'] else None
            inter_df = xl.parse(sheets['Inter'], index_col=0) if sheets['Inter'] else None

            if intra_df is not None and inter_df is not None:
                combined_df = intra_df + inter_df
                combined_sheets[base_name] = combined_df
            elif intra_df is not None:
                combined_sheets[base_name] = intra_df
            elif inter_df is not None:
                combined_sheets[base_name] = inter_df

    # Write the combined matrices and other sheets to the output Excel file
    with pd.ExcelWriter(output_excel_file, engine='openpyxl') as writer:
        # Write combined intra-inter matrices
        for base_name, combined_df in combined_sheets.items():
            combined_df.to_excel(writer, sheet_name=base_name)
        
        # Transfer the remaining sheets
        with pd.ExcelFile(input_excel_file) as xl:  # Reopen for reading remaining sheets
            for sheet_name in sheets_to_transfer:
                df = xl.parse(sheet_name, index_col=0)
                df.to_excel(writer, sheet_name=sheet_name)
    
    normalized_output_file = normalize_path(output_excel_file)
    print(f"Intra and Inter matrices were combined and written to '{normalized_output_file}'")


def compute_all_standard_led_int_en_matrices(system_labels, conversion_factor, method, LEDAW_output_path, main_subsystem_matching_labels, main_filenames, alternative_filenames, bsse_found):
    """Process LED matrices of supersystem and its subsystems to compute LED interaction energy map and its components, and write results to a new Excel file."""

    normalized_LEDAW_output_path = normalize_path(LEDAW_output_path)

    input_excel_file = os.path.join(normalized_LEDAW_output_path, 'Unprocessed_LED_matrices3_combined_INTRA-INTER.xlsx')
    output_excel_file = os.path.join(normalized_LEDAW_output_path, 'Unrelabeled_All_Standard_LED_matrices4.xlsx')

    supersystem_label = system_labels[0]
    subsystem_labels = system_labels[1:]
    matrices = {}

    special_hfld_sheets = ['Disp SP', 'WP', 'Disp HFLD']

    with pd.ExcelFile(input_excel_file) as xl:
        for sheet_name in xl.sheet_names:
            if supersystem_label in sheet_name:
                new_sheet_name = sheet_name.replace(supersystem_label, '').strip()
                df_supersystem = xl.parse(sheet_name, index_col=0)

                df_sum_subsystems = pd.DataFrame(0, index=df_supersystem.index, columns=df_supersystem.columns)

                for subsystem_label in subsystem_labels:
                    subsystem_sheet_name = sheet_name.replace(str(supersystem_label), str(subsystem_label))
                    if subsystem_sheet_name in xl.sheet_names:
                        df_subsystem = xl.parse(subsystem_sheet_name, index_col=0)

                        df_subsystem = df_subsystem[~df_subsystem.index.duplicated(keep='first')]
                        df_subsystem = df_subsystem.loc[:, ~df_subsystem.columns.duplicated(keep='first')]

                        df_sum_subsystems += df_subsystem

                is_hfld = method.lower() == "hfld"
                is_disp_related = any(sheet_name.startswith(prefix + ' ') or sheet_name == prefix for prefix in special_hfld_sheets)

                if is_hfld and is_disp_related:
                    # Always use main_subsystem_matching_labels for zeroing intra-subsystem terms,
                    # regardless of whether BSSE is present or not.
                    for sublist in main_subsystem_matching_labels:
                        for j in range(len(sublist)):
                            for k in range(j + 1, len(sublist)):
                                row = sublist[j] - 1
                                col = sublist[k] - 1
                                df_supersystem.iat[row, col] = 0
                                df_supersystem.iat[col, row] = 0
                
                    df_result = df_supersystem * conversion_factor
                else:
                    df_result = (df_supersystem - df_sum_subsystems) * conversion_factor
                
                # Keep only upper triangle
                df_result = df_result.where(np.triu(np.ones(df_result.shape), k=0).astype(bool))
                matrices[new_sheet_name] = df_result

    # Handle Dispersion matrices for DLPNO-CCSD(T), DLPNO-CCSD, and HFLD methods
    if 'Disp SP' in matrices:
        df_disp_wp = matrices.get('WP', pd.DataFrame(np.nan, index=matrices['SP'].index, columns=matrices['SP'].columns)).copy()
        np.fill_diagonal(df_disp_wp.values, np.nan)
        matrices['Disp WP'] = df_disp_wp

        df_disp_sp = matrices['Disp SP']
        df_disp_ccsd = df_disp_sp + df_disp_wp
        
        if method.lower() == "hfld":
            matrices['Disp HFLD'] = df_disp_ccsd
            matrices['C-HFLD'] = df_disp_ccsd

        if method.lower() == "dlpno-ccsd":
            matrices['Disp CCSD'] = df_disp_ccsd

    if method.lower() == "dlpno-ccsd(t)" and 'T' in matrices:
        df_disp_t = matrices['T']
        df_disp_sp = matrices['Disp SP']
        df_disp_wp = matrices['Disp WP']
        
        with np.errstate(divide='ignore', invalid='ignore'):
            df_disp_t = np.divide(df_disp_sp * df_disp_t, matrices['SP'])
            df_disp_t[np.isnan(df_disp_t)] = 0

        df_disp_t = df_disp_t.where(np.triu(np.ones(df_disp_t.shape), k=0).astype(bool))
        matrices['Disp T'] = df_disp_t
        
        df_disp_ccsd_t = df_disp_sp + df_disp_wp + df_disp_t
        matrices['Disp CCSD(T)'] = df_disp_ccsd_t

        if 'Singles' in matrices:
            df_singles = matrices['Singles']
            df_c_ccsd_t = matrices['SP'] + matrices['WP'] + matrices['T'] + df_singles
            matrices['C-CCSD(T)'] = df_c_ccsd_t
        else:
            df_c_ccsd_t = matrices['SP'] + matrices['WP'] + matrices['T']
            matrices['C-CCSD(T)'] = df_c_ccsd_t

        df_inter_non_disp_ccsd_t = matrices['C-CCSD(T)'] - matrices['Disp CCSD(T)']
        matrices['Inter-NonDisp-C-CCSD(T)'] = df_inter_non_disp_ccsd_t

        if 'REF' in matrices:
            df_total = matrices['REF'] + matrices['C-CCSD(T)']
            matrices['TOTAL'] = df_total

    if method.lower() == "dlpno-ccsd" and 'Singles' in matrices:
        df_singles = matrices['Singles']
        df_c_ccsd = matrices['SP'] + matrices['WP'] + df_singles
        matrices['C-CCSD'] = df_c_ccsd

        df_inter_non_disp_ccsd = matrices['C-CCSD'] - matrices['Disp CCSD']
        matrices['Inter-NonDisp-C-CCSD'] = df_inter_non_disp_ccsd

        if 'REF' in matrices:
            df_total = matrices['REF'] + matrices['C-CCSD']
            matrices['TOTAL'] = df_total

    if method.lower() == "hfld" and 'C-HFLD' in matrices:
        if 'REF' in matrices:
            df_total = matrices['REF'].copy()
            df_c_hfld = matrices['C-HFLD']
            
            for i in range(df_c_hfld.shape[0]):
                for j in range(df_c_hfld.shape[1]):
                    if not np.isnan(df_c_hfld.iat[i, j]):
                        df_total.iat[i, j] += df_c_hfld.iat[i, j]
            
            matrices['TOTAL'] = df_total

    # Compute the total SOLV interaction energy
    solv_values, solv_int_energy, _, _, _ = compute_bulk_solvation_contribution(main_filenames, alternative_filenames, conversion_factor, method)

    # Add the SOLV matrix and update TOTAL if applicable
    if 'TOTAL' in matrices and solv_int_energy != 0:
        total_matrix = matrices['TOTAL']
        total_sum_wo_solv = total_matrix.sum().sum()

        if total_sum_wo_solv != 0:
            solv_matrix = total_matrix * (solv_int_energy / total_sum_wo_solv)
            matrices['SOLV'] = solv_matrix
            matrices['TOTAL'] = total_matrix + solv_matrix

    # Write the matrices from the dictionary to the output Excel file
    with pd.ExcelWriter(output_excel_file, engine='openpyxl') as writer:
        if 'TOTAL' in matrices:
            matrices['TOTAL'].to_excel(writer, sheet_name='TOTAL')
        if 'SOLV' in matrices:
            matrices['SOLV'].to_excel(writer, sheet_name='SOLV')

        # Write all other matrices
        for sheet_name, df in matrices.items():
            if sheet_name not in ['TOTAL', 'SOLV']:
                df.to_excel(writer, sheet_name=sheet_name)

    normalized_output_file = normalize_path(output_excel_file)
    print(f"All standard '{method}/LED' interaction energy matrices were written to '{normalized_output_file}'")


def relabel_and_sort_fragments(relabel_mapping, LEDAW_output_path, method):
    """Relabels and sorts the LED matrices to change fragment labels and ensures symmetry with diagonal elements of 'Electrostat' and 'Exchange' set to None."""

    normalized_LEDAW_output_path = normalize_path(LEDAW_output_path)

    input_excel_file = os.path.join(normalized_LEDAW_output_path, 'Unrelabeled_All_Standard_LED_matrices4.xlsx')
    output_excel_file = os.path.join(normalized_LEDAW_output_path, 'Uncleaned_All_Standard_LED_matrices5.xlsx')

    # If relabel_mapping is empty or None, just copy the file
    if not relabel_mapping:
        shutil.copy(input_excel_file, output_excel_file)
        print(f"No relabeling applied. The file was copied as '{normalize_path(output_excel_file)}' without changes.")
        return

    # Safely open Excel file for reading and writing using context managers
    with pd.ExcelFile(input_excel_file) as xl:
        with pd.ExcelWriter(output_excel_file, engine='openpyxl', mode='w') as writer:
            for sheet_name in xl.sheet_names:
                df = xl.parse(sheet_name, index_col=0)

                # Convert to integers
                df.columns = [int(c) for c in df.columns]
                df.index = [int(i) for i in df.index]

                # Apply relabeling using mapping
                df.columns = [df.columns[i - 1] for i in relabel_mapping]
                df.index = [df.index[i - 1] for i in relabel_mapping]

                # Sort rows and columns
                df.sort_index(axis=0, inplace=True)
                df.sort_index(axis=1, inplace=True)

                # Set diagonal of Electrostat and Exchange to None
                if sheet_name in ['Electrostat', 'Exchange']:
                    np.fill_diagonal(df.values, None)

                # Ensure symmetry: copy upper to lower triangle
                for i in range(df.shape[0]):
                    for j in range(i):
                        if pd.isna(df.iloc[j, i]) and not pd.isna(df.iloc[i, j]):
                            df.iloc[j, i] = df.iloc[i, j]
                            df.iloc[i, j] = np.nan

                # Write updated sheet
                df.to_excel(writer, sheet_name=sheet_name)

    print(f"All standard '{method}/LED' interaction energy matrices after relabeling and sorting fragments were written to '{normalize_path(output_excel_file)}' without cleaning redundant fragments")


def clean_redundant_frags(LEDAW_output_path, method):
    """Removes zero-valued rows and columns from all standard LED sheets taking the upper triangle of the TOTAL sheet as reference.
       If all diagonal elements in a sheet are exactly zero, they are also replaced with NaN.
    """
    
    normalized_LEDAW_output_path = normalize_path(LEDAW_output_path)
    input_excel_file = os.path.join(normalized_LEDAW_output_path, 'Uncleaned_All_Standard_LED_matrices5.xlsx')
    output_excel_file = os.path.join(normalized_LEDAW_output_path, 'All_Standard_LED_matrices.xlsx')

    xl = pd.ExcelFile(input_excel_file)
    if 'TOTAL' not in xl.sheet_names:
        raise ValueError("Sheet 'TOTAL' not found in input Excel file.")

    df_total = xl.parse('TOTAL', index_col=0)

    # Ensure symmetry and float dtype
    df_total = df_total.astype(float)
    df_total = df_total.where(np.triu(np.ones(df_total.shape)).astype(bool))

    # Identify redundant rows/columns (all zeros in upper triangle including diagonal)
    mask_upper = np.triu(np.ones(df_total.shape), k=0).astype(bool)
    redundant_labels = df_total.columns[(df_total.where(mask_upper).fillna(0) == 0).all()].tolist()

    print(f"Redundant labels corresponding to no actual fragment: {redundant_labels}")

    # Save cleaned sheets to new Excel
    with pd.ExcelWriter(output_excel_file, engine='openpyxl') as writer:
        for sheet_name in xl.sheet_names:
            df = xl.parse(sheet_name, index_col=0)
            if not df.empty:
                # Drop redundant fragments
                df_cleaned = df.drop(index=redundant_labels, columns=redundant_labels, errors='ignore')

                # Replace full zero diagonal with NaN
                diag_values = [df_cleaned.iat[i, i] for i in range(min(df_cleaned.shape[0], df_cleaned.shape[1]))]
                if all(v == 0 for v in diag_values):
                    for i in range(len(diag_values)):
                        df_cleaned.iat[i, i] = np.nan

                df_cleaned.to_excel(writer, sheet_name=sheet_name)
            else:
                df.to_excel(writer, sheet_name=sheet_name)

    normalized_output_file = normalize_path(output_excel_file)
    print(f"All standard '{method}/LED' interaction energy matrices after cleaning redundant labels were written to '{normalized_output_file}'")


def extract_final_fragment_labels_from_summary(LEDAW_output_path_nbody=None):
    """Returns the final fragment labels from the 'REF' sheet in the Summary_Standard_LED_matrices.xlsx file."""

    if not LEDAW_output_path_nbody:
        return None

    normalized_path = normalize_path(LEDAW_output_path_nbody)
    summary_file = os.path.join(normalized_path, 'Summary_Standard_LED_matrices.xlsx')

    if not os.path.isfile(summary_file):
        return None

    try:
        df_ref = pd.read_excel(summary_file, sheet_name='REF', index_col=0)
        return [int(label) for label in df_ref.index]
    except Exception:
        return None


def write_standard_LED_summary_int_en_matrices(method, LEDAW_output_path):
    """Write summary standard LED interaction energy maps to an Excel file, preserving index formatting."""

    # Normalize the LEDAW_output_path
    normalized_LEDAW_output_path = normalize_path(LEDAW_output_path)

    input_excel_file = os.path.join(normalized_LEDAW_output_path, 'All_Standard_LED_matrices.xlsx')
    output_excel_file = os.path.join(normalized_LEDAW_output_path, 'Summary_Standard_LED_matrices.xlsx')

    # Define the sheets to be written
    if method.lower() == "dlpno-ccsd(t)":
        sheets_to_transfer = [
            'TOTAL', 'SOLV', 'REF', 'Electrostat', 'Exchange',
            'C-CCSD(T)', 'Disp CCSD(T)', 'Inter-NonDisp-C-CCSD(T)'
        ]
    elif method.lower() == "dlpno-ccsd":
        sheets_to_transfer = [
            'TOTAL', 'SOLV', 'REF', 'Electrostat', 'Exchange',
            'C-CCSD', 'Disp CCSD', 'Inter-NonDisp-C-CCSD'
        ]
    elif method.lower() == "hfld":
        sheets_to_transfer = [
            'TOTAL', 'SOLV', 'REF', 'Electrostat', 'Exchange',
            'Disp HFLD'
        ]
    else:
        raise ValueError(f"Unsupported method: {method}. Please specify 'DLPNO-CCSD(T)', 'DLPNO-CCSD', or 'HFLD'.")

    # Write selected sheets using pandas to preserve index appearance in Excel
    with pd.ExcelFile(input_excel_file) as xl, pd.ExcelWriter(output_excel_file, engine='openpyxl') as writer:
        for sheet_name in sheets_to_transfer:
            if sheet_name in xl.sheet_names:
                df = xl.parse(sheet_name, index_col=0)
                df.to_excel(writer, sheet_name=sheet_name)
            else:
                print(f"Warning: Sheet '{sheet_name}' not found in {input_excel_file}.")

    print(f"The summary standard '{method}/LED' matrices were written to '{normalize_path(output_excel_file)}'")


def compute_fp_el_prep(df_ref):
    """Compute fragment pairwise electronic preparation matrices."""

    diagonal_elements = np.diag(df_ref.values)
    distributed_matrix = np.zeros(df_ref.shape)

    # Calculate the non-diagonal sums using the updated denominator function
    non_diagonal_sums = np.array([compute_denominator_for_fp_el_prep(df_ref, i) for i in range(df_ref.shape[0])])

    for i in range(df_ref.shape[0]):
        for j in range(i + 1, df_ref.shape[1]):  # Only consider the upper triangle (i < j)
            if non_diagonal_sums[i] != 0 and non_diagonal_sums[j] != 0:
                term_1 = (diagonal_elements[i] * abs(df_ref.iloc[i, j])) / non_diagonal_sums[i]
                term_2 = (diagonal_elements[j] * abs(df_ref.iloc[i, j])) / non_diagonal_sums[j]
                distributed_value = term_1 + term_2
                distributed_matrix[i, j] = distributed_value

    # Convert to DataFrame and set diagonal and below-diagonal elements to NaN
    distributed_df = pd.DataFrame(distributed_matrix, index=df_ref.index, columns=df_ref.columns)
    np.fill_diagonal(distributed_df.values, np.nan)
    distributed_df = distributed_df.where(np.triu(np.ones(distributed_df.shape), k=1).astype(bool))

    return distributed_df


def compute_denominator_for_fp_el_prep(df, index):
    """Compute the sum of absolute non-diagonal elements involving a particular index (both row and column)."""
    row_sum = abs(df.iloc[index, :]).sum() - abs(df.iloc[index, index])
    col_sum = abs(df.iloc[:, index]).sum() - abs(df.iloc[index, index])
    return row_sum + col_sum


def process_fp_LED_matrices(method, main_filenames, alternative_filenames, conversion_factor, LEDAW_output_path):
    """Process the matrices as per the given method and write to the output Excel file."""

    # Normalize the LEDAW_output_path
    normalized_LEDAW_output_path = normalize_path(LEDAW_output_path)

    # Define file paths with normalized paths
    input_excel_file = os.path.join(normalized_LEDAW_output_path, 'Summary_Standard_LED_matrices.xlsx')
    output_excel_file = os.path.join(normalized_LEDAW_output_path, 'Summary_fp-LED_matrices.xlsx')

    # Load the input Excel file
    xl = pd.ExcelFile(input_excel_file)

    # Process REF sheet
    df_ref = pd.read_excel(xl, sheet_name='REF', index_col=0)
    df_ref.columns = [int(c) for c in df_ref.columns]
    df_ref.index = [int(i) for i in df_ref.index]

    df_ref_distributed = compute_fp_el_prep(df_ref)
    df_ref_distributed.columns.name = 'REF-EL-PREP'

    # Get SOLV values and calculate SOLV interaction energy
    solv_values, e_solv_int_en, _, _, _ = compute_bulk_solvation_contribution(main_filenames, alternative_filenames, conversion_factor, method)

    # Initialize a new Excel writer object
    with pd.ExcelWriter(output_excel_file, engine='openpyxl') as writer:

        # Load Electrostat and Exchange sheets for all methods
        df_electrostat = pd.read_excel(xl, sheet_name='Electrostat', index_col=0)
        df_exchange = pd.read_excel(xl, sheet_name='Exchange', index_col=0)
        df_electrostat.columns = [int(c) for c in df_electrostat.columns]
        df_electrostat.index = [int(i) for i in df_electrostat.index]
        df_exchange.columns = [int(c) for c in df_exchange.columns]
        df_exchange.index = [int(i) for i in df_exchange.index]

        if method.lower() in ["dlpno-ccsd(t)", "dlpno-ccsd"]:
            # Read TOTAL
            df_total = pd.read_excel(xl, sheet_name='TOTAL', index_col=0)
            df_total.columns = [int(c) for c in df_total.columns]
            df_total.index = [int(i) for i in df_total.index]

            # Subtract SOLV to get standard LED total without SOLV
            try:
                df_solv = pd.read_excel(xl, sheet_name='SOLV', index_col=0)
                df_solv.columns = [int(c) for c in df_solv.columns]
                df_solv.index = [int(i) for i in df_solv.index]
                std_df_total_no_solv = df_total - df_solv
            except ValueError:
                std_df_total_no_solv = df_total

            # Correlation EL-PREP on standard total without SOLV
            df_total_distributed = compute_fp_el_prep(std_df_total_no_solv)
            df_total_distributed.columns.name = 'C-CCSD(T)-EL-PREP' if method.lower() == "dlpno-ccsd(t)" else 'C-CCSD-EL-PREP'

            df_c_ccsd = df_total_distributed - df_ref_distributed
            df_c_ccsd.columns.name = 'CCSD(T)-EL-PREP' if method.lower() == "dlpno-ccsd(t)" else 'CCSD-EL-PREP'

            # DISP - INTER-NODISP decomposition
            df_disp_ccsd = pd.read_excel(xl, sheet_name='Disp CCSD' if method.lower() == "dlpno-ccsd" else 'Disp CCSD(T)', index_col=0)
            df_inter_nondisp = pd.read_excel(xl, sheet_name='Inter-NonDisp-C-CCSD' if method.lower() == "dlpno-ccsd" else 'Inter-NonDisp-C-CCSD(T)', index_col=0)
            df_disp_ccsd.columns = [int(c) for c in df_disp_ccsd.columns]
            df_disp_ccsd.index = [int(i) for i in df_disp_ccsd.index]
            df_inter_nondisp.columns = [int(c) for c in df_inter_nondisp.columns]
            df_inter_nondisp.index = [int(i) for i in df_inter_nondisp.index]

            # Reconstruct total without SOLV from fp-LED logic
            df_ref_final = df_electrostat + df_exchange + df_ref_distributed
            df_c_ccsd_final = df_c_ccsd + df_disp_ccsd + df_inter_nondisp
            df_total_no_solv = df_ref_final + df_c_ccsd_final

            # Apply fp-LED SOLV if needed
            if e_solv_int_en != 0:
                solv_matrix = df_total_no_solv * (e_solv_int_en / df_total_no_solv.sum().sum())
                df_total_final = df_total_no_solv + solv_matrix
                df_total_final.to_excel(writer, sheet_name='TOTAL')
                solv_matrix.to_excel(writer, sheet_name='SOLV')
            else:
                df_total_final = df_total_no_solv
                df_total_final.to_excel(writer, sheet_name='TOTAL')

            # Write all matrices
            df_ref_final.to_excel(writer, sheet_name='REF')
            df_electrostat.to_excel(writer, sheet_name='Electrostat')
            df_exchange.to_excel(writer, sheet_name='Exchange')
            df_ref_distributed.to_excel(writer, sheet_name='REF-EL-PREP')
            df_c_ccsd_final.to_excel(writer, sheet_name='C-CCSD(T)' if method.lower() == "dlpno-ccsd(t)" else 'C-CCSD')
            df_disp_ccsd.to_excel(writer, sheet_name='Disp CCSD(T)' if method.lower() == "dlpno-ccsd(t)" else 'Disp CCSD')
            df_inter_nondisp.to_excel(writer, sheet_name='Inter-NonDisp-C-CCSD(T)' if method.lower() == "dlpno-ccsd(t)" else 'Inter-NonDisp-C-CCSD')
            df_c_ccsd.to_excel(writer, sheet_name='C-CCSD(T)-EL-PREP' if method.lower() == "dlpno-ccsd(t)" else 'C-CCSD-EL-PREP')
            df_total_distributed.to_excel(writer, sheet_name='CCSD(T)-EL-PREP' if method.lower() == "dlpno-ccsd(t)" else 'CCSD-EL-PREP')

        elif method.lower() == "hfld":
            # Calculate REF = Electrostat + Exchange + REF-EL-PREP
            df_ref_final = df_electrostat + df_exchange + df_ref_distributed

            # Calculate TOTAL = REF + Disp HFLD
            df_disp_hfld = pd.read_excel(xl, sheet_name='Disp HFLD', index_col=0)
            df_disp_hfld.columns = [int(c) for c in df_disp_hfld.columns]
            df_disp_hfld.index = [int(i) for i in df_disp_hfld.index]

            df_total_hfld = df_ref_final + df_disp_hfld
            df_total_hfld.columns.name = 'TOTAL'

            # If SOLV interaction energy is not zero, calculate and add the SOLV matrix
            if e_solv_int_en != 0:
                solv_matrix = df_total_hfld * (e_solv_int_en / df_total_hfld.sum().sum())
                df_total_hfld += solv_matrix
                df_total_hfld.to_excel(writer, sheet_name='TOTAL')
                solv_matrix.to_excel(writer, sheet_name='SOLV')
            else:
                df_total_hfld.to_excel(writer, sheet_name='TOTAL')

            # Write the remaining HFLD matrices
            df_ref_final.to_excel(writer, sheet_name='REF')
            df_electrostat.to_excel(writer, sheet_name='Electrostat')
            df_exchange.to_excel(writer, sheet_name='Exchange')
            df_ref_distributed.to_excel(writer, sheet_name='REF-EL-PREP')
            df_disp_hfld.to_excel(writer, sheet_name='Disp HFLD')

        normalized_output_file = normalize_path(output_excel_file)
        print(f"fp-LED interaction energy matrices were written to '{normalized_output_file}'")


def write_nbody_solv_files(LEDAW_output_path, total_solv_int_energy, ref_diel_int_energy, ref_cds_int_energy, corr_diel_int_energy):
    """
    Writes SOLV-fp.xlsx and SOLV-STD.xlsx for the individual componennts of SOLV contribution.
    Reads base SOLV matrices from Summary_Standard_LED_matrices.xlsx and Summary_fp-LED_matrices.xlsx
    and scales them with appropriate factors for computing individual components.
    Applies specific NaN masking based on FP (diagonal and below NaN) or STD (below diagonal NaN) schemes.
    """
    normalized_path = normalize_path(LEDAW_output_path)
    solv_fp_output_file = os.path.join(normalized_path, 'SOLV-fp.xlsx')
    solv_std_output_file = os.path.join(normalized_path, 'SOLV-STD.xlsx')
    summary_std_input_file = os.path.join(normalized_path, 'Summary_Standard_LED_matrices.xlsx')
    summary_fp_input_file = os.path.join(normalized_path, 'Summary_fp-LED_matrices.xlsx')

    os.makedirs(normalized_path, exist_ok=True)

    # Helper function to calculate scaled matrix, handles division by zero and empty matrices
    def get_scaled_matrix(base_matrix, component_int_energy, total_int_energy):
        if base_matrix is None or base_matrix.empty:
            return pd.DataFrame() # Return an empty DataFrame if base is not valid
        if total_int_energy == 0:
            # If total_int_energy is 0, all components are 0.
            # Create a zero matrix of the same shape as base_matrix.
            return pd.DataFrame(np.zeros(base_matrix.shape), index=base_matrix.index, columns=base_matrix.columns)
        
        scale_factor = component_int_energy / total_int_energy
        return base_matrix * scale_factor

    # Helper to apply the masking (below diagonal NaN)
    def mask_below_diagonal(df):
        if df.empty:
            return df
        # Create a boolean mask for the upper triangle (k=0 includes diagonal)
        mask = np.triu(np.ones(df.shape), k=0).astype(bool)
        return df.where(mask)

    # Helper to apply the masking (diagonal and below NaN)
    def mask_diagonal_and_below(df):
        if df.empty:
            return df
        # Create a boolean mask for the upper triangle (k=1 excludes diagonal)
        mask = np.triu(np.ones(df.shape), k=1).astype(bool)
        return df.where(mask)

    # --- Process and Write SOLV-fp.xlsx ---
    try:
        # Read the base SOLV matrix from Summary_fp-LED_matrices.xlsx
        # Use .fillna(0.0) to treat missing values as zero for scaling
        base_solv_fp_df = pd.read_excel(summary_fp_input_file, sheet_name='SOLV', index_col=0).fillna(0.0)
        # Ensure numeric column/index names if they were read as strings, keeping original if not numeric
        base_solv_fp_df.columns = [int(c) if str(c).isdigit() else c for c in base_solv_fp_df.columns]
        base_solv_fp_df.index = [int(i) if str(i).isdigit() else i for i in base_solv_fp_df.index]

        with pd.ExcelWriter(solv_fp_output_file, engine='openpyxl', mode='w') as writer:
            # Write the base SOLV matrix for FP scheme, with diagonal and below as NaN
            masked_solv_fp_df = mask_diagonal_and_below(base_solv_fp_df)
            masked_solv_fp_df.to_excel(writer, sheet_name='SOLV')

            # Calculate and write individual components for FP scheme, with diagonal and below as NaN
            ref_diel_fp_df = get_scaled_matrix(base_solv_fp_df, ref_diel_int_energy, total_solv_int_energy)
            masked_ref_diel_fp_df = mask_diagonal_and_below(ref_diel_fp_df)
            masked_ref_diel_fp_df.to_excel(writer, sheet_name='REF-DIEL')

            ref_cds_fp_df = get_scaled_matrix(base_solv_fp_df, ref_cds_int_energy, total_solv_int_energy)
            masked_ref_cds_fp_df = mask_diagonal_and_below(ref_cds_fp_df)
            masked_ref_cds_fp_df.to_excel(writer, sheet_name='REF-CDS')

            corr_diel_fp_df = get_scaled_matrix(base_solv_fp_df, corr_diel_int_energy, total_solv_int_energy)
            masked_corr_diel_fp_df = mask_diagonal_and_below(corr_diel_fp_df)
            masked_corr_diel_fp_df.to_excel(writer, sheet_name='CORR-DIEL')

        print(f"SOLV-fp.xlsx created successfully at '{normalize_path(solv_fp_output_file)}'")

    except FileNotFoundError:
        print(f"Warning: '{os.path.basename(summary_fp_input_file)}' not found. Skipping SOLV-fp.xlsx creation.")
    except KeyError:
        print(f"Warning: 'SOLV' sheet not found in '{os.path.basename(summary_fp_input_file)}'. Skipping SOLV-fp.xlsx creation.")
    except Exception as e:
        print(f"An error occurred while creating SOLV-fp.xlsx: {e}")


    # --- Process and Write SOLV-STD.xlsx ---
    try:
        # Read the base SOLV matrix from Summary_Standard_LED_matrices.xlsx
        # Use .fillna(0.0) to treat missing values as zero for scaling
        base_solv_std_df = pd.read_excel(summary_std_input_file, sheet_name='SOLV', index_col=0).fillna(0.0)
        # Ensure numeric column/index names if they were read as strings
        base_solv_std_df.columns = [int(c) if str(c).isdigit() else c for c in base_solv_std_df.columns]
        base_solv_std_df.index = [int(i) if str(i).isdigit() else i for i in base_solv_std_df.index]

        with pd.ExcelWriter(solv_std_output_file, engine='openpyxl', mode='w') as writer:
            # Write the base SOLV matrix for STD scheme, with below diagonal as NaN
            masked_solv_std_df = mask_below_diagonal(base_solv_std_df)
            masked_solv_std_df.to_excel(writer, sheet_name='SOLV')

            # Calculate and write individual components for STD scheme, with below diagonal as NaN
            ref_diel_std_df = get_scaled_matrix(base_solv_std_df, ref_diel_int_energy, total_solv_int_energy)
            masked_ref_diel_std_df = mask_below_diagonal(ref_diel_std_df)
            masked_ref_diel_std_df.to_excel(writer, sheet_name='REF-DIEL')

            ref_cds_std_df = get_scaled_matrix(base_solv_std_df, ref_cds_int_energy, total_solv_int_energy)
            masked_ref_cds_std_df = mask_below_diagonal(ref_cds_std_df)
            masked_ref_cds_std_df.to_excel(writer, sheet_name='REF-CDS')

            corr_diel_std_df = get_scaled_matrix(base_solv_std_df, corr_diel_int_energy, total_solv_int_energy)
            masked_corr_diel_std_df = mask_below_diagonal(corr_diel_std_df)
            masked_corr_diel_std_df.to_excel(writer, sheet_name='CORR-DIEL')

        print(f"SOLV-STD.xlsx created successfully at '{normalize_path(solv_std_output_file)}'")

    except FileNotFoundError:
        print(f"Warning: '{os.path.basename(summary_std_input_file)}' not found. Skipping SOLV-STD.xlsx creation.")
    except KeyError:
        print(f"Warning: 'SOLV' sheet not found in '{os.path.basename(summary_std_input_file)}'. Skipping SOLV-STD.xlsx creation.")
    except Exception as e:
        print(f"An error occurred while creating SOLV-STD.xlsx: {e}")


def delete_unprocessed_files(LEDAW_output_path):
    """Delete temporary excel files"""

    normalized_LEDAW_output_path = normalize_path(LEDAW_output_path)
    pattern = os.path.join(normalized_LEDAW_output_path, 'Un*')
    pattern = normalize_path(pattern)
    files_to_delete = glob.glob(pattern)

    print("Temporary files written will be deleted:")
    for file_path in files_to_delete:
        normalized_file_path = normalize_path(file_path)
        try:
            gc.collect()
            os.remove(normalized_file_path)
            print(f"Deleted: {normalized_file_path}")
        except PermissionError as e:
            print(f"PermissionError: Could not delete {normalized_file_path} — file may be open.")
            if platform.system() == 'Windows':
                try:
                    for proc in psutil.process_iter(['pid', 'open_files']):
                        for file in proc.info.get('open_files') or []:
                            if file.path == os.path.abspath(normalized_file_path):
                                print(f"File is open in process {proc.pid}")
                except Exception as ex:
                    print(f"Could not check open processes: {ex}")
            else:
                print("Skipping process check: OS allows deleting open files.")
        except Exception as ex:
            print(f"Unexpected error deleting {normalized_file_path}: {ex}")


def engine_LED_N_body(main_filenames, alternative_filenames, conversion_factor, method, LEDAW_output_path, relabel_mapping=None, use_ref_as_rhf_in_hfld=None):
    """Engine function to process LED N-body interaction energy matrices, standardize and reorder labels, compute final matrices, and provide summaries."""

    normalized_LEDAW_output_path = normalize_path(LEDAW_output_path)

    # Step 1: Generate labeled filenames
    labeled_main_filenames = label_systems(main_filenames)
    labeled_alt_filenames = label_systems(alternative_filenames)

    system_labels = list(labeled_main_filenames.values())
    alternative_labels = list(labeled_alt_filenames.values())

    # Step 2: Process multi-fragment files (generates tmp1.xlsx)
    multifrag_system_processing(main_filenames, alternative_filenames, normalized_LEDAW_output_path, system_labels)

    # Step 3: Process single-fragment files (generates tmp2.xlsx)
    singlefrag_system_processing(labeled_main_filenames=labeled_main_filenames, labeled_alt_filenames=labeled_alt_filenames,
        LEDAW_output_path=normalized_LEDAW_output_path, method=method, use_ref_as_rhf_in_hfld=use_ref_as_rhf_in_hfld)

    # Step 4: Combine tmp1.xlsx and tmp2.xlsx into Unprocessed_LED_matrices.xlsx
    combine_unprocessed_LED_data_fies(LEDAW_output_path=normalized_LEDAW_output_path)

    # Step 5: Construct label mappings
    main_label_mappings_ghost_free, alt_label_mappings_ghost_free, bsse_found = construct_label_mappings(main_filenames, alternative_filenames, LEDAW_output_path)

    # Step 6: Construct label list of each subsystem alligning them to the main SUPERSYS labels, excluding ghost fragments.
    main_subsystem_matching_labels = subsystem_label_lists_alligned_to_main_supersystem(system_labels, main_label_mappings_ghost_free)

    # Step 7: Unify fragment labels
    unify_labels(system_labels=system_labels, alternative_labels=alternative_labels, main_label_mappings_ghost_free=main_label_mappings_ghost_free,
        alternative_label_mappings_ghost_free=alt_label_mappings_ghost_free, LEDAW_output_path=normalized_LEDAW_output_path)

    # Step 8: Remove redundant ALT sheets
    combine_main_ALT(LEDAW_output_path=normalized_LEDAW_output_path)

    # Step 9: Combine Intra and Inter matrices
    combine_intra_inter_matrices(LEDAW_output_path=normalized_LEDAW_output_path)

    # Step 10: Compute standard LED interaction matrices
    compute_all_standard_led_int_en_matrices(system_labels=system_labels, conversion_factor=conversion_factor, method=method,
        LEDAW_output_path=normalized_LEDAW_output_path, main_subsystem_matching_labels=main_subsystem_matching_labels,
        main_filenames=main_filenames, alternative_filenames=alternative_filenames, bsse_found=bsse_found)

    # Step 11: Relabel and sort fragments (if mapping is provided)
    relabel_and_sort_fragments(relabel_mapping=relabel_mapping, LEDAW_output_path=normalized_LEDAW_output_path, method=method)

    # Step 12: Clean redundant fragments from the interaction energy maps
    clean_redundant_frags(LEDAW_output_path=normalized_LEDAW_output_path, method=method)

    # Step 13: Write standard summary LED matrices
    write_standard_LED_summary_int_en_matrices(method=method, LEDAW_output_path=normalized_LEDAW_output_path)

    # Step 14: Process and write fp-LED matrices
    process_fp_LED_matrices(method=method, main_filenames=main_filenames, alternative_filenames=alternative_filenames,
        conversion_factor=conversion_factor, LEDAW_output_path=normalized_LEDAW_output_path)

    # Step 15: Compute SOLV components and write them to SOLV-fp.xlsx and SOLV-STD.xlsx
    solv_values, total_solv_int_energy, ref_diel_int_energy, ref_cds_int_energy, corr_diel_int_energy = \
        compute_bulk_solvation_contribution(main_filenames, alternative_filenames, conversion_factor, method)
    write_nbody_solv_files(normalized_LEDAW_output_path, total_solv_int_energy,
                           ref_diel_int_energy, ref_cds_int_energy, corr_diel_int_energy)

    # Step 16: Delete intermediate excel files
    delete_unprocessed_files(LEDAW_output_path)

    print('\n')
    print('*'*120)
    print(f"  N-body LED analyses were terminated NORMALLY. Standard and fp-LED matrices are at {normalized_LEDAW_output_path}")
    print('*'*120)
