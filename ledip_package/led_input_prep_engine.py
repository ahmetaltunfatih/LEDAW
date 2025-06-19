import os
import shutil
import re
from collections import Counter


def apply_conditional_heading_removals(heading_content):
    """
    Removes keywords and blocks from an ORCA input heading involving a single real fragment.
    For single fragment cases:
    - The 'LED' keyword is removed.
    - 'DoLEDHF' directive is removed from the mdci block.
    - If the method is 'HFLD':
        - 'HFLD', 'LoosePNO', 'NormalPNO', and 'TightPNO' keywords are removed.
        - The %mdci block is removed.
	Some stylistic modifications are also done after removal of these keywords and blocks.
    """

    lines = heading_content.splitlines()

    base_removal_keywords = re.compile(r'\b(LED)\b', re.IGNORECASE)
    hfld_related_keywords = re.compile(r'\b(HFLD|LoosePNO|NormalPNO|TightPNO)\b', re.IGNORECASE)
    doledhf_regex = re.compile(r'\bDoLEDHF\s+(True|False)\b', re.IGNORECASE)

    hfld_present = bool(re.search(r'\bHFLD\b', heading_content, re.IGNORECASE))

    processed_lines = []
    i = 0
    n = len(lines)

    while i < n:
        line = lines[i]

        # -- Handle %mdci block --
        if re.match(r'^\s*%\s*mdci\b', line, re.IGNORECASE):
            block_lines = []
            block_start_index = len(processed_lines)
            while i < n:
                block_lines.append(lines[i])
                if re.match(r'^\s*end\s*$', lines[i], re.IGNORECASE):
                    break
                i += 1

            if hfld_present:
                i += 1  # skip 'end'
                # clean up extra empty lines above
                while block_start_index > 0 and processed_lines and not processed_lines[-1].strip():
                    processed_lines.pop()
                    block_start_index -= 1
                # add a single empty line if needed
                if processed_lines and (i < n and lines[i].strip()):
                    processed_lines.append('')
                continue

            cleaned_block = []
            has_substance = False
            doledhf_removed = False

            for bl in block_lines:
                bl_mod = doledhf_regex.sub('', bl)
                if bl_mod != bl:
                    doledhf_removed = True

                stripped = re.sub(r'\s{2,}', ' ', bl_mod).strip()
                if not stripped or stripped == '!':
                    continue

                if not re.match(r'^\s*%\s*mdci\b', bl, re.IGNORECASE) and not re.match(r'^\s*end\s*$', bl, re.IGNORECASE):
                    if not stripped.startswith(('!', '*')):
                        has_substance = True

                if stripped.startswith(('!', '%', '*')):
                    cleaned_block.append(stripped)
                else:
                    leading = bl[:len(bl) - len(bl.lstrip())]
                    cleaned_block.append(leading + stripped)

            if doledhf_removed and not has_substance:
                i += 1  # skip past 'end'
                while block_start_index > 0 and processed_lines and not processed_lines[-1].strip():
                    processed_lines.pop()
                    block_start_index -= 1
                if processed_lines and (i < n and lines[i].strip()):
                    processed_lines.append('')
                continue
            else:
                processed_lines.extend(cleaned_block)
                i += 1
                continue

        # -- Handle outside %mdci block --
        line_mod = doledhf_regex.sub('', line)
        if hfld_present:
            line_mod = hfld_related_keywords.sub('', line_mod)
        line_mod = base_removal_keywords.sub('', line_mod)

        stripped = re.sub(r'\s{2,}', ' ', line_mod).strip()

        # Remove if line is just '!'
        if stripped == '!':
            i += 1
            continue

        # Keep user-inserted empty lines
        if not stripped:
            processed_lines.append('')
        else:
            leading = line[:len(line) - len(line.lstrip())]
            processed_lines.append(leading + stripped)

        i += 1

    return "\n".join(processed_lines)


def led_input_prep_engine(xyzfile, subsystems, orca_inp_heading, output_dir):
    """
    Generate NBODY and TWOBODY input files for LED analysis based on a labeled XYZ file.

    Parameters:
    - xyzfile: path to labeled XYZ file (with atom labels like C(2) for fragment 2)
    - subsystems: list of subsystems, each a list of fragment indices (e.g. [[1,9],[2,3]]).
                  If empty ([]), each fragment from the XYZ file will be treated as a separate subsystem.
    - orca_inp_heading: ORCA input heading for all systems (string)
    - output_dir: target directory to write output files
    """
    # Clean headings: strip overall leading/trailing whitespace including empty lines
    orca_inp_heading_stripped = orca_inp_heading.strip()

    # Parse the labeled XYZ file (skip the first two lines)
    atoms = []
    label_pattern = re.compile(r'([A-Za-z]+)(?:\((\d+)\))?') # Regex to parse C(3) or C

    with open(xyzfile, 'r') as f:
        lines = f.readlines()
    if len(lines) < 3:
        raise ValueError("XYZ file must have at least 3 lines (atom count, comment, and coordinates).")
    
    for line in lines[2:]:
        parts = line.strip().split()
        if not parts:
            continue
        
        symbol_with_label = parts[0]
        match = label_pattern.fullmatch(symbol_with_label)
        if not match:
            raise ValueError(f"Could not parse element label: '{symbol_with_label}' in {xyzfile}")
        
        element = match.group(1)
        frag_index = int(match.group(2)) if match.group(2) else None 

        coords_str = parts[1:]
        if len(coords_str) != 3:
            raise ValueError(f"Invalid coordinate line in XYZ file: '{line.strip()}' - Expected 3 coordinates.")
        
        x, y, z = coords_str
        atoms.append({'element': element, 'frag': frag_index, 'x': x, 'y': y, 'z': z})

    # Group atoms by fragment index
    frags = {}
    for atom in atoms:
        frag = atom['frag']
        if frag is None:
            continue
        frags.setdefault(frag, []).append(atom)
    all_fragments = sorted(frags.keys())

    # --- Handle subsystems=[] and Validate subsystems ---
    if not subsystems: # If subsystems list is empty, treat each fragment as a separate subsystem
        subsystems_to_process = [[frag] for frag in all_fragments]
    else:
        subsystems_to_process = subsystems # Use the provided subsystems

    # Validate the subsystems_to_process list
    flattened_subsystems_frags = []
    for sublist in subsystems_to_process:
        for frag in sublist:
            flattened_subsystems_frags.append(frag)

    # Check 1: All fragment labels in supersystem are included
    missing_frags = set(all_fragments) - set(flattened_subsystems_frags)
    if missing_frags:
        raise ValueError(
            f"Error: The following fragments from the XYZ file are not included in the provided subsystems list: "
            f"{sorted(list(missing_frags))}. Please correct subsystem definitions and run again."
        )

    # Check 2: No fragments are listed in subsystems_to_process that are not in all_fragments (extra fragments)
    extra_frags = set(flattened_subsystems_frags) - set(all_fragments)
    if extra_frags:
        raise ValueError(
            f"Error: The following fragments are listed in subsystems but not found in the XYZ file: "
            f"{sorted(list(extra_frags))}. Please correct subsystem definitions and run again."
        )

    # Check 3: All fragments are included only once (no duplicates)
    if len(flattened_subsystems_frags) != len(set(flattened_subsystems_frags)):
        duplicates = [item for item, count in Counter(flattened_subsystems_frags).items() if count > 1]
        raise ValueError(
            f"Error: The following fragments are entered more than once in the subsystems list: "
            f"{sorted(duplicates)}. Please correct subsystem definitions and run again."
        )
    
    # --- End Subsystems Handling ---

    # Prepare NBODY and TWOBODY directories (remove existing if present)
    nbody_dir = os.path.join(output_dir, 'NBODY')
    twobody_dir = os.path.join(output_dir, 'TWOBODY')
    for d in [nbody_dir, twobody_dir]:
        if os.path.exists(d):
            shutil.rmtree(d)
        os.makedirs(d)

    # NBODY subdirectories
    subsys_woBSSE_dir = os.path.join(nbody_dir, 'SUBSYS_woBSSE')
    subsys_withBSSE_dir = os.path.join(nbody_dir, 'SUBSYS_withBSSE')
    os.makedirs(subsys_woBSSE_dir)
    os.makedirs(subsys_withBSSE_dir)

    # TWOBODY subdirectories
    dimer_dir = os.path.join(twobody_dir, 'DIMER')
    onebody_woBSSE_dir = os.path.join(twobody_dir, 'ONEBODY_woBSSE')
    onebody_withBSSE_dir = os.path.join(twobody_dir, 'ONEBODY_withBSSE') 
    os.makedirs(dimer_dir)
    os.makedirs(onebody_woBSSE_dir)
    os.makedirs(onebody_withBSSE_dir)

    # ---------------- NBODY: supersystem ----------------
    supersys_file = os.path.join(nbody_dir, 'supersys.inp')
    with open(supersys_file, 'w') as f:
        f.write(orca_inp_heading_stripped + '\n')
        for atom in atoms:
            el = atom['element']
            frag = atom['frag']
            x, y, z = atom['x'], atom['y'], atom['z']
            if frag is not None:
                f.write(f"{el}({frag}) {x} {y} {z}\n")
            else:
                f.write(f"{el} {x} {y} {z}\n")
        f.write('*\n')

    # ---------------- NBODY: SUBSYS_woBSSE ----------------
    for i, subs in enumerate(subsystems_to_process, start=1):
        filepath = os.path.join(subsys_woBSSE_dir, f"subsys{i}.inp")
        with open(filepath, 'w') as f:
            head_to_write = orca_inp_heading_stripped
            # Apply removals only if this subsystem consists of a single fragment
            if len(subs) == 1:
                head_to_write = apply_conditional_heading_removals(head_to_write)
            f.write(head_to_write + '\n')
            
            frag_map = {old: new for new, old in enumerate(sorted(subs), start=1)}
            single_frag_in_subsys = (len(subs) == 1)
            
            for old_frag in sorted(subs):
                new_id = frag_map[old_frag]
                for atom in frags.get(old_frag, []):
                    el = atom['element']
                    x, y, z = atom['x'], atom['y'], atom['z']
                    if single_frag_in_subsys:
                        f.write(f"{el} {x} {y} {z}\n")
                    else:
                        f.write(f"{el}({new_id}) {x} {y} {z}\n")
            f.write('*\n')

    # ---------------- NBODY: SUBSYS_withBSSE ----------------
    for i, subs in enumerate(subsystems_to_process, start=1):
        filepath = os.path.join(subsys_withBSSE_dir, f"subsys{i}.inp")
        with open(filepath, 'w') as f:
            head_to_write = orca_inp_heading_stripped
            # Apply removals only if this subsystem has a single real fragment (many ghost fragments may exist)
            if len(subs) == 1: # 'subs' here directly represents the real fragments
                head_to_write = apply_conditional_heading_removals(head_to_write)
            f.write(head_to_write + '\n')
            
            real_frags = sorted(subs)
            frag_map = {old: new for new, old in enumerate(real_frags, start=1)}
            
            dummy_frags = [frag for frag in all_fragments if frag not in real_frags]
            
            next_id_for_dummy = len(real_frags) + 1
            for old_frag in sorted(dummy_frags):
                frag_map[old_frag] = next_id_for_dummy
                next_id_for_dummy += 1
            
            for old_frag in real_frags:
                new_id = frag_map[old_frag]
                for atom in frags.get(old_frag, []):
                    el = atom['element']
                    x, y, z = atom['x'], atom['y'], atom['z']
                    f.write(f"{el}({new_id}) {x} {y} {z}\n")
            
            for old_frag in dummy_frags:
                new_id = frag_map[old_frag]
                for atom in frags.get(old_frag, []):
                    el = atom['element']
                    x, y, z = atom['x'], atom['y'], atom['z']
                    f.write(f"{el}:({new_id}) {x} {y} {z}\n")
            f.write('*\n')

    # ---------------- TWOBODY: DIMER ----------------
    num_sub = len(subsystems_to_process)
    for idx1 in range(num_sub):
        for idx2 in range(idx1+1, num_sub):
            # Extract all individual fragments from the two subsystems involved in this dimer
            current_subsystem_frags1 = subsystems_to_process[idx1]
            current_subsystem_frags2 = subsystems_to_process[idx2]

            # Iterate through all combinations of individual fragments from these two subsystems
            for frag1_id in current_subsystem_frags1:
                for frag2_id in current_subsystem_frags2:
                    fname = f"dimer{frag1_id}-{frag2_id}.inp"
                    filepath = os.path.join(dimer_dir, fname)
                    with open(filepath, 'w') as f:
                        f.write(orca_inp_heading_stripped + '\n')
                        # Write atoms for frag1_id (labeled 1)
                        for atom in frags.get(frag1_id, []):
                            el = atom['element']
                            x, y, z = atom['x'], atom['y'], atom['z']
                            f.write(f"{el}(1) {x} {y} {z}\n")
                        # Write atoms for frag2_id (labeled 2)
                        for atom in frags.get(frag2_id, []):
                            el = atom['element']
                            x, y, z = atom['x'], atom['y'], atom['z']
                            f.write(f"{el}(2) {x} {y} {z}\n")
                        f.write('*\n')

    # ---------------- TWOBODY: ONEBODY_woBSSE ----------------
    for frag_index in all_fragments:
        filepath = os.path.join(onebody_woBSSE_dir, f"mono{frag_index}.inp")
        with open(filepath, 'w') as f:
            head_to_write = apply_conditional_heading_removals(orca_inp_heading_stripped)
            f.write(head_to_write + '\n')

            for atom in frags.get(frag_index, []):
                el = atom['element']
                x, y, z = atom['x'], atom['y'], atom['z']
                f.write(f"{el} {x} {y} {z}\n")
            f.write('*\n')

    # ---------------- TWOBODY: ONEBODY_withBSSE ----------------
    monomer_heading_with_bsse_removed = apply_conditional_heading_removals(orca_inp_heading_stripped)

    num_sub = len(subsystems_to_process)
    for idx1 in range(num_sub):
        for idx2 in range(idx1+1, num_sub):
            # Extract all individual fragments from the two subsystems involved
            current_subsystem_frags1 = subsystems_to_process[idx1]
            current_subsystem_frags2 = subsystems_to_process[idx2]

            # Iterate through all combinations of individual fragments from these two subsystems
            for frag1_id in current_subsystem_frags1:
                for frag2_id in current_subsystem_frags2:
                    # File for frag1_id as real (labeled 1), frag2_id as dummy (labeled 2)
                    fname1 = f"dimer{frag1_id}-{frag2_id}_{frag1_id}.inp"
                    filepath1 = os.path.join(onebody_withBSSE_dir, fname1)
                    with open(filepath1, 'w') as f:
                        f.write(monomer_heading_with_bsse_removed + '\n')
                        for atom in frags.get(frag1_id, []):
                            el = atom['element']
                            x, y, z = atom['x'], atom['y'], atom['z']
                            f.write(f"{el}(1) {x} {y} {z}\n")
                        for atom in frags.get(frag2_id, []): # This is the dummy frag for frag2_id
                            el = atom['element']
                            x, y, z = atom['x'], atom['y'], atom['z']
                            f.write(f"{el}:(2) {x} {y} {z}\n") # Dummy fragment label is always (2)
                        f.write('*\n')

                    # File for frag2_id as real (labeled 1), frag1_id as dummy (labeled 2)
                    fname2 = f"dimer{frag1_id}-{frag2_id}_{frag2_id}.inp"
                    filepath2 = os.path.join(onebody_withBSSE_dir, fname2)
                    with open(filepath2, 'w') as f:
                        f.write(monomer_heading_with_bsse_removed + '\n')
                        for atom in frags.get(frag2_id, []):
                            el = atom['element']
                            x, y, z = atom['x'], atom['y'], atom['z']
                            f.write(f"{el}(1) {x} {y} {z}\n")
                        for atom in frags.get(frag1_id, []): # This is the dummy frag for frag1_id
                            el = atom['element']
                            x, y, z = atom['x'], atom['y'], atom['z']
                            f.write(f"{el}:(2) {x} {y} {z}\n") # Dummy fragment label is always (2)
                        f.write('*\n')

    print(f"LED input preparation complete. Files written to '{output_dir}'.\n")
