import math
import os
import re


def read_xyz(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    num_atoms = int(lines[0])
    comment = lines[1].strip()
    atoms = []
    for i in range(num_atoms):
        parts = lines[i + 2].split()
        symbol = parts[0]
        x, y, z = map(float, parts[1:4])
        atoms.append((symbol, x, y, z))
    return atoms, comment


def write_xyz(filename, atoms, comment):
    with open(filename, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{comment}\n")
        for (symbol, frag_id, x, y, z) in atoms:
            f.write(f"{symbol}({frag_id}) {x:.9f} {y:.9f} {z:.9f}\n")


def distance(a1, a2):
    return math.sqrt(sum((c1 - c2) ** 2 for c1, c2 in zip(a1, a2)))


def detect_fragments(atoms, cutoff=1.8, cut_bonds=None):
    bonded = [[] for _ in atoms]
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            if distance(atoms[i][1:], atoms[j][1:]) <= cutoff:
                bonded[i].append(j)
                bonded[j].append(i)

    if cut_bonds:
        for (a, b) in cut_bonds:
            if b in bonded[a]:
                bonded[a].remove(b)
            if a in bonded[b]:
                bonded[b].remove(a)

    visited = [False] * len(atoms)
    fragments = []

    def dfs(i, current_frag):
        visited[i] = True
        current_frag.append(i)
        for j in bonded[i]:
            if not visited[j]:
                dfs(j, current_frag)

    for i in range(len(atoms)):
        if not visited[i]:
            frag = []
            dfs(i, frag)
            fragments.append(frag)

    return fragments


def centroid(fragment, atoms):
    coords = [atoms[i][1:] for i in fragment]
    x = sum(c[0] for c in coords) / len(coords)
    y = sum(c[1] for c in coords) / len(coords)
    z = sum(c[2] for c in coords) / len(coords)
    return (x, y, z)


class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, i):
        if self.parent[i] != i:
            self.parent[i] = self.find(self.parent[i])
        return self.parent[i]

    def union(self, i, j):
        pi = self.find(i)
        pj = self.find(j)
        if pi != pj:
            self.parent[pi] = pj


def parse_index_keys(d):
    index_based = {}
    symbol_based = {}

    for key, val in d.items():
        # Convert val to int if possible, otherwise keep as is
        try:
            val_int = int(val)
        except Exception:
            val_int = val  # fallback: keep original if cannot convert

        if any(c.isdigit() or c in "-," for c in str(key)):
            indices = set()
            parts = str(key).split(',')
            for part in parts:
                part = part.strip()
                if '-' in part:
                    start, end = part.split('-')
                    indices.update(range(int(start), int(end) + 1))
                else:
                    indices.add(int(part))
            for idx in indices:
                index_based[idx] = val_int
        else:
            symbol_based[key] = val_int

    return index_based, symbol_based


def get_n_target(atom_index_0based, symbol, index_based, symbol_based):
    if atom_index_0based in index_based:
        return index_based[atom_index_0based]
    if symbol in symbol_based:
        return symbol_based[symbol]
    return None


def merge_fragments(fragments, atoms, nonstandard_coordination_numbers):
    num_frags = len(fragments)
    frag_centroids = [centroid(f, atoms) for f in fragments]

    index_based, symbol_based = parse_index_keys(nonstandard_coordination_numbers or {})

    uf = UnionFind(num_frags)
    special_atoms = []

    for atom_idx0, (sym, x, y, z) in enumerate(atoms):
        n_target = get_n_target(atom_idx0, sym, index_based, symbol_based)
        if n_target is None:
            continue

        if n_target == 0:
            special_atoms.append([atom_idx0])
            continue

        n_target = min(n_target, num_frags)
        dists = [(i, distance((x, y, z), c)) for i, c in enumerate(frag_centroids)]
        dists.sort(key=lambda t: t[1])
        closest_frag_ids = [i for i, _ in dists[:n_target]]

        for i in range(len(closest_frag_ids) - 1):
            uf.union(closest_frag_ids[i], closest_frag_ids[i + 1])

        special_atoms.append([atom_idx0] + closest_frag_ids)

    root_to_frag = {}
    for i in range(num_frags):
        root = uf.find(i)
        if root not in root_to_frag:
            root_to_frag[root] = []
        root_to_frag[root].extend(fragments[i])

    for group in special_atoms:
        atom_idx0 = group[0]
        if len(group) == 1:
            new_root = max(root_to_frag.keys(), default=-1) + 1
            root_to_frag[new_root] = [atom_idx0]
        else:
            target_root = uf.find(group[1])
            root_to_frag[target_root].append(atom_idx0)

    return list(root_to_frag.values())


def fragmentation_engine(xyzfile, cutoff=1.80, cut_bonds=None, nonstandard_coordination_numbers=None):

    atoms_raw, original_comment = read_xyz(xyzfile) # atoms_raw is [(symbol_str_from_file, x,y,z), ...]

    if not atoms_raw:
        print("Input XYZ file is empty.")
        return

    # --- Part 1: Parse input, check if all labeled, gather original label info ---
    # Stores dicts: {'elem': str, 'x':f, 'y':f, 'z':f, 'original_label': int or None, 'original_index': int}
    parsed_atom_details = []
    
    # original_label_value -> list of atom_indices (original_index from atoms_raw)
    atoms_by_original_label = {}
    # Original labels in order of their first appearance in the input file
    ordered_unique_original_labels = []
    all_atoms_fully_pre_labeled = True # Assume true initially

    for i, (symbol_in_file, x, y, z) in enumerate(atoms_raw):
        match = re.match(r"([A-Za-z]+)\s*(?:\(\s*(\d+)\s*\))?", symbol_in_file)
        current_atom_original_label = None
        element_symbol = ""

        if match:
            element_symbol = match.group(1)
            label_str = match.group(2)
            if label_str is not None and label_str.strip(): # Ensure label string is not empty
                try:
                    current_atom_original_label = int(label_str)
                    if current_atom_original_label not in atoms_by_original_label:
                        atoms_by_original_label[current_atom_original_label] = []
                        ordered_unique_original_labels.append(current_atom_original_label)
                    atoms_by_original_label[current_atom_original_label].append(i)
                except ValueError: # Handle non-integer labels if necessary, or count as unlabeled
                    print(f"Warning: Atom {i} has non-integer label '{label_str}'. Treating as unlabeled.")
                    all_atoms_fully_pre_labeled = False
            else: # No label group, or label group is empty
                all_atoms_fully_pre_labeled = False
        else:
            # If regex doesn't match, it might be just an element symbol or unparsable
            element_symbol = symbol_in_file # Assume it's an element for now
            all_atoms_fully_pre_labeled = False
            # Consider raising an error if element_symbol is not valid, e.g. by checking against a periodic table

        parsed_atom_details.append({
            'elem': element_symbol, 'x': x, 'y': y, 'z': z,
            'original_label': current_atom_original_label, 'original_index': i
        })
    
    # Further check: if all_atoms_fully_pre_labeled is true, atoms_by_original_label must not be empty
    # (unless atoms_raw was empty, which is checked above)
    if all_atoms_fully_pre_labeled and not atoms_by_original_label:
        # This could happen if all atoms had "Elem()" with empty parens, effectively no valid labels
        all_atoms_fully_pre_labeled = False 
        print("Note: Atoms appeared to have label syntax but no valid labels found; using standard fragmentation.")


    # --- Part 2: Conditional logic based on pre-labeling state ---
    if all_atoms_fully_pre_labeled:
        print("All atoms are pre-labeled. Applying sequential re-labeling based on input order.")

        # Map original labels to new sequential labels (1, 2, 3...)
        # based on the order of appearance of original labels.
        original_to_new_label_map = {
            original_label_val: new_idx + 1
            for new_idx, original_label_val in enumerate(ordered_unique_original_labels)
        }

        # Prepare atoms for output: list of (elem, new_label, x, y, z, original_index)
        output_atom_tuples = []
        for atom_data in parsed_atom_details:
            # All atoms have an original_label in this branch
            original_label = atom_data['original_label']
            new_sequential_label = original_to_new_label_map[original_label]
            output_atom_tuples.append((
                atom_data['elem'], new_sequential_label,
                atom_data['x'], atom_data['y'], atom_data['z'],
                atom_data['original_index'] # Keep original_index for secondary sort
            ))

        # Sort by the new sequential label, then by original index to maintain intra-fragment order
        output_atom_tuples.sort(key=lambda atom_tuple: (atom_tuple[1], atom_tuple[5]))
        
        # Final list of atoms for writing: (elem, new_label, x, y, z)
        final_atoms_for_writing = [
            (elem, lbl, x, y, z) for elem, lbl, x, y, z, _orig_idx in output_atom_tuples
        ]

        outname = os.path.splitext(xyzfile)[0] + "_relabeled.xyz"
        num_final_fragments = len(ordered_unique_original_labels)
        
        write_xyz(outname, final_atoms_for_writing, 
                  f"{num_final_fragments} fragments (sequentially relabeled from input) in {xyzfile}")

        print(f"{len(parsed_atom_details)} atoms processed.")
        print(f"{num_final_fragments} fragments identified from input and relabeled sequentially: {list(original_to_new_label_map.values())}")
        print(f"Output written to {outname}")
        print("Note: `cutoff`, `cut_bonds`, `nonstandard_coordination_numbers` were ignored as all atoms were pre-labeled.")
        return # End of this specific path

    # --- Else branch: Fallback to original logic for partially or un-labeled files ---
    else:
        if not all_atoms_fully_pre_labeled: # This will be true if we are in this else block due to partial/no labeling
            print("Unlabeled atoms will be automatically fragmented. Existing fragment labels (if any) will be preserved.\n")

        # 1. Segregate atoms: identify pre-labeled fragments and collect unlabeled atoms.
        #    `parsed_atom_details` is already populated from the initial parsing step:
        #    Each element is: {'elem': str, 'x':f, 'y':f, 'z':f, 'original_label': int or None, 'original_index': int}
        
        pre_labeled_fragments_data = {} # original_label -> list of atom_data_dicts for that fragment
        unlabeled_atom_global_indices = []    # list of original_atom_indices that were not labeled in input
        
        existing_original_label_values = set() # To avoid label collisions later

        for atom_data in parsed_atom_details:
            original_label = atom_data['original_label']
            atom_idx = atom_data['original_index']

            if original_label is not None:
                if original_label not in pre_labeled_fragments_data:
                    pre_labeled_fragments_data[original_label] = []
                pre_labeled_fragments_data[original_label].append(atom_data)
                existing_original_label_values.add(original_label)
            else:
                unlabeled_atom_global_indices.append(atom_idx)

        # This list will hold all atoms with their final labels for sorting and output
        # Each element: (elem_symbol, final_label, x, y, z, original_global_index)
        all_atoms_final_tuples = []

        # Add atoms from pre-labeled fragments directly to the final list.
        # These fragments keep their original labels and atom compositions.
        for original_label, atoms_in_fragment_list in pre_labeled_fragments_data.items():
            for atom_data in atoms_in_fragment_list:
                all_atoms_final_tuples.append((
                    atom_data['elem'],
                    original_label, # Use the original label
                    atom_data['x'], atom_data['y'], atom_data['z'],
                    atom_data['original_index']
                ))
            print(f"Preserved pre-labeled fragment {original_label} with {len(atoms_in_fragment_list)} atoms.")

        # 2. Process unlabeled atoms if any exist.
        if unlabeled_atom_global_indices:
            print(f"Processing {len(unlabeled_atom_global_indices)} unlabeled atoms...\n")
            
            # Create a temporary list of atom data specifically for the unlabeled subset.
            # This list will have new, 0-based indices relative to itself.
            unlabeled_subset_atom_data = [] # list of (elem_symbol, x, y, z)
            # Map: index_in_unlabeled_subset -> original_global_atom_index
            unlabeled_subset_idx_to_global_idx = {} 
            # Map: original_global_atom_index -> index_in_unlabeled_subset
            global_idx_to_unlabeled_subset_idx = {}
            
            for new_subset_idx, global_idx in enumerate(unlabeled_atom_global_indices):
                atom_data = parsed_atom_details[global_idx] # Get full data from original parsing
                unlabeled_subset_atom_data.append((atom_data['elem'], atom_data['x'], atom_data['y'], atom_data['z']))
                unlabeled_subset_idx_to_global_idx[new_subset_idx] = global_idx
                global_idx_to_unlabeled_subset_idx[global_idx] = new_subset_idx

            # --- Remap `cut_bonds` and `nonstandard_coordination_numbers` for the unlabeled subset ---
            # Only consider bonds/NCNs that solely involve atoms within the current unlabeled subset.
            
            remapped_cut_bonds_for_subset = None
            if cut_bonds:
                temp_remapped_bonds = []
                for a_glob, b_glob in cut_bonds:
                    if a_glob in global_idx_to_unlabeled_subset_idx and b_glob in global_idx_to_unlabeled_subset_idx:
                        a_sub = global_idx_to_unlabeled_subset_idx[a_glob]
                        b_sub = global_idx_to_unlabeled_subset_idx[b_glob]
                        temp_remapped_bonds.append(sorted((a_sub, b_sub)))
                if temp_remapped_bonds:
                    remapped_cut_bonds_for_subset = temp_remapped_bonds
            # print(f"Remapped cut_bonds for unlabeled subset: {remapped_cut_bonds_for_subset}")

            remapped_ncn_for_subset = None
            if nonstandard_coordination_numbers:
                temp_remapped_ncn = {}
                idx_based_orig_nc, sym_based_orig_nc = parse_index_keys(nonstandard_coordination_numbers)
                # Add symbol-based NCNs directly (they apply universally by element type)
                for sym, val in sym_based_orig_nc.items():
                    temp_remapped_ncn[sym] = val
                # Remap index-based NCNs if the global index is in the unlabeled set
                for global_idx, n_target in idx_based_orig_nc.items(): # global_idx is an int
                    if global_idx in global_idx_to_unlabeled_subset_idx:
                        subset_idx = global_idx_to_unlabeled_subset_idx[global_idx]
                        temp_remapped_ncn[subset_idx] = n_target # Key is now subset_idx
                if temp_remapped_ncn:
                    remapped_ncn_for_subset = temp_remapped_ncn
            # print(f"Remapped NCN for unlabeled subset: {remapped_ncn_for_subset}")

            # Perform geometric fragmentation ONLY on the `unlabeled_subset_atom_data`.
            base_frags_from_unlabeled = detect_fragments(
                unlabeled_subset_atom_data, 
                cutoff=cutoff, 
                cut_bonds=remapped_cut_bonds_for_subset
            )
            
            final_frags_from_unlabeled_subset_indices = merge_fragments(
                base_frags_from_unlabeled, 
                unlabeled_subset_atom_data, 
                remapped_ncn_for_subset
            )
            
            print(f"{len(final_frags_from_unlabeled_subset_indices)} new fragments generated from unlabeled atoms.\n")

            # Assign new, unique labels to these newly generated fragments.
            # New labels start from 1 and avoid collision with `existing_original_label_values`.
            _next_new_label_for_unlabeled_frags = 1 
            for frag_made_of_subset_indices in final_frags_from_unlabeled_subset_indices:
                # Find an available new label
                while _next_new_label_for_unlabeled_frags in existing_original_label_values:
                    _next_new_label_for_unlabeled_frags += 1
                current_new_label = _next_new_label_for_unlabeled_frags
                _next_new_label_for_unlabeled_frags += 1 # Prepare for the next new fragment

                # Add atoms of this new fragment to the global list `all_atoms_final_tuples`
                for subset_idx in frag_made_of_subset_indices:
                    global_idx = unlabeled_subset_idx_to_global_idx[subset_idx]
                    atom_data = parsed_atom_details[global_idx] # Get original full atom data
                    all_atoms_final_tuples.append((
                        atom_data['elem'],
                        current_new_label, # Assign the new label
                        atom_data['x'], atom_data['y'], atom_data['z'],
                        global_idx
                    ))
        else:
            print("No unlabeled atoms to process for new fragmentation.")

        # 3. Sort all atoms (pre-labeled and newly labeled from unlabeled) for final output.
        all_atoms_final_tuples.sort(key=lambda atom_tuple: (atom_tuple[1], atom_tuple[5])) # Sort by label, then by original_idx

        # 4. Update original `cut_bonds` and `nonstandard_coordination_numbers` to reflect the new global atom ordering in `all_atoms_final_tuples`.
        old_global_idx_to_new_sorted_idx = {
            atom_tuple[5]: new_idx for new_idx, atom_tuple in enumerate(all_atoms_final_tuples)
        }

        updated_cut_bonds_for_output = []
        if cut_bonds: # Use original `cut_bonds`
            for a_glob, b_glob in cut_bonds:
                if a_glob in old_global_idx_to_new_sorted_idx and b_glob in old_global_idx_to_new_sorted_idx:
                    new_a = old_global_idx_to_new_sorted_idx[a_glob]
                    new_b = old_global_idx_to_new_sorted_idx[b_glob]
                    updated_cut_bonds_for_output.append(sorted((new_a, new_b)))
        print("Updated cut_bonds (mapped to new final atom indices):")
        print(updated_cut_bonds_for_output,"\n")

        updated_ncn_for_output = {}
        if nonstandard_coordination_numbers: # Use original `nonstandard_coordination_numbers`
            idx_based_orig_nc, sym_based_orig_nc = parse_index_keys(nonstandard_coordination_numbers)
            for global_idx, n_target in idx_based_orig_nc.items():
                if global_idx in old_global_idx_to_new_sorted_idx:
                    new_sorted_idx = old_global_idx_to_new_sorted_idx[global_idx]
                    updated_ncn_for_output[str(new_sorted_idx)] = n_target # Store with string key for consistency
            for sym, n_target in sym_based_orig_nc.items():
                updated_ncn_for_output[sym] = n_target
        print("Updated nonstandard_coordination_numbers (keys mapped to new final atom indices):")
        print(updated_ncn_for_output,"\n")

        # 5. Write the final XYZ file.
        outname_partial_preserved = os.path.splitext(xyzfile)[0] + "_labeled.xyz"
        
        final_atom_count = len(all_atoms_final_tuples)
        num_output_fragments = len(set(t[1] for t in all_atoms_final_tuples)) # Count unique final labels
        
        final_comment_text = f"{num_output_fragments} fragments from {xyzfile}"
        
        # Prepare list of (symbol, label, x, y, z) for write_xyz
        atoms_for_writing = [(s,l,x,y,z) for s,l,x,y,z,_orig_idx in all_atoms_final_tuples]
        write_xyz(outname_partial_preserved, atoms_for_writing, final_comment_text)
        
        print(f"{final_atom_count} atoms processed.\n")
        print(f"Output written to {outname_partial_preserved}\n")
