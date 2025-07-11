class Patterns:
    """A class to hold all regex patterns used in the script."""
    def __init__(self):
        self.local_energy_decomp_pattern = r"LOCAL ENERGY DECOMPOSITION"
        self.pattern_e0 = r"E\(0\)\s+\.+\s+([-]?\d*\.\d+)"
        self.pattern_strong_corr = r"E\(CORR\)\(strong-pairs\)\s+\.+\s+([-]?\d*\.\d+)"
        self.pattern_weak_corr = r"E\(CORR\)\(weak-pairs\)\s+\.+\s+([-]?\d*\.\d+)"
        self.pattern_triples_corr = r"Triples Correction \(T\)\s+\.+\s+([-]?\d*\.\d+)"
        self.PATTERNS = {
            "intra_ref": r"Intra REF\. energy\s+([-]?\d*\.\d+)(\s+[-]?\d*\.\d+)*\s+sum=",
            "intra_ref_alt": (
                r"INTRA-FRAGMENT REF\. ENERGY FOR FRAGMENT\s+(\d+)" # Group 1: Fragment number
                r"[\s\S]*?" # Non-greedy match for anything between the header and "Total energy"
                r"Total energy\s+=\s+([-]?\d*\.\d+)" # Group 2: Total energy
            ),
            "intra_corr": r"Intra Correlation energy\s+([-]?\d*\.\d+)(\s+[-]?\d*\.\d+)*\s+sum=",
            "intra_strong_pairs": r"Intra strong pairs\s+([-]?\d*\.\d+)(\s+[-]?\d*\.\d+)*\s+sum=",
            "intra_triples": r"Intra triples\s+([-]?\d*\.\d+)(\s+[-]?\d*\.\d+)*\s+sum=",
            "intra_weak_pairs": r"Intra weak pairs\s+([-]?\d*\.\d+)(\s+[-]?\d*\.\d+)*\s+sum=",
            "singles_contribution": r"Singles contribution\s+([-]?\d*\.\d+)(\s+[-]?\d*\.\d+)*\s+sum=",
            "dispersion_strong_pairs": r"Dispersion\s+(\d+),(\d+)\s+([-]?\d*\.\d+)",
            "ref_inter": (
                r"INTER-FRAGMENT REF. ENERGY FOR FRAGMENTs\s+(\d+)\s+AND\s+(\d+)\s+[-]+\s+"
                r"Nuclear repulsion\s+=\s+[-]?\d*\.\d+\s+"
                r"Nuclear attraction\s+=\s+[-]?\d*\.\d+\s+"
                r"Coulomb repulsion\s+=\s+[-]?\d*\.\d+\s+"
                r"[-]+\s+"
                r"Sum of electrostatics\s+=\s+([-]?\d*\.\d+)\s+\(.*?\)\s+"
                r"Two electron exchange\s+=\s+([-]?\d*\.\d+)\s+\(.*?\)\s+"
                r"[-]+\s+"
                r"Total REF. interaction\s+=\s+([-]?\d*\.\d+)\s+\(.*?\)"
            ),
            "ref_corr_inter": (
                r"Interaction of Fragments\s+(\d+)\s+and\s+(\d+):\s+[-]+\s+"
                r"Interfragment reference\s+([-]?\d*\.\d+)\s+\(.*?\)\s+"
                r"Interfragment correlation\s+([-]?\d*\.\d+)\s+\(.*?\)"
            ),
            "corr_inter_full": (
                r"Interaction correlation for Fragments\s+(\d+)\s+and\s+(\d+):\s+[-]+\s+"
                r"Inter strong pairs\s+([-]?\d*\.\d+)\s+\(.*?\)\s+"
                r"Inter triples\s+([-]?\d*\.\d+)\s+\(.*?\)\s+"
                r"Inter weak pairs\s+([-]?\d*\.\d+)\s+\(.*?\)\s+"
                r"[-]+\s+"
            ),
            "corr_inter_partial": (
                r"Interaction correlation for Fragments\s+(\d+)\s+and\s+(\d+):\s+[-]+\s+"
                r"Inter strong pairs\s+([-]?\d*\.\d+)\s+\(.*?\)\s+"
                r"Inter weak pairs\s+([-]?\d*\.\d+)\s+\(.*?\)\s+"
                r"[-]+\s+"
            )
        }


class MatrixDifferenceError(Exception):
    """Custom exception for handling significant matrix differences."""
    pass
