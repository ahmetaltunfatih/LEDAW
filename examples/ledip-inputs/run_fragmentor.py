import sys
sys.path.insert(0,"../..")
from ledip_package import *


xyzfile = "wt.xyz"

cutoff = 1.62

cut_bonds = [
    [0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [10, 11], 
    [12, 13], [14, 15], [16, 17], [18, 19], [20, 21], 
    [22, 23], [24, 25], [26, 27]
]

nonstandard_coordination_numbers = {
#    "Bi": 4,
#    "Fe": 0,
#    "7-9,11,13": 3
}


try:
    fragmentation_engine(xyzfile, cutoff=cutoff, cut_bonds=cut_bonds,
                         nonstandard_coordination_numbers=nonstandard_coordination_numbers)
    print("Fragmentation engine completed successfully!\n")
except Exception as e:
    print(f"Error: {e}")
    print("Fragmentation engine failed to complete.\n")


if sys.platform.startswith('win'):
    os.system("pause")
else:
    input("Press Enter to exit...") # For Linux/macOS