import os
import glob
import importlib

# Get the directory where this __init__.py is located
module_dir = os.path.dirname(__file__)

# List all the Python files in this directory (except __init__.py)
module_files = glob.glob(os.path.join(module_dir, "*.py"))

# Import all functions from all modules
for module_file in module_files:
    module_name = os.path.basename(module_file)[:-3]  # remove ".py" from the name
    if module_name != "__init__":
        module = importlib.import_module(f'.{module_name}', package=__name__)
        globals().update({name: getattr(module, name) for name in dir(module) if not name.startswith('_')})