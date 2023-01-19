import os

# path to python
path_to_python = "python"

# path to uq_pc.py
path_to_uq_pc = os.path.join(os.environ.get("UQTK_INS"), "examples", "uqpc", "uq_pc.py")

# path to the dfi app
path_to_dfi = os.path.join(os.environ.get("UQTK_INS"), "bin", "dfi")

# settings to save figures if desired
save = False
save_location = "."