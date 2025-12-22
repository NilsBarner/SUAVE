# Minimal example to test ComputeCFDMesh (see trunk\SUAVE\Input_Output\OpenVSP\write_vsp_mesh.py)

# NILS: add OpenVSP python binding to path (add to system path long-term)
import sys
sys.path.insert(0, r"C:\Users\nmb48\Documents\GitHub\SUAVE\OpenVSP-3.46.0-win64-Python3.9\OpenVSP-3.46.0-win64\python\openvsp")

import openvsp as vsp

vsp.ClearVSPModel()
vsp.ReadVSPFile("base.vsp3")

# NILS: fails silently when below lines are uncommented
# Requires further debugging!
# vsp.SetCFDMeshVal(vsp.CFD_HALF_MESH_FLAG,1)
# vsp.SetCFDMeshVal(vsp.CFD_FAR_FIELD_FLAG,1)

vsp.SetComputationFileName(vsp.CFD_STL_TYPE, "test.stl")

vsp.ComputeCFDMesh(
    vsp.SET_ALL,      # surfaces
    vsp.SET_NONE,     # degenset
    vsp.CFD_STL_TYPE  # output
)

print("DONE")

"""
# NILS: with two lines commented
(suave) C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\BWB_CFD>python "C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\BWB_CFD\test.py"
Invalid reason 17
Invalid reason 17
DONE

# NILS: with two lines uncommented (silent failure)
(suave) C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\BWB_CFD>python "C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\BWB_CFD\test.py"

# NILS: current state of running "C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\BWB_CFD\BWB.py"
(suave) C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\BWB_CFD>python "C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\BWB_CFD\BWB.py"
C:\Users\nmb48\Documents\GitHub\SUAVE\trunk\SUAVE\Plugins\pint\__init__.py:17: UserWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html. The pkg_resources package is slated for removal as early as 2025-11-30. Refrain from using this package or pin to Setuptools<81.
  import pkg_resources
Reseting OpenVSP Model in Memory
Writing main_wing to OpenVSP Model
Writing nacelle to OpenVSP Model
Writing nacelle_2 to OpenVSP Model
Writing nacelle_3 to OpenVSP Model
Saving OpenVSP File at C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\BWB_CFD/base.vsp3
Starting mesh for base (This may take several minutes)
"""


