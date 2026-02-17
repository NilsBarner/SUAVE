"""
This scipt demos the use of the ways to visualise the geometry and the solution in OptVL.
It must be run from within an \examples folder, which is at the same directory level as
an \airfoils folder. Check out \OptVL\examples\plot_aircraft.py,
\OptVL\examples\plot_airfoils.py, and \OptVL\examples\plot_sectional_data.py.
"""

import numpy as np
import matplotlib.pyplot as plt
from optvl import OVLSolver

# from matplotlib_custom_settings import *

# Load .avl file
# ovl_solver = OVLSolver(geo_file=r"C:\Users\nmb48\Documents\GitHub\SUAVE\avl3.52\AVL3.52rel09032025\runs\b737.avl", debug=False)  # test case
# ovl_solver = OVLSolver(geo_file=r"C:\Users\nmb48\avl_files_6\vehicle.avl", debug=False)  # my designs
ovl_solver = OVLSolver(geo_file=r"C:\Users\nmb48\avl_files_6\vehicle.jvl", debug=False)  # my designs

# Uncomment if want to plot lift distribution
# ovl_solver.set_variable("alpha", 5.0)  # set AOA
# ovl_solver.set_variable("beta", 0.0)  # set AOS
ovl_solver.execute_run()

# Produce figure
# ovl_solver.plot_geom_nils(colors=[colors[0], colors[1]])
ovl_solver.plot_geom()  # original version

ovl_solver.set_variable("alpha", 5.00)
ovl_solver.execute_run()

# ovl_solver.plot_cp()
# ovl_solver.plot_cp_nils()

# =============================================================================
import sys
sys.exit()
# =============================================================================

#%% ..\..\..\OptVL\examples\plot_airfoils.py

# """This scipt demos looking at the airfoil data that is used by AVL"""
r"""
from optvl import OVLSolver
import matplotlib.pyplot as plt

# ovl = OVLSolver(geo_file="aircraft.avl", debug=True)
ovl = OVLSolver(geo_file=r"C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\B737_AVL_Tutorial\aircraft.avl", debug=True)
# ovl = OVLSolver(geo_file=r"C:\Users\nmb48\avl_files_6\vehicle.avl", debug=True)
surf_data = ovl.get_surface_params(include_geom=True, include_airfoils=True)

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
idx_color = 0

for surf_key in ["Wing", "Horizontal Tail"]:
# for surf_key in ["main_wing"]:
    # x coorindates for airfoil
    airfoil_coords = surf_data[surf_key]["airfoils"]
    
    # leading edge points
    yles = surf_data[surf_key]["yles"]
    xles = surf_data[surf_key]["xles"]

    # airfoil files
    afiles = surf_data[surf_key]["afiles"]

    for idx_airfoil in range(len(airfoil_coords)):
        x_offset = xles[idx_airfoil]
        y_offset = yles[idx_airfoil]
        coords = airfoil_coords[idx_airfoil]
        
        label = afiles[idx_airfoil]
        plt.plot(x_offset + coords[0,:], y_offset + coords[1,:], color=colors[idx_color], label=label)
        idx_color += 1

plt.axis("equal")
plt.legend()
plt.show()

# geo_file="Boeing_737-800.avl"
# ["main_wing", "horizontal_stabilizer"]
"""
#%% ..\..\..\OptVL\examples\plot_sectional_data.py

from optvl import OVLSolver
import numpy as np
import matplotlib.pyplot as plt

# ovl = OVLSolver(geo_file="aircraft.avl", debug=False)
# ovl = OVLSolver(geo_file=r"C:\Users\nmb48\avl_files_6\vehicle.avl", debug=False)
ovl = OVLSolver(geo_file=r"C:\Users\nmb48\avl_files_6\vehicle.jvl", debug=False)
ovl.set_variable("alpha", 5.0)
ovl.set_variable("beta", 0.0)
ovl.execute_run()

# keys-start
strip_data = ovl.get_strip_forces()
first_surf = list(strip_data.keys())[0]
print(strip_data[first_surf].keys())
# keys-end

for surf_key in strip_data:
    span_distance = strip_data[surf_key]["Y LE"]
    # plt.plot(span_distance, strip_data[surf_key]["lift dist"], color="red")
    plt.plot(span_distance, strip_data[surf_key]["CL"], color="red")
    # plt.plot(span_distance, strip_data[surf_key]["CL perp"], color="firebrick", linestyle="--")
    # plt.plot(span_distance, strip_data[surf_key]["spanloading"], color="orange")
    
    # plt.plot(span_distance, strip_data[surf_key]["CD"], color="blue")
    # plt.plot(span_distance, strip_data[surf_key]["drag dist"], color="blue")
    
    # plt.plot(span_distance, strip_data[surf_key]["Cm"], color="green")
    # plt.plot(span_distance, strip_data[surf_key]["Cl"], color="green")
    # plt.plot(span_distance, strip_data[surf_key]["Cn"], color="green")
    
# plt.legend(["roll distribution", "yaw distribution"])
plt.legend(["lift dist", "CL", "CL perp."])
# plt.title("roll and yaw spanwise data")
plt.title("lift spanwise data")
plt.xlabel("spanwise position")
plt.show()

