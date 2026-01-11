import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt

mesh = pv.read("base_surface_flow.vtu")

# Ensure triangles
mesh = mesh.triangulate()

# Convert point → cell data (as you already did)
mesh = mesh.point_data_to_cell_data()

# Extract surface: UnstructuredGrid → PolyData
surf = mesh.extract_surface()

# Now compute normals (this WILL work)
surf = surf.compute_normals(
    cell_normals=True,
    point_normals=False,
    auto_orient_normals=True,
    inplace=False,
)

# Use surf from here on
centers = surf.cell_centers().points
areas   = surf.compute_cell_sizes(length=False, area=True).cell_data["Area"]
normals = surf.cell_data["Normals"]
Cp      = surf.cell_data["Pressure_Coefficient"]

# Lift direction (global Z)
lift_dir = np.array([0.0, 0.0, 1.0])

# Differential force per cell (non-dimensional)
dF = -Cp[:, None] * areas[:, None] * normals
dL = dF @ lift_dir

# Spanwise coordinate (Y)
y = centers[:, 1]

# Bin spanwise
nbins = 50
y_bins = np.linspace(y.min(), y.max(), nbins + 1)
y_mid  = 0.5 * (y_bins[:-1] + y_bins[1:])

L_span = np.zeros(nbins)
for i in range(nbins):
    mask = (y >= y_bins[i]) & (y < y_bins[i+1])
    L_span[i] = dL[mask].sum()

# Normalize (AVL-like)
L_total = L_span.sum()
cl_span = L_span / L_total

# Plot
plt.figure()
plt.plot(y_mid, cl_span, "-o")
plt.xlabel("Spanwise coordinate y [m]")
plt.ylabel("Normalized lift distribution")
plt.title("Spanwise Lift Distribution (SU2)")
plt.grid(True)
plt.show()
