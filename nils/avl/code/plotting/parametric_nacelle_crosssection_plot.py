"""
This script demonstrates the parametrisation of a circular engine
cross-section.
"""

__all__ = []

import numpy as np
import matplotlib.pyplot as plt

radius = 0.7430612571316905

origin = [[10.591785098336818, 4.323656451350929, 1.0134693714341547]]

X0 = origin[0][0]     # use segment 0 origin X
Yc = origin[0][1]
Zc = origin[0][2]

# Angular locations (in percent span)
n_segments = 13
segment_percent_span = np.linspace(0, 1, num=n_segments,endpoint=False)

# Convert to angle around circle
thetas = 2.0 * np.pi * np.array(segment_percent_span)

fig, ax = plt.subplots()

for i_segs in range(n_segments):
    
    theta = thetas[i_segs]
    
    # Circular coordinates
    dy = radius * np.cos(theta)
    dz = radius * np.sin(theta)
    
    ax.scatter(Yc + dy, Zc + dz)
    
ax.set_aspect('equal')
plt.show()

#%%

panel_coords = np.array([
    [10.5918, 3.4218,  1.2358],
    [10.5918, 3.4218,  0.7912],
    [10.5918, 3.6284,  1.6294],
    [10.5918, 3.6284,  0.3975],
    [10.5918, 3.9943,  1.8819],
    [10.5918, 3.9943,  0.145],
    [10.5918, 4.4356,  1.9355],
    [10.5918, 4.4356,  0.0914],
    [10.5918, 4.8513,  1.7779],
    [10.5918, 4.8513,  0.2491],
    [10.5918, 5.1461,  1.4451],
    [10.5918, 5.1461,  0.5818],
    [10.5918, 5.2525,  1.0135],
])


panel_coords_ref = np.array([
     [0.00,  0.0,    1.0],
     [0.00,  0.5,    0.866],
     [0.00,  0.866,  0.5],
     [0.00,  1.0,    0.0],
     [0.00,  0.866, -0.5],
     [0.00,  0.5,   -0.866],
     [0.00,  0.0,   -1.0],
     [0.00, -0.5,   -0.866],
     [0.00, -0.866, -0.5],
     [0.00, -1.0,    0.0],
     [0.00, -0.866,  0.5],
     [0.00, -0.5,    0.866],
     [0.00,  0.0,    1.0],
])

fig, ax = plt.subplots()

ax.scatter(panel_coords[:, 1], panel_coords[:, 2])
ax.scatter(panel_coords_ref[:, 1], panel_coords_ref[:, 2])

ax.set_aspect('equal')
plt.show()




