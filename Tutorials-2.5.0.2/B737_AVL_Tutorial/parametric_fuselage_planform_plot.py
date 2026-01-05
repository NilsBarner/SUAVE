import sys
import copy
import numpy as np
from scipy.interpolate import Akima1DInterpolator
import matplotlib.pyplot as plt

#%% Vertical fuselage cross-section

N_points_long = 500
N_points_vert = 50
N_points_horz = 50
r_fuse = 2
x = np.array([0, 5, 15, 25, 35])
z = np.array([-1, 0, 0, 0, 1])
x_nose_start = x[0]
x_nose_end = x[1]
x_tail_start = x[-2]
x_tail_end = x[-1]
z_centreline = z[2]
xs = np.linspace(min(x), max(x), num=N_points_long)
xs_nose = copy.deepcopy(xs)[xs < x_nose_end]
xs_tail = copy.deepcopy(xs)[xs > x_tail_start]
nnose = len(xs_nose)
ntail = len(xs_tail)

# z_akima = Akima1DInterpolator(x, z, method="akima")(xs)
z_makima = Akima1DInterpolator(x, z, method="makima")(xs)
z_offset_nose = z_makima[0]
z_offset_tail = z_makima[-1]

# Nose top

a_nose_top = 1.5
r_nose_top = r_fuse + z[1] - z[0]

x_top_nose = []
z_top_nose = []
for i, _x in enumerate(xs_nose, start=1):
    fraci = (i - 1) / (nnose - 1)
    fracx = np.cos(0.5 * np.pi * fraci)
    x_top_nose.append(x_nose_end + (x_nose_start - x_nose_end) * fracx)
    z_top_nose.append(z_offset_nose + (r_nose_top * (1.0 - fracx**a_nose_top)**(1.0 / a_nose_top)))

# Nose bottom

a_nose_bottom = 2
r_nose_bottom = r_fuse - (z[1] - z[0])

x_bottom_nose = []
z_bottom_nose = []
for i, _x in enumerate(xs_nose, start=1):
    fraci = (i - 1) / (nnose - 1)
    fracx = np.cos(0.5 * np.pi * fraci)
    x_bottom_nose.append(x_nose_end + (x_nose_start - x_nose_end) * fracx)
    z_bottom_nose.append(z_offset_nose - (r_nose_bottom * (1.0 - fracx**a_nose_bottom)**(1.0 / a_nose_bottom)))

z_makima_nose = np.linspace(-r_fuse, r_fuse, num=N_points_vert,endpoint=True)
x_makima_vert_nose = Akima1DInterpolator(
    np.concatenate((z_bottom_nose[::-1][:-1], z_top_nose)),
    np.concatenate((x_bottom_nose[::-1][:-1], x_top_nose)),
    method="makima",
)(z_makima_nose)
# sys.exit()

# Tail top

b_tail_top = 3
r_tail_top = r_fuse - (z[-1] - z[-2])

x_top_tail = []
z_top_tail = []
for i, _x in enumerate(xs_tail, start=1):
    fraci = (i - 1) / (ntail - 1)
    fracx = np.cos(0.5 * np.pi * fraci)
    x_top_tail.append(x_tail_start + (x_tail_end - x_tail_start) * fracx)
    z_top_tail.append(z_offset_tail + (r_tail_top + (-0.4*r_tail_top) * fracx**b_tail_top))

# Tail bottom

b_tail_bottom = 2
r_tail_bottom = r_fuse + (z[-1] - z[-2])

x_bottom_tail = []
z_bottom_tail = []
for i, _x in enumerate(xs_tail, start=1):
    fraci = (i - 1) / (ntail - 1)
    fracx = np.cos(0.5 * np.pi * fraci)
    x_bottom_tail.append(x_tail_start + (x_tail_end - x_tail_start) * fracx)
    z_bottom_tail.append(z_offset_tail - (r_tail_bottom + (-0.9 * r_tail_bottom) * fracx**b_tail_bottom))
    
z_makima_tail = np.linspace(-r_fuse, r_fuse, num=N_points_horz,endpoint=True)
x_makima_vert_tail = Akima1DInterpolator(
    np.concatenate((
        z_bottom_tail[::-1][:-1],
        [z_bottom_tail[0], (z_bottom_tail[0] + z_top_tail[0]) / 2, z_top_tail[0]],
        z_top_tail[1:],
    )),
    np.concatenate((
        x_bottom_tail[::-1][:-1],
        [x_bottom_tail[0], (x_bottom_tail[0] + x_top_tail[0]) / 2, x_top_tail[0]],
        x_top_tail[1:],
    )),
    method="makima",
)(z_makima_tail)
# sys.exit()

# Plot curves

fig, ax = plt.subplots()
ax.plot(x, z, "o", label="data")
# ax.plot(xs, z_akima, label="akima")
ax.plot(xs, z_makima, label="makima")

ax.plot(x_top_nose, z_top_nose)
ax.plot(x_bottom_nose, z_bottom_nose)

ax.plot(x_top_tail, z_top_tail)
ax.plot(x_bottom_tail, z_bottom_tail)

ax.plot(
    [x_top_nose[-1], x_top_tail[-1]],
    [z_top_nose[-1], z_top_tail[-1]]
)
ax.plot(
    [x_bottom_nose[-1], x_bottom_tail[-1]],
    [z_bottom_nose[-1], z_bottom_tail[-1]]
)

ax.plot(
    [x_top_tail[0], x_bottom_tail[0]],
    [z_top_tail[0], z_bottom_tail[0]]
)

ax.scatter(x_makima_vert_nose, z_makima_nose, marker='.', color='red', zorder=100)
ax.scatter(x_makima_vert_tail, z_makima_tail, marker='.', color='red', zorder=100)

ax.set_aspect('equal')
ax.legend()
plt.show()

# import sys
# sys.exit()

#%% Horizontal fuselage cross-section

# Nose left

a_nose_left = 1.5

x_left_nose = []
y_left_nose = []
for i, _x in enumerate(xs_nose, start=1):
    fraci = (i - 1) / (nnose - 1)
    fracx = np.cos(0.5 * np.pi * fraci)
    x_left_nose.append(x_nose_end + (x_nose_start - x_nose_end) * fracx)
    y_left_nose.append(r_fuse * (1.0 - fracx**a_nose_left)**(1.0 / a_nose_left))

y_makima_nose = np.linspace(-r_fuse, r_fuse, num=N_points_vert,endpoint=True)
x_makima_horz_nose = Akima1DInterpolator(
    np.concatenate((-np.array(y_left_nose[::-1][:-1]), y_left_nose)),
    np.concatenate((x_left_nose[::-1][:-1], x_left_nose)),
    method="makima",
)(y_makima_nose)

# Tail left

btail = 3

x_left_tail = []
y_left_tail = []
for i, _x in enumerate(xs_tail, start=1):
    fraci = (i - 1) / (ntail - 1)
    fracx = np.cos(0.5 * np.pi * fraci)
    x_left_tail.append(x_tail_start + (x_tail_end - x_tail_start) * fracx)
    y_left_tail.append(r_fuse + (-0.7 * r_fuse) * fracx**btail)
    
y_makima_tail = np.linspace(-r_fuse, r_fuse, num=N_points_vert,endpoint=True)
x_makima_horz_tail = Akima1DInterpolator(
    np.concatenate((
        -np.array(y_left_tail[::-1][:-1]),
        [-y_left_tail[0], 0, y_left_tail[0]],
        y_left_tail[1:],
    )),
    np.concatenate((
        x_left_tail[::-1][:-1],
        [x_left_tail[0], x_left_tail[0], x_left_tail[0]],
        x_left_tail[1:],
    )),
    method="makima",
)(y_makima_tail)

# Plot curves

fig, ax = plt.subplots()

ax.plot(x_left_nose, np.array(y_left_nose))
ax.plot(x_left_nose, -np.array(y_left_nose))

ax.plot(x_left_tail, np.array(y_left_tail))
ax.plot(x_left_tail, -np.array(y_left_tail))

ax.plot(
    [x_left_nose[-1], x_left_tail[-1]],
    [y_left_nose[-1], y_left_tail[-1]]
)
ax.plot(
    [x_left_nose[-1], x_left_tail[-1]],
    [-y_left_nose[-1], -y_left_tail[-1]]
)

ax.plot(
    [x_left_tail[0], x_left_tail[0]],
    [y_left_tail[0], -y_left_tail[0]]
)

ax.scatter(x_makima_horz_nose, y_makima_nose, marker='.', color='red', zorder=100)
ax.scatter(x_makima_horz_tail, y_makima_tail, marker='.', color='red', zorder=100)

ax.set_aspect('equal')
ax.legend()
plt.show()

