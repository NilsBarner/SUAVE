import numpy as np
from scipy.interpolate import Akima1DInterpolator, CubicHermiteSpline
import matplotlib.pyplot as plt

# Coordinates
x = np.array([0, 5, 15, 25, 35])
y = np.array([1.25, 2, 2, 2, 3.25]) + 3
xs = np.linspace(min(x), max(x), num=100)
r_fuse = 2
r_nose = 0.25
y_centreline = y[2]

# Centerline
y_makima = Akima1DInterpolator(x, y, method="makima")(xs)

# Get key points from the actual data
y_tip = y[0]  # First point: y = 4
y_tail = y[-1]  # Last point: y = 6.5
x_nose_junction = x[1]  # x = 5
x_tail_start = x[-2]  # x = 35

# Upper and lower horizontal line y-coordinates
y_upper = y_centreline + r_fuse  # y = 7
y_lower = y_centreline - r_fuse  # y = 3

# Tail extension
x_tail_end = x[-1]  # Extend 5 units beyond

# NOSE SECTION: Monotonic splines with bounded derivatives
# For monotonicity without overshoot: |d0| ≤ 3*|Δy|/Δx and |d1| ≤ 3*|Δy|/Δx
# But we also want vertical at tip, so we use the maximum allowable derivative

dx_nose = x_nose_junction - 0
dy_upper_nose = y_upper - y_tip  # 3
dy_lower_nose = y_lower - y_tip  # -1

# Upper nose: curve upward from tip to upper line
# Max derivative for monotonic Hermite: 3*Δy/Δx
max_deriv_upper_nose = 3 * abs(dy_upper_nose) / dx_nose  # 3*3/5 = 1.8
x_nose = np.array([0, x_nose_junction])
y_upper_nose_vals = np.array([y_tip, y_upper])
dydx_upper_nose = np.array([max_deriv_upper_nose, 0])  # Steep at tip, horizontal at junction

upper_nose_spline = CubicHermiteSpline(x_nose, y_upper_nose_vals, dydx_upper_nose)

# Lower nose: curve downward from tip to lower line
max_deriv_lower_nose = 3 * abs(dy_lower_nose) / dx_nose  # 3*1/5 = 0.6
y_lower_nose_vals = np.array([y_tip, y_lower])
dydx_lower_nose = np.array([-max_deriv_lower_nose, 0])  # Steep downward at tip, horizontal at junction

lower_nose_spline = CubicHermiteSpline(x_nose, y_lower_nose_vals, dydx_lower_nose)

# TAIL SECTION: Monotonic splines with bounded derivatives
dx_tail = x_tail_end - x_tail_start
dy_upper_tail = y_tail - y_upper  # -0.5
dy_lower_tail = y_tail - y_lower  # 3.5

# Upper tail: curve downward from upper line to tail
max_deriv_upper_tail = 3 * abs(dy_upper_tail) / dx_tail  # 3*0.5/5 = 0.3
x_tail = np.array([x_tail_start, x_tail_end])
y_upper_tail_vals = np.array([y_upper, y_tail])
dydx_upper_tail = np.array([0, -max_deriv_upper_tail])  # Horizontal at start, steep downward at end

upper_tail_spline = CubicHermiteSpline(x_tail, y_upper_tail_vals, dydx_upper_tail)

# Lower tail: curve upward from lower line to tail
max_deriv_lower_tail = 3 * abs(dy_lower_tail) / dx_tail  # 3*3.5/5 = 2.1
y_lower_tail_vals = np.array([y_lower, y_tail])
dydx_lower_tail = np.array([0, max_deriv_lower_tail])  # Horizontal at start, steep upward at end

lower_tail_spline = CubicHermiteSpline(x_tail, y_lower_tail_vals, dydx_lower_tail)

# Generate fine points for plotting
x_nose_fine = np.linspace(0, x_nose_junction, 100)
x_straight_fine = np.linspace(x_nose_junction, x_tail_start, 100)
x_tail_fine = np.linspace(x_tail_start, x_tail_end, 100)

y_upper_nose_fine = upper_nose_spline(x_nose_fine)
y_lower_nose_fine = lower_nose_spline(x_nose_fine)
y_upper_tail_fine = upper_tail_spline(x_tail_fine)
y_lower_tail_fine = lower_tail_spline(x_tail_fine)

# Plot
fig, ax = plt.subplots(figsize=(14, 7))

# Centerline
ax.plot(x, y, "o", label="centerline data", markersize=8, zorder=5, color='blue')
ax.plot(xs, y_makima, label="centerline (makima)", linewidth=2, zorder=4, color='blue')

# Upper surface
ax.plot(x_nose_fine, y_upper_nose_fine, 'r-', linewidth=2.5, label="upper surface")
ax.plot(x_straight_fine, np.ones_like(x_straight_fine) * y_upper, 'r-', linewidth=2.5)
ax.plot(x_tail_fine, y_upper_tail_fine, 'r-', linewidth=2.5)

# Lower surface
ax.plot(x_nose_fine, y_lower_nose_fine, 'g-', linewidth=2.5, label="lower surface")
ax.plot(x_straight_fine, np.ones_like(x_straight_fine) * y_lower, 'g-', linewidth=2.5)
ax.plot(x_tail_fine, y_lower_tail_fine, 'g-', linewidth=2.5)

# Reference horizontal lines (dashed)
ax.axhline(y_centreline - r_fuse, linestyle='--', color='gray', alpha=0.5, linewidth=1)
ax.axhline(y_centreline + r_fuse, linestyle='--', color='gray', alpha=0.5, linewidth=1)

# Mark junction points
ax.plot([x_nose_junction, x_nose_junction], [y_lower, y_upper], 'ko', markersize=5, zorder=6)
ax.plot([x_tail_start, x_tail_start], [y_lower, y_upper], 'ko', markersize=5, zorder=6)

# Mark end points
ax.plot(0, y_tip, 'ko', markersize=8, zorder=6)
ax.plot(x_tail_end, y_tail, 'ko', markersize=8, zorder=6)

ax.set_aspect('equal')
ax.legend(loc='upper left')
ax.grid(True, alpha=0.3)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Fuselage Outline with Monotonic Hermite Splines')
plt.tight_layout()
plt.show()

# Print diagnostics
print(f"Nose section: x = {x_nose[0]:.1f} to {x_nose[1]:.1f}")
print(f"  Upper nose: Δy = {dy_upper_nose:.1f}, max derivative = {max_deriv_upper_nose:.2f}")
print(f"  Lower nose: Δy = {dy_lower_nose:.1f}, max derivative = {max_deriv_lower_nose:.2f}")
print(f"\nStraight section: x = {x_nose_junction:.1f} to {x_tail_start:.1f}")
print(f"\nTail section: x = {x_tail_start:.1f} to {x_tail_end:.1f}")
print(f"  Upper tail: Δy = {dy_upper_tail:.1f}, max derivative = {max_deriv_upper_tail:.2f}")
print(f"  Lower tail: Δy = {dy_lower_tail:.1f}, max derivative = {max_deriv_lower_tail:.2f}")