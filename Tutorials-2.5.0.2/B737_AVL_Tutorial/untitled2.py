import numpy as np
from scipy.interpolate import CubicHermiteSpline
import matplotlib.pyplot as plt

# Coordinates
x = np.array([0, 4, 25, 35])
y = np.array([1, 2, 2, 4])

# Set derivatives: zero at middle points for horizontal segment
derivatives = np.zeros(4)
derivatives[1] = 0.0  # horizontal at point 2
derivatives[2] = 0.0  # horizontal at point 3

# For endpoints, we can use natural spline conditions (zero second derivative)
# or choose values that avoid inflection points

# Option 1: Use natural spline approach - estimate from neighboring segments
# A simple heuristic: use a fraction of the average slope
derivatives[0] = 0.5 * (y[1] - y[0]) / (x[1] - x[0])  # gentler slope
derivatives[3] = 0.5 * (y[3] - y[2]) / (x[3] - x[2])  # gentler slope

spline = CubicHermiteSpline(x, y, derivatives)

# Check for inflection points by examining second derivative
x_fine = np.linspace(0, 35, 1000)
y_fine = spline(x_fine)
y_second_deriv = spline.derivative(2)(x_fine)

# Find sign changes in second derivative (inflection points)
sign_changes = np.where(np.diff(np.sign(y_second_deriv)))[0]

print(f"Number of inflection points: {len(sign_changes)}")
print(f"Derivative at x=0: {derivatives[0]:.6f}")
print(f"Derivative at x=4: {derivatives[1]:.6f}")
print(f"Derivative at x=25: {derivatives[2]:.6f}")
print(f"Derivative at x=35: {derivatives[3]:.6f}")

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

ax1.plot(x_fine, y_fine, 'b-', label='Spline', linewidth=2)
ax1.plot(x, y, 'ro', markersize=8, label='Data points')
ax1.grid(True, alpha=0.3)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title('Cubic Hermite Spline with Horizontal Segment')
ax1.legend()

ax2.plot(x_fine, y_second_deriv, 'g-', linewidth=2)
ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax2.grid(True, alpha=0.3)
ax2.set_xlabel('x')
ax2.set_ylabel("y''")
ax2.set_title('Second Derivative (inflection points where this crosses zero)')

plt.tight_layout()
plt.show()