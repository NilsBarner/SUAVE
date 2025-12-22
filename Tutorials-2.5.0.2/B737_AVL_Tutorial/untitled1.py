import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

# Coordinates
x = np.array([0, 4, 25, 35])
y = np.array([1, 2, 2, 4])

# To ensure horizontal segment between points 1 and 2,
# we need to set the derivatives at those points to zero
bc_type = ((1, 0.0),  # derivative at first point (will be computed naturally)
           (1, 0.0))   # derivative at last point (will be computed naturally)

# Actually, we need a different approach - specify derivatives at middle points
# Using CubicSpline with custom boundary conditions won't work directly
# Let's manually set derivatives

# Calculate derivatives: set to 0 at indices 1 and 2 for horizontal segment
derivatives = np.zeros(4)
derivatives[0] = (y[1] - y[0]) / (x[1] - x[0])  # slope at first point
derivatives[1] = 0.0  # horizontal at point 2
derivatives[2] = 0.0  # horizontal at point 3  
derivatives[3] = (y[3] - y[2]) / (x[3] - x[2])  # slope at last point

# Create cubic Hermite spline with specified derivatives
from scipy.interpolate import CubicHermiteSpline

spline = CubicHermiteSpline(x, y, derivatives)

# Generate points for plotting
x_fine = np.linspace(0, 35, 500)
y_fine = spline(x_fine)

# Plot
plt.figure(figsize=(10, 6))
plt.plot(x_fine, y_fine, 'b-', label='Spline', linewidth=2)
plt.plot(x, y, 'ro', markersize=8, label='Data points')
plt.grid(True, alpha=0.3)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Cubic Hermite Spline with Horizontal Segment')
plt.legend()
plt.show()

# Print the spline object for use
print("Spline created successfully!")
print(f"Derivative at point 2 (x=4): {spline.derivative()(4):.6f}")
print(f"Derivative at point 3 (x=25): {spline.derivative()(25):.6f}")