import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere

# --- user-provided constants & functions assumed present:
# cbar = 4.2350
# V = 0.45 * Atmosphere(11e3).speed_of_sound[0]
# zeta(lamda), T0p5(lamda, V, cbar), T2(...) (we'll call T2 = -T0p5(...)),
# omega0(lamda, V, cbar), CAP(lamda, V, cbar, n, alpha), tau(lamda, V, cbar)

cbar = 4.2350  # from "C:\Users\nmb48\avl_files\body_axis_derivatives_case_01_01.txt"
V = 0.45 * Atmosphere(11e3).speed_of_sound[0]

def zeta(lamda):
    return lamda.real / np.sqrt(lamda.real**2 + lamda.imag**2)  # page 129 in AE3202 Flight Dynamics Lecture Notes

def T0p5(lamda, V, cbar):
    return np.log(0.5) / lamda.real * cbar / V  # pages 125 and 127 in AE3202 Flight Dynamics Lecture Notes

def T2(lamda):
    return -T0p5(lamda)  # page 127 in AE3202 Flight Dynamics Lecture Notes

def omega0(lamda, V, cbar):
    return np.sqrt(lamda.real**2 + lamda.imag**2) * V / cbar  # page 129 in AE3202 Flight Dynamics Lecture Notes

def CAP(lamda, V, cbar, n, alpha):  # (3.111) in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris
    return omega0(lamda, V, cbar)**2 / (n / alpha)

def tau(lamda, V, cbar):
    return -np.log(0.5) * T0p5(lamda, V, cbar)
#
# If T2 was defined without V,cbar in your snippet, we will compute it here as -T0p5(...).

def compute_lambdas_from_zeta_omega(zeta0, omega0, V, cbar):
    if not (-1.0 <= zeta0 <= 1.0):
        raise ValueError("zeta must be in [-1, 1].")
    if omega0 <= 0:
        raise ValueError("omega0 must be > 0.")

    r0 = omega0 * cbar / V
    x = zeta0 * r0
    y_mag = r0 * np.sqrt(max(0.0, 1.0 - zeta0**2))

    lam_plus  = x + 1j * y_mag
    lam_minus = x - 1j * y_mag
    return lam_plus, lam_minus, r0

# Example inputs you can change:
zeta0  = 0.6     # target zeta
_omega0 = 2.5     # target omega0 (rad/s or whatever units consistent with V,cbar)
# (V and cbar come from your workspace)

# Compute lambdas
lam_p, lam_m, r0 = compute_lambdas_from_zeta_omega(zeta0, _omega0, V, cbar)

# Compute derived quantities using your functions
# T0p5 accepts (lamda, V, cbar)
T0p5_p = T0p5(lam_p, V, cbar)
T0p5_m = T0p5(lam_m, V, cbar)

T2_p = -T0p5_p   # use the intended definition T2 = -T0p5(...)
T2_m = -T0p5_m

omega0_p = omega0(lam_p, V, cbar)   # should equal omega0 (up to numerical roundoff)
omega0_m = omega0(lam_m, V, cbar)

# example CAP parameters (set to whatever n and alpha you need)
n = 1.0
alpha = 1.0
CAP_p = CAP(lam_p, V, cbar, n, alpha)
CAP_m = CAP(lam_m, V, cbar, n, alpha)

tau_p = tau(lam_p, V, cbar)
tau_m = tau(lam_m, V, cbar)

# print numbers
print("r0 (radius) =", r0)
print("lambda + =", lam_p)
print("lambda - =", lam_m)
print("")
print("T0.5 + =", T0p5_p, "  T2 + =", T2_p)
print("T0.5 - =", T0p5_m, "  T2 - =", T2_m)
print("")
print("omega0(lam +) =", omega0_p, "  omega0(lam -) =", omega0_m)
print("CAP + =", CAP_p, "  CAP - =", CAP_m)
print("tau + =", tau_p, "  tau - =", tau_m)

# --- optional plot of the locus (circle for omega0 and rays for zeta)
theta = np.linspace(0, 2*np.pi, 400)
circle_x = r0 * np.cos(theta)
circle_y = r0 * np.sin(theta)

# ray line (plot as a line through origin with angle arccos(zeta0))
tmax = 1.5 * r0
t = np.linspace(-tmax, tmax, 400)
ray_x = t * zeta0
ray_y = t * np.sqrt(max(0.0, 1 - zeta0**2))

fig, ax = plt.subplots()
ax.plot(circle_x, circle_y, label=f"omega0 = {omega0} -> r={r0:.3g}")
ax.plot(ray_x, ray_y, '--', label=f"zeta = {zeta0} (ray)")
ax.scatter([lam_p.real, lam_m.real], [lam_p.imag, lam_m.imag],
           color='C2', zorder=5, label='intersection')
ax.set_aspect('equal', 'box')
ax.set_xlabel('Re')
ax.set_ylabel('Im')
ax.legend()
ax.grid(True)
plt.show()
