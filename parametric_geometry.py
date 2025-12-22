"""
This script defines airframe and engine geometrical elements parametrically.
NOTE: this is only the tip of the iceberg of what Kulfan's shape transforms
(low_paper_2008_univparamgeomreprmeth_kulfan) can do - I can also model
cross-sectional shapes, as described later in the above paper and evaluate
these using CFD.
"""

import sys
import numpy as np
from scipy.integrate import simpson
from typing import Annotated, Tuple
from scipy.optimize import minimize

from matplotlib_custom_settings import *

def generate_streamlined_body_geometry(
    R_le_over_c: Annotated[float, "[-]"],
    beta_tail: Annotated[float, "[deg]"],
    Psi_zeta_max: Annotated[float, "[-]"], # this value does not fully work as expected - to be investigated once less pressing stuff has been dealt with
    zeta_max: Annotated[float, "[-]"],
    zeta_te: Annotated[float, "[-]"],
    dimensional_known_dict: Annotated[dict, "[m] or [m^3]"] # length or volume
) -> Tuple[float, float, float, float, float, float]:
    """
    Generate nacelle geometry from airfoil-shaped body of revolution based on low_paper_2008_univparamgeomreprmeth_kulfan
    by either specifying the volume and fineness ratio or the length and fineness ratio (diameter).
    """
    
    def _define_quadratic(x0, y0, x1, y1):
        """
        Generate second-order polynomial tangent to point (x0, y0)
        and secant through point (x1, y1).
        """
        a = (y1 - y0) / (x1 - x0)**2
        b = -2 * a * x0
        c = y0 + a * x0**2
        return np.poly1d([a, b, c])
    
    # Fixed parameters
    N_1, N_2 = 0.5, 1 # for round-nose airfoil
    Psi_le = 0
    Psi_te = 1
    
    l_over_d = 1/zeta_max/2 # concerns HALF of airfoil (see definition of zeta_max on page 2 and in Fig. 1 of low_paper_2008_univparamgeomreprmeth_kulfan)
    
    Vol = None
    x_list_sorted = None
    y_list_sorted = None
    Psi_list = None
    zeta_list = None
    zeta_list_sorted = None
    
    def _calculate_volume_from_length(l):
        nonlocal Vol, y_list_sorted, x_list_sorted, Psi_list, zeta_list, zeta_list_sorted
        
        d_max = l[0]/zeta_max
        
        S_le = np.sqrt(2 * R_le_over_c) # (4)
        S_te = np.tan(beta_tail * np.pi/180) + zeta_te # (5)
        S_zeta_max = (zeta_max - Psi_zeta_max * zeta_te) / (np.sqrt(Psi_zeta_max) * (1 - Psi_zeta_max)) # (2)

        quadratic_front = _define_quadratic(Psi_zeta_max, S_zeta_max, Psi_le, S_le)
        quadratic_aft = _define_quadratic(Psi_zeta_max, S_zeta_max, Psi_te, S_te)
        
        zeta_list = []
        Psi_list = np.linspace(1e-4, 1, 100)
        
        for Psi in Psi_list:
            if Psi >= Psi_zeta_max:
                S = quadratic_aft(Psi)
            elif Psi < Psi_zeta_max:
                S = quadratic_front(Psi)
            C = Psi**N_1 * (1 - Psi)**N_2 # (6)
            zeta = C * S + Psi * zeta_te # (7)
            zeta_list.append(zeta)
        
        # Sort points in ascending x-order for integration
        sorted_indices = np.argsort(Psi_list)
        Psi_list_sorted = Psi_list[sorted_indices]
        zeta_list_sorted = np.array(zeta_list)[sorted_indices]
        
        # Go from nondimensional coordinates to dimensional lengths
        x_list_sorted = Psi_list_sorted * l[0]
        y_list_sorted = zeta_list_sorted * l[0]  # zeta_max refers to max. ndim. radius, hence x2
        
        # Calculate the volume of revolution around the x-axis using the disk method: V = pi * integral( y^2 dx )
        Vol = np.pi * simpson(y_list_sorted**2, x=x_list_sorted)
        
        f = abs(Vol - _Vol)
        
        return f
    
    if list(dimensional_known_dict.keys())[0] == 'length':
        _Vol = np.inf # dummy value
        l = list(dimensional_known_dict.values())[0]
        _calculate_volume_from_length([l]) # input must be of type list
    elif list(dimensional_known_dict.keys())[0] == 'volume':
        _Vol = list(dimensional_known_dict.values())[0]
        bnds = [(0, np.inf)]
        x0 = (_Vol * 4 * l_over_d**2 / np.pi)**(1/3) # use cylinder as initial guess
        outputs = minimize(
            _calculate_volume_from_length, x0=x0, method='SLSQP', bounds=bnds
        )
        l = outputs.x[0]
    
    # Calculate the surface area of revolution around the x-axis: S = 2 * pi * integral(f * sqrt(1 + fprime**2))
    dy_dx_list = np.gradient(y_list_sorted, x_list_sorted)
    integrand = 2 * np.pi * y_list_sorted * np.sqrt(1 + dy_dx_list**2)
    surface_area = simpson(integrand, x=x_list_sorted)
    
    # Calculate the surface area up to the point of maximum thickness
    idx_zeta_max = np.argmin(abs(zeta_list_sorted - zeta_max))
    y_list_masked = y_list_sorted[:idx_zeta_max + 1]
    x_list_masked = x_list_sorted[:idx_zeta_max + 1]
    dy_dx_list_masked = np.gradient(y_list_masked, x_list_masked)
    integrand_masked = 2 * np.pi * y_list_masked * np.sqrt(1 + dy_dx_list_masked**2)
    surface_area_le_to_zeta_max = simpson(integrand_masked, x=x_list_masked)
    
    # Calculate the frontal area
    frontal_area = np.pi * (zeta_max * l)**2
    
    # Calculate the circumference (Euclidean distance between consecutive points) in 2D
    circumference = np.sum(np.sqrt(np.diff(x_list_sorted)**2 + np.diff(y_list_sorted)**2)) * 2 # x2 since there is an upper and a lower half of the airfoil
    
    # Create closed area from half-airfoil contour
    Psi_list_closed = np.append(Psi_list, np.flip(Psi_list))
    zeta_list_closed = np.append(np.array(zeta_list), -1 * np.flip(np.array(zeta_list)))
    
    # Check: plot side-view of nacelle contour
    # fig, ax = plt.subplots()
    # ax.plot(Psi_list_closed, zeta_list_closed, color = 'black')
    # ax.set_aspect('equal')
    # ax.set_xlabel(r'$\Psi=x/c$')
    # ax.set_ylabel(r'$\zeta=z/c$')
    # ax.set_xlim(0, 1)
    # plt.show()
    
    return Vol, l, circumference, frontal_area, surface_area, Psi_list_closed, zeta_list_closed, surface_area_le_to_zeta_max
    # return Vol, l, circumference, frontal_area, surface_area, Psi_list, zeta_list, surface_area_le_to_zeta_max
    
#%%

if __name__ == '__main__':
    
    Vol, l, _, _, surface_area, Psi_list_closed, zeta_list_closed, _ = \
        generate_streamlined_body_geometry(
            R_le_over_c=0.06, beta_tail=30, Psi_zeta_max=0.35, zeta_max=1/5/2, zeta_te=0,
            dimensional_known_dict={'volume': 6}
        )
    
    # Plot airfoils
    fig, ax = plt.subplots()
    ax.plot(Psi_list_closed * l, zeta_list_closed * l, color='blue', label='Manual result', linestyle='dashed')
    ax.set_xlabel('Chord (m)')
    ax.set_ylabel('Thickness (m)')
    ax.set_aspect('equal')
    ax.spines[['right', 'top']].set_visible(False)
    ax.tick_params(axis='x', which='both', top=False)
    ax.tick_params(axis='y', which='both', right=False)
    ax.tick_params(bottom=False)
    # plt.savefig(f'{CONFIG.extHEX_type}.png')
    ax.set_aspect('equal')
    plt.show()
    
    # Test whether length- and volume-input yield the same result
    Vol, l, _, _, surface_area, Psi_list_closed, zeta_list_closed, _ = \
        generate_streamlined_body_geometry(
            R_le_over_c=0.06, beta_tail=30, Psi_zeta_max=0.35, zeta_max=1/5/2, zeta_te=0,
            dimensional_known_dict={'length': 3}
        )
    print('surface_area =', surface_area)
    Vol, l, _, _, surface_area, Psi_list_closed, zeta_list_closed, _ = \
        generate_streamlined_body_geometry(
            R_le_over_c=0.06, beta_tail=30, Psi_zeta_max=0.35, zeta_max=1/5/2, zeta_te=0,
            dimensional_known_dict={'volume': Vol})
    print('surface_area =', surface_area)
    
    sys.exit()
    
    #%% Reproduce Fig. 2
    
    def define_quadratic(x0, y0, x1, y1):
        """
        Generate second-order polynomial tangent to point (x0, y0)
        and secant through point (x1, y1).
        """
        a = (y1 - y0) / (x1 - x0)**2
        b = -2 * a * x0
        c = y0 + a * x0**2
        return np.poly1d([a, b, c])
    
    N_1, N_2 = 0.5, 1 # for round-nose airfoil
    zeta_T = 0
    Psi_Z_max = 0.35
    S_Z_max = 0.15 # not provided but matches Fig. 2 reasonably well
    Psi_le = 0
    Psi_te = 1
    
    S_le_range = [0.11, 0.08, 0.05]
    S_te_range = [0.14, 0.08, 0.04]
    Psi_range = np.linspace(1e-4, 1, 1000)
    
    # Vary LE shape function
    
    S_list = []
    zeta_list = []
    Psi_list = []
    
    for S_le in S_le_range:
        S_sublist = []
        zeta_sublist = []
        quadratic_front = define_quadratic(Psi_Z_max, S_Z_max, Psi_le, S_le)
        for Psi in Psi_range:       
            if Psi >= Psi_Z_max:
                S = S_Z_max
            elif Psi < Psi_Z_max:
                S = quadratic_front(Psi)
            
            C = Psi**N_1 * (1 - Psi)**N_2
            zeta = C * S + Psi * zeta_T
            
            S_sublist.append(S)
            zeta_sublist.append(zeta)
            
        S_list.append(S_sublist)
        zeta_subarray = np.append(np.array(zeta_sublist), -1 * np.flip(np.array(zeta_sublist)))
        Psi_subarray = np.append(Psi_range, np.flip(Psi_range))
        zeta_list.append(zeta_subarray)
        Psi_list.append(Psi_subarray)
        
    # Plot LHS of Fig. 2
        
    fig, ax = plt.subplots()
    for i, S_sublist in enumerate(S_list):
        ax.plot(Psi_range, S_sublist)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.2)
    ax.set_xlabel(r'$\Psi=x/c$')
    ax.set_ylabel('$S$')
    plt.show()
    
    fig, ax = plt.subplots()
    for i, zeta_sublist in enumerate(zeta_list):
        ax.plot(Psi_list[i], zeta_sublist - 0.2 * i, color = 'black')
    ax.set_aspect('equal')
    ax.set_xlim(0, 1)
    ax.set_xlabel(r'$\Psi=x/c$')
    ax.set_ylabel(r'$\zeta=z/c$')
    plt.show()
    
    # Vary TE shape function
    
    S_list = []
    zeta_list = []
    Psi_list = []
    
    for S_te in S_te_range:
        S_sublist = []
        zeta_sublist = []
        quadratic_aft = define_quadratic(Psi_Z_max, S_Z_max, Psi_te, S_te)
        for Psi in Psi_range:       
            if Psi <= Psi_Z_max:
                S = S_Z_max
            elif Psi > Psi_Z_max:
                S = quadratic_aft(Psi)
                
            C = Psi**N_1 * (1 - Psi)**N_2
            zeta = C * S + Psi * zeta_T
                
            S_sublist.append(S)
            zeta_sublist.append(zeta)
            
        S_list.append(S_sublist)
        zeta_subarray = np.append(np.array(zeta_sublist), -1 * np.flip(np.array(zeta_sublist)))
        Psi_subarray = np.append(Psi_range, np.flip(Psi_range))
        zeta_list.append(zeta_subarray)
        Psi_list.append(Psi_subarray)
        
    # Plot RHS of Fig. 2
        
    fig, ax = plt.subplots()
    for i, S_sublist in enumerate(S_list):
        ax.plot(Psi_range, S_sublist)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.2)
    ax.set_xlabel(r'$\Psi=x/c$')
    ax.set_ylabel('$S$')
    plt.show()
    
    fig, ax = plt.subplots()
    for i, zeta_sublist in enumerate(zeta_list):
        ax.plot(Psi_list[i], zeta_sublist - 0.2 * i, color = 'black')
    ax.set_aspect('equal')
    ax.set_xlim(0, 1)
    ax.set_xlabel(r'$\Psi=x/c$')
    ax.set_ylabel(r'$\zeta=z/c$')
    plt.show()

    #%% Reproduce S_Z_max in Fig. 1
    
    Psi_range = np.linspace(1e-4, 1 - 1e-4, 100)

    Z_max = 0.02
    zeta_T = 0

    Z_max_transf_list = []

    for Psi in Psi_range:
        zeta = Z_max
        Z_max_transf = (zeta - Psi * zeta_T) / (np.sqrt(Psi) * (1 - Psi))
        Z_max_transf_list.append(Z_max_transf)
        
    fig, ax = plt.subplots()
    ax.plot(Psi_range, Z_max_transf_list)
    ax.set_ylim(0, 0.2)
    ax.set_xlabel(r'$\Psi=x/c$')
    ax.set_ylabel('$S$')
    plt.show()
    


