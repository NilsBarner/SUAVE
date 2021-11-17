## @ingroup Methods-Aerodynamics-Common-Fidelity_Zero-Lift
# compute_HFW_inflow_velocities.py
#
# Created:  Sep 2021, R. Erhard
# Modified:    

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
from SUAVE.Core import Data
from SUAVE.Methods.Aerodynamics.Common.Fidelity_Zero.Lift.generate_propeller_wake_distribution import generate_propeller_wake_distribution
from SUAVE.Methods.Aerodynamics.Common.Fidelity_Zero.Lift.compute_wake_induced_velocity import compute_wake_induced_velocity

# package imports
import numpy as np
from scipy.interpolate import interp1d

def compute_HFW_inflow_velocities( prop ):
    """
    Assumptions:
        None

    Source:
        N/A
    Inputs:
        prop - rotor instance
    Outputs:
        Va   - axial velocity array of shape (ctrl_pts, Nr, Na)        [m/s]
        Vt   - tangential velocity array of shape (ctrl_pts, Nr, Na)   [m/s]
    """     
    # use results from prior bemt iteration
    prop_outputs  = prop.outputs
    cpts          = len(prop_outputs.velocity)
    Na            = prop.number_azimuthal_stations
    Nr            = len(prop.chord_distribution)
    r             = prop.radius_distribution


    # check if VD for rotor system has already been computed
    if prop.system_vortex_distribution:# is not None and prop.self.wake_settings.converge_HFW:
        # use pre-computed system vortex distribution
        VD = prop.system_vortex_distribution
        WD = VD.Wake_collapsed
    else:
        print("Vortex distribution not yet computed. Generating distribution for single rotor...")
    
        props=Data()
        props.propeller = prop
        
        # generate wake distribution using initial circulation from BEMT
        VD = Data()
        WD, _, _  = generate_propeller_wake_distribution(props,cpts,VD )   

    # set shape of velocitie arrays
    Va = np.zeros((cpts,Nr,Na))
    Vt = np.zeros((cpts,Nr,Na))
    
    
    prop.wake_skew_angle = WD.wake_skew_angle
    
    # ----------------------------------------------------------------
    # Compute the wake-induced velocities at propeller blade
    # ----------------------------------------------------------------
    # set the evaluation points in the vortex distribution: (Na, nblades, Nr)
    r = prop.radius_distribution 
    Yb   = prop.Wake_VD.Yblades_cp[:,0,:] 
    Zb   = prop.Wake_VD.Zblades_cp[:,0,:] 
    Xb   = prop.Wake_VD.Xblades_cp[:,0,:] 
    

    VD.YC = (Yb[:,1:] + Yb[:,:-1])/2
    VD.ZC = (Zb[:,1:] + Zb[:,:-1])/2
    VD.XC = (Xb[:,1:] + Xb[:,:-1])/2
     
    
    VD.n_cp = np.size(VD.YC)

    # Compute induced velocities at blade from the system of prescribed wakes
    V_ind   = compute_wake_induced_velocity(WD, VD, cpts)
    u_2d       = V_ind[0,:,:,0]   # velocity in vehicle x-frame
    v_2d       = V_ind[0,:,:,1]   # velocity in vehicle y-frame
    w_2d       = V_ind[0,:,:,2]   # velocity in vehicle z-frame
    
    # interpolate to get values at rotor radial stations
    r_midpts = (r[1:] + r[:-1])/2
    for i in range(Na):
        blade_angle  = (i*(2*np.pi/(Na))) * prop.rotation
        
        u = u_2d[i,:]
        v = v_2d[i,:]
        w = w_2d[i,:]
        
        u_r = interp1d(r_midpts, u, fill_value="extrapolate")
        v_r = interp1d(r_midpts, v, fill_value="extrapolate")
        w_r = interp1d(r_midpts, w, fill_value="extrapolate")
        
        up = u_r(r)
        vp = v_r(r)
        wp = w_r(r)       

        # Update velocities at the disc
        Va[:,:,i]  = -up
        Vt[:,:,i]  = (-vp*np.cos(blade_angle) + wp*np.sin(blade_angle)) 


    prop.system_vortex_distribution = VD

    return Va, Vt