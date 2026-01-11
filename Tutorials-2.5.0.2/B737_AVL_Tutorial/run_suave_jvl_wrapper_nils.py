"""
This file is adapted from tut_mission_B737_AVL.py
to contain the minimum set of SUAVE methods needed to
perform a static and dynamic stability analysis on
aircraft geometry data imported as a .csv from TASOPT.jl.

NOTE: the definition of a wing segment always requires at
least two instances of SUAVE.Components.Wings.Segment() -
one at segment.percent_span_location = 0.0 and
one at segment.percent_span_location = 1.0, so that the
wing itself can be interpolated in between.
"""

__all__ = []

import sys
import ast
import numpy as np
import pandas as pd
import pylab as plt
from copy import deepcopy
from datetime import date

import SUAVE
assert SUAVE.__version__=='2.5.0', 'These tutorials only work with the SUAVE 2.5.0 release'
from SUAVE.Core import Units
from SUAVE.Methods.Geometry.Two_Dimensional.Planform import segment_properties

# ----------------------------------------------------------------------
#   Analysis Setup
# ----------------------------------------------------------------------

def full_setup(vehicle=None):

    # vehicle data
    if vehicle == None:
        vehicle  = vehicle_setup()
        # print(id(vehicle))
    configs  = configs_setup(vehicle)

    # vehicle analyses
    configs_analyses = analyses_setup(configs)

    analyses = SUAVE.Analyses.Analysis.Container()
    analyses.configs  = configs_analyses
    
    return configs, analyses

# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ----------------------------------------------------------------------

def analyses_setup(configs):

    analyses = SUAVE.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

def base_analysis(vehicle):

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = SUAVE.Analyses.Vehicle()

    # ------------------------------------------------------------------
    #  Stability Analysis
    stability = SUAVE.Analyses.Stability.AVL()
    # stability.settings.filenames.avl_bin_name = r"C:\Users\nmb48\Documents\GitHub\SUAVE\avl3.52\avl.exe"  # NILS: set AVL executable name
    stability.settings.filenames.avl_bin_name = r"C:\Users\nmb48\Documents\GitHub\SUAVE\jvl2.16\jvl.exe"
    #stability.settings.spanwise_vortex_density                  = 3
    stability.geometry = vehicle
    analyses.append(stability)
    
    return analyses

def vehicle_setup():
    
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------    
    
    vehicle = SUAVE.Vehicle()
    vehicle.tag = 'Boeing_737-800'
    
    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------    

    # mass properties
    # NILS: note here that I EITHER have to
    # set aircraft.mass_properties.mass OR
    # aircraft.mass_properties.max_takeoff
    # but NOT both, else `mass` in the xxx.mass
    # file might be incorrect!
    # NOTE: the following mass properties are
    # STRICTLY REQUIRED to run successfully!
    # vehicle.mass_properties.max_takeoff               = 79015.8 * Units.kilogram
    vehicle.mass_properties.max_takeoff               = 0.0  # NILS
    vehicle.mass_properties.mass                      = 79015.8 * Units.kilogram  # NILS
    vehicle.mass_properties.takeoff                   = 79015.8 * Units.kilogram   
    vehicle.mass_properties.max_zero_fuel             = 62732.0 * Units.kilogram 
    # NILS: copied from regression\scripts\Vehicles\Boeing_737.py - needed for DYNAMIC stability analysis
    vehicle.mass_properties.center_of_gravity         = [[ 15.30987849,   0.        ,  -0.48023939]]
    vehicle.mass_properties.moments_of_inertia.tensor = [[3173074.17, 0 , 28752.77565],[0 , 3019041.443, 0],[0, 0, 5730017.433]] # estimated, not correct
    
    # ------------------------------------------------------------------        
    #   Main Wing
    # ------------------------------------------------------------------        
    
    # NILS: SINGLE-SECTION wing (from Tutorials-2.5.0.2\B737_AVL_Tutorial\tut_mission_B737_AVL.py)
    
    # wing = SUAVE.Components.Wings.Main_Wing()
    # wing.tag = 'main_wing'
    
    # wing.aspect_ratio            = 10.18
    # wing.sweeps.quarter_chord    = 25 * Units.deg
    # wing.thickness_to_chord      = 0.1
    # wing.taper                   = 0.1
    # wing.spans.projected         = 34.32 * Units.meter
    # wing.chords.root             = 7.760 * Units.meter
    # wing.chords.tip              = 0.782 * Units.meter
    # wing.chords.mean_aerodynamic = 4.235 * Units.meter
    # wing.areas.reference         = 124.862 * Units['meters**2']  
    # wing.twists.root             = 4.0 * Units.degrees
    # wing.twists.tip              = 0.0 * Units.degrees
    # wing.origin                  = [[13.61 * Units.meter, 0, -1.27 * Units.meter]]
    # wing.vertical                = False
    # wing.symmetric               = True
    # wing.high_lift               = False
    # wing.dynamic_pressure_ratio  = 1.0

    # # add to vehicle
    # vehicle.append_component(wing)
    
    # NILS: MULTI-SECTION wing (from regression\scripts\Vehicles\Boeing_737.py)
    
    wing = SUAVE.Components.Wings.Main_Wing()
    wing.tag = 'main_wing'

    wing.aspect_ratio            = 10.18
    wing.sweeps.quarter_chord    = 25 * Units.deg
    wing.thickness_to_chord      = 0.1
    wing.taper                   = 0.1

    wing.spans.projected         = 34.32

    wing.chords.root             = 7.760 * Units.meter
    wing.chords.tip              = 0.782 * Units.meter
    wing.chords.mean_aerodynamic = 4.235 * Units.meter

    wing.areas.reference         = 124.862
    wing.areas.wetted            = 225.08
    
    wing.twists.root             = 4.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees

    wing.origin                  = [[13.61,0,-0.93]]
    wing.aerodynamic_center      = [0,0,0]   

    wing.vertical                = False
    wing.symmetric               = True
    wing.high_lift               = True

    wing.dynamic_pressure_ratio  = 1.0

    # Wing Segments
    root_airfoil                          = SUAVE.Components.Airfoils.Airfoil()
    root_airfoil.coordinate_file          = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737a.txt'
    segment                               = SUAVE.Components.Wings.Segment()
    segment.tag                           = 'Root'
    segment.percent_span_location         = 0.0
    segment.twist                         = 4. * Units.deg
    segment.root_chord_percent            = 1.
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 2.5 * Units.degrees
    segment.sweeps.quarter_chord          = 28.225 * Units.degrees
    segment.append_airfoil(root_airfoil)
    wing.append_segment(segment)

    yehudi_airfoil                        = SUAVE.Components.Airfoils.Airfoil()
    yehudi_airfoil.coordinate_file        = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737b.txt'
    segment                               = SUAVE.Components.Wings.Segment()
    segment.tag                           = 'Yehudi'
    segment.percent_span_location         = 0.324
    segment.twist                         = 0.047193 * Units.deg
    segment.root_chord_percent            = 0.5
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 5.5 * Units.degrees
    segment.sweeps.quarter_chord          = 25. * Units.degrees
    segment.append_airfoil(yehudi_airfoil)
    wing.append_segment(segment)

    mid_airfoil                           = SUAVE.Components.Airfoils.Airfoil()
    mid_airfoil.coordinate_file           = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737c.txt'
    segment                               = SUAVE.Components.Wings.Segment()
    segment.tag                           = 'Section_2'
    segment.percent_span_location         = 0.963
    segment.twist                         = 0.00258 * Units.deg
    segment.root_chord_percent            = 0.220
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 5.5 * Units.degrees
    segment.sweeps.quarter_chord          = 56.75 * Units.degrees
    segment.append_airfoil(mid_airfoil)
    wing.append_segment(segment)

    tip_airfoil                           =  SUAVE.Components.Airfoils.Airfoil()
    tip_airfoil.coordinate_file           = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737d.txt'
    segment                               = SUAVE.Components.Wings.Segment()
    segment.tag                           = 'Tip'
    segment.percent_span_location         = 1.
    segment.twist                         = 0. * Units.degrees
    segment.root_chord_percent            = 0.10077
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 0.
    segment.sweeps.quarter_chord          = 0.
    segment.append_airfoil(tip_airfoil)
    wing.append_segment(segment)
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing)    
    
    # NILS: the below ctrl surfaces will invoke the `if num_ctrl != 0: `
    # block in trunk\SUAVE\Methods\Aerodynamics\AVL\read_results.py
    
    # control surfaces -------------------------------------------
    slat                          = SUAVE.Components.Wings.Control_Surfaces.Slat()
    slat.tag                      = 'slat'
    slat.span_fraction_start      = 0.2
    slat.span_fraction_end        = 0.963
    slat.deflection               = 0.0 * Units.degrees
    slat.chord_fraction           = 0.075
    wing.append_control_surface(slat)

    flap                          = SUAVE.Components.Wings.Control_Surfaces.Flap()
    flap.tag                      = 'flap'
    flap.span_fraction_start      = 0.2
    flap.span_fraction_end        = 0.7
    flap.deflection               = 0.0 * Units.degrees
    flap.configuration_type       = 'double_slotted'
    flap.chord_fraction           = 0.30
    wing.append_control_surface(flap)

    aileron                       = SUAVE.Components.Wings.Control_Surfaces.Aileron()
    aileron.tag                   = 'aileron'
    aileron.span_fraction_start   = 0.7
    aileron.span_fraction_end     = 0.963
    aileron.deflection            = 0.0 * Units.degrees
    aileron.chord_fraction        = 0.16
    wing.append_control_surface(aileron)
    
    # add to vehicle
    vehicle.append_component(wing)

    # # ------------------------------------------------------------------        
    # #  Horizontal Stabilizer
    # # ------------------------------------------------------------------        
    
    # NILS: SINGLE-SECTION h-tail (from Tutorials-2.5.0.2\B737_AVL_Tutorial\tut_mission_B737_AVL.py)
    
    # wing = SUAVE.Components.Wings.Horizontal_Tail()
    # wing.tag = 'horizontal_stabilizer'
    
    # wing.aspect_ratio            = 6.16     
    # wing.sweeps.quarter_chord    = 40 * Units.deg
    # wing.thickness_to_chord      = 0.08
    # wing.taper                   = 0.2
    # wing.spans.projected         = 14.2 * Units.meter
    # wing.chords.root             = 4.7  * Units.meter
    # wing.chords.tip              = .955 * Units.meter
    # wing.chords.mean_aerodynamic = 3.0  * Units.meter
    # wing.areas.reference         = 32.488   * Units['meters**2']  
    # wing.twists.root             = 3.0 * Units.degrees
    # wing.twists.tip              = 3.0 * Units.degrees  
    # wing.origin                  = [[32.83 * Units.meter, 0 , 1.14 * Units.meter]]
    # wing.vertical                = False 
    # wing.symmetric               = True
    # wing.dynamic_pressure_ratio  = 0.9  
    
    # # add to vehicle
    # vehicle.append_component(wing)
    
    # NILS: MULTI-SECTION h-tail (from regression\scripts\Vehicles\Boeing_737.py)
    # NOTE: for the `trim_aircraft = True` in trunk\SUAVE\Analyses\Stability\AVL.py
    # to work, ALSO THESE MUST BE DEFINED AS MULTI-SEGMENT WINGS, otherwise
    # 'elevator' won't feature in `wing.control_surfaces` in
    # trunk\SUAVE\Methods\Flight_Dynamics\Dynamic_Stability\compute_dynamic_flight_modes.py
    
    wing = SUAVE.Components.Wings.Horizontal_Tail()
    wing.tag = 'horizontal_stabilizer'

    wing.aspect_ratio            = 4.99
    wing.sweeps.quarter_chord    = 28.2250 * Units.deg  
    wing.thickness_to_chord      = 0.08
    wing.taper                   = 0.3333 

    wing.spans.projected         = 14.4

    wing.chords.root             = 4.2731 
    wing.chords.tip              = 1.4243 
    wing.chords.mean_aerodynamic = 8.0

    wing.areas.reference         = 41.49
    wing.areas.exposed           = 59.354    # Exposed area of the horizontal tail
    wing.areas.wetted            = 71.81     # Wetted area of the horizontal tail
    wing.twists.root             = 3.0 * Units.degrees
    wing.twists.tip              = 3.0 * Units.degrees

    wing.origin                  = [[33.02,0,1.466]]
    wing.aerodynamic_center      = [0,0,0]

    wing.vertical                = False
    wing.symmetric               = True

    wing.dynamic_pressure_ratio  = 0.9

    # Wing Segments
    segment                        = SUAVE.Components.Wings.Segment()
    segment.tag                    = 'root_segment'
    segment.percent_span_location  = 0.0
    segment.twist                  = 0. * Units.deg
    segment.root_chord_percent     = 1.0
    segment.dihedral_outboard      = 8.63 * Units.degrees
    segment.sweeps.quarter_chord   = 28.2250  * Units.degrees 
    segment.thickness_to_chord     = .1
    wing.append_segment(segment)

    segment                        = SUAVE.Components.Wings.Segment()
    segment.tag                    = 'tip_segment'
    segment.percent_span_location  = 1.
    segment.twist                  = 0. * Units.deg
    segment.root_chord_percent     = 0.3333               
    segment.dihedral_outboard      = 0 * Units.degrees
    segment.sweeps.quarter_chord   = 0 * Units.degrees  
    segment.thickness_to_chord     = .1
    wing.append_segment(segment)
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing)        

    # control surfaces -------------------------------------------
    elevator                       = SUAVE.Components.Wings.Control_Surfaces.Elevator()
    elevator.tag                   = 'elevator'
    elevator.span_fraction_start   = 0.09
    elevator.span_fraction_end     = 0.92
    elevator.deflection            = 0.0  * Units.deg
    elevator.chord_fraction        = 0.3
    wing.append_control_surface(elevator)

    # add to vehicle
    vehicle.append_component(wing)
    
    # ------------------------------------------------------------------
    #   Vertical Stabilizer - NILS: treat as SINGLE-SECTION wing
    # ------------------------------------------------------------------
    
    wing = SUAVE.Components.Wings.Vertical_Tail()
    wing.tag = 'vertical_stabilizer'    

    wing.aspect_ratio            = 1.91
    wing.sweeps.quarter_chord    = 25. * Units.deg
    wing.thickness_to_chord      = 0.08
    wing.taper                   = 0.25
    wing.spans.projected         = 7.777 * Units.meter
    wing.chords.root             = 8.19  * Units.meter
    wing.chords.tip              = 0.95  * Units.meter
    wing.chords.mean_aerodynamic = 4.0   * Units.meter
    wing.areas.reference         = 27.316 * Units['meters**2']  
    wing.twists.root             = 0.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees  
    wing.origin                  = [[28.79 * Units.meter, 0, 1.54 * Units.meter]] # meters
    wing.vertical                = True 
    wing.symmetric               = False
    wing.t_tail                  = False
    wing.dynamic_pressure_ratio  = 1.0
        
    # add to vehicle
    vehicle.append_component(wing)

    # ------------------------------------------------------------------
    #  Fuselage
    # ------------------------------------------------------------------
    
    fuselage = SUAVE.Components.Fuselages.Fuselage()
    fuselage.tag = 'fuselage'
    
    # fuselage.number_coach_seats    = vehicle.passengers
    # fuselage.seats_abreast         = 6
    # fuselage.seat_pitch            = 1     * Units.meter
    # fuselage.fineness.nose         = 1.6
    # fuselage.fineness.tail         = 2.
    fuselage.lengths.nose          = 6.4   * Units.meter
    fuselage.lengths.tail          = 8.0   * Units.meter
    # fuselage.lengths.cabin         = 28.85 * Units.meter
    fuselage.lengths.total         = 38.02 * Units.meter
    # fuselage.lengths.fore_space    = 6.    * Units.meter
    # fuselage.lengths.aft_space     = 5.    * Units.meter
    fuselage.width                 = 3.74  * Units.meter
    fuselage.heights.maximum       = 3.74  * Units.meter
    # fuselage.effective_diameter    = 3.74     * Units.meter
    # fuselage.areas.side_projected  = 142.1948 * Units['meters**2'] 
    # fuselage.areas.wetted          = 446.718  * Units['meters**2'] 
    # fuselage.areas.front_projected = 12.57    * Units['meters**2'] 
    # fuselage.differential_pressure = 5.0e4 * Units.pascal # Maximum differential pressure
    
    # fuselage.heights.at_quarter_length          = 3.74 * Units.meter
    # fuselage.heights.at_three_quarters_length   = 3.65 * Units.meter
    # fuselage.heights.at_wing_root_quarter_chord = 3.74 * Units.meter
    
    # add to vehicle
    vehicle.append_component(fuselage)
    
    # =============================================================================
    # ------------------------------------------------------------------
    #   Nacelles
    # ------------------------------------------------------------------ 
    nacelle                       = SUAVE.Components.Nacelles.Nacelle()
    nacelle.tag                   = 'nacelle_1'
    nacelle.length                = 2.71
    nacelle.inlet_diameter        = 1.90
    nacelle.diameter              = 2.05
    nacelle.areas.wetted          = 1.1*np.pi*nacelle.diameter*nacelle.length
    nacelle.origin                = [[13.72, -4.86,-1.9]]
    nacelle.flow_through          = True  
    nacelle_airfoil               = SUAVE.Components.Airfoils.Airfoil() 
    nacelle_airfoil.naca_4_series_airfoil = '2410'
    nacelle.append_airfoil(nacelle_airfoil)

    nacelle_2                     = deepcopy(nacelle)
    nacelle_2.tag                 = 'nacelle_2'
    nacelle_2.origin              = [[13.72, 4.86,-1.9]]
    
    vehicle.append_component(nacelle)  
    vehicle.append_component(nacelle_2)
    # =============================================================================
    
    return vehicle

# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def configs_setup(vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------
    configs = SUAVE.Components.Configs.Config.Container()

    base_config = SUAVE.Components.Configs.Config(vehicle)
    base_config.tag = 'base'
    configs.append(base_config)

    return configs

#%%

if __name__ == '__main__':
    """
    This block replaces the default B737 geometry created in full_setup() with a
    clean-slate SUAVE.Vehicle() suitable for custom AVL geometry generation.
    
    SUAVE’s Data objects (including Vehicle and Config) use an internal diffing
    mechanism based on the attributes `_base` and `_diff`.  The AVL interface
    relies on this structure and will fail (e.g. AttributeError: 'Vehicle' object
    has no attribute '_base') if `_base` is missing or inconsistent.
    
    The default vehicle created in full_setup() is stored in
        configs.base            (a Config wrapping a Vehicle)
        analyses.configs.base.stability.geometry
    
    To wipe the B737 defaults correctly, the existing Vehicle must be replaced by
    a *new* Vehicle instance that also has its own `_base` Vehicle:
    
        aircraft = SUAVE.Vehicle()
        aircraft._base = SUAVE.Vehicle()
    
    Both configs.base and the AVL geometry reference must then point to this same
    new Vehicle instance.  This preserves SUAVE’s required structure while giving
    a completely empty aircraft definition.  The id() checks verify that the old
    Vehicle has been replaced everywhere AVL expects it.
    
    In short: SUAVE cannot operate without a valid `_base` attribute, so resetting
    to a clean aircraft must be done by constructing a new Vehicle + new _base and
    replacing both the config and analysis geometry references accordingly.
    
    
    print-statements to highlight this:
    
    id(analyses.configs.base.stability.geometry)
    Out[12]: 1532983067328

    id(analyses.configs.base.stability.geometry._base)
    Out[13]: 1532446661984

    id(configs.base)
    Out[14]: 1532983067328

    id(configs.base._base)
    Out[15]: 1532446661984
    """

    """
    configs.base._base.keys()
    Out[19]: dict_keys(['tag', 'fuselages', 'wings', 'networks', 'nacelles', 'systems', 'mass_properties', 'payload', 'costs', 'envelope', 'landing_gear', 'reference_area', 'passengers', 'performance'])

    configs.base.keys()
    Out[20]: dict_keys(['tag', 'fuselages', 'wings', 'networks', 'nacelles', 'systems', 'mass_properties', 'payload', 'costs', 'envelope', 'landing_gear', 'reference_area', 'passengers', 'performance', '_base', '_diff'])
    """

    """
    id(avl_object.geometry._base)
    Out[29]: 1532937625664

    id(configs.base._base)
    Out[30]: 1532937625664
    """
    
    # =============================================================================
    # ----------------------------
    # Command-line: study index
    # ----------------------------
    # Usage: python run_suave_avl_wrapper_nils.py study_idx
    # where study_idx is an integer 1..7
    if len(sys.argv) > 1:
        try:
            study_idx = int(sys.argv[1])
        except ValueError:
            raise SystemExit("Usage: python run_suave_avl_wrapper_nils.py study_idx (integer 1..7)")
    else:
        # default (choose sensible default) — change if you want no default
        study_idx = 6
    
    if not (1 <= study_idx <= 7):
        raise SystemExit("study_idx must be an integer between 1 and 7")
    # =============================================================================
    
    configs, analyses = full_setup()
    avl_object = analyses.configs.base.stability
    
    # NOTE: `._base` required to avoid `AttributeError: 'NoneType' object has no attribute 'items'`
    # To understand why I have to use `aircraft = avl_object.geometry._base` instead of
    # `aircraft = SUAVE.Vehicle()`, compare `print(id(aircraft))` with `print(id(vehicle))`
    # in `full_setup()` above - they are not the same object instance!
    
    aircraft_is_b737 = False  # False
    keep_b737_defaults = False  # True
    wing_has_segments = True  # False
    htail_has_segments = True  # False
    ac_segment = "regional"  # "narrowbody" or "regional"
    
    if not keep_b737_defaults:
        """
        NOTE: configs.base is a Config, not a Vehicle -> see configs_setup()
        """
        
        assert id(configs.base) == id(avl_object.geometry)
        id_old = id(configs.base)
        
        # Define clean-slate aircraft

        aircraft = SUAVE.Vehicle()   # new clean vehicle
        aircraft_base = SUAVE.Vehicle()   # base version for diffing
    
        # SUAVE requirement: `_base` holds the baseline configuration
        aircraft._base = aircraft_base
    
        # Attach this clean-slate "aircraft" to configs AND avl_object
        configs.base = aircraft
        avl_object.geometry = aircraft
        
        id_new = id(configs.base)
        assert id_new != id_old
        assert id(configs.base) == id(avl_object.geometry)
        
        # Check that configs.base has been wiped indeed
        # print('aircraft.mass_properties.max_zero_fuel =', aircraft.mass_properties.max_zero_fuel)
        # print('aircraft.wings.keys() =', aircraft.wings.keys())
        # sys.exit()
        
        # Read geometry data from .csv file
        if aircraft_is_b737:
            df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\B737_AVL_Tutorial\suave_avl_wrapper_tasopt_inputs_b737.csv')
            aircraft.tag = 'Boeing_737800'
            t_tail_bool = False
        elif not aircraft_is_b737:
            date_str = date.today().strftime("%d%m%y")
            if ac_segment == "narrowbody":
                df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_narrowbody_kerosene_090126.csv')
                # df_geom = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_narrowbody_kerosene_{date_str}.csv')
                df_mass = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_narrowbody_kerosene_090126_6.csv')
                # df_mass = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_narrowbody_kerosene_{date_str}_{study_idx}.csv')
                aircraft.tag = 'Airbus_A220-100'
                t_tail_bool = False
            elif ac_segment == "regional":
                df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_kerosene_090126.csv')
                # df_geom = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_kerosene_{date_str}.csv')
                df_mass = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_regional_kerosene_090126_6.csv')
                # df_mass = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_regional_kerosene_{date_str}_{study_idx}.csv')
                aircraft.tag = 'ATR_72-600'
                t_tail_bool = True
            
        # sys.exit('Stop here.')
            
        # Aircraft mass properties
        
        for counter, row in df_mass.iterrows():
            print('counter =', counter)
            
            # =============================================================================
            configs, analyses = full_setup()
            avl_object = analyses.configs.base.stability
            
            # Define clean-slate aircraft

            aircraft = SUAVE.Vehicle()   # new clean vehicle
            aircraft_base = SUAVE.Vehicle()   # base version for diffing
        
            # SUAVE requirement: `_base` holds the baseline configuration
            aircraft._base = aircraft_base
        
            # Attach this clean-slate "aircraft" to configs AND avl_object
            configs.base = aircraft
            avl_object.geometry = aircraft
            
            if ac_segment == "narrowbody":
                aircraft.tag = 'Airbus_A220-100'
            elif ac_segment == "regional":
                aircraft.tag = 'ATR_72-600'
            # =============================================================================
            
            sigma_fcs = row['sigma_fcs']
            span_loc = row['span_loc']
            fcs_loc = row['fcs_loc']
            wing_frac = row['wing_frac']
            nacelle_frac = row['nacelle_frac']
            
            aircraft.mass_properties.center_of_gravity[0][0] = row['x_cg']
            aircraft.mass_properties.center_of_gravity[0][1] = row['y_cg']
            aircraft.mass_properties.center_of_gravity[0][2] = row['z_cg']
            aircraft.mass_properties.mass = row['mass']
            aircraft.mass_properties.max_takeoff = row['max_takeoff']
            aircraft.mass_properties.takeoff = row['takeoff']
            aircraft.mass_properties.max_zero_fuel = row['max_zero_fuel']
            moments_of_inertia = aircraft.mass_properties.moments_of_inertia.tensor
            moments_of_inertia[0][0] = row['Ixx']
            moments_of_inertia[1][1] = row['Iyy']
            moments_of_inertia[2][2] = row['Izz']
            moments_of_inertia[0][1] = row['Ixy']
            moments_of_inertia[1][2] = row['Iyz']
            moments_of_inertia[2][0] = row['Izx']
        
            # Fuselage
            
            fuselage = SUAVE.Components.Fuselages.Fuselage()
            fuselage.tag = 'fuselage'
            
            # for suave_body in aircraft.fuselages:  # write_geometry.py
            fuselage.lengths.total = df_geom['body_lengths_total'][0]
            fuselage.lengths.nose = df_geom['body_lengths_nose'][0]
            fuselage.lengths.tail = df_geom['body_lengths_tail'][0]
            fuselage.width = df_geom['body_widths_maximum'][0]
            fuselage.heights.maximum = df_geom['body_heights_maximum'][0]
            
            aircraft.append_component(fuselage)
                
            # Main wing
            
            wing = SUAVE.Components.Wings.Main_Wing()
            wing.tag = 'main_wing'
            
            if not wing_has_segments:
                
                wing.aspect_ratio = df_geom['wing_aspect_ratio'][0]
                wing.sweeps.quarter_chord = df_geom['wing_sweeps_quarter_chord'][0]
                wing.thickness_to_chord = df_geom['wing_thickness_to_chord'][0]
                wing.taper = df_geom['wing_taper'][0]
                wing.spans.projected = df_geom['wing_spans_projected'][0]
                wing.chords.root = df_geom['wing_chords_root'][0]
                wing.chords.tip = df_geom['wing_chords_tip'][0]
                wing.chords.mean_aerodynamic = df_geom['wing_chords_mean_aerodynamic'][0]
                wing.areas.reference = df_geom['wing_areas_reference'][0]
                wing.twists.root = df_geom['wing_twists_root'][0]
                wing.twists.tip = df_geom['wing_twists_tip'][0]
                wing.origin = [[
                    df_geom['wing_origin_x'][0],
                    df_geom['wing_origin_y'][0],
                    df_geom['wing_origin_z'][0],
                ]]
                wing.vertical = df_geom['wing_vertical'][0]
                wing.symmetric = df_geom['wing_symmetric'][0]
                wing.high_lift = df_geom['wing_high_lift'][0]
                wing.dihedral = df_geom['wing_dihedral'][0]
                
            elif wing_has_segments:  # <<< @NILS: CHECK UNITS! >>>
                
                wing.aspect_ratio = df_geom['wing_aspect_ratio'][0]
                wing.sweeps.quarter_chord = df_geom['wing_sweeps_quarter_chord'][0]
                # wing.thickness_to_chord = df_geom['wing_thickness_to_chord'][0]  # defined per segment below
                # wing.taper = df_geom['wing_taper'][0]  # defined per segment below
                wing.spans.projected = df_geom['wing_spans_projected'][0]
                wing.chords.root = df_geom['wing_chords_root'][0]
                wing.chords.tip = df_geom['wing_chords_tip'][0]
                wing.chords.mean_aerodynamic = df_geom['wing_chords_mean_aerodynamic'][0]
                wing.areas.reference = df_geom['wing_areas_reference'][0]
                wing.twists.root = df_geom['wing_twists_root'][0]
                wing.twists.tip = df_geom['wing_twists_tip'][0]
                wing.origin = [[
                    df_geom['wing_origin_x'][0],
                    df_geom['wing_origin_y'][0],
                    df_geom['wing_origin_z'][0],
                ]]
                wing.vertical = df_geom['wing_vertical'][0]
                wing.symmetric = df_geom['wing_symmetric'][0]
                wing.high_lift = df_geom['wing_high_lift'][0]
                # wing.dihedral = df_geom['wing_dihedral'][0]  # defined per segment below
                
                if aircraft_is_b737:  # 4 segments (see ..\..\regression\scripts\Vehicles\Boeing_737.py)
    
                    root_airfoil                          = SUAVE.Components.Airfoils.Airfoil()
                    root_airfoil.coordinate_file          = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737a.txt'
                    segment                               = SUAVE.Components.Wings.Segment()
                    segment.tag                           = 'Root'
                    segment.percent_span_location         = df_geom['wing_root_percent_span_location'][0]
                    segment.twist                         = df_geom['wing_root_twist'][0]
                    segment.root_chord_percent            = df_geom['wing_root_root_chord_percent'][0]
                    segment.thickness_to_chord            = df_geom['wing_root_thickness_to_chord'][0]
                    segment.dihedral_outboard             = df_geom['wing_root_dihedral_outboard'][0]
                    segment.sweeps.quarter_chord          = df_geom['wing_root_sweeps_quarter_chord'][0]
                    segment.append_airfoil(root_airfoil)
                    wing.append_segment(segment)
        
                    yehudi_airfoil                       = SUAVE.Components.Airfoils.Airfoil()
                    yehudi_airfoil.coordinate_file       = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737b.txt'
                    segment                               = SUAVE.Components.Wings.Segment()
                    segment.tag                           = 'Yehudi'
                    segment.percent_span_location         = df_geom['wing_yehudi_percent_span_location'][0]
                    segment.twist                         = df_geom['wing_yehudi_twist'][0]
                    segment.root_chord_percent            = df_geom['wing_yehudi_root_chord_percent'][0]
                    segment.thickness_to_chord            = df_geom['wing_yehudi_thickness_to_chord'][0]
                    segment.dihedral_outboard             = df_geom['wing_yehudi_dihedral_outboard'][0]
                    segment.sweeps.quarter_chord          = df_geom['wing_yehudi_sweeps_quarter_chord'][0]
                    segment.append_airfoil(yehudi_airfoil)
                    wing.append_segment(segment)
        
                    section2_airfoil                      =  SUAVE.Components.Airfoils.Airfoil()
                    section2_airfoil.coordinate_file      = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737c.txt'
                    segment                               = SUAVE.Components.Wings.Segment()
                    segment.tag                           = 'Section 2'
                    segment.percent_span_location         = df_geom['wing_section2_percent_span_location'][0]
                    segment.twist                         = df_geom['wing_section2_twist'][0]
                    segment.root_chord_percent            = df_geom['wing_section2_root_chord_percent'][0]
                    segment.thickness_to_chord            = df_geom['wing_section2_thickness_to_chord'][0]
                    segment.dihedral_outboard             = df_geom['wing_section2_dihedral_outboard'][0]
                    segment.sweeps.quarter_chord          = df_geom['wing_section2_sweeps_quarter_chord'][0]
                    segment.append_airfoil(section2_airfoil)
                    wing.append_segment(segment)
                    
                    tip_airfoil                      =  SUAVE.Components.Airfoils.Airfoil()
                    tip_airfoil.coordinate_file      = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737d.txt'
                    segment                               = SUAVE.Components.Wings.Segment()
                    segment.tag                           = 'Tip'
                    segment.percent_span_location         = df_geom['wing_tip_percent_span_location'][0]
                    segment.twist                         = df_geom['wing_tip_twist'][0]
                    segment.root_chord_percent            = df_geom['wing_tip_root_chord_percent'][0]
                    segment.thickness_to_chord            = df_geom['wing_tip_thickness_to_chord'][0]
                    segment.dihedral_outboard             = df_geom['wing_tip_dihedral_outboard'][0]
                    segment.sweeps.quarter_chord          = df_geom['wing_tip_sweeps_quarter_chord'][0]
                    segment.append_airfoil(tip_airfoil)
                    wing.append_segment(segment)
                    
                elif not aircraft_is_b737:  # 3 segments (see TASOPT.jl)
                    
                    ###
                    import re

                    section_indices = sorted(
                        {int(re.search(r"sect_(\d+)_", c).group(1))
                         for c in df_geom.columns
                         if c.startswith("sect_")}
                    )
                    
                    for i in section_indices:
                        segment_airfoil                          = SUAVE.Components.Airfoils.Airfoil()
                        segment_airfoil.coordinate_file          = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737a.txt'
                        segment                               = SUAVE.Components.Wings.Segment()
                        segment.tag                           = df_geom[f'sect_{i}_label'][0]
                        print('segment.tag =', segment.tag)
                        segment.percent_span_location         = df_geom[f'sect_{i}_percent_span_location'][0]
                        segment.twist                         = 0.0  # not modelled in TASOPT.jl
                        segment.root_chord_percent            = df_geom[f'sect_{i}_root_chord_percent'][0]
                        segment.thickness_to_chord            = df_geom[f'sect_{i}_thickness_to_chord'][0]
                        segment.dihedral_outboard             = df_geom[f'sect_{i}_dihedral_outboard'][0]
                        segment.sweeps.quarter_chord          = df_geom[f'sect_{i}_sweeps_quarter_chord'][0]
                        segment.append_airfoil(segment_airfoil)
                        wing.append_segment(segment)
                    # sys.exit('Done.')
                    ###
                    
                    # center_airfoil                          = SUAVE.Components.Airfoils.Airfoil()
                    # center_airfoil.coordinate_file          = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737a.txt'
                    # segment                               = SUAVE.Components.Wings.Segment()
                    # segment.tag                           = 'Root'
                    # segment.percent_span_location         = df_geom['wing_center_percent_span_location'][0]
                    # segment.twist                         = df_geom['wing_center_twist'][0]
                    # segment.root_chord_percent            = df_geom['wing_center_root_chord_percent'][0]
                    # segment.thickness_to_chord            = df_geom['wing_center_thickness_to_chord'][0]
                    # segment.dihedral_outboard             = df_geom['wing_center_dihedral_outboard'][0]
                    # segment.sweeps.quarter_chord          = df_geom['wing_center_sweeps_quarter_chord'][0]
                    # segment.append_airfoil(center_airfoil)
                    # wing.append_segment(segment)
        
                    # inboard_airfoil                       = SUAVE.Components.Airfoils.Airfoil()
                    # inboard_airfoil.coordinate_file       = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737b.txt'
                    # segment                               = SUAVE.Components.Wings.Segment()
                    # segment.tag                           = 'Yehudi'
                    # segment.percent_span_location         = df_geom['wing_inboard_percent_span_location'][0]
                    # segment.twist                         = df_geom['wing_inboard_twist'][0]
                    # segment.root_chord_percent            = df_geom['wing_inboard_root_chord_percent'][0]
                    # segment.thickness_to_chord            = df_geom['wing_inboard_thickness_to_chord'][0]
                    # segment.dihedral_outboard             = df_geom['wing_inboard_dihedral_outboard'][0]
                    # segment.sweeps.quarter_chord          = df_geom['wing_inboard_sweeps_quarter_chord'][0]
                    # segment.append_airfoil(inboard_airfoil)
                    # wing.append_segment(segment)
        
                    # outboard_airfoil                      =  SUAVE.Components.Airfoils.Airfoil()
                    # outboard_airfoil.coordinate_file      = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737c.txt'
                    # segment                               = SUAVE.Components.Wings.Segment()
                    # segment.tag                           = 'Section 2'
                    # segment.percent_span_location         = df_geom['wing_outboard_percent_span_location'][0]
                    # segment.twist                         = df_geom['wing_outboard_twist'][0]
                    # segment.root_chord_percent            = df_geom['wing_outboard_root_chord_percent'][0]
                    # segment.thickness_to_chord            = df_geom['wing_outboard_thickness_to_chord'][0]
                    # segment.dihedral_outboard             = df_geom['wing_outboard_dihedral_outboard'][0]
                    # segment.sweeps.quarter_chord          = df_geom['wing_outboard_sweeps_quarter_chord'][0]
                    # segment.append_airfoil(outboard_airfoil)
                    # wing.append_segment(segment)
                    
                    # tip_airfoil                      =  SUAVE.Components.Airfoils.Airfoil()
                    # tip_airfoil.coordinate_file      = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737d.txt'
                    # segment                               = SUAVE.Components.Wings.Segment()
                    # segment.tag                           = 'Tip'
                    # segment.percent_span_location         = df_geom['wing_tip_percent_span_location'][0]
                    # segment.twist                         = df_geom['wing_tip_twist'][0]
                    # segment.root_chord_percent            = df_geom['wing_tip_root_chord_percent'][0]
                    # segment.thickness_to_chord            = df_geom['wing_tip_thickness_to_chord'][0]
                    # segment.dihedral_tip             = df_geom['wing_tip_dihedral_tip'][0]
                    # segment.sweeps.quarter_chord          = df_geom['wing_tip_sweeps_quarter_chord'][0]
                    # segment.append_airfoil(tip_airfoil)
                    # wing.append_segment(segment)
                    
                # # control surfaces -------------------------------------------
                # slat                          = SUAVE.Components.Wings.Control_Surfaces.Slat()
                # slat.tag                      = 'slat'
                # slat.span_fraction_start      = df_geom['wing_slat_span_fraction_start'][0]
                # slat.span_fraction_end        = df_geom['wing_slat_span_fraction_end'][0]
                # slat.deflection               = df_geom['wing_slat_deflection'][0]
                # slat.chord_fraction           = df_geom['wing_slat_chord_fraction'][0]
                # wing.append_control_surface(slat)
    
                # flap                          = SUAVE.Components.Wings.Control_Surfaces.Flap()
                # flap.tag                      = 'flap'
                # flap.span_fraction_start      = df_geom['wing_flap_span_fraction_start'][0]
                # flap.span_fraction_end        = df_geom['wing_flap_span_fraction_end'][0]
                # flap.deflection               = df_geom['wing_flap_deflection'][0]
                # flap.configuration_type       = df_geom['wing_flap_configuration_type'][0]
                # flap.chord_fraction           = df_geom['wing_flap_chord_fraction'][0]
                # wing.append_control_surface(flap)
    
                # aileron                       = SUAVE.Components.Wings.Control_Surfaces.Aileron()
                # aileron.tag                   = 'aileron'
                # aileron.span_fraction_start   = df_geom['wing_aileron_span_fraction_start'][0]
                # aileron.span_fraction_end     = df_geom['wing_aileron_span_fraction_end'][0]
                # aileron.deflection            = df_geom['wing_aileron_deflection'][0]
                # aileron.chord_fraction        = df_geom['wing_aileron_chord_fraction'][0]
                # wing.append_control_surface(aileron)
                
            aircraft.append_component(wing)
            
            # Horizontal stabiliser
            
            wing = SUAVE.Components.Wings.Horizontal_Tail()
            wing.tag = 'horizontal_stabilizer'
            
            if not htail_has_segments:
            
                wing.aspect_ratio = df_geom['htail_aspect_ratio'][0]
                wing.sweeps.quarter_chord = df_geom['htail_sweeps_quarter_chord'][0]
                wing.thickness_to_chord = df_geom['htail_thickness_to_chord'][0]
                wing.taper = df_geom['htail_taper'][0]
                wing.spans.projected = df_geom['htail_spans_projected'][0]
                wing.chords.root = df_geom['htail_chords_root'][0]
                wing.chords.tip = df_geom['htail_chords_tip'][0]
                # wing.chords.mean_aerodynamic = df_geom['htail_chords_mean_aerodynamic'][0]
                wing.areas.reference = df_geom['htail_areas_reference'][0]
                wing.twists.root = df_geom['htail_twists_root'][0]
                wing.twists.tip = df_geom['htail_twists_tip'][0]
                wing.origin = [[
                    df_geom['htail_origin_x'][0],
                    df_geom['htail_origin_y'][0],
                    df_geom['htail_origin_z'][0],
                ]]
                wing.vertical = df_geom['htail_vertical'][0]
                wing.symmetric = df_geom['htail_symmetric'][0]
                wing.high_lift = df_geom['htail_high_lift'][0]
                wing.dihedral = df_geom['htail_dihedral'][0]
                
                
            elif htail_has_segments:
                
                wing.aspect_ratio = df_geom['htail_aspect_ratio'][0]
                wing.sweeps.quarter_chord = df_geom['htail_sweeps_quarter_chord'][0]
                # wing.thickness_to_chord = df_geom['htail_thickness_to_chord'][0]  # defined per segment below
                # wing.taper = df_geom['htail_taper'][0]  # defined per segment below
                wing.spans.projected = df_geom['htail_spans_projected'][0]
                wing.chords.root = df_geom['htail_chords_root'][0]
                wing.chords.tip = df_geom['htail_chords_tip'][0]
                # wing.chords.mean_aerodynamic = df_geom['htail_chords_mean_aerodynamic'][0]
                wing.areas.reference = df_geom['htail_areas_reference'][0]
                wing.twists.root = df_geom['htail_twists_root'][0]
                wing.twists.tip = df_geom['htail_twists_tip'][0]
                wing.origin = [[
                    df_geom['htail_origin_x'][0],
                    df_geom['htail_origin_y'][0],
                    df_geom['htail_origin_z'][0],
                ]]
                wing.vertical = df_geom['htail_vertical'][0]
                wing.symmetric = df_geom['htail_symmetric'][0]
                wing.high_lift = df_geom['htail_high_lift'][0]
                # wing.dihedral = df_geom['htail_dihedral'][0]  # defined per segment below
                
                segment                               = SUAVE.Components.Wings.Segment()
                segment.tag                           = 'root_segment'
                segment.percent_span_location         = df_geom['htail_root_percent_span_location'][0]
                segment.twist                         = df_geom['htail_root_twist'][0]
                segment.root_chord_percent            = df_geom['htail_root_root_chord_percent'][0]
                segment.thickness_to_chord            = df_geom['htail_root_thickness_to_chord'][0]
                segment.dihedral_outboard             = df_geom['htail_root_dihedral_outboard'][0]
                segment.sweeps.quarter_chord          = df_geom['htail_root_sweeps_quarter_chord'][0]
                wing.append_segment(segment)
                
                segment                               = SUAVE.Components.Wings.Segment()
                segment.tag                           = 'tip_segment'
                segment.percent_span_location         = df_geom['htail_tip_percent_span_location'][0]
                segment.twist                         = df_geom['htail_tip_twist'][0]
                segment.root_chord_percent            = df_geom['htail_tip_root_chord_percent'][0]
                segment.thickness_to_chord            = df_geom['htail_tip_thickness_to_chord'][0]
                segment.dihedral_outboard             = df_geom['htail_tip_dihedral_outboard'][0]
                segment.sweeps.quarter_chord          = df_geom['htail_tip_sweeps_quarter_chord'][0]
                wing.append_segment(segment)
                
                # control surfaces -------------------------------------------
                elevator                       = SUAVE.Components.Wings.Control_Surfaces.Elevator()
                elevator.tag                   = 'elevator'
                elevator.span_fraction_start   = df_geom['htail_elevator_span_fraction_start'][0]
                elevator.span_fraction_end     = df_geom['htail_elevator_span_fraction_end'][0]
                elevator.deflection            = df_geom['htail_elevator_deflection'][0]
                elevator.chord_fraction        = df_geom['htail_elevator_chord_fraction'][0]
                wing.append_control_surface(elevator)
                
            aircraft.append_component(wing)
            
            # Vertical stabiliser
            wing = SUAVE.Components.Wings.Vertical_Tail()
            wing.tag = 'vertical_stabilizer'
            wing.aspect_ratio = df_geom['vtail_aspect_ratio'][0]
            wing.sweeps.quarter_chord = df_geom['vtail_sweeps_quarter_chord'][0]
            wing.thickness_to_chord = df_geom['vtail_thickness_to_chord'][0]
            wing.taper = df_geom['vtail_taper'][0]
            wing.spans.projected = df_geom['vtail_spans_projected'][0]
            wing.chords.root = df_geom['vtail_chords_root'][0]
            wing.chords.tip = df_geom['vtail_chords_tip'][0]
            # wing.chords.mean_aerodynamic = df_geom['vtail_chords_mean_aerodynamic'][0]
            wing.areas.reference = df_geom['vtail_areas_reference'][0]
            wing.twists.root = df_geom['vtail_twists_root'][0]
            wing.twists.tip = df_geom['vtail_twists_tip'][0]
            wing.origin = [[
                df_geom['vtail_origin_x'][0],
                df_geom['vtail_origin_y'][0],
                df_geom['vtail_origin_z'][0],
            ]]
            wing.vertical = df_geom['vtail_vertical'][0]
            wing.symmetric = df_geom['vtail_symmetric'][0]
            wing.high_lift = df_geom['vtail_high_lift'][0]
            wing.dihedral = df_geom['vtail_dihedral'][0]
            wing.t_tail = t_tail_bool
            aircraft.append_component(wing)
            
            # Nacelles
            '''
            nacelle = SUAVE.Components.Nacelles.Nacelle()
            # nacelle.tag = 'nacelle_1'
            nacelle.tag = 'nacelle'
            nacelle.length = df_geom['nacelle_length'][0]
            nacelle.inlet_diameter = df_geom['nacelle_inlet_diameter'][0]
            nacelle.diameter = df_geom['nacelle_diameter'][0]
            nacelle.areas.wetted = df_geom['nacelle_areas_wetted'][0]
            # nacelle.origin = df_geom['nacelle_origin'].to_numpy()
            nacelle.origin = np.array([
                ast.literal_eval(row)[0]
                for row in df_geom['nacelle_origin']
            ])  # NILS
            nacelle.flow_through = df_geom['nacelle_flow_through'][0]
            nacelle_airfoil = SUAVE.Components.Airfoils.Airfoil() 
            nacelle_airfoil.naca_4_series_airfoil = '2410'
            nacelle.append_airfoil(nacelle_airfoil)
            
            ###
            nacelle_2 = deepcopy(nacelle)
            nacelle_2.tag = 'nacelle_2'
            nacelle_2.origin = np.array([
                ast.literal_eval(row)[1]
                for row in df_geom['nacelle_origin']
            ])  # NILS
            
            # nacelle_3 = deepcopy(nacelle)
            # nacelle_3.tag = 'nacelle_3'
            # nacelle_3.origin = np.array([
            #     ast.literal_eval(row)[2]
            #     for row in df_geom['nacelle_origin']
            # ])  # NILS
            
            # nacelle_4 = deepcopy(nacelle)
            # nacelle_4.tag = 'nacelle_4'
            # nacelle_4.origin = np.array([
            #     ast.literal_eval(row)[3]
            #     for row in df_geom['nacelle_origin']
            # ])  # NILS
            
            nacelle_5 = deepcopy(nacelle_2)
            nacelle_5.tag = 'nacelle_5'
            nacelle_5.origin[1] *= -1
            
            # nacelle_6 = deepcopy(nacelle_3)
            # nacelle_6.tag = 'nacelle_6'
            # nacelle_6.origin[1] *= -1
            
            # nacelle_7 = deepcopy(nacelle_4)
            # nacelle_7.tag = 'nacelle_7'
            # nacelle_7.origin[1] *= -1
            ###
            
            aircraft.append_component(nacelle)
            aircraft.append_component(nacelle_2)
            # aircraft.append_component(nacelle_3)
            # aircraft.append_component(nacelle_4)
            aircraft.append_component(nacelle_5)
            # aircraft.append_component(nacelle_6)
            # aircraft.append_component(nacelle_7)
            
            # print(aircraft.nacelles.nacelle.Airfoil)
            # sys.exit('Stop.')
            '''
            #####
            # Number of nacelles (must be even)
            N_nacelles = len(ast.literal_eval(df_geom['nacelle_origin'].iloc[0]))
            assert N_nacelles % 2 == 0, "Number of nacelles must be even"
            
            N_half = N_nacelles // 2
            
            # --- Base nacelle ---
            base_nacelle = SUAVE.Components.Nacelles.Nacelle()
            base_nacelle.tag = 'nacelle'
            base_nacelle.length = df_geom['nacelle_length'][0]
            base_nacelle.inlet_diameter = df_geom['nacelle_inlet_diameter'][0]
            base_nacelle.diameter = df_geom['nacelle_diameter'][0]
            base_nacelle.areas.wetted = df_geom['nacelle_areas_wetted'][0]
            base_nacelle.flow_through = df_geom['nacelle_flow_through'][0]
            
            nacelle_airfoil = SUAVE.Components.Airfoils.Airfoil()
            nacelle_airfoil.naca_4_series_airfoil = '2410'
            base_nacelle.append_airfoil(nacelle_airfoil)
            
            # --- Create nacelles ---
            nacelles = []
            
            # Positive-Y side
            for i in range(N_half):
                nac = deepcopy(base_nacelle)
                nac.tag = 'nacelle' if i == 0 else f'nacelle_{i+1}'
                nac.origin = np.array([
                    ast.literal_eval(row)[i]
                    for row in df_geom['nacelle_origin']
                ])
                nacelles.append(nac)
            
            # Mirrored negative-Y side
            for i in range(N_half):
                nac = deepcopy(nacelles[i])
                nac.tag = f'nacelle_{i + 1 + N_half}'
                nac.origin[1] *= -1
                nacelles.append(nac)
            
            # --- Append to aircraft ---
            for nac in nacelles:
                aircraft.append_component(nac)
            #####
            
            # ###
            # hdisk = np.pi / 4 * Dprop
            # ###
            
            # =============================================================================
            avl_object.sample_training(
                study_idx=study_idx, counter=counter,
                sigma_fcs = sigma_fcs,
                span_loc = span_loc,
                fcs_loc = fcs_loc,
                wing_frac = wing_frac,
                nacelle_frac = nacelle_frac,
            )
            # =============================================================================
            
        # elif keep_b737_defaults:
            
        #     # aircraft = avl_object.geometry._base  # instance of SUAVE.Vehicle()
        #     # print('id(aircraft) =', id(aircraft))
        #     # # OR
        #     # aircraft = configs.base._base  # instance of SUAVE.Vehicle()
        #     # print('id(aircraft) =', id(aircraft))
            
        #     pass
            
        # sys.exit('Stop.')
        
        #%% Run sample_training() only
        
        # # avl_object.sample_training()
        # # =============================================================================
        # avl_object.sample_training(study_idx=study_idx, counter=counter)
        # # =============================================================================


    #%% OUTDATED (TO BE DELETED EVENTUALLY)
    
    # for tag, suave_wing in aircraft.wings.items():  # write_geometry.py
    #     # suave_wing.spans.projected = df['wing_spans_projected'][0]
    #     # suave_wing.origin = [[
    #     #     df['wing_origin_x'][0],
    #     #     df['wing_origin_y'][0],
    #     #     df['wing_origin_z'][0],
    #     # ]]
    #     # suave_wing.dihedral = df['wing_dihedral'][0]
        
    #     # for i in range(3):
    #     #     suave_wing.Segments.append(SUAVE.Components.Wings.Segment())
    #     # center_segment = suave_wing.Segments[0]
    #     # inboard_segment = suave_wing.Segments[1]
    #     # outboard_segment = suave_wing.Segments[2]
                
    #     # center_segment.sweeps.leading_edge = df['wing_segment_sweep_leading_edge_center'][0]
    #     # center_segment.root_chord_percent = df['wing_segment_root_chord_percent_center'][0]
    #     # center_segment.percent_span_location = df['wing_segment_percent_span_location_center'][0]
    #     # center_segment.sweeps.quarter_chord = df['wing_segment_sweep_quarter_chord_center'][0]
    #     # center_segment.twist = df['wing_segment_twist_center'][0]
        
    #     # inboard_segment.sweeps.leading_edge = df['wing_segment_sweep_leading_edge_inboard'][0]
    #     # inboard_segment.root_chord_percent = df['wing_segment_root_chord_percent_inboard'][0]
    #     # inboard_segment.percent_span_location = df['wing_segment_percent_span_location_inboard'][0]
    #     # inboard_segment.sweeps.quarter_chord = df['wing_segment_sweep_quarter_chord_inboard'][0]
    #     # inboard_segment.twist = df['wing_segment_twist_inboard'][0]
        
    #     # outboard_segment.sweeps.leading_edge = df['wing_segment_sweep_leading_edge_outboard'][0]
    #     # outboard_segment.root_chord_percent = df['wing_segment_root_chord_percent_outboard'][0]
    #     # outboard_segment.percent_span_location = df['wing_segment_percent_span_location_outboard'][0]
    #     # outboard_segment.sweeps.quarter_chord = df['wing_segment_sweep_quarter_chord_outboard'][0]
    #     # outboard_segment.twist = df['wing_segment_twist_outboard'][0]
    
    # =============================================================================
    
    # if tag == 'main_wing':
        
    # elif tag == 'horizontal_stabilizer':
        
    # elif tag == 'vertical_stabilizer':
        
    # root_airfoil                          = SUAVE.Components.Airfoils.Airfoil()
    # root_airfoil.coordinate_file          = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737a.txt'
    # segment                               = SUAVE.Components.Wings.Segment()
    # segment.tag                           = 'root'
    # segment.percent_span_location         = df['wing_root_percent_span_location'][0]
    # segment.twist                         = df['wing_root_twist'][0]
    # segment.root_chord_percent            = df['wing_root_root_chord_percent'][0]
    # segment.thickness_to_chord            = df['wing_root_thickness_to_chord'][0]
    # segment.dihedral_outboard             = df['wing_root_dihedral_outboard'][0]
    # segment.sweeps.quarter_chord          = df['wing_root_sweeps_quarter_chord'][0]
    # segment.append_airfoil(root_airfoil)
    # wing.append_segment(segment)

    # inboard_airfoil                       = SUAVE.Components.Airfoils.Airfoil()
    # inboard_airfoil.coordinate_file       = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737c.txt'
    # segment                               = SUAVE.Components.Wings.Segment()
    # segment.tag                           = 'inboard'
    # segment.percent_span_location         = df['wing_inboard_percent_span_location'][0]
    # segment.twist                         = df['wing_inboard_twist'][0]
    # segment.root_chord_percent            = df['wing_inboard_root_chord_percent'][0]
    # segment.thickness_to_chord            = df['wing_inboard_thickness_to_chord'][0]
    # segment.dihedral_outboard             = df['wing_inboard_dihedral_outboard'][0]
    # segment.sweeps.quarter_chord          = df['wing_inboard_sweeps_quarter_chord'][0]
    # segment.append_airfoil(inboard_airfoil)
    # wing.append_segment(segment)

    # outboard_airfoil                      =  SUAVE.Components.Airfoils.Airfoil()
    # outboard_airfoil.coordinate_file      = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737d.txt'
    # segment                               = SUAVE.Components.Wings.Segment()
    # segment.tag                           = 'outboard'
    # segment.percent_span_location         = df['wing_outboard_percent_span_location'][0]
    # segment.twist                         = df['wing_outboard_twist'][0]
    # segment.root_chord_percent            = df['wing_outboard_root_chord_percent'][0]
    # segment.thickness_to_chord            = df['wing_outboard_thickness_to_chord'][0]
    # segment.dihedral_outboard             = df['wing_outboard_dihedral_outboard'][0]
    # segment.sweeps.quarter_chord          = df['wing_outboard_sweeps_quarter_chord'][0]
    # segment.append_airfoil(outboard_airfoil)
    # wing.append_segment(segment)
    
    
    
    