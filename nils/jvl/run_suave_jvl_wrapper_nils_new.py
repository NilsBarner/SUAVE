r"""
This file is adapted from tut_mission_B737_AVL.py
to contain the minimum set of SUAVE methods needed to
perform a static and dynamic stability analysis on
aircraft geometry data imported as a .csv from TASOPT.jl.

NOTE: the definition of a wing segment always requires at
least two instances of SUAVE.Components.Wings.Segment() -
one at segment.percent_span_location = 0.0 and
one at segment.percent_span_location = 1.0, so that the
wing itself can be interpolated in between.


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

configs.base._base.keys()
Out[19]: dict_keys(['tag', 'fuselages', 'wings', 'networks', 'nacelles', 'systems', 'mass_properties', 'payload', 'costs', 'envelope', 'landing_gear', 'reference_area', 'passengers', 'performance'])

configs.base.keys()
Out[20]: dict_keys(['tag', 'fuselages', 'wings', 'networks', 'nacelles', 'systems', 'mass_properties', 'payload', 'costs', 'envelope', 'landing_gear', 'reference_area', 'passengers', 'performance', '_base', '_diff'])

id(avl_object.geometry._base)
Out[29]: 1532937625664

id(configs.base._base)
Out[30]: 1532937625664


# NOTE: `._base` required to avoid `AttributeError: 'NoneType' object has no attribute 'items'`
# To understand why I have to use `aircraft = avl_object.geometry._base` instead of
# `aircraft = SUAVE.Vehicle()`, compare `print(id(aircraft))` with `print(id(vehicle))`
# in `full_setup()` above - they are not the same object instance!

NOTE: configs.base is a Config, not a Vehicle -> see configs_setup()


To be run as module from Anaconda Prompt via:
(suave) C:\Users\nmb48\Documents\GitHub\SUAVE>python -m nils.jvl.run_suave_jvl_wrapper_nils_new
"""

__all__ = []

import re
import sys
import ast
import numpy as np
import pandas as pd
from copy import deepcopy
from datetime import date

import SUAVE
assert SUAVE.__version__=='2.5.0', 'These tutorials only work with the SUAVE 2.5.0 release'
from SUAVE.Core import Units
from SUAVE.Methods.Geometry.Two_Dimensional.Planform import segment_properties
from SUAVE.Methods.Propulsion import propeller_design
from SUAVE.Components.Energy.Networks.Battery_Propeller import Battery_Propeller
from nils.suave_base_ac_object import full_setup

#%% High-level user inputs

ac_segment = "regional"  # "narrowbody" or "regional"

#%% Command-line parsing

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

#%% Aircraft object

configs, analyses = full_setup()
avl_object = analyses.configs.base.stability

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

# Read geometry data from .csv file
date_str = date.today().strftime("%d%m%y")
if ac_segment == "narrowbody":
    df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\mass\data\mass_distr_geometry_narrowbody_kerosene_090126.csv')
    # df_geom = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_narrowbody_kerosene_{date_str}.csv')
    df_mass = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\mass\data\mass_distr_results_narrowbody_kerosene_090126_6.csv')
    # df_mass = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_narrowbody_kerosene_{date_str}_{study_idx}.csv')
    aircraft.tag = 'Airbus_A220-100'
    t_tail_bool = False
elif ac_segment == "regional":
    df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\mass\data\mass_distr_geometry_regional_kerosene_230126.csv')
    # df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_kerosene_090126.csv')
    # df_geom = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\mass\data\mass_distr_geometry_regional_kerosene_{date_str}.csv')
    df_mass = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\mass\data\mass_distr_results_regional_kerosene_230126_6.csv')
    # df_mass = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_regional_kerosene_090126_6.csv')
    # df_mass = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\mass\data\mass_distr_results_regional_kerosene_{date_str}_{study_idx}.csv')
    aircraft.tag = 'ATR_72-600'
    t_tail_bool = True

# Aircraft mass properties

for counter, row in df_mass.iterrows():
    print('counter =', counter)
    
    configs, analyses = full_setup()
    avl_object = analyses.configs.base.stability
    
    avl_object.settings.filenames.avl_bin_name = r"C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\jvl2.16\jvl.exe"
    
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
    
    ###
    section_indices = sorted(
        {int(re.search(r"sect_(\d+)_", c).group(1))
         for c in df_geom.columns
         if c.startswith("sect_")}
    )
    
    for i in section_indices:
        segment_airfoil = SUAVE.Components.Airfoils.Airfoil()
        segment_airfoil.coordinate_file = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737a.txt'
        segment = SUAVE.Components.Wings.Segment()
        segment.tag = df_geom[f'sect_{i}_label'][0]
        segment.percent_span_location = df_geom[f'sect_{i}_percent_span_location'][0]
        segment.twist = 0.0  # not modelled in TASOPT.jl
        segment.root_chord_percent = df_geom[f'sect_{i}_root_chord_percent'][0]
        segment.thickness_to_chord = df_geom[f'sect_{i}_thickness_to_chord'][0]
        segment.dihedral_outboard = df_geom[f'sect_{i}_dihedral_outboard'][0]
        segment.sweeps.quarter_chord = df_geom[f'sect_{i}_sweeps_quarter_chord'][0]
        segment.append_airfoil(segment_airfoil)
        wing.append_segment(segment)
    
    aircraft.append_component(wing)
    
    # Horizontal stabiliser
    
    wing = SUAVE.Components.Wings.Horizontal_Tail()
    wing.tag = 'horizontal_stabilizer'
    
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
    
    segment = SUAVE.Components.Wings.Segment()
    segment.tag = 'root_segment'
    segment.percent_span_location = df_geom['htail_root_percent_span_location'][0]
    segment.twist = df_geom['htail_root_twist'][0]
    segment.root_chord_percent = df_geom['htail_root_root_chord_percent'][0]
    segment.thickness_to_chord = df_geom['htail_root_thickness_to_chord'][0]
    segment.dihedral_outboard = df_geom['htail_root_dihedral_outboard'][0]
    segment.sweeps.quarter_chord = df_geom['htail_root_sweeps_quarter_chord'][0]
    wing.append_segment(segment)
    
    segment = SUAVE.Components.Wings.Segment()
    segment.tag = 'tip_segment'
    segment.percent_span_location = df_geom['htail_tip_percent_span_location'][0]
    segment.twist = df_geom['htail_tip_twist'][0]
    segment.root_chord_percent = df_geom['htail_tip_root_chord_percent'][0]
    segment.thickness_to_chord = df_geom['htail_tip_thickness_to_chord'][0]
    segment.dihedral_outboard = df_geom['htail_tip_dihedral_outboard'][0]
    segment.sweeps.quarter_chord = df_geom['htail_tip_sweeps_quarter_chord'][0]
    wing.append_segment(segment)
    
    """ Commented on 23.01.2026 as do not want to trim aircraft for blown-wing design space exploration
    # control surfaces -------------------------------------------
    elevator                       = SUAVE.Components.Wings.Control_Surfaces.Elevator()
    elevator.tag                   = 'elevator'
    elevator.span_fraction_start   = df_geom['htail_elevator_span_fraction_start'][0]
    elevator.span_fraction_end     = df_geom['htail_elevator_span_fraction_end'][0]
    elevator.deflection            = df_geom['htail_elevator_deflection'][0]
    elevator.chord_fraction        = df_geom['htail_elevator_chord_fraction'][0]
    wing.append_control_surface(elevator)
    """
        
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
    #####
    # Number of nacelles (must be even)
    N_nacelles = len(ast.literal_eval(df_geom['nacelle_origin'].iloc[0])) * 2  # .csv file only contains nacelle positions for one wing half
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
    
    # =============================================================================
    ##### NILS: from X57_Maxwell_Mod2.py
    #---------------------------------------------------------------------------------------------
    # DEFINE PROPELLER
    #---------------------------------------------------------------------------------------------
    # build network
    net = Battery_Propeller()
    net.number_of_propeller_engines  = N_nacelles
    net.identical_propellers         = True
    
    for nacelle in nacelles:
    
        # Component 2 the Propeller 
        prop = SUAVE.Components.Energy.Converters.Propeller()
        prop.tag = 'propeller_1'
        prop.tip_radius = df_geom['Dprop'][0] / 2
        prop.Jgain = df_geom['Jgain'][0]
        net.propellers.append(prop)
    
    # add the network to the vehicle
    aircraft.append_component(net)
    #####
    
    # IMPORTANT inputs to AVL class
    tag = 'jvl'
    settings_trim_aircraft = False
    training_angle_of_attack = np.array([3.0]) * Units.degrees
    training_Mach = np.array([0.15])
    training_side_slip_angle = np.zeros_like(training_Mach) * Units.degrees
    training_altitude = np.array([0.0])
    training_load_factor = np.array([1.0])
    training_mass = np.array([aircraft.mass_properties.takeoff])
    backend = 'JVL'
    run_modal = False
    settings_number_spanwise_vortices = 30
    
    # avl_object.sample_training(
    #     study_idx=study_idx, counter=counter,
    #     sigma_fcs = sigma_fcs,
    #     span_loc = span_loc,
    #     fcs_loc = fcs_loc,
    #     wing_frac = wing_frac,
    #     nacelle_frac = nacelle_frac,
    #     # NILS: inputs added on 17.03.2026 to distinguish AVL from JVL calls
    #     N_eng=N_nacelles,
    #     Prop_PR_des=1.005,
    #     tag=tag,
    #     settings_trim_aircraft=settings_trim_aircraft,
    #     training_angle_of_attack=training_angle_of_attack,
    #     training_Mach=training_Mach,
    #     training_side_slip_angle=training_side_slip_angle,
    #     training_altitude=training_altitude,
    #     training_load_factor=training_load_factor,
    #     training_mass=training_mass,
    #     backend=backend,
    #     run_modal=run_modal,
    #     settings_number_spanwise_vortices=settings_number_spanwise_vortices,
    # )
    
    # =============================================================================
    # 26.03.2026
    from SUAVE.Methods.Aerodynamics.AVL.read_results import read_results
    results = read_results(avl_object, backend='JVL')
    # =============================================================================
        
    
    