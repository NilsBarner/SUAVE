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


To be run as module from Anaconda Prompt via:
(suave) C:\Users\nmb48\Documents\GitHub\SUAVE>python -m nils.avl.run_suave_avl_wrapper_nils_new
"""

__all__ = []

import sys
import numpy as np
import pandas as pd
import pylab as plt
from copy import deepcopy
from datetime import date

import SUAVE
assert SUAVE.__version__=='2.5.0', 'These tutorials only work with the SUAVE 2.5.0 release'
from SUAVE.Core import Units
from SUAVE.Methods.Geometry.Two_Dimensional.Planform import segment_properties

from nils.suave_base_ac_object import full_setup

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

configs, analyses = full_setup()
avl_object = analyses.configs.base.stability

# NOTE: `._base` required to avoid `AttributeError: 'NoneType' object has no attribute 'items'`
# To understand why I have to use `aircraft = avl_object.geometry._base` instead of
# `aircraft = SUAVE.Vehicle()`, compare `print(id(aircraft))` with `print(id(vehicle))`
# in `full_setup()` above - they are not the same object instance!

wing_has_segments = True  # False
htail_has_segments = True  # False
ac_segment = "narrowbody"  # "narrowbody" or "regional"
fuel = "kerosene"  # "kerosene" or "LH2"
include_nacelles = False  # NILS: added on 21.03.2026 to exclude simple OpenVSP nacelles from geometry in favour of my own 3D CST models

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
date_str = date.today().strftime("%d%m%y")
if ac_segment == "narrowbody":
    if fuel == "kerosene":
        # # df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\mass\data\mass_distr_geometry_narrowbody_kerosene_211225.csv')
        # df_geom = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_narrowbody_kerosene_{date_str}.csv')
        # # df_mass = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\mass\data\mass_distr_results_narrowbody_kerosene_211225_6.csv')
        # df_mass = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_narrowbody_kerosene_{date_str}_{study_idx}.csv')
        # df_geom = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_narrowbody_kerosene_{date_str}.csv')
        # df_mass = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_narrowbody_kerosene_{date_str}_{study_idx}.csv')
        df_geom = pd.read_csv(r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_narrowbody_kerosene_260326.csv")
        df_mass = pd.read_csv(r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_narrowbody_kerosene_260326.csv")
    elif fuel == "LH2":
        # =============================================================================
        # =============================================================================
        # =============================================================================
        # =============================================================================
        # 25.03.2026: something is wrong with the narrowbody LH2 geometry file
        # TEMPORARY FIX: use KEROSENE GEOMETRY but LH2 MASS file!
        # df_geom = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_narrowbody_LH2_{date_str}.csv')
        df_geom = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_narrowbody_kerosene_250326.csv')
        # =============================================================================
        # =============================================================================
        # =============================================================================
        # =============================================================================
        # df_mass = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_narrowbody_LH2_{date_str}_{study_idx}.csv')
        df_mass = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_narrowbody_LH2_250326_{study_idx}.csv')
    aircraft.tag = 'Airbus_A220-100'
    t_tail_bool = False
elif ac_segment == "regional":
    if fuel == "kerosene":
        # # df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\mass\data\mass_distr_geometry_regional_kerosene_211225.csv')
        # df_geom = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_kerosene_{date_str}.csv')
        # # df_mass = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\mass\data\mass_distr_results_regional_kerosene_211225_6.csv')
        # df_mass = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_regional_kerosene_{date_str}_{study_idx}.csv')
        # df_geom = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_kerosene_{date_str}.csv')
        # df_mass = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_regional_kerosene_{date_str}_{study_idx}.csv')
        df_geom = pd.read_csv(r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_kerosene_260326.csv")
        df_mass = pd.read_csv(r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_regional_kerosene_260326.csv")
    elif fuel == "LH2":
        # df_geom = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_LH2_{date_str}.csv')
        # df_mass = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_regional_LH2_{date_str}_{study_idx}.csv')
        df_geom = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_LH2_250326.csv')
        df_mass = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_regional_LH2_250326_{study_idx}.csv')
    aircraft.tag = 'ATR_72-600'
    t_tail_bool = True
    
# sys.exit('Stop here.')
    
# Aircraft mass properties

for counter, row in df_mass.iterrows():
    print('counter =', counter)
    
    configs, analyses = full_setup()
    avl_object = analyses.configs.base.stability
    
    avl_object.settings.filenames.avl_bin_name = r"C:\Users\nmb48\Documents\GitHub\SUAVE\nils\avl\avl3.52\avl.exe"  # NILS: set AVL executable name
    
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
        
        center_airfoil = SUAVE.Components.Airfoils.Airfoil()
        center_airfoil.coordinate_file = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737a.txt'
        segment = SUAVE.Components.Wings.Segment()
        segment.tag = 'Root'
        segment.percent_span_location = df_geom['wing_center_percent_span_location'][0]
        segment.twist = df_geom['wing_center_twist'][0]
        segment.root_chord_percent = df_geom['wing_center_root_chord_percent'][0]
        segment.thickness_to_chord = df_geom['wing_center_thickness_to_chord'][0]
        segment.dihedral_outboard = df_geom['wing_center_dihedral_outboard'][0]
        segment.sweeps.quarter_chord = df_geom['wing_center_sweeps_quarter_chord'][0]
        segment.append_airfoil(center_airfoil)
        wing.append_segment(segment)

        inboard_airfoil = SUAVE.Components.Airfoils.Airfoil()
        inboard_airfoil.coordinate_file = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737b.txt'
        segment = SUAVE.Components.Wings.Segment()
        segment.tag = 'Yehudi'
        segment.percent_span_location = df_geom['wing_inboard_percent_span_location'][0]
        segment.twist = df_geom['wing_inboard_twist'][0]
        segment.root_chord_percent = df_geom['wing_inboard_root_chord_percent'][0]
        segment.thickness_to_chord = df_geom['wing_inboard_thickness_to_chord'][0]
        segment.dihedral_outboard = df_geom['wing_inboard_dihedral_outboard'][0]
        segment.sweeps.quarter_chord = df_geom['wing_inboard_sweeps_quarter_chord'][0]
        segment.append_airfoil(inboard_airfoil)
        wing.append_segment(segment)

        outboard_airfoil = SUAVE.Components.Airfoils.Airfoil()
        outboard_airfoil.coordinate_file = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737c.txt'
        segment = SUAVE.Components.Wings.Segment()
        segment.tag = 'Section 2'
        segment.percent_span_location = df_geom['wing_outboard_percent_span_location'][0]
        segment.twist = df_geom['wing_outboard_twist'][0]
        segment.root_chord_percent = df_geom['wing_outboard_root_chord_percent'][0]
        segment.thickness_to_chord = df_geom['wing_outboard_thickness_to_chord'][0]
        segment.dihedral_outboard = df_geom['wing_outboard_dihedral_outboard'][0]
        segment.sweeps.quarter_chord = df_geom['wing_outboard_sweeps_quarter_chord'][0]
        segment.append_airfoil(outboard_airfoil)
        wing.append_segment(segment)
        
        tip_airfoil = SUAVE.Components.Airfoils.Airfoil()
        tip_airfoil.coordinate_file = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737d.txt'
        segment = SUAVE.Components.Wings.Segment()
        segment.tag = 'Tip'
        segment.percent_span_location = df_geom['wing_tip_percent_span_location'][0]
        segment.twist = df_geom['wing_tip_twist'][0]
        segment.root_chord_percent = df_geom['wing_tip_root_chord_percent'][0]
        segment.thickness_to_chord = df_geom['wing_tip_thickness_to_chord'][0]
        segment.dihedral_tip = df_geom['wing_tip_dihedral_tip'][0]
        segment.sweeps.quarter_chord = df_geom['wing_tip_sweeps_quarter_chord'][0]
        segment.append_airfoil(tip_airfoil)
        wing.append_segment(segment)
            
        # control surfaces -------------------------------------------
        slat = SUAVE.Components.Wings.Control_Surfaces.Slat()
        slat.tag = 'slat'
        slat.span_fraction_start = df_geom['wing_slat_span_fraction_start'][0]
        slat.span_fraction_end = df_geom['wing_slat_span_fraction_end'][0]
        slat.deflection = df_geom['wing_slat_deflection'][0]
        slat.chord_fraction = df_geom['wing_slat_chord_fraction'][0]
        wing.append_control_surface(slat)

        flap = SUAVE.Components.Wings.Control_Surfaces.Flap()
        flap.tag = 'flap'
        flap.span_fraction_start = df_geom['wing_flap_span_fraction_start'][0]
        flap.span_fraction_end = df_geom['wing_flap_span_fraction_end'][0]
        flap.deflection = df_geom['wing_flap_deflection'][0]
        flap.configuration_type = df_geom['wing_flap_configuration_type'][0]
        flap.chord_fraction = df_geom['wing_flap_chord_fraction'][0]
        wing.append_control_surface(flap)

        aileron = SUAVE.Components.Wings.Control_Surfaces.Aileron()
        aileron.tag = 'aileron'
        aileron.span_fraction_start = df_geom['wing_aileron_span_fraction_start'][0]
        aileron.span_fraction_end = df_geom['wing_aileron_span_fraction_end'][0]
        aileron.deflection = df_geom['wing_aileron_deflection'][0]
        aileron.chord_fraction = df_geom['wing_aileron_chord_fraction'][0]
        wing.append_control_surface(aileron)
        
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
        
        # control surfaces -------------------------------------------
        elevator = SUAVE.Components.Wings.Control_Surfaces.Elevator()
        elevator.tag = 'elevator'
        elevator.span_fraction_start = df_geom['htail_elevator_span_fraction_start'][0]
        elevator.span_fraction_end = df_geom['htail_elevator_span_fraction_end'][0]
        elevator.deflection = df_geom['htail_elevator_deflection'][0]
        elevator.chord_fraction = df_geom['htail_elevator_chord_fraction'][0]
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
    nacelle = SUAVE.Components.Nacelles.Nacelle()
    # nacelle.tag = 'nacelle_1'
    nacelle.tag = 'nacelle'
    nacelle.length = df_geom['nacelle_length'][0]
    nacelle.inlet_diameter = df_geom['nacelle_inlet_diameter'][0]
    nacelle.diameter = df_geom['nacelle_diameter'][0]
    nacelle.areas.wetted = df_geom['nacelle_areas_wetted'][0]
    # nacelle.origin = df_geom['nacelle_origin'].to_numpy()
    # =============================================================================
    import ast
    nacelle.origin = df_geom['nacelle_origin'].apply(lambda x: ast.literal_eval(x)[0]).to_numpy()
    # =============================================================================
    nacelle.flow_through = df_geom['nacelle_flow_through'][0]
    nacelle_airfoil = SUAVE.Components.Airfoils.Airfoil() 
    nacelle_airfoil.naca_4_series_airfoil = '2410'
    nacelle.append_airfoil(nacelle_airfoil)
    
    # n_segments = 18
    
    # # Wing Segments
    # for i_segs in range(n_segments):
    #     root_airfoil                          = SUAVE.Components.Airfoils.Airfoil()
    #     root_airfoil.coordinate_file          = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737a.txt'
    #     segment                               = SUAVE.Components.Wings.Segment()
    #     # segment                               = SUAVE.Components.Nacelles.Segment()
    #     segment.tag                           = 'Root'
    #     segment.percent_span_location         = 0.0
    #     segment.twist                         = 0.0 * Units.deg
    #     segment.root_chord_percent            = 1.0
    #     segment.thickness_to_chord            = 0.1
    #     segment.dihedral_outboard             = 0.0 * Units.degrees
    #     segment.sweeps.quarter_chord          = 0.0 * Units.degrees
    #     segment.thickness_to_chord            = 0.1
    #     segment.append_airfoil(root_airfoil)
    #     nacelle.append_segment(segment)
    
    # # # Fill out more segment properties automatically
    # # wing = segment_properties(wing)  
    
    # nacelle_2 = deepcopy(nacelle)
    # nacelle_2.tag = 'nacelle_2'
    # nacelle_2_origin = deepcopy(nacelle.origin)
    # nacelle_2_origin[1] *= -1
    # nacelle_2.origin = nacelle_2_origin
    
    if include_nacelles:  # NILS: added on 26.03.2026 to exclude simple OpenVSP nacelles from AVL analysis for comparison with SU2
        aircraft.append_component(nacelle)  
        # aircraft.append_component(nacelle_2)
    
    # print(aircraft.nacelles.nacelle_1.Airfoil)
    
    # IMPORTANT inputs to AVL class
    tag = 'avl'
    settings_trim_aircraft = True
    backend = 'AVL'
    run_modal = True
    settings_number_spanwise_vortices = 30
    
    # NILS: longitudinal test cases from Table 6.2 in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris
    # NOTE: training inputs always have to be at least 1D, otherwise get
    # `TypeError: object of type 'float' has no len()` in
    # Documents\GitHub\SUAVE\trunk\SUAVE\Methods\Aerodynamics\AVL\translate_data.py
    # training_load_factor = np.array([1.0, 2.5, 1.0, 1.0, 2.5, 1.0])
    training_load_factor = np.array([1.0])
    # training_altitude = np.array([0.0, 10e3, 35e3, 35e3, 10e3, 0.0]) * 0.3048
    # training_Mach = np.array([0.2, 0.5, 0.7, 0.7, 0.5, 0.2])
    if ac_segment == 'regional':
        # training_altitude = np.array([0.0, 10e3 * 20/35, 20e3, 20e3, 10e3 * 20/35, 0.0]) * 0.3048
        # training_Mach = np.array([0.2 * 0.43/0.7, 0.5 * 0.43/0.7, 0.43, 0.43, 0.5 * 0.43/0.7, 0.2 * 0.43/0.7])
        # =============================================================================
        training_altitude = np.array([20e3]) * 0.3048
        training_Mach = np.array([0.43])
        # =============================================================================
    elif ac_segment == 'narrowbody':
        # training_altitude = np.array([0.0, 10e3 * 39/35, 39e3, 39e3, 10e3 * 39/35, 0.0]) * 0.3048
        # training_Mach = np.array([0.2 * 0.78/0.7, 0.5 * 0.78/0.7, 0.78, 0.78, 0.5 * 0.78/0.7, 0.2 * 0.78/0.7])
        # =============================================================================
        training_altitude = np.array([39e3]) * 0.3048
        training_Mach = np.array([0.78])
        # =============================================================================
    training_side_slip_angle = np.zeros_like(training_Mach) * Units.degrees
    # NOTE: 6x faster if use `np.array([0])` instead of `np.zeros_like(self.training.Mach)` (6x duplication)
    training_angle_of_attack = np.array([0])  # to be trimmed
    # training_mass = np.array([
    #     aircraft.mass_properties.takeoff,
    #     aircraft.mass_properties.takeoff,
    #     aircraft.mass_properties.takeoff,
    #     aircraft.mass_properties.max_zero_fuel,
    #     aircraft.mass_properties.max_zero_fuel,
    #     aircraft.mass_properties.max_zero_fuel,
    # ])
    # =============================================================================
    training_mass = np.array([aircraft.mass_properties.max_zero_fuel])
    # =============================================================================
    
    # Run sample_training() only
    avl_object.sample_training(
        study_idx=study_idx, counter=counter,
        sigma_fcs = sigma_fcs,
        span_loc = span_loc,
        fcs_loc = fcs_loc,
        wing_frac = wing_frac,
        nacelle_frac = nacelle_frac,
        # NILS: inputs added on 17.03.2026 to distinguish AVL from JVL calls
        tag=tag,
        settings_trim_aircraft=settings_trim_aircraft,
        training_angle_of_attack=training_angle_of_attack,
        training_Mach=training_Mach,
        training_side_slip_angle=training_side_slip_angle,
        training_altitude=training_altitude,
        training_load_factor=training_load_factor,
        training_mass=training_mass,
        backend=backend,
        run_modal=run_modal,
        settings_number_spanwise_vortices=settings_number_spanwise_vortices,
        aircraft_tag=aircraft.tag,
    )
    

    