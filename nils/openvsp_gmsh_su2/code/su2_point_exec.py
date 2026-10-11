r"""
This script is based off \SUAVE\Tutorials-2.5.0.2\BWB_CFD\bwb\BWB.py
It is used both for running CFD simulations using SU2, as well as
to generate 3D surface meshes of the external airframe geometry used
for visualising packaging and stability studies in Chapter 7 of my thesis.
In the latter case, the process aborts after the surface mesh generation step.

Note that this script MUST be run from the Anaconda Prompt as

(suave-vsp3190) C:\Users\nmb48\Documents\GitHub\SUAVE>python "C:\Users\nmb48\Documents\GitHub\SUAVE\nils\openvsp_gmsh_su2\code\su2_point_exec.py"

since Python 3.6.15 required for the OpenVSP 3.19 is too old for any
spyder-kernels compatible with the Spyder 6 GUI.
"""

__all__ = []

# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# NILS: add OpenVSP python binding to path (add to system path long-term)
import sys
# sys.path.insert(0, r"C:\Users\nmb48\Documents\GitHub\SUAVE\nils\openvsp_gmsh_su2\OpenVSP-3.46.0-win64-Python3.9\OpenVSP-3.46.0-win64\python\openvsp")  # NILS: use with OpenVSP 3.46.0
sys.path.insert(0, r'C:\Users\nmb48\Documents\GitHub\SUAVE\nils\openvsp_gmsh_su2\OpenVSP-3.19.0-win64\python\openvsp')  # NILS: use with OpenVSP 3.19.0

import SUAVE
assert SUAVE.__version__=='2.5.0', 'These tutorials only work with the SUAVE 2.5.0 release'

import numpy as np
import pylab as plt

from SUAVE.Core import Data, Units
from SUAVE.Input_Output.OpenVSP import write
from SUAVE.Input_Output.OpenVSP import get_vsp_measurements
from SUAVE.Methods.Geometry.Two_Dimensional.Planform import segment_properties
from SUAVE.Plots.Performance.Mission_Plots import *

# NILS
import ast
import pandas as pd
from copy import deepcopy
from datetime import date
from SUAVE.Methods.Flight_Dynamics.Static_Stability.Approximations.Supporting_Functions.convert_sweep import convert_sweep

AC_SEGMENT = "narrowbody"  # NILS: "narrowbody" or "regional"
HALF_MESH_FLAG = False  # NILS: set `False` for visualisation, set `True` to use symmetry to halve domain size
CFD_FAR_FIELD_FLAG = 0  # NILS: set `1` for CFD simulation, set `0` for visualisation to generate surface .stl mesh only
MAXIMUM_ITERATIONS = 3000  # NILS: increase from original value of 10 to 1500, then to 3000
INCLUDE_NACELLES = False  # NILS: added on 21.03.2026 to exclude simple OpenVSP nacelles from geometry in favour of my own 3D CST models

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main(vehicle=None):  # NILS: added default argument

    configs, analyses = full_setup(vehicle)  # NILS: added argument

    simple_sizing(configs)

    configs.finalize()
    analyses.finalize()
    
    return

# ----------------------------------------------------------------------
#   Analysis Setup
# ----------------------------------------------------------------------

def full_setup(vehicle=None):  # NILS: from run_suave_avl_wrapper_nils.py

    # vehicle data
    vehicle  = vehicle_setup()
    configs  = configs_setup(vehicle)

    # vehicle analyses
    configs_analyses = analyses_setup(configs, vehicle)  # NILS: added second argument
    
    analyses = SUAVE.Analyses.Analysis.Container()
    analyses.configs  = configs_analyses

    return configs, analyses

# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ----------------------------------------------------------------------

def analyses_setup(configs, vehicle):  # NILS: added second argument

    analyses = SUAVE.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config, vehicle)  # NILS: added second argument
        analyses[tag] = analysis
    
    return analyses

def base_analysis(config, vehicle):  # NILS: added second argument
    
    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = SUAVE.Analyses.Vehicle()

    # ------------------------------------------------------------------
    #  Basic Geometry Relations
    sizing = SUAVE.Analyses.Sizing.Sizing()
    sizing.features.vehicle = config
    analyses.append(sizing)
    
    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = SUAVE.Analyses.Aerodynamics.SU2_Euler()
    aerodynamics.geometry = config  # NILS: replaced by below line
    # aerodynamics.geometry = vehicle  # NILS: use vehicle instead of config
    # print(id(aerodynamics.geometry))  # NILS: compare with `print(id(self.geometry))` in SU2_Euler.py (should match)
    # print('BEFORE: aerodynamics.geometry.tag =', aerodynamics.geometry.tag)
    # aerodynamics.geometry.tag = vehicle.tag  # NILS: set geometry tag to vehicle tag, NOT config tag
    # This would be sensible for aerodynamics.geometry = config above, but somehow `tag = self.geometry.tag`
    # in SU2_Euler.py still returns 'base' rather than e.g. Airbus_A220-100.
    # print('AFTER: aerodynamics.geometry.tag =', aerodynamics.geometry.tag)
    
    # NILS: require MS MPI for parallel execution of SU2 (implement in future if have time)
    # aerodynamics.process.compute.lift.inviscid.settings.parallel          = True
    # aerodynamics.process.compute.lift.inviscid.settings.processors        = 12
    # aerodynamics.process.compute.lift.inviscid.training_file              = 'base_data_1500.txt'  # NILS: commented, othwerwise surrogate model based on precomputed CFD results will be used
    
    # NILS: increased to 3000 to obtain smooth lift distribution
    aerodynamics.process.compute.lift.inviscid.settings.maximum_iterations = MAXIMUM_ITERATIONS
    aerodynamics.settings.drag_coefficient_increment = 0.0000  # NILS: I could technically use this to model a nacelle/radiator
    aerodynamics.settings.half_mesh_flag             = HALF_MESH_FLAG
    aerodynamics.settings.cfd_far_field_flag         = CFD_FAR_FIELD_FLAG
    aerodynamics.settings.span_efficiency            = 0.85  # NILS: kept as in BWB.py for now, as not present in tut_mission_B737_AVL.py
    
    if vehicle.tag == 'Airbus_A220-100':
        aerodynamics.process.compute.lift.inviscid.training.Mach               = np.array([.78])  # NILS: previously `np.array([.3, .5, .7, .85])`
        aerodynamics.process.compute.lift.inviscid.training.angle_of_attack    = np.array([3.]) * Units.deg  # NILS: previously `np.array([0.,3.,6.]) * Units.deg`
        # aerodynamics.process.compute.lift.inviscid.training.freestream_pressure = np.array([22699.93683700412])  # NILS: added to support analysis at different altitudes (ISA at 11 km)
        # aerodynamics.process.compute.lift.inviscid.training.freestream_temperature = np.array([216.77351270445553])  # NILS: added to support analysis at different altitudes (ISA at 11 km)
        aerodynamics.process.compute.lift.inviscid.training.freestream_pressure = np.array([19746.224626609113])  # NILS: added to support analysis at different altitudes (ISA at 39000 ft)
        aerodynamics.process.compute.lift.inviscid.training.freestream_temperature = np.array([216.65])  # NILS: added to support analysis at different altitudes (ISA at 39000 ft)
    elif vehicle.tag == 'ATR_72-600':
        aerodynamics.process.compute.lift.inviscid.training.Mach = np.array([.43])
        aerodynamics.process.compute.lift.inviscid.training.angle_of_attack = np.array([3.]) * Units.deg
        aerodynamics.process.compute.lift.inviscid.training.freestream_pressure = np.array([46600.63382092578])
        aerodynamics.process.compute.lift.inviscid.training.freestream_temperature = np.array([248.5639621472224])
    
    analyses.append(aerodynamics)
    
    # done!
    return analyses    
    
# ----------------------------------------------------------------------
#   Define the Vehicle
# ----------------------------------------------------------------------

def vehicle_setup(
    aircraft_is_b737 = False,  # False
    keep_b737_defaults = False,  # True
    wing_has_segments = True,  # False
    htail_has_segments = True,  # False
    ac_segment = AC_SEGMENT,
):
    
    # Read geometry data from .csv file
    if aircraft_is_b737:
        df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\B737_AVL_Tutorial\suave_avl_wrapper_tasopt_inputs_b737.csv')
        t_tail_bool = False
    elif not aircraft_is_b737:
        date_str = date.today().strftime("%d%m%y")
        if ac_segment == "narrowbody":
            # df_geom = pd.read_csv(r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\airbus\data\mass_distr_geometry_narrowbody_kerosene_260326.csv")
            # df_mass = pd.read_csv(r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\airbus\data\mass_distr_results_narrowbody_kerosene_260326.csv")
            df_geom = pd.read_csv(r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_narrowbody_LH2_031026.csv")
            df_mass = pd.read_csv(r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_narrowbody_LH2_031026.csv")
            t_tail_bool = False
        elif ac_segment == "regional":
            # df_geom = pd.read_csv(r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\airbus\data\mass_distr_geometry_regional_kerosene_260326.csv")
            # df_mass = pd.read_csv(r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\airbus\data\mass_distr_results_regional_kerosene_260326.csv")
            # df_geom = pd.read_csv(r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_LH2_021026.csv")
            # df_mass = pd.read_csv(r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_regional_LH2_021026.csv")
            df_geom = pd.read_csv(r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_LH2_031026.csv")
            df_mass = pd.read_csv(r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_regional_LH2_031026.csv")
            t_tail_bool = True
    
    # Define clean-slate aircraft

    aircraft = SUAVE.Vehicle()   # new clean vehicle
    
    if ac_segment == "narrowbody":
        aircraft.tag = 'Airbus_A220-100'
    elif ac_segment == "regional":
        aircraft.tag = 'ATR_72-600'
        
    for counter, row in df_mass.iterrows():
    
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
        
        # NILS: specify wing area so SU2 can correctly scale physical dimensions and forces
        aircraft.reference_area = df_geom['wing_areas_reference'][0]  # NILS: see `f.write(f'REF_AREA = {float(ref_area)}\n\n')` in write_SU2_cfg.py

        # NILS: fuselage geometry definition adapted from b737_su2_nils.py (AVL definition insufficient for SU2)
        # See C:\Users\nmb48\Documents\GitHub\SUAVE\regression\scripts\Vehicles\Concorde.py and
        # C:\Users\nmb48\Documents\GitHub\SUAVE\trunk\SUAVE\Components\Fuselages\Fuselage.py for other
        # fuselage parametrisation options.
        
        fuselage = SUAVE.Components.Fuselages.Fuselage()
        fuselage.tag = 'fuselage'
        
        fuselage.number_coach_seats    = aircraft.passengers
        fuselage.seats_abreast         = 6  # NILS: internal geometry irrelevant for CFD analysis
        fuselage.seat_pitch            = 1     * Units.meter  # NILS: internal geometry irrelevant for CFD analysis
        fuselage.lengths.nose          = df_geom['body_lengths_nose'][0]
        fuselage.lengths.tail          = df_geom['body_lengths_tail'][0]
        #fuselage.fineness.nose         = fuselage.lengths.nose / fuselage.heights.maximum   # NILS: calculated for consistency
        #fuselage.fineness.tail         = fuselage.lengths.tail / fuselage.heights.maximum   # NILS: calculated for consistency
        fuselage.fineness.nose         = 1.6  # NILS: above value too low, resulting in forward bulge
        fuselage.fineness.tail         = 2.  # NILS: above value too low, resulting in rearward bulge
        fuselage.lengths.cabin         = df_geom['body_lengths_total'][0]  # NILS: internal geometry irrelevant for CFD analysis
        fuselage.lengths.total         = df_geom['body_lengths_total'][0]
        fuselage.lengths.fore_space    = 0.0  # NILS: internal geometry irrelevant for CFD analysis
        fuselage.lengths.aft_space     = 0.0  # NILS: internal geometry irrelevant for CFD analysis
        fuselage.width                 = df_geom['body_widths_maximum'][0]
        fuselage.heights.maximum       = df_geom['body_heights_maximum'][0]
        fuselage.effective_diameter    = df_geom['body_heights_maximum'][0]
        fuselage.heights.at_quarter_length          = df_geom['body_heights_maximum'][0]
        fuselage.heights.at_three_quarters_length   = df_geom['body_heights_maximum'][0] * 3.65 / 3.74  # NILS: might have to be reduced in accordance with 
        fuselage.heights.at_wing_root_quarter_chord = df_geom['body_heights_maximum'][0]
        
        ### NILS: fuselage shape customisation inspired by regression\scripts\Vehicles\Concorde.py
        # See vsp_fuselage.py for details.
        fuselage.OpenVSP_values = Data() # VSP uses degrees directly
    
        fuselage.OpenVSP_values.nose = Data()
        fuselage.OpenVSP_values.nose.top = Data()
        fuselage.OpenVSP_values.nose.side = Data()
        fuselage.OpenVSP_values.nose.top.angle = 55.0  # NILS: increase angle (deg) to make nose blunter
        fuselage.OpenVSP_values.nose.top.strength = 0.75  # NILS: increase to enforce angle more
        fuselage.OpenVSP_values.nose.side.angle = 55.0  # NILS: increase angle (deg) to make nose blunter
        fuselage.OpenVSP_values.nose.side.strength = 0.75  # NILS: increase to enforce angle more
        fuselage.OpenVSP_values.nose.TB_Sym = True
        fuselage.OpenVSP_values.nose.z_pos = -.01  # NILS: lower z-position of nose (for tail default is +0.02 of fuselage length, measured from centre line)
        
        fuselage.OpenVSP_values.tail = Data()
        fuselage.OpenVSP_values.tail.top = Data()
        fuselage.OpenVSP_values.tail.side = Data()    
        fuselage.OpenVSP_values.tail.bottom = Data()
        fuselage.OpenVSP_values.tail.top.angle = 0.0
        fuselage.OpenVSP_values.tail.top.strength = 0.0
        ###
        
        aircraft.append_component(fuselage)
            
        # Main wing
        
        wing = SUAVE.Components.Wings.Main_Wing()
        wing.tag = 'main_wing'
        
        wing.aspect_ratio = df_geom['wing_aspect_ratio'][0]
        wing.sweeps.quarter_chord = df_geom['wing_sweeps_quarter_chord'][0]
        # NILS: defined per segment below, but MUST be specified nonetheless, else get zero thickness at symmetry plane!
        wing.thickness_to_chord = df_geom['wing_center_thickness_to_chord'][0]
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
        
        if aircraft_is_b737:  # 4 segments (see ..\..\regression\scripts\Vehicles\Boeing_737.py)

            root_airfoil                          = SUAVE.Components.Airfoils.Airfoil()
            root_airfoil.coordinate_file          = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737a.txt'
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
            yehudi_airfoil.coordinate_file       = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737b.txt'
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
            section2_airfoil.coordinate_file      = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737c.txt'
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
            tip_airfoil.coordinate_file      = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737d.txt'
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
        
            # NILS: convert sweep from 35% of chord (TASOPT convention - CHECK IN MANUAL) to quarter-chord value 
            quarter_chord_sweep = convert_sweep(wing, 0.35, 0.25)
            # print('wing.sweeps.quarter_chord =', np.rad2deg(wing.sweeps.quarter_chord))
            # print('quarter_chord_sweep =', np.rad2deg(quarter_chord_sweep))
            
            center_airfoil                          = SUAVE.Components.Airfoils.Airfoil()
            center_airfoil.coordinate_file          = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737a.txt'
            segment                               = SUAVE.Components.Wings.Segment()
            segment.tag                           = 'Root'
            segment.percent_span_location         = df_geom['wing_center_percent_span_location'][0]
            segment.twist                         = df_geom['wing_center_twist'][0]
            segment.root_chord_percent            = df_geom['wing_center_root_chord_percent'][0]
            segment.thickness_to_chord            = df_geom['wing_center_thickness_to_chord'][0]
            segment.dihedral_outboard             = df_geom['wing_center_dihedral_outboard'][0]
            segment.sweeps.quarter_chord          = df_geom['wing_center_sweeps_quarter_chord'][0]  # NILS: center segment is unswept
            segment.append_airfoil(center_airfoil)
            wing.append_segment(segment)

            inboard_airfoil                       = SUAVE.Components.Airfoils.Airfoil()
            inboard_airfoil.coordinate_file       = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737b.txt'
            segment                               = SUAVE.Components.Wings.Segment()
            segment.tag                           = 'Yehudi'
            segment.percent_span_location         = df_geom['wing_inboard_percent_span_location'][0]
            segment.twist                         = df_geom['wing_inboard_twist'][0]
            segment.root_chord_percent            = df_geom['wing_inboard_root_chord_percent'][0]
            segment.thickness_to_chord            = df_geom['wing_inboard_thickness_to_chord'][0]
            segment.dihedral_outboard             = df_geom['wing_inboard_dihedral_outboard'][0]
            # segment.sweeps.quarter_chord          = df_geom['wing_inboard_sweeps_quarter_chord'][0]
            segment.sweeps.quarter_chord          = quarter_chord_sweep
            segment.append_airfoil(inboard_airfoil)
            wing.append_segment(segment)

            outboard_airfoil                      =  SUAVE.Components.Airfoils.Airfoil()
            outboard_airfoil.coordinate_file      = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737c.txt'
            segment                               = SUAVE.Components.Wings.Segment()
            segment.tag                           = 'Section 2'
            segment.percent_span_location         = df_geom['wing_outboard_percent_span_location'][0]
            segment.twist                         = df_geom['wing_outboard_twist'][0]
            segment.root_chord_percent            = df_geom['wing_outboard_root_chord_percent'][0]
            segment.thickness_to_chord            = df_geom['wing_outboard_thickness_to_chord'][0]
            segment.dihedral_outboard             = df_geom['wing_outboard_dihedral_outboard'][0]
            # segment.sweeps.quarter_chord          = df_geom['wing_outboard_sweeps_quarter_chord'][0]
            segment.sweeps.quarter_chord          = quarter_chord_sweep
            segment.append_airfoil(outboard_airfoil)
            wing.append_segment(segment)
            
            tip_airfoil                      =  SUAVE.Components.Airfoils.Airfoil()
            tip_airfoil.coordinate_file      = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737d.txt'
            segment                               = SUAVE.Components.Wings.Segment()
            segment.tag                           = 'Tip'
            segment.percent_span_location         = df_geom['wing_tip_percent_span_location'][0]
            segment.twist                         = df_geom['wing_tip_twist'][0]
            segment.root_chord_percent            = df_geom['wing_tip_root_chord_percent'][0]
            segment.thickness_to_chord            = df_geom['wing_tip_thickness_to_chord'][0]
            segment.dihedral_tip             = df_geom['wing_tip_dihedral_tip'][0]
            # segment.sweeps.quarter_chord          = df_geom['wing_tip_sweeps_quarter_chord'][0]
            segment.sweeps.quarter_chord          = quarter_chord_sweep
            segment.append_airfoil(tip_airfoil)
            wing.append_segment(segment)
            
        # Fill out more segment properties automatically
        wing = segment_properties(wing)
        
        aircraft.append_component(wing)
        
        # Horizontal stabiliser
        
        wing = SUAVE.Components.Wings.Horizontal_Tail()
        wing.tag = 'horizontal_stabilizer'
        
        wing.aspect_ratio = df_geom['htail_aspect_ratio'][0]
        wing.sweeps.quarter_chord = df_geom['htail_sweeps_quarter_chord'][0]
        # NILS: defined per segment below, but MUST be specified nonetheless, else get zero thickness at symmetry plane!
        wing.thickness_to_chord = df_geom['htail_root_thickness_to_chord'][0]
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
            
        # Fill out more segment properties automatically
        wing = segment_properties(wing)
            
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
            df_geom['vtail_origin_z'][0] + 0.2,
        ]]
        wing.vertical = df_geom['vtail_vertical'][0]
        wing.symmetric = df_geom['vtail_symmetric'][0]
        wing.high_lift = df_geom['vtail_high_lift'][0]
        wing.dihedral = df_geom['vtail_dihedral'][0]
        wing.t_tail = t_tail_bool
        
        # Fill out more segment properties automatically
        wing = segment_properties(wing)
        
        aircraft.append_component(wing)
        
        # Nacelles
        
        # Number of nacelles (must be even)
        N_half = len(ast.literal_eval(df_geom['nacelle_origin'].iloc[0]))
        # assert N_half % 2 == 0, "Number of nacelles must be even"
        
        
        # --- Base nacelle ---
        base_nacelle = SUAVE.Components.Nacelles.Nacelle()
        base_nacelle.tag = 'nacelle_1'
        base_nacelle.length = df_geom['nacelle_length'][0]
        # base_nacelle.inlet_diameter = df_geom['nacelle_inlet_diameter'][0]  # NILS: original line (TASOPT may not model difference between inlet and maximum nacelle diameter)
        base_nacelle.inlet_diameter = df_geom['nacelle_diameter'][0] * 1.90 / 2.05  # NILS: based on b737_su2_nils.py (for OpenVSP/Gmsh to mesh nacelles, require non-zero height;
        # see `height = nacelle.diameter - nacelle.inlet_diameter` in vsp_nacelle.py)
        base_nacelle.diameter = df_geom['nacelle_diameter'][0]
        # base_nacelle.areas.wetted = df_geom['nacelle_areas_wetted'][0]
        base_nacelle.areas.wetted = 1.1*np.pi*base_nacelle.diameter*base_nacelle.length  # NILS
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
            nac.origin = np.array([[
                ast.literal_eval(row)[i]
                for row in df_geom['nacelle_origin']
            ]])  # NILS: added extra dimension relative to run_suave_jvl_wrapper_nils.py
            nacelles.append(nac)
        
        # Mirrored negative-Y side
        for i in range(N_half):
            nac = deepcopy(nacelles[i])
            nac.tag = f'nacelle_{i + 1 + N_half}'
            nac.origin[0][1] *= -1  # NILS: added [0] relative to run_suave_jvl_wrapper_nils.py
            nacelles.append(nac)
        
        # --- Append to aircraft ---
        if INCLUDE_NACELLES == True:  # NILS: added on 21.03.2026 to exclude simple OpenVSP nacelles from geometry in favour of my own 3D CST models
            for nac in nacelles:
                aircraft.append_component(nac)
        
    # ------------------------------------------------------------------
    #   Vehicle Definition Complete
    # ------------------------------------------------------------------

    return aircraft

# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def configs_setup(vehicle):

    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------
    configs = SUAVE.Components.Configs.Config.Container()

    base_config = SUAVE.Components.Configs.Config(vehicle)
    base_config.tag = 'base'  # NILS: leave as 'base' - a vehicle (e.g. Airbus_A220-100) can have multiple configs (here only 1 though)
    configs.append(base_config)
    
    # NILS: use vehicle name as opposed
    # to config name to allow running
    # multiple instances of OpenVSP/Gmsh/SU2
    # in parallel!
    write(vehicle,base_config.tag)  # NILS: original line
    # write(vehicle,vehicle.tag)  # NILS: modified line

    return configs

def simple_sizing(configs):

    base = configs.base
    base.pull_base()

    # zero fuel weight
    base.mass_properties.max_zero_fuel = 0.9 * base.mass_properties.max_takeoff 

    # Areas
    wetted_areas = get_vsp_measurements(base.tag)  # NILS - change tag to vehicle.tag rather than config.tag too?

    for wing in base.wings:
        wing.areas.wetted   = wetted_areas[wing.tag]
        wing.areas.exposed  = wetted_areas[wing.tag]
        wing.areas.affected = 0.6 * wing.areas.wetted

    # diff the new data
    base.store_diff()

    return

#%%

if __name__ == '__main__':
    
    main()
            

