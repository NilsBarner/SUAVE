# BWB.py
# 
# Created:  Jan 2017, E. Botero
# Modified: Mar 2018, T. MacDonald

# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# NILS: add OpenVSP python binding to path (add to system path long-term)
import sys
# sys.path.insert(0, r"C:\Users\nmb48\Documents\GitHub\SUAVE\OpenVSP-3.46.0-win64-Python3.9\OpenVSP-3.46.0-win64\python\openvsp")  # NILS: use with OpenVSP 3.46.0
sys.path.insert(0, r'C:\Users\nmb48\Documents\GitHub\SUAVE\OpenVSP-3.19.0-win64\python\openvsp')  # NILS: use with OpenVSP 3.19.0

import SUAVE
assert SUAVE.__version__=='2.5.0', 'These tutorials only work with the SUAVE 2.5.0 release'

import numpy as np
import pylab as plt

from SUAVE.Core import Data, Units

from SUAVE.Input_Output.OpenVSP import write
from SUAVE.Input_Output.OpenVSP import get_vsp_measurements

from SUAVE.Methods.Propulsion.turbofan_sizing import turbofan_sizing
from SUAVE.Methods.Geometry.Two_Dimensional.Cross_Section.Propulsion import compute_turbofan_geometry
from SUAVE.Methods.Geometry.Two_Dimensional.Planform import segment_properties

from SUAVE.Plots.Performance.Mission_Plots import *

from copy import deepcopy

# NILS
from datetime import date
import pandas as pd
from SUAVE.Methods.Propulsion import propeller_design
from SUAVE.Components.Energy.Networks.Battery_Propeller import Battery_Propeller

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main(vehicle=None):  # NILS: added default argument

    configs, analyses = full_setup(vehicle)  # NILS: added argument

    simple_sizing(configs)
    # sys.exit('Stop after simple_sizing.')

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
    aerodynamics.process.compute.lift.inviscid.settings.maximum_iterations = 1500  # NILS: increase from original value of 10
    
    aerodynamics.settings.drag_coefficient_increment = 0.0000
    # aerodynamics.settings.half_mesh_flag             = False  # NILS: original setting
    aerodynamics.settings.half_mesh_flag             = True  # NILS: use symmetry to halve domain size
    aerodynamics.settings.span_efficiency            = 0.85  # NILS: kept as in BWB.py for now, as not present in tut_mission_B737_AVL.py
    
    aerodynamics.process.compute.lift.inviscid.training.Mach               = np.array([.4])  # NILS: previously `np.array([.3, .5, .7, .85])`
    aerodynamics.process.compute.lift.inviscid.training.angle_of_attack    = np.array([3.]) * Units.deg  # NILS: previously `np.array([0.,3.,6.]) * Units.deg`
    aerodynamics.process.compute.lift.inviscid.training.freestream_pressure = np.array([54048.26223756018])  # NILS: added to support analysis at different altitudes (ISA at 11 km)
    aerodynamics.process.compute.lift.inviscid.training.freestream_temperature = np.array([255.67554322180348])  # NILS: added to support analysis at different altitudes (ISA at 11 km)
    
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
    ac_segment = "regional",  # "narrowbody" or "regional"
):
    
    # Read geometry data from .csv file
    if aircraft_is_b737:
        df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\B737_AVL_Tutorial\suave_avl_wrapper_tasopt_inputs_b737.csv')
        t_tail_bool = False
    elif not aircraft_is_b737:
        date_str = date.today().strftime("%d%m%y")
        if ac_segment == "narrowbody":
            # df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_narrowbody_kerosene_211225.csv')
            # df_mass = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_narrowbody_kerosene_211225_6.csv')
            df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_narrowbody_kerosene_090126.csv')
            df_mass = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_narrowbody_kerosene_211225_6.csv')
            t_tail_bool = False
        elif ac_segment == "regional":
            # df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_kerosene_211225.csv')
            # df_mass = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_regional_kerosene_211225_6.csv')
            df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_kerosene_090126.csv')
            df_mass = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_regional_kerosene_090126_6.csv')
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
            
            center_airfoil                          = SUAVE.Components.Airfoils.Airfoil()
            center_airfoil.coordinate_file          = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737a.txt'
            segment                               = SUAVE.Components.Wings.Segment()
            segment.tag                           = 'Root'
            segment.percent_span_location         = df_geom['wing_center_percent_span_location'][0]
            segment.twist                         = df_geom['wing_center_twist'][0]
            segment.root_chord_percent            = df_geom['wing_center_root_chord_percent'][0]
            segment.thickness_to_chord            = df_geom['wing_center_thickness_to_chord'][0]
            segment.dihedral_outboard             = df_geom['wing_center_dihedral_outboard'][0]
            segment.sweeps.quarter_chord          = df_geom['wing_center_sweeps_quarter_chord'][0]
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
            segment.sweeps.quarter_chord          = df_geom['wing_inboard_sweeps_quarter_chord'][0]
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
            segment.sweeps.quarter_chord          = df_geom['wing_outboard_sweeps_quarter_chord'][0]
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
            segment.sweeps.quarter_chord          = df_geom['wing_tip_sweeps_quarter_chord'][0]
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
        '''
        nacelle = SUAVE.Components.Nacelles.Nacelle()
        nacelle.tag = 'nacelle_1'
        nacelle.length = df_geom['nacelle_length'][0]
        # nacelle.inlet_diameter = df_geom['nacelle_inlet_diameter'][0]  # NILS: original line (TASOPT may not model difference between inlet and maximum nacelle diameter)
        nacelle.inlet_diameter = df_geom['nacelle_diameter'][0] * 1.90 / 2.05  # NILS: based on b737_su2_nils.py (for OpenVSP/Gmsh to mesh nacelles, require non-zero height;
        # see `height = nacelle.diameter - nacelle.inlet_diameter` in vsp_nacelle.py)
        nacelle.diameter = df_geom['nacelle_diameter'][0]
        # nacelle.areas.wetted = df_geom['nacelle_areas_wetted'][0]
        nacelle.areas.wetted = 1.1*np.pi*nacelle.diameter*nacelle.length  # NILS
        nacelle.origin = [list(array) for array in np.expand_dims(df_geom['nacelle_origin'].to_numpy(), axis=0)]  # NILS: added np.expand_dims() to avoid "IndexError: invalid index to scalar variable." in trunk\SUAVE\Input_Output\OpenVSP\vsp_nacelle.py
        print('nacelle.origin =', nacelle.origin)
        # NILS: could it be that this should be a list() rather than a numpy array?
        nacelle.flow_through = df_geom['nacelle_flow_through'][0]
        nacelle_airfoil = SUAVE.Components.Airfoils.Airfoil() 
        nacelle_airfoil.naca_4_series_airfoil = '2410'
        nacelle.append_airfoil(nacelle_airfoil)
        
        nacelle_2 = deepcopy(nacelle)
        nacelle_2.tag = 'nacelle_2'
        nacelle_2_origin = deepcopy(nacelle.origin)
        nacelle_2_origin[0][1] *= -1  # NILS: added [0] relative to run_suave_avl_wrapper_nils.py
        nacelle_2.origin = nacelle_2_origin
        # NILS: could it be that this should be a list() rather than a numpy array?
        
        aircraft.append_component(nacelle)  
        aircraft.append_component(nacelle_2)
        '''
        ###
        import ast
        
        # Number of nacelles (must be even)
        N_half = len(ast.literal_eval(df_geom['nacelle_origin'].iloc[0]))
        print('N_half =', N_half)
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
            print('Hello world')
        
        # Mirrored negative-Y side
        for i in range(N_half):
            nac = deepcopy(nacelles[i])
            nac.tag = f'nacelle_{i + 1 + N_half}'
            nac.origin[0][1] *= -1  # NILS: added [0] relative to run_suave_jvl_wrapper_nils.py
            nacelles.append(nac)
            print('Hello world 2')
        
        # --- Append to aircraft ---
        for nac in nacelles:
            aircraft.append_component(nac)
        ###
        r'''
        # NILS: \regression\scripts files containing 'propeller'
        # C:\Users\nmb48\Documents\GitHub\SUAVE\regression\scripts\Vehicles\Cessna_172.py
        # C:\Users\nmb48\Documents\GitHub\SUAVE\regression\scripts\Vehicles\Electric_Multicopter.py
        # C:\Users\nmb48\Documents\GitHub\SUAVE\regression\scripts\Vehicles\Solar_UAV.py
        # C:\Users\nmb48\Documents\GitHub\SUAVE\regression\scripts\Vehicles\Stopped_Rotor.py
        # C:\Users\nmb48\Documents\GitHub\SUAVE\regression\scripts\Vehicles\Tiltwing.py
        # C:\Users\nmb48\Documents\GitHub\SUAVE\regression\scripts\Vehicles\X57_Maxwell_Mod2.py
        # C:\Users\nmb48\Documents\GitHub\SUAVE\regression\scripts\Vehicles\Propellers\APC_10x7_thin_electric.py
        
        ##### NILS: from X57_Maxwell_Mod2.py
        #---------------------------------------------------------------------------------------------
        # DEFINE PROPELLER
        #---------------------------------------------------------------------------------------------
        # build network
        net = Battery_Propeller()
        net.number_of_propeller_engines  = 2. 
        net.identical_propellers         = True
        
        # Component 2 the Propeller 
        prop = SUAVE.Components.Energy.Converters.Propeller()
        prop.tag = 'propeller_1'
        prop.number_of_blades       = 2.0
        prop.freestream_velocity    = 135.*Units['mph']
        prop.angular_velocity       = 1300.  * Units.rpm
        prop.tip_radius             = 76./2. * Units.inches
        prop.hub_radius             = 8.     * Units.inches
        prop.design_Cl              = 0.8
        prop.design_altitude        = 12000. * Units.feet
        prop.design_altitude        = 12000. * Units.feet
        prop.design_thrust          = 1200.
        prop.origin                 = nacelle.origin  # [[2.,2.5,0.784]]
        prop.rotation               = -1
        prop.symmetry               = True
        prop.variable_pitch         = True 
        prop.airfoil_geometry       =  [r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/NACA_4412.txt']
        prop.airfoil_polars         = [[r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/Polars/NACA_4412_polar_Re_50000.txt' ,
                                        r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/Polars/NACA_4412_polar_Re_100000.txt' ,
                                        r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/Polars/NACA_4412_polar_Re_200000.txt' ,
                                        r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/Polars/NACA_4412_polar_Re_500000.txt' ,
                                        r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/Polars/NACA_4412_polar_Re_1000000.txt' ]]

        prop.airfoil_polar_stations = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        prop                        = propeller_design(prop)

        prop_left = deepcopy(prop)
        prop_left.tag = 'propeller_2' 
        prop_left.origin   = nacelle_2.origin  # [[2.,-2.5,0.784]]
        prop_left.rotation = 1
        
        net.propellers.append(prop)
        net.propellers.append(prop_left)
        """#####
        
        ####### NILS: from Cessna_172.py
        # ------------------------------------------------------------------
        #   Piston Propeller Network
        # ------------------------------------------------------------------    
        
        # build network
        net                                         = SUAVE.Components.Energy.Networks.Internal_Combustion_Propeller()
        net.tag                                     = 'internal_combustion'
        net.number_of_engines                       = 1.
        net.identical_propellers                    = True
                                                    
        # the engine                    
        engine                                  = SUAVE.Components.Energy.Converters.Internal_Combustion_Engine()
        engine.sea_level_power                  = 180. * Units.horsepower
        engine.flat_rate_altitude               = 0.0
        engine.rated_speed                      = 2700. * Units.rpm
        engine.power_specific_fuel_consumption  = 0.52 
        net.engines.append(engine)
        
        # the prop
        prop = SUAVE.Components.Energy.Converters.Propeller()
        prop.number_of_blades        = 2.0
        prop.freestream_velocity     = 119.   * Units.knots
        prop.angular_velocity        = 2650.  * Units.rpm
        prop.tip_radius              = 5 * 76./2. * Units.inches
        prop.hub_radius              = 8.     * Units.inches
        prop.design_Cl               = 0.8
        prop.design_altitude         = 12000. * Units.feet
        prop.design_power            = .64 * 180. * Units.horsepower
        prop.variable_pitch          = True

        prop.airfoil_geometry        =  [r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/NACA_4412.txt'] 
        prop.airfoil_polars          = [[r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/Polars/NACA_4412_polar_Re_50000.txt' ,
                                         r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/Polars/NACA_4412_polar_Re_100000.txt' ,
                                         r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/Polars/NACA_4412_polar_Re_200000.txt' ,
                                         r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/Polars/NACA_4412_polar_Re_500000.txt' ,
                                         r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/Polars/NACA_4412_polar_Re_1000000.txt']]

        prop.airfoil_polar_stations  = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]      
        prop                         = propeller_design(prop)   
        
        ## NILS
        #prop.rotation = -1
        #prop.origin = nacelle.origin
        ## prop.origin[0][1] *= -1  # NILS: added [0] relative to run_suave_avl_wrapper_nils.py
        
        #prop_left = deepcopy(prop)
        #prop_left.tag = 'propeller_2' 
        #prop_left.origin[0][1] *= -1  # NILS: added [0] relative to run_suave_avl_wrapper_nils.py
        #prop_left.rotation = 1
        
        net.propellers.append(prop)
        #net.propellers.append(prop_left)  # NILS
        """#######
        
        # add the network to the vehicle
        aircraft.append_component(net)
        '''
        ###
        
    # ------------------------------------------------------------------
    #   Vehicle Definition Complete
    # ------------------------------------------------------------------

    return aircraft

r"""
# NILS: from b737_su2_nils.py
def vehicle_setup():
    
    # NILS: vehicle description mirroring tut_mission_B737_AVL.py,
    # with modifications allowing CFD analysis
    
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------    
    vehicle = SUAVE.Vehicle()
    vehicle.tag = 'Boeing_737-800'    

    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------    
    # mass properties
    # vehicle.mass_properties.max_takeoff               = 79015.8 * Units.kilogram 
    vehicle.mass_properties.max_takeoff               = 0.0  # NILS
    vehicle.mass_properties.takeoff                   = 79015.8 * Units.kilogram   
    vehicle.mass_properties.operating_empty           = 62746.4 * Units.kilogram 
    vehicle.mass_properties.takeoff                   = 79015.8 * Units.kilogram 
    vehicle.mass_properties.max_zero_fuel             = 62732.0 * Units.kilogram 
    vehicle.mass_properties.cargo                     = 10000.  * Units.kilogram   
    # =============================================================================
    vehicle.mass_properties.mass = 79015.8 * Units.kilogram  # NILS
    # =============================================================================
    # NILS: copied from regression\scripts\Vehicles\Boeing_737.py (inertias required for dynamic stability analysis)
    vehicle.mass_properties.center_of_gravity         = [[ 15.30987849,   0.        ,  -0.48023939]]
    vehicle.mass_properties.moments_of_inertia.tensor = [[3173074.17, 0 , 28752.77565],[0 , 3019041.443, 0],[0, 0, 5730017.433]] # estimated, not correct
    
    # envelope properties
    vehicle.envelope.ultimate_load = 2.5
    vehicle.envelope.limit_load    = 1.5

    # basic parameters
    vehicle.reference_area         = 124.862 * Units['meters**2']  
    vehicle.passengers             = 170
    vehicle.systems.control        = "fully powered" 
    vehicle.systems.accessories    = "medium range"

    # ------------------------------------------------------------------        
    #   Main Wing
    # ------------------------------------------------------------------        
    # NILS: regression\scripts\Vehicles\Boeing_737.py defines a multi-segment wing (as I would get from TASOPT.jl)
    
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
    root_airfoil.coordinate_file          = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737a.txt'
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
    yehudi_airfoil.coordinate_file        = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737b.txt'
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
    mid_airfoil.coordinate_file           = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737c.txt'
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
    tip_airfoil.coordinate_file           = r'C:/Users/nmb48/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737d.txt'
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

    # add to vehicle
    vehicle.append_component(wing)
    
    # NILS: for the `trim_aircraft = True` in trunk\SUAVE\Analyses\Stability\AVL.py
    # to work, ALSO THESE MUST BE DEFINED AS MULTI-SEGMENT WINGS, otherwise
    # 'elevator' won't feature in `wing.control_surfaces` in
    # trunk\SUAVE\Methods\Flight_Dynamics\Dynamic_Stability\compute_dynamic_flight_modes.py
    
    # ------------------------------------------------------------------
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------

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
    
    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------

    wing = SUAVE.Components.Wings.Vertical_Tail()
    wing.tag = 'vertical_stabilizer'

    wing.aspect_ratio            = 1.98865
    wing.sweeps.quarter_chord    = 31.2  * Units.deg   
    wing.thickness_to_chord      = 0.08
    wing.taper                   = 0.1183

    wing.spans.projected         = 8.33
    wing.total_length            = wing.spans.projected 
    
    wing.chords.root             = 10.1 
    wing.chords.tip              = 1.20 
    wing.chords.mean_aerodynamic = 4.0

    wing.areas.reference         = 34.89
    wing.areas.wetted            = 57.25 
    
    wing.twists.root             = 0.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees

    wing.origin                  = [[26.944,0,1.54]]
    wing.aerodynamic_center      = [0,0,0]

    wing.vertical                = True
    wing.symmetric               = False
    wing.t_tail                  = False

    wing.dynamic_pressure_ratio  = 1.0


    # Wing Segments
    segment                               = SUAVE.Components.Wings.Segment()
    segment.tag                           = 'root'
    segment.percent_span_location         = 0.0
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 1.
    segment.dihedral_outboard             = 0 * Units.degrees
    segment.sweeps.quarter_chord          = 61.485 * Units.degrees  
    segment.thickness_to_chord            = .1
    wing.append_segment(segment)

    segment                               = SUAVE.Components.Wings.Segment()
    segment.tag                           = 'segment_1'
    segment.percent_span_location         = 0.2962
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 0.45
    segment.dihedral_outboard             = 0. * Units.degrees
    segment.sweeps.quarter_chord          = 31.2 * Units.degrees   
    segment.thickness_to_chord            = .1
    wing.append_segment(segment)

    segment                               = SUAVE.Components.Wings.Segment()
    segment.tag                           = 'segment_2'
    segment.percent_span_location         = 1.0
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 0.1183 
    segment.dihedral_outboard             = 0.0 * Units.degrees
    segment.sweeps.quarter_chord          = 0.0    
    segment.thickness_to_chord            = .1  
    wing.append_segment(segment)
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing)        

    # add to vehicle
    vehicle.append_component(wing)
    
    # ------------------------------------------------------------------
    #  Fuselage
    # ------------------------------------------------------------------
    
    fuselage = SUAVE.Components.Fuselages.Fuselage()
    fuselage.tag = 'fuselage'
    
    fuselage.number_coach_seats    = vehicle.passengers
    fuselage.seats_abreast         = 6
    fuselage.seat_pitch            = 1     * Units.meter
    fuselage.fineness.nose         = 1.6
    fuselage.fineness.tail         = 2.
    fuselage.lengths.nose          = 6.4   * Units.meter
    fuselage.lengths.tail          = 8.0   * Units.meter
    fuselage.lengths.cabin         = 28.85 * Units.meter
    fuselage.lengths.total         = 38.02 * Units.meter
    fuselage.lengths.fore_space    = 6.    * Units.meter
    fuselage.lengths.aft_space     = 5.    * Units.meter
    fuselage.width                 = 3.74  * Units.meter
    fuselage.heights.maximum       = 3.74  * Units.meter
    fuselage.effective_diameter    = 3.74     * Units.meter
    fuselage.areas.side_projected  = 142.1948 * Units['meters**2'] 
    fuselage.areas.wetted          = 446.718  * Units['meters**2'] 
    fuselage.areas.front_projected = 12.57    * Units['meters**2'] 
    fuselage.differential_pressure = 5.0e4 * Units.pascal # Maximum differential pressure
    
    fuselage.heights.at_quarter_length          = 3.74 * Units.meter
    fuselage.heights.at_three_quarters_length   = 3.65 * Units.meter
    fuselage.heights.at_wing_root_quarter_chord = 3.74 * Units.meter
    
    # add to vehicle
    vehicle.append_component(fuselage)
    
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

    # ------------------------------------------------------------------
    #   Turbofan Network
    # ------------------------------------------------------------------    
    
    #instantiate the gas turbine network
    turbofan = SUAVE.Components.Energy.Networks.Turbofan()
    turbofan.tag = 'turbofan'
    
    # setup
    turbofan.number_of_engines = 2
    turbofan.bypass_ratio      = 5.4
    turbofan.origin            = [[13.72, 4.86,-1.9],[13.72, -4.86,-1.9]]
    
    # working fluid
    turbofan.working_fluid = SUAVE.Attributes.Gases.Air()
    
    # ------------------------------------------------------------------
    #   Component 1 - Ram
    
    # to convert freestream static to stagnation quantities
    # instantiate
    ram = SUAVE.Components.Energy.Converters.Ram()
    ram.tag = 'ram'
    
    # add to the network
    turbofan.append(ram)

    # ------------------------------------------------------------------
    #  Component 2 - Inlet Nozzle
    
    # instantiate
    inlet_nozzle = SUAVE.Components.Energy.Converters.Compression_Nozzle()
    inlet_nozzle.tag = 'inlet_nozzle'
    
    # setup
    inlet_nozzle.polytropic_efficiency = 0.98
    inlet_nozzle.pressure_ratio        = 0.98
    
    # add to network
    turbofan.append(inlet_nozzle)
    
    # ------------------------------------------------------------------
    #  Component 3 - Low Pressure Compressor
    
    # instantiate 
    compressor = SUAVE.Components.Energy.Converters.Compressor()    
    compressor.tag = 'low_pressure_compressor'

    # setup
    compressor.polytropic_efficiency = 0.91
    compressor.pressure_ratio        = 1.14    
    
    # add to network
    turbofan.append(compressor)
    
    # ------------------------------------------------------------------
    #  Component 4 - High Pressure Compressor
    
    # instantiate
    compressor = SUAVE.Components.Energy.Converters.Compressor()    
    compressor.tag = 'high_pressure_compressor'
    
    # setup
    compressor.polytropic_efficiency = 0.91
    compressor.pressure_ratio        = 13.415    
    
    # add to network
    turbofan.append(compressor)

    # ------------------------------------------------------------------
    #  Component 5 - Low Pressure Turbine
    
    # instantiate
    turbine = SUAVE.Components.Energy.Converters.Turbine()   
    turbine.tag='low_pressure_turbine'
    
    # setup
    turbine.mechanical_efficiency = 0.99
    turbine.polytropic_efficiency = 0.93     
    
    # add to network
    turbofan.append(turbine)
      
    # ------------------------------------------------------------------
    #  Component 6 - High Pressure Turbine
    
    # instantiate
    turbine = SUAVE.Components.Energy.Converters.Turbine()   
    turbine.tag='high_pressure_turbine'

    # setup
    turbine.mechanical_efficiency = 0.99
    turbine.polytropic_efficiency = 0.93     
    
    # add to network
    turbofan.append(turbine)  
    
    # ------------------------------------------------------------------
    #  Component 7 - Combustor
    
    # instantiate    
    combustor = SUAVE.Components.Energy.Converters.Combustor()   
    combustor.tag = 'combustor'
    
    # setup
    combustor.efficiency                = 0.99 
    combustor.alphac                    = 1.0   
    combustor.turbine_inlet_temperature = 1450 # K
    combustor.pressure_ratio            = 0.95
    combustor.fuel_data                 = SUAVE.Attributes.Propellants.Jet_A()    
    
    # add to network
    turbofan.append(combustor)

    # ------------------------------------------------------------------
    #  Component 8 - Core Nozzle
    
    # instantiate
    nozzle = SUAVE.Components.Energy.Converters.Expansion_Nozzle()   
    nozzle.tag = 'core_nozzle'
    
    # setup
    nozzle.polytropic_efficiency = 0.95
    nozzle.pressure_ratio        = 0.99    
    
    # add to network
    turbofan.append(nozzle)

    # ------------------------------------------------------------------
    #  Component 9 - Fan Nozzle
    
    # instantiate
    nozzle = SUAVE.Components.Energy.Converters.Expansion_Nozzle()   
    nozzle.tag = 'fan_nozzle'

    # setup
    nozzle.polytropic_efficiency = 0.95
    nozzle.pressure_ratio        = 0.99    
    
    # add to network
    turbofan.append(nozzle)
    
    # ------------------------------------------------------------------
    #  Component 10 - Fan
    
    # instantiate
    fan = SUAVE.Components.Energy.Converters.Fan()   
    fan.tag = 'fan'

    # setup
    fan.polytropic_efficiency = 0.93
    fan.pressure_ratio        = 1.7    
    
    # add to network
    turbofan.append(fan)
    
    # ------------------------------------------------------------------
    #Component 10 : thrust (to compute the thrust)
    thrust = SUAVE.Components.Energy.Processes.Thrust()       
    thrust.tag ='compute_thrust'
 
    #total design thrust (includes all the engines)
    thrust.total_design             = 2*24000. * Units.N #Newtons
 
    #design sizing conditions
    altitude      = 35000.0*Units.ft
    mach_number   = 0.78 
    isa_deviation = 0.
    
    #Engine setup for noise module    
    # add to network
    turbofan.thrust = thrust

    #size the turbofan
    turbofan_sizing(turbofan,mach_number,altitude)   
    
    # add  gas turbine network turbofan to the vehicle 
    vehicle.append_component(turbofan)      

    # ------------------------------------------------------------------
    #   Vehicle Definition Complete
    # ------------------------------------------------------------------

    return vehicle
"""
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
            

