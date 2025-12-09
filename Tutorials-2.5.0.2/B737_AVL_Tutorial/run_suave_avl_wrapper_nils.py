"""
This file is adapted from tut_mission_B737_AVL.py
to contain the minimum set of SUAVE methods needed to
perform a static and dynamic stability analysis on
aircraft geometry data imported as a .csv from TASOPT.jl.
"""

__all__ = []

import numpy as np
import pylab as plt
from copy import deepcopy

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
    stability.settings.filenames.avl_bin_name = r"C:\Users\nmb48\Documents\GitHub\SUAVE\avl3.52\avl.exe"  # NILS: set AVL executable name
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
    vehicle.mass_properties.max_takeoff               = 79015.8 * Units.kilogram
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

    # ------------------------------------------------------------------        
    #  Horizontal Stabilizer - NILS: treat as SINGLE-SECTION wing
    # ------------------------------------------------------------------        
    
    wing = SUAVE.Components.Wings.Horizontal_Tail()
    wing.tag = 'horizontal_stabilizer'
    
    wing.aspect_ratio            = 6.16     
    wing.sweeps.quarter_chord    = 40 * Units.deg
    wing.thickness_to_chord      = 0.08
    wing.taper                   = 0.2
    wing.spans.projected         = 14.2 * Units.meter
    wing.chords.root             = 4.7  * Units.meter
    wing.chords.tip              = .955 * Units.meter
    wing.chords.mean_aerodynamic = 3.0  * Units.meter
    wing.areas.reference         = 32.488   * Units['meters**2']  
    wing.twists.root             = 3.0 * Units.degrees
    wing.twists.tip              = 3.0 * Units.degrees  
    wing.origin                  = [[32.83 * Units.meter, 0 , 1.14 * Units.meter]]
    wing.vertical                = False 
    wing.symmetric               = True
    wing.dynamic_pressure_ratio  = 0.9  
    
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
    
    import pandas as pd

    configs, analyses = full_setup()
    avl_object = analyses.configs.base.stability
    
    keep_b737_defaults = False  # True
    wing_has_segments = True  # False
    
    if keep_b737_defaults == False:
        
        # Read geometry data from .csv file
        df = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl\suave_avl_wrapper_tasopt_inputs.csv')
        # df = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl\suave_avl_wrapper_tasopt_inputs_b737.csv')
        
        # # aircraft = avl_object.geometry._base  # correspons to `vehicle`
        # aircraft = avl_object.geometry  # correspons to `vehicle`
        aircraft = SUAVE.Vehicle()
        aircraft.tag = 'Airbus_A220-100'
    
        # r'''
        aircraft.mass_properties.center_of_gravity[0][0] = df['x_cg'][0]
        aircraft.mass_properties.center_of_gravity[0][1] = df['y_cg'][0]
        aircraft.mass_properties.center_of_gravity[0][2] = df['z_cg'][0]
        aircraft.mass_properties.max_takeoff = df['mass'][0]
        moments_of_inertia = aircraft.mass_properties.moments_of_inertia.tensor
        moments_of_inertia[0][0] = df['Ixx'][0]
        moments_of_inertia[1][1] = df['Iyy'][0]
        moments_of_inertia[2][2] = df['Izz'][0]
        moments_of_inertia[0][1] = df['Ixy'][0]
        moments_of_inertia[1][2] = df['Iyz'][0]
        moments_of_inertia[2][0] = df['Izx'][0]
        
        fuselage = SUAVE.Components.Fuselages.Fuselage()
        fuselage.tag = 'fuselage'
        
        # for suave_body in aircraft.fuselages:  # write_geometry.py
        fuselage.lengths.total = df['body_lengths_total'][0]
        fuselage.lengths.nose = df['body_lengths_nose'][0]
        fuselage.lengths.tail = df['body_lengths_tail'][0]
        fuselage.width = df['body_widths_maximum'][0]
        fuselage.heights.maximum = df['body_heights_maximum'][0]
        
        aircraft.append_component(fuselage)
        
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
            
        # Main wing
        
        wing = SUAVE.Components.Wings.Main_Wing()
        wing.tag = 'main_wing'
        
        if wing_has_segments == False:
            
            wing.aspect_ratio = df['wing_aspect_ratio'][0]
            wing.sweeps.quarter_chord = df['wing_sweeps_quarter_chord'][0]
            wing.thickness_to_chord = df['wing_thickness_to_chord'][0]
            wing.taper = df['wing_taper'][0]
            wing.spans.projected = df['wing_spans_projected'][0]
            wing.chords.root = df['wing_chords_root'][0]
            wing.chords.tip = df['wing_chords_tip'][0]
            wing.chords.mean_aerodynamic = df['wing_chords_mean_aerodynamic'][0]
            wing.areas.reference = df['wing_areas_reference'][0]
            wing.twists.root = df['wing_twists_root'][0]
            wing.twists.tip = df['wing_twists_tip'][0]
            wing.origin = [[
                df['wing_origin_x'][0],
                df['wing_origin_y'][0],
                df['wing_origin_z'][0],
            ]]
            wing.vertical = df['wing_vertical'][0]
            wing.symmetric = df['wing_symmetric'][0]
            wing.high_lift = df['wing_high_lift'][0]
            wing.dihedral = df['wing_dihedral'][0]
            
        elif wing_has_segments == True:
            
            # <<< @NILS: CHECK UNITS! >>>
            
            wing.aspect_ratio = df['wing_aspect_ratio'][0]
            wing.sweeps.quarter_chord = df['wing_sweeps_quarter_chord'][0]
            # wing.thickness_to_chord = df['wing_thickness_to_chord'][0]
            # wing.taper = df['wing_taper'][0]
            wing.spans.projected = df['wing_spans_projected'][0]
            wing.chords.root = df['wing_chords_root'][0]
            wing.chords.tip = df['wing_chords_tip'][0]
            wing.chords.mean_aerodynamic = df['wing_chords_mean_aerodynamic'][0]
            wing.areas.reference = df['wing_areas_reference'][0]
            wing.twists.root = df['wing_twists_root'][0]
            wing.twists.tip = df['wing_twists_tip'][0]
            wing.origin = [[
                df['wing_origin_x'][0],
                df['wing_origin_y'][0],
                df['wing_origin_z'][0],
            ]]
            wing.vertical = df['wing_vertical'][0]
            wing.symmetric = df['wing_symmetric'][0]
            wing.high_lift = df['wing_high_lift'][0]
            # wing.dihedral = df['wing_dihedral'][0]
            
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
            
            # # =============================================================================
            # root_airfoil                          = SUAVE.Components.Airfoils.Airfoil()
            # root_airfoil.coordinate_file          = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737a.txt'
            # segment                               = SUAVE.Components.Wings.Segment()
            # segment.tag                           = 'Root'
            # segment.percent_span_location         = df['wing_root_percent_span_location'][0]
            # segment.twist                         = df['wing_root_twist'][0]
            # segment.root_chord_percent            = df['wing_root_root_chord_percent'][0]
            # segment.thickness_to_chord            = df['wing_root_thickness_to_chord'][0]
            # segment.dihedral_outboard             = df['wing_root_dihedral_outboard'][0]
            # segment.sweeps.quarter_chord          = df['wing_root_sweeps_quarter_chord'][0]
            # segment.append_airfoil(root_airfoil)
            # wing.append_segment(segment)

            # yehudi_airfoil                       = SUAVE.Components.Airfoils.Airfoil()
            # yehudi_airfoil.coordinate_file       = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737b.txt'
            # segment                               = SUAVE.Components.Wings.Segment()
            # segment.tag                           = 'Yehudi'
            # segment.percent_span_location         = df['wing_yehudi_percent_span_location'][0]
            # segment.twist                         = df['wing_yehudi_twist'][0]
            # segment.root_chord_percent            = df['wing_yehudi_root_chord_percent'][0]
            # segment.thickness_to_chord            = df['wing_yehudi_thickness_to_chord'][0]
            # segment.dihedral_outboard             = df['wing_yehudi_dihedral_outboard'][0]
            # segment.sweeps.quarter_chord          = df['wing_yehudi_sweeps_quarter_chord'][0]
            # segment.append_airfoil(yehudi_airfoil)
            # wing.append_segment(segment)

            # section2_airfoil                      =  SUAVE.Components.Airfoils.Airfoil()
            # section2_airfoil.coordinate_file      = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737c.txt'
            # segment                               = SUAVE.Components.Wings.Segment()
            # segment.tag                           = 'Section 2'
            # segment.percent_span_location         = df['wing_section2_percent_span_location'][0]
            # segment.twist                         = df['wing_section2_twist'][0]
            # segment.root_chord_percent            = df['wing_section2_root_chord_percent'][0]
            # segment.thickness_to_chord            = df['wing_section2_thickness_to_chord'][0]
            # segment.dihedral_outboard             = df['wing_section2_dihedral_outboard'][0]
            # segment.sweeps.quarter_chord          = df['wing_section2_sweeps_quarter_chord'][0]
            # segment.append_airfoil(section2_airfoil)
            # wing.append_segment(segment)
            
            # tip_airfoil                      =  SUAVE.Components.Airfoils.Airfoil()
            # tip_airfoil.coordinate_file      = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737c.txt'
            # segment                               = SUAVE.Components.Wings.Segment()
            # segment.tag                           = 'Tip'
            # segment.percent_span_location         = df['wing_tip_percent_span_location'][0]
            # segment.twist                         = df['wing_tip_twist'][0]
            # segment.root_chord_percent            = df['wing_tip_root_chord_percent'][0]
            # segment.thickness_to_chord            = df['wing_tip_thickness_to_chord'][0]
            # segment.dihedral_outboard             = df['wing_tip_dihedral_outboard'][0]
            # segment.sweeps.quarter_chord          = df['wing_tip_sweeps_quarter_chord'][0]
            # segment.append_airfoil(tip_airfoil)
            # wing.append_segment(segment)
            # # =============================================================================
            
            # =============================================================================
            center_airfoil                          = SUAVE.Components.Airfoils.Airfoil()
            center_airfoil.coordinate_file          = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737a.txt'
            segment                               = SUAVE.Components.Wings.Segment()
            segment.tag                           = 'Root'
            segment.percent_span_location         = df['wing_center_percent_span_location'][0]
            segment.twist                         = df['wing_center_twist'][0]
            segment.center_chord_percent            = df['wing_center_root_chord_percent'][0]
            segment.thickness_to_chord            = df['wing_center_thickness_to_chord'][0]
            segment.dihedral_outboard             = df['wing_center_dihedral_outboard'][0]
            segment.sweeps.quarter_chord          = df['wing_center_sweeps_quarter_chord'][0]
            segment.append_airfoil(center_airfoil)
            wing.append_segment(segment)

            inboard_airfoil                       = SUAVE.Components.Airfoils.Airfoil()
            inboard_airfoil.coordinate_file       = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737b.txt'
            segment                               = SUAVE.Components.Wings.Segment()
            segment.tag                           = 'Yehudi'
            segment.percent_span_location         = df['wing_inboard_percent_span_location'][0]
            segment.twist                         = df['wing_inboard_twist'][0]
            segment.root_chord_percent            = df['wing_inboard_root_chord_percent'][0]
            segment.thickness_to_chord            = df['wing_inboard_thickness_to_chord'][0]
            segment.dihedral_outboard             = df['wing_inboard_dihedral_outboard'][0]
            segment.sweeps.quarter_chord          = df['wing_inboard_sweeps_quarter_chord'][0]
            segment.append_airfoil(inboard_airfoil)
            wing.append_segment(segment)

            outboard_airfoil                      =  SUAVE.Components.Airfoils.Airfoil()
            outboard_airfoil.coordinate_file      = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737c.txt'
            segment                               = SUAVE.Components.Wings.Segment()
            segment.tag                           = 'Section 2'
            segment.percent_span_location         = df['wing_outboard_percent_span_location'][0]
            segment.twist                         = df['wing_outboard_twist'][0]
            segment.root_chord_percent            = df['wing_outboard_root_chord_percent'][0]
            segment.thickness_to_chord            = df['wing_outboard_thickness_to_chord'][0]
            segment.dihedral_outboard             = df['wing_outboard_dihedral_outboard'][0]
            segment.sweeps.quarter_chord          = df['wing_outboard_sweeps_quarter_chord'][0]
            segment.append_airfoil(outboard_airfoil)
            wing.append_segment(segment)
            # =============================================================================
            
        aircraft.append_component(wing)
        
        # elif tag == 'horizontal_stabilizer':
        
        # Horizontal stabiliser
        
        wing = SUAVE.Components.Wings.Horizontal_Tail()
        wing.tag = 'horizontal_stabilizer'
        
        wing.aspect_ratio = df['htail_aspect_ratio'][0]
        wing.sweeps.quarter_chord = df['htail_sweeps_quarter_chord'][0]
        wing.thickness_to_chord = df['htail_thickness_to_chord'][0]
        wing.taper = df['htail_taper'][0]
        wing.spans.projected = df['htail_spans_projected'][0]
        wing.chords.root = df['htail_chords_root'][0]
        wing.chords.tip = df['htail_chords_tip'][0]
        # wing.chords.mean_aerodynamic = df['htail_chords_mean_aerodynamic'][0]
        wing.areas.reference = df['htail_areas_reference'][0]
        wing.twists.root = df['htail_twists_root'][0]
        wing.twists.tip = df['htail_twists_tip'][0]
        wing.origin = [[
            df['htail_origin_x'][0],
            df['htail_origin_y'][0],
            df['htail_origin_z'][0],
        ]]
        wing.vertical = df['htail_vertical'][0]
        wing.symmetric = df['htail_symmetric'][0]
        wing.high_lift = df['htail_high_lift'][0]
        wing.dihedral = df['htail_dihedral'][0]
        
        aircraft.append_component(wing)
        
        # elif tag == 'vertical_stabilizer':
        
        # Vertical stabiliser
            
        wing.aspect_ratio = df['vtail_aspect_ratio'][0]
        wing.sweeps.quarter_chord = df['vtail_sweeps_quarter_chord'][0]
        wing.thickness_to_chord = df['vtail_thickness_to_chord'][0]
        wing.taper = df['vtail_taper'][0]
        wing.spans.projected = df['vtail_spans_projected'][0]
        wing.chords.root = df['vtail_chords_root'][0]
        wing.chords.tip = df['vtail_chords_tip'][0]
        # wing.chords.mean_aerodynamic = df['vtail_chords_mean_aerodynamic'][0]
        wing.areas.reference = df['vtail_areas_reference'][0]
        wing.twists.root = df['vtail_twists_root'][0]
        wing.twists.tip = df['vtail_twists_tip'][0]
        wing.origin = [[
            df['vtail_origin_x'][0],
            df['vtail_origin_y'][0],
            df['vtail_origin_z'][0],
        ]]
        wing.vertical = df['vtail_vertical'][0]
        wing.symmetric = df['vtail_symmetric'][0]
        wing.high_lift = df['vtail_high_lift'][0]
        wing.dihedral = df['vtail_dihedral'][0]
        
        # =============================================================================
        # '''
        #%%
        
        # print(aircraft.mass_properties.center_of_gravity[0][0])
        # print(aircraft.mass_properties.center_of_gravity[0][1])
        # print(aircraft.mass_properties.center_of_gravity[0][2])
        # print(aircraft.mass_properties.mass)
        # moments_of_inertia = aircraft.mass_properties.moments_of_inertia.tensor
        # print(moments_of_inertia)
        # print(moments_of_inertia[0][0])
        # print(moments_of_inertia[1][1])
        # print(moments_of_inertia[2][2])
        # print(moments_of_inertia[0][1])
        # print(moments_of_inertia[1][2])
        # print(moments_of_inertia[2][0])
        
        # for suave_body in aircraft.fuselages:  # write_geometry.py
        #     print(suave_body.lengths.total)
        #     print(suave_body.lengths.nose)
        #     print(suave_body.lengths.tail)
        #     print(suave_body.width)
        #     print(suave_body.heights.maximum)
        
        # for tag, suave_wing in aircraft.wings.items():  # write_geometry.py
                    
        #     print(suave_wing.spans.projected)
        #     print(suave_wing.origin)
        #     print(suave_wing.dihedral)
            
        #     # print(suave_wing.sweeps.leading_edge)
            
        #     # define root section 
        #     print(suave_wing.chords.root)
        #     print(suave_wing.twists.root)
    
        #     # define tip section
        #     print(suave_wing.chords.tip)
        #     print(suave_wing.twists.tip)
            
        #     print(suave_wing.sweeps.quarter_chord)
        #     print(suave_wing.thickness_to_chord)
        #     print(suave_wing.taper)
        #     print(suave_wing.chords.mean_aerodynamic)
        #     print(suave_wing.areas.reference)
        #     print(suave_wing.vertical)
        #     print(suave_wing.symmetric)
        #     print(suave_wing.high_lift)
            
    # sys.exit('Stop.')
    
    #%%    
    
    avl_object.sample_training()



    
    
    
    