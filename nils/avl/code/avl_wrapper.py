__all__ = ["wrap_avl"]

import ast
import numpy as np
from copy import deepcopy

import SUAVE
from SUAVE.Core import Units

from nils.suave_base_ac_object import full_setup


def wrap_avl(
    ac_segment, wing_has_segments, htail_has_segments, include_nacelles, study_idx,
    df_geom, counter, row,
) -> list:

    # Define clean-slate aircraft
    configs, analyses = full_setup()
    avl_object = analyses.configs.base.stability
    avl_object.settings.filenames.avl_bin_name = r"C:\Users\nmb48\Documents\GitHub\SUAVE\nils\avl\avl3.52\avl.exe"  # NILS: set AVL executable name
    assert id(configs.base) == id(avl_object.geometry)  # NOTE: configs.base is a Config, not a Vehicle -> see configs_setup()
    id_old = id(configs.base)
    aircraft = SUAVE.Vehicle()   # new clean vehicle
    aircraft_base = SUAVE.Vehicle()   # base version for diffing
    aircraft._base = aircraft_base  # SUAVE requirement: `_base` holds the baseline configuration
    configs.base = aircraft  # attach clean-slate "aircraft" to configs
    avl_object.geometry = aircraft  # attach clean-slate "aircraft" to avl_object
    assert id(configs.base) == id(avl_object.geometry)
    # print('aircraft.mass_properties.max_zero_fuel, aircraft.wings.keys() =', aircraft.mass_properties.max_zero_fuel, aircraft.wings.keys())  # check that configs.base has been wiped indeed
    
    # Set aircraft-dependent parameters
    if ac_segment == "narrowbody":
        aircraft.tag = 'Airbus_A220-100'
        t_tail_bool = False
    elif ac_segment == "regional":
        aircraft.tag = 'ATR_72-600'
        t_tail_bool = True
    
    # Extract mass configuration parameters
    sigma_fcs = row['sigma_fcs']
    span_loc = row['span_loc']
    fcs_loc = row['fcs_loc']
    wing_frac = row['wing_frac']
    nacelle_frac = row['nacelle_frac']
    
    # Set mass properties
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
    
    if not wing_has_segments:  # insufficient for stability analysis in AVL
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
        
    elif wing_has_segments:  # NILS: CHECK UNITS
        
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

        # Wing segments
        
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
            
        # Control surfaces

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

        # Wing segments
        
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
        
        # Control surfaces
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
    nacelle.origin = df_geom['nacelle_origin'].apply(
        lambda x: ast.literal_eval(x)[0]
    ).to_numpy()
    nacelle.flow_through = df_geom['nacelle_flow_through'][0]
    nacelle_airfoil = SUAVE.Components.Airfoils.Airfoil() 
    nacelle_airfoil.naca_4_series_airfoil = '2410'
    nacelle.append_airfoil(nacelle_airfoil)
    
    # nacelle_2 = deepcopy(nacelle)
    # nacelle_2.tag = 'nacelle_2'
    # nacelle_2_origin = deepcopy(nacelle.origin)
    # nacelle_2_origin[1] *= -1
    # nacelle_2.origin = nacelle_2_origin
    
    if include_nacelles:  # NILS: added on 26.03.2026 to exclude simple OpenVSP nacelles from AVL analysis for comparison with SU2
        aircraft.append_component(nacelle)  
        # aircraft.append_component(nacelle_2)
    # print(aircraft.nacelles.nacelle_1.Airfoil)
    
    #%% AVL-specific inputs to AVL class (different for use with JVL)
    
    tag = 'avl'
    settings_trim_aircraft = True
    backend = 'AVL'
    run_modal = True
    settings_number_spanwise_vortices = 30
    
    # NILS: longitudinal test cases from Table 6.2 in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris
    # NOTE: training inputs always have to be at least 1D, otherwise get
    # `TypeError: object of type 'float' has no len()` in
    # Documents\GitHub\SUAVE\trunk\SUAVE\Methods\Aerodynamics\AVL\translate_data.py
    
    # 6 mission points
    # training_load_factor = np.array([1.0, 2.5, 1.0, 1.0, 2.5, 1.0])
    # 1 mission point (ToC)
    training_load_factor = np.array([1.0])
    if ac_segment == 'regional':
        # 6 mission points
        # training_altitude = np.array([0.0, 10e3 * 20/35, 20e3, 20e3, 10e3 * 20/35, 0.0]) * 0.3048  # NILS: correction factors to account for different cruise altitude of ATR 72-600 than B737 in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris
        # training_Mach = np.array([0.2 * 0.43/0.7, 0.5 * 0.43/0.7, 0.43, 0.43, 0.5 * 0.43/0.7, 0.2 * 0.43/0.7])  # NILS: correction factors to account for different cruise Mach number of ATR 72-600 than B737 in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris
        # 1 mission point (ToC)
        training_altitude = np.array([20e3]) * 0.3048
        training_Mach = np.array([0.43])
    elif ac_segment == 'narrowbody':
        # 6 mission points
        # training_altitude = np.array([0.0, 10e3 * 39/35, 39e3, 39e3, 10e3 * 39/35, 0.0]) * 0.3048  # NILS: correction factors to account for different cruise altitude of A220-100 than B737 in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris
        # training_Mach = np.array([0.2 * 0.78/0.7, 0.5 * 0.78/0.7, 0.78, 0.78, 0.5 * 0.78/0.7, 0.2 * 0.78/0.7])  # NILS: correction factors to account for different cruise Mach number of A220-100 than B737 in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris
        # 1 mission point (ToC)
        training_altitude = np.array([39e3]) * 0.3048
        training_Mach = np.array([0.78])
    training_side_slip_angle = np.zeros_like(training_Mach) * Units.degrees
    training_angle_of_attack = np.array([0])  # to be trimmed  # NOTE: 6x faster if use `np.array([0])` instead of `np.zeros_like(self.training.Mach)` (6x duplication)
    # 6 mission points
    # training_mass = np.array([
    #     aircraft.mass_properties.takeoff,
    #     aircraft.mass_properties.takeoff,
    #     aircraft.mass_properties.takeoff,
    #     aircraft.mass_properties.max_zero_fuel,
    #     aircraft.mass_properties.max_zero_fuel,
    #     aircraft.mass_properties.max_zero_fuel,
    # ])
    # 1 mission point (ToC)
    training_mass = np.array([aircraft.mass_properties.max_zero_fuel])
    
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
    
    return
    
    
    