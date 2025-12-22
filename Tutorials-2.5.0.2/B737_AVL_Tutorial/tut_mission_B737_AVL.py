# tut_mission_B737_AVL.py
# 
# Created:  Mar 2018, SUAVE Team

# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# SUAVE Imports
import SUAVE
assert SUAVE.__version__=='2.5.0', 'These tutorials only work with the SUAVE 2.5.0 release'

from SUAVE.Core import Data, Units
from SUAVE.Plots.Performance.Mission_Plots import *
from SUAVE.Methods.Propulsion.turbofan_sizing import turbofan_sizing
from SUAVE.Methods.Geometry.Two_Dimensional.Cross_Section.Propulsion import compute_turbofan_geometry
from SUAVE.Input_Output.Results import  print_parasite_drag,  \
     print_compress_drag, \
     print_engine_data,   \
     print_mission_breakdown, \
     print_weight_breakdown
from SUAVE.Methods.Geometry.Two_Dimensional.Planform import segment_properties  # added by NILS

# Python Imports
import numpy as np
import pandas as pd
import pylab as plt

from copy import deepcopy

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():

    configs, analyses = full_setup()

    simple_sizing(configs)
    configs.finalize()
    analyses.finalize()

    # weight analysis
    weights = analyses.configs.base.weights
    breakdown = weights.evaluate()      

    # mission analysis
    mission = analyses.missions.base
    results = mission.evaluate()

    # NILS: commented to suppress mission plots for time being
    # plt the old results
    plot_mission(results)
    
    return

# ----------------------------------------------------------------------
#   Analysis Setup
# ----------------------------------------------------------------------

def full_setup():

    # vehicle data
    vehicle  = vehicle_setup()
    configs  = configs_setup(vehicle)

    # vehicle analyses
    configs_analyses = analyses_setup(configs)

    # mission analyses
    mission  = mission_setup(configs_analyses)
    missions_analyses = missions_setup(mission)

    analyses = SUAVE.Analyses.Analysis.Container()
    analyses.configs  = configs_analyses
    analyses.missions = missions_analyses

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
    #  Basic Geometry Relations
    sizing = SUAVE.Analyses.Sizing.Sizing()
    sizing.features.vehicle = vehicle
    analyses.append(sizing)

    # ------------------------------------------------------------------
    #  Weights
    weights = SUAVE.Analyses.Weights.Weights_Transport()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = SUAVE.Analyses.Aerodynamics.AVL()
    aerodynamics.process.compute.lift.inviscid.settings.filenames.avl_bin_name = r"C:\Users\nmb48\Documents\GitHub\SUAVE\avl3.52\avl.exe"  # NILS: set AVL executable name
    #aerodynamics.process.compute.lift.inviscid.settings.spanwise_vortex_density    = 3 
    aerodynamics.geometry = vehicle
    analyses.append(aerodynamics)

    # ------------------------------------------------------------------
    #  Stability Analysis
    stability = SUAVE.Analyses.Stability.AVL()
    stability.settings.filenames.avl_bin_name = r"C:\Users\nmb48\Documents\GitHub\SUAVE\avl3.52\avl.exe"  # NILS: set AVL executable name
    #stability.settings.spanwise_vortex_density                  = 3
    stability.geometry = vehicle
    analyses.append(stability)

    # ------------------------------------------------------------------
    #  Energy
    energy= SUAVE.Analyses.Energy.Energy()
    energy.network = vehicle.networks
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = SUAVE.Analyses.Planets.Planet()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = SUAVE.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   

    return analyses    

# ----------------------------------------------------------------------
#   Define the Vehicle
# ----------------------------------------------------------------------

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
    #  Landing Gear
    # ------------------------------------------------------------------        
    # used for noise calculations
    landing_gear = SUAVE.Components.Landing_Gear.Landing_Gear()
    landing_gear.tag = "main_landing_gear"
    
    landing_gear.main_tire_diameter = 1.12000 * Units.m
    landing_gear.nose_tire_diameter = 0.6858 * Units.m
    landing_gear.main_strut_length  = 1.8 * Units.m
    landing_gear.nose_strut_length  = 1.3 * Units.m
    landing_gear.main_units  = 2    #number of main landing gear units
    landing_gear.nose_units  = 1    #number of nose landing gear
    landing_gear.main_wheels = 2    #number of wheels on the main landing gear
    landing_gear.nose_wheels = 2    #number of wheels on the nose landing gear      
    vehicle.landing_gear = landing_gear

    # ------------------------------------------------------------------        
    #   Main Wing
    # ------------------------------------------------------------------        
    
    # NILS: Tutorials-2.5.0.2\B737_AVL_Tutorial\tut_mission_B737_AVL.py originally used a single-segment simple wing
    
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
    
    """ NILS: the below horizontal and vertical stabilisers are defined as single-segment wings
    # ------------------------------------------------------------------        
    #  Horizontal Stabilizer
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
    #   Vertical Stabilizer
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
    """
    
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

def simple_sizing(configs):

    base = configs.base
    base.pull_base()

    # zero fuel weight
    base.mass_properties.max_zero_fuel = 0.9 * base.mass_properties.max_takeoff 

    # wing areas
    for wing in base.wings:
        wing.areas.wetted   = 2.0 * wing.areas.reference
        wing.areas.exposed  = 0.8 * wing.areas.wetted
        wing.areas.affected = 0.6 * wing.areas.wetted

    # diff the new data
    base.store_diff()

    return

# ----------------------------------------------------------------------
#   Define the Mission
# ----------------------------------------------------------------------

def mission_setup(analyses):

    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = SUAVE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'the_mission'

    #airport
    airport = SUAVE.Attributes.Airports.Airport()
    airport.altitude   =  0.0  * Units.ft
    airport.delta_isa  =  0.0
    airport.atmosphere = SUAVE.Attributes.Atmospheres.Earth.US_Standard_1976()

    mission.airport = airport    

    # unpack Segments module
    Segments = SUAVE.Analyses.Mission.Segments

    # base segment
    base_segment = Segments.Segment()

    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_1"
    segment.analyses.extend( analyses.base)  
    ones_row = segment.state.ones_row
    segment.state.unknowns.body_angle = ones_row(1) * 7. * Units.deg    
    segment.altitude_start = 0.0   * Units.km
    segment.altitude_end   = 3.0   * Units.km
    segment.air_speed      = 125.0 * Units['m/s']
    segment.climb_rate     = 6.0   * Units['m/s']

    # add to misison
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Second Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_2"
    segment.analyses.extend( analyses.base)
    ones_row = segment.state.ones_row
    segment.state.unknowns.body_angle = ones_row(1) * 5. * Units.deg  
    segment.altitude_end   = 8.0   * Units.km
    segment.air_speed      = 190.0 * Units['m/s']
    segment.climb_rate     = 6.0   * Units['m/s']

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Third Climb Segment: constant Speed, Constant Rate
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_3"
    segment.analyses.extend( analyses.base)
    ones_row = segment.state.ones_row
    segment.state.unknowns.body_angle = ones_row(1) * 5. * Units.deg  
    segment.altitude_end = 10.668 * Units.km
    segment.air_speed    = 226.0  * Units['m/s']
    segment.climb_rate   = 3.0    * Units['m/s']

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------    
    #   Cruise Segment: Constant Speed, Constant Altitude
    # ------------------------------------------------------------------    

    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "cruise"
    segment.analyses.extend( analyses.base)
    segment.air_speed  = 230.412 * Units['m/s']
    segment.distance   = 2490. * Units.nautical_miles

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   First Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_1"
    segment.analyses.extend( analyses.base)
    ones_row = segment.state.ones_row
    segment.state.unknowns.body_angle = ones_row(1) * 5. * Units.deg      
    segment.altitude_end = 8.0   * Units.km
    segment.air_speed    = 220.0 * Units['m/s']
    segment.descent_rate = 4.5   * Units['m/s']

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Second Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_2"
    segment.analyses.extend( analyses.base)
    segment.altitude_end = 6.0   * Units.km
    segment.air_speed    = 195.0 * Units['m/s']
    segment.descent_rate = 5.0   * Units['m/s']

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Third Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_3"
    segment.analyses.extend( analyses.base)
    segment.altitude_end = 4.0   * Units.km
    segment.air_speed    = 170.0 * Units['m/s']
    segment.descent_rate = 5.0   * Units['m/s']

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Fourth Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_4"
    segment.analyses.extend( analyses.base)
    segment.altitude_end = 2.0   * Units.km
    segment.air_speed    = 150.0 * Units['m/s']
    segment.descent_rate = 5.0   * Units['m/s']

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Fifth Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_5"
    segment.analyses.extend( analyses.base)
    ones_row = segment.state.ones_row
    segment.state.unknowns.body_angle = ones_row(1) * 5. * Units.deg       
    segment.altitude_end = 0.0   * Units.km
    segment.air_speed    = 145.0 * Units['m/s']
    segment.descent_rate = 3.0   * Units['m/s']

    # append to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Mission definition complete    
    # ------------------------------------------------------------------

    return mission

def missions_setup(base_mission):

    # the mission container
    missions = SUAVE.Analyses.Mission.Mission.Container()

    # ------------------------------------------------------------------
    #   Base Mission
    # ------------------------------------------------------------------

    missions.base = base_mission

    return missions  

# ----------------------------------------------------------------------
#   Plot Mission
# ----------------------------------------------------------------------

def plot_mission(results,line_style='bo-'):
    
    # NILS: plot_stability_coefficients()
    plot_stability_coefficients(results, line_style)
    # import sys
    # sys.exit('Plot only stability coefficients.')

    # Plot Aerodynamic Forces 
    plot_aerodynamic_forces(results, line_style)
    
    # Plot Aerodynamic Coefficients 
    plot_aerodynamic_coefficients(results, line_style)
    
    # Drag Components
    plot_drag_components(results, line_style)    
    
    # Plot Velocities 
    plot_aircraft_velocities(results, line_style)          
        
    return

#%%

if __name__ == '__main__':
    
    # NILS: original call (runs various analyses, including stability)
    # main()    
    # plt.show()
    
    # NILS: perform stability analysis only
    configs, analyses = full_setup()
    avl_object = analyses.configs.base.stability
    aircraft = avl_object.geometry  # corresponds to `vehicle`
    
    # To modify STATIC stability
    # aircraft.mass_properties.center_of_gravity[0][0] += 2
    # To modify DYNAMIC stability
    # aircraft.mass_properties.moments_of_inertia.tensor[0][0] *= 2
    # aircraft.mass_properties.moments_of_inertia.tensor[1][1] *= 2
    # aircraft.mass_properties.moments_of_inertia.tensor[2][2] *= 2
    
    # Run sweep over flight condition parameters only
    # avl_object.sample_training()
    
    #%% Extract data for verification of TASOPT-SUAVE-AVL interface
    # NOTE: comment `avl_object.sample_training()` line above, else will get
    # `UnsupportedOperation: fileno` error (SUAVE-AVL WRAPPER CAN ONLY BE
    # EXECUTED FROM A COMMAND PROMPT, NOT FROM AN IPYTHON CONSOLE LIKE SPYDER)!
    
    avl_object = analyses.configs.base.stability
    # aircraft = avl_object.geometry._base  # corresponds to `vehicle`
    aircraft = avl_object.geometry  # corresponds to `vehicle`
    
    x_cg = aircraft.mass_properties.center_of_gravity[0][0]
    y_cg = aircraft.mass_properties.center_of_gravity[0][1]
    z_cg = aircraft.mass_properties.center_of_gravity[0][2]
    mass = aircraft.mass_properties.mass
    max_takeoff = aircraft.mass_properties.max_takeoff
    takeoff = aircraft.mass_properties.takeoff
    max_zero_fuel = aircraft.mass_properties.max_zero_fuel
    moments_of_inertia = aircraft.mass_properties.moments_of_inertia.tensor
    Ixx = moments_of_inertia[0][0]
    Iyy = moments_of_inertia[1][1]
    Izz = moments_of_inertia[2][2]
    Ixy = moments_of_inertia[0][1]
    Iyz = moments_of_inertia[1][2]
    Izx = moments_of_inertia[2][0]
    
    for suave_body in aircraft.fuselages:  # write_geometry.py
        body_lengths_total = suave_body.lengths.total
        body_lengths_nose = suave_body.lengths.nose
        body_lengths_tail = suave_body.lengths.tail
        body_widths_maximum = suave_body.width
        body_heights_maximum = suave_body.heights.maximum
    
    wing_dict = {name: {} for name in aircraft.wings.keys()}
    for tag, suave_wing in aircraft.wings.items():  # write_geometry.py
    
        wing_dict[tag]['spans_projected'] = suave_wing.spans.projected
        wing_dict[tag]['origin_x'] = suave_wing.origin[0][0]
        wing_dict[tag]['origin_y'] = suave_wing.origin[0][1]
        wing_dict[tag]['origin_z'] = suave_wing.origin[0][2]
        wing_dict[tag]['dihedral'] = suave_wing.dihedral
        
        # define root section 
        wing_dict[tag]['chords_root'] = suave_wing.chords.root
        wing_dict[tag]['twists_root'] = suave_wing.twists.root
    
        # define tip section
        wing_dict[tag]['chords_tip'] = suave_wing.chords.tip
        wing_dict[tag]['twists_tip'] = suave_wing.twists.tip
        
        
        wing_dict[tag]['aspect_ratio'] = suave_wing.aspect_ratio
        wing_dict[tag]['sweeps_quarter_chord'] = suave_wing.sweeps.quarter_chord
        wing_dict[tag]['thickness_to_chord'] = suave_wing.thickness_to_chord
        wing_dict[tag]['taper'] = suave_wing.taper
        wing_dict[tag]['chords_mean_aerodynamic'] = suave_wing.chords.mean_aerodynamic
        wing_dict[tag]['areas_reference'] = suave_wing.areas.reference
        wing_dict[tag]['vertical'] = suave_wing.vertical
        wing_dict[tag]['symmetric'] = suave_wing.symmetric
        wing_dict[tag]['high_lift'] = suave_wing.high_lift
        
    #%%
    
    main_wing_dict = wing_dict['main_wing']
    horizontal_stabilizer_dict = wing_dict['horizontal_stabilizer']
    vertical_stabilizer_dict = wing_dict['vertical_stabilizer']
    
    main_wing = aircraft.wings['main_wing']
    htail = aircraft.wings['horizontal_stabilizer']
    
    df = pd.DataFrame({
        "x_cg": [x_cg],
        "y_cg": [y_cg],
        "z_cg": [z_cg],
        "mass": [mass],
        "max_takeoff": [max_takeoff],
        "takeoff": [takeoff],
        "max_zero_fuel": [max_zero_fuel],
    
        "Ixx": Ixx,
        "Iyy": Iyy,
        "Izz": Izz,
        "Ixy": Ixy,
        "Iyz": Iyz,
        "Izx": Izx,
    
        "body_lengths_total": body_lengths_total,
        "body_lengths_nose": body_lengths_nose,
        "body_lengths_tail": body_lengths_tail,
        "body_widths_maximum": body_widths_maximum,
        "body_heights_maximum": body_heights_maximum,
        
        # Main wing
        
        # "wing_aspect_ratio": main_wing_dict['aspect_ratio'],
        # "wing_sweeps_quarter_chord": main_wing_dict['sweeps_quarter_chord'],
        # "wing_thickness_to_chord": main_wing_dict['thickness_to_chord'],
        # "wing_taper": main_wing_dict['taper'],
        # "wing_spans_projected": main_wing_dict['spans_projected'],
        # "wing_chords_root": main_wing_dict['chords_root'],
        # "wing_chords_tip": main_wing_dict['chords_tip'],
        # "wing_chords_mean_aerodynamic": main_wing_dict['chords_mean_aerodynamic'],
        # "wing_areas_reference": main_wing_dict['areas_reference'],
        # "wing_twists_root": main_wing_dict['twists_root'],
        # "wing_twists_tip": main_wing_dict['twists_tip'],
        # "wing_origin_x": main_wing_dict['origin_x'],
        # "wing_origin_y": main_wing_dict['origin_y'],
        # "wing_origin_z": main_wing_dict['origin_z'],
        # "wing_vertical": main_wing_dict['vertical'],
        # "wing_symmetric": main_wing_dict['symmetric'],
        # "wing_high_lift": main_wing_dict['high_lift'],
        # "wing_dihedral": main_wing_dict['dihedral'],
        
        # =============================================================================
        'wing_aspect_ratio': [main_wing.aspect_ratio],
        'wing_sweeps_quarter_chord': [main_wing.sweeps.quarter_chord],
        'wing_thickness_to_chord': [main_wing.thickness_to_chord],
        'wing_taper': [main_wing.taper],
        'wing_spans_projected': [main_wing.spans.projected],
        'wing_chords_root': [main_wing.chords.root],
        'wing_chords_tip': [main_wing.chords.tip],
        'wing_chords_mean_aerodynamic': [main_wing.chords.mean_aerodynamic],
        'wing_areas_reference': [main_wing.areas.reference],
        'wing_twists_root': [main_wing.twists.root],
        'wing_twists_tip': [main_wing.twists.tip],
        'wing_origin_x': [main_wing.origin[0][0]],
        'wing_origin_y': [main_wing.origin[0][1]],
        'wing_origin_z': [main_wing.origin[0][2]],
        'wing_vertical': [main_wing.vertical],
        'wing_symmetric': [main_wing.symmetric],
        'wing_high_lift': [main_wing.high_lift],
        'wing_dihedral': [main_wing.dihedral],
        
        'wing_root_percent_span_location': [main_wing.Segments['root'].percent_span_location],
        'wing_root_twist': [main_wing.Segments['root'].twist],
        'wing_root_root_chord_percent': [main_wing.Segments['root'].root_chord_percent],
        'wing_root_thickness_to_chord': [main_wing.Segments['root'].thickness_to_chord],
        'wing_root_dihedral_outboard': [main_wing.Segments['root'].dihedral_outboard],
        'wing_root_sweeps_quarter_chord': [main_wing.Segments['root'].sweeps.quarter_chord],
        
        'wing_yehudi_percent_span_location': [main_wing.Segments['yehudi'].percent_span_location],
        'wing_yehudi_twist': [main_wing.Segments['yehudi'].twist],
        'wing_yehudi_root_chord_percent': [main_wing.Segments['yehudi'].root_chord_percent],
        'wing_yehudi_thickness_to_chord': [main_wing.Segments['yehudi'].thickness_to_chord],
        'wing_yehudi_dihedral_outboard': [main_wing.Segments['yehudi'].dihedral_outboard],
        'wing_yehudi_sweeps_quarter_chord': [main_wing.Segments['yehudi'].sweeps.quarter_chord],

        'wing_section2_percent_span_location': [main_wing.Segments['section_2'].percent_span_location],
        'wing_section2_twist': [main_wing.Segments['section_2'].twist],
        'wing_section2_root_chord_percent': [main_wing.Segments['section_2'].root_chord_percent],
        'wing_section2_thickness_to_chord': [main_wing.Segments['section_2'].thickness_to_chord],
        'wing_section2_dihedral_outboard': [main_wing.Segments['section_2'].dihedral_outboard],
        'wing_section2_sweeps_quarter_chord': [main_wing.Segments['section_2'].sweeps.quarter_chord],
        
        'wing_tip_percent_span_location': [main_wing.Segments['tip'].percent_span_location],
        'wing_tip_twist': [main_wing.Segments['tip'].twist],
        'wing_tip_root_chord_percent': [main_wing.Segments['tip'].root_chord_percent],
        'wing_tip_thickness_to_chord': [main_wing.Segments['tip'].thickness_to_chord],
        'wing_tip_dihedral_outboard': [main_wing.Segments['tip'].dihedral_outboard],
        'wing_tip_sweeps_quarter_chord': [main_wing.Segments['tip'].sweeps.quarter_chord],
        
        'wing_slat_span_fraction_start': [main_wing.control_surfaces['slat'].span_fraction_start],
        'wing_slat_span_fraction_end': [main_wing.control_surfaces['slat'].span_fraction_end],
        'wing_slat_deflection': [main_wing.control_surfaces['slat'].deflection],
        'wing_slat_chord_fraction': [main_wing.control_surfaces['slat'].chord_fraction],
        
        'wing_flap_span_fraction_start': [main_wing.control_surfaces['flap'].span_fraction_start],
        'wing_flap_span_fraction_end': [main_wing.control_surfaces['flap'].span_fraction_end],
        'wing_flap_deflection': [main_wing.control_surfaces['flap'].deflection],
        'wing_flap_configuration_type': [main_wing.control_surfaces['flap'].configuration_type],
        'wing_flap_chord_fraction': [main_wing.control_surfaces['flap'].chord_fraction],
        
        'wing_aileron_span_fraction_start': [main_wing.control_surfaces['aileron'].span_fraction_start],
        'wing_aileron_span_fraction_end': [main_wing.control_surfaces['aileron'].span_fraction_end],
        'wing_aileron_deflection': [main_wing.control_surfaces['aileron'].deflection],
        'wing_aileron_chord_fraction': [main_wing.control_surfaces['aileron'].chord_fraction],
        # =============================================================================
        
        # Horizontal stabiliser
        # "htail_aspect_ratio": horizontal_stabilizer_dict['aspect_ratio'],
        # "htail_sweeps_quarter_chord": horizontal_stabilizer_dict['sweeps_quarter_chord'],
        # "htail_thickness_to_chord": horizontal_stabilizer_dict['thickness_to_chord'],
        # "htail_taper": horizontal_stabilizer_dict['taper'],
        # "htail_spans_projected": horizontal_stabilizer_dict['spans_projected'],
        # "htail_chords_root": horizontal_stabilizer_dict['chords_root'],
        # "htail_chords_tip": horizontal_stabilizer_dict['chords_tip'],
        # "htail_chords_mean_aerodynamic": horizontal_stabilizer_dict['chords_mean_aerodynamic'],
        # "htail_areas_reference": horizontal_stabilizer_dict['areas_reference'],
        # "htail_twists_root": horizontal_stabilizer_dict['twists_root'],
        # "htail_twists_tip": horizontal_stabilizer_dict['twists_tip'],
        # "htail_origin_x": horizontal_stabilizer_dict['origin_x'],
        # "htail_origin_y": horizontal_stabilizer_dict['origin_y'],
        # "htail_origin_z": horizontal_stabilizer_dict['origin_z'],
        # "htail_vertical": horizontal_stabilizer_dict['vertical'],
        # "htail_symmetric": horizontal_stabilizer_dict['symmetric'],
        # "htail_high_lift": horizontal_stabilizer_dict['high_lift'],
        # "htail_dihedral": horizontal_stabilizer_dict['dihedral'],
        
        'htail_aspect_ratio': [htail.aspect_ratio],
        'htail_sweeps_quarter_chord': [htail.sweeps.quarter_chord],
        'htail_thickness_to_chord': [htail.thickness_to_chord],
        'htail_taper': [htail.taper],
        'htail_spans_projected': [htail.spans.projected],
        'htail_chords_root': [htail.chords.root],
        'htail_chords_tip': [htail.chords.tip],
        'htail_chords_mean_aerodynamic': [htail.chords.mean_aerodynamic],
        'htail_areas_reference': [htail.areas.reference],
        'htail_twists_root': [htail.twists.root],
        'htail_twists_tip': [htail.twists.tip],
        'htail_origin_x': [htail.origin[0][0]],
        'htail_origin_y': [htail.origin[0][1]],
        'htail_origin_z': [htail.origin[0][2]],
        'htail_vertical': [htail.vertical],
        'htail_symmetric': [htail.symmetric],
        'htail_high_lift': [htail.high_lift],
        'htail_dihedral': [htail.dihedral],
        
        'htail_root_percent_span_location': [htail.Segments['root_segment'].percent_span_location],
        'htail_root_twist': [htail.Segments['root_segment'].twist],
        'htail_root_root_chord_percent': [htail.Segments['root_segment'].root_chord_percent],
        'htail_root_thickness_to_chord': [htail.Segments['root_segment'].thickness_to_chord],
        'htail_root_dihedral_outboard': [htail.Segments['root_segment'].dihedral_outboard],
        'htail_root_sweeps_quarter_chord': [htail.Segments['root_segment'].sweeps.quarter_chord],
        
        'htail_tip_percent_span_location': [htail.Segments['tip_segment'].percent_span_location],
        'htail_tip_twist': [htail.Segments['tip_segment'].twist],
        'htail_tip_root_chord_percent': [htail.Segments['tip_segment'].root_chord_percent],
        'htail_tip_thickness_to_chord': [htail.Segments['tip_segment'].thickness_to_chord],
        'htail_tip_dihedral_outboard': [htail.Segments['tip_segment'].dihedral_outboard],
        'htail_tip_sweeps_quarter_chord': [htail.Segments['tip_segment'].sweeps.quarter_chord],
        
        'htail_elevator_span_fraction_start': [htail.control_surfaces['elevator'].span_fraction_start],
        'htail_elevator_span_fraction_end': [htail.control_surfaces['elevator'].span_fraction_end],
        'htail_elevator_deflection': [htail.control_surfaces['elevator'].deflection],
        'htail_elevator_chord_fraction': [htail.control_surfaces['elevator'].chord_fraction],
        
        # Vertical stabiliser
        "vtail_aspect_ratio": vertical_stabilizer_dict['aspect_ratio'],
        "vtail_sweeps_quarter_chord": vertical_stabilizer_dict['sweeps_quarter_chord'],
        "vtail_thickness_to_chord": vertical_stabilizer_dict['thickness_to_chord'],
        "vtail_taper": vertical_stabilizer_dict['taper'],
        "vtail_spans_projected": vertical_stabilizer_dict['spans_projected'],
        "vtail_chords_root": vertical_stabilizer_dict['chords_root'],
        "vtail_chords_tip": vertical_stabilizer_dict['chords_tip'],
        "vtail_chords_mean_aerodynamic": vertical_stabilizer_dict['chords_mean_aerodynamic'],
        "vtail_areas_reference": vertical_stabilizer_dict['areas_reference'],
        "vtail_twists_root": vertical_stabilizer_dict['twists_root'],
        "vtail_twists_tip": vertical_stabilizer_dict['twists_tip'],
        "vtail_origin_x": vertical_stabilizer_dict['origin_x'],
        "vtail_origin_y": vertical_stabilizer_dict['origin_y'],
        "vtail_origin_z": vertical_stabilizer_dict['origin_z'],
        "vtail_vertical": vertical_stabilizer_dict['vertical'],
        "vtail_symmetric": vertical_stabilizer_dict['symmetric'],
        "vtail_high_lift": vertical_stabilizer_dict['high_lift'],
        "vtail_dihedral": vertical_stabilizer_dict['dihedral'],
    })
    
    df.to_csv(r"C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\B737_AVL_Tutorial\suave_avl_wrapper_tasopt_inputs_b737.csv", index=False)

    
    