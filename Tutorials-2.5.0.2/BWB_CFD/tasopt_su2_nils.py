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
from datetime import date
import pandas as pd

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main(vehicle=None):  # NILS: added default argument

    configs, analyses = full_setup(vehicle)  # NILS: added argument

    simple_sizing(configs)

    configs.finalize()
    analyses.finalize()
    '''
    # weight analysis
    weights = analyses.configs.base.weights   

    # mission analysis
    mission = analyses.missions.base
    results = mission.evaluate()

    # plt the old results
    plot_mission(results)
    '''
    return

# ----------------------------------------------------------------------
#   Analysis Setup
# ----------------------------------------------------------------------

def full_setup(vehicle=None):  # NILS: from run_suave_avl_wrapper_nils.py

    # vehicle data
    # vehicle  = vehicle_setup()
    if vehicle == None:
        vehicle  = vehicle_setup()  # NILS: from run_suave_avl_wrapper_nils.py
    configs  = configs_setup(vehicle)

    # vehicle analyses
    configs_analyses = analyses_setup(configs)
    '''
    # mission analyses
    mission  = mission_setup(configs_analyses)
    missions_analyses = missions_setup(mission)
    '''
    analyses = SUAVE.Analyses.Analysis.Container()
    analyses.configs  = configs_analyses
    # analyses.missions = missions_analyses  # NILS

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
    '''
    # ------------------------------------------------------------------
    #  Weights
    weights = SUAVE.Analyses.Weights.Weights_BWB()  # NILS: inapplicable to tube-and-wing aircraft, but not used anyway
    weights.vehicle = vehicle
    analyses.append(weights)
    '''
    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = SUAVE.Analyses.Aerodynamics.SU2_Euler()
    aerodynamics.geometry = vehicle
    
    #aerodynamics.process.compute.lift.inviscid.settings.parallel          = True
    #aerodynamics.process.compute.lift.inviscid.settings.processors        = 12  
    # aerodynamics.process.compute.lift.inviscid.training_file              = 'base_data_1500.txt'  # NILS: commented, othwerwise surrogate model based on precomputed CFD results will be used
    ###
    # aerodynamics.process.compute.lift.inviscid.training_file              = 'base_data.txt'  # NILS: uncomment if want to perform mission analysis using surrogate model
    ###
    aerodynamics.process.compute.lift.inviscid.settings.maximum_iterations = 500  # NILS: increase from original value of 10
    
    aerodynamics.settings.drag_coefficient_increment = 0.0000
    aerodynamics.settings.half_mesh_flag             = False
    aerodynamics.settings.span_efficiency            = 0.85  # NILS: kept as in BWB.py for now, as not present in tut_mission_B737_AVL.py
    
    aerodynamics.process.compute.lift.inviscid.training.Mach               = np.array([.3, .5, .7, .85]) 
    aerodynamics.process.compute.lift.inviscid.training.angle_of_attack    = np.array([0.,3.,6.]) * Units.deg    
    '''
    wing_segments = vehicle.wings.main_wing.Segments
    wing_segments.section_1.vsp_mesh = Data()
    wing_segments.section_1.vsp_mesh.inner_radius  = 4.
    wing_segments.section_1.vsp_mesh.outer_radius  = 4.
    wing_segments.section_1.vsp_mesh.inner_length  = .14
    wing_segments.section_1.vsp_mesh.outer_length  = .14
    
    wing_segments.section_2.vsp_mesh = Data()
    wing_segments.section_2.vsp_mesh.inner_radius  = 4.
    wing_segments.section_2.vsp_mesh.outer_radius  = 4.
    wing_segments.section_2.vsp_mesh.inner_length  = .14
    wing_segments.section_2.vsp_mesh.outer_length  = .14
    
    wing_segments.section_3.vsp_mesh = Data()
    wing_segments.section_3.vsp_mesh.inner_radius  = 4.
    wing_segments.section_3.vsp_mesh.outer_radius  = 4.
    wing_segments.section_3.vsp_mesh.inner_length  = .14
    wing_segments.section_3.vsp_mesh.outer_length  = .14
    
    wing_segments.section_4.vsp_mesh = Data()
    wing_segments.section_4.vsp_mesh.inner_radius  = 4.
    wing_segments.section_4.vsp_mesh.outer_radius  = 2.8
    wing_segments.section_4.vsp_mesh.inner_length  = .14
    wing_segments.section_4.vsp_mesh.outer_length  = .14      
    '''
    analyses.append(aerodynamics)
    '''
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
    '''
    # done!
    return analyses    

# ----------------------------------------------------------------------
#   Define the Vehicle
# ----------------------------------------------------------------------

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
    # base_config.tag = 'Airbus_A220-100'  # NILS
    configs.append(base_config)
    
    write(vehicle,base_config.tag) 

    return configs

# ----------------------------------------------------------------------
#   Plot Mission
# ----------------------------------------------------------------------

def plot_mission(results,line_style='bo-'):

    # Plot Aerodynamic Forces 
    plot_aerodynamic_forces(results, line_style)

    # Plot Aerodynamic Coefficients 
    plot_aerodynamic_coefficients(results, line_style)

    # Drag Components
    plot_drag_components(results, line_style)

    # Plot Altitude, sfc, vehicle weight 
    plot_altitude_sfc_weight(results, line_style)

    # Plot Velocities 
    plot_aircraft_velocities(results, line_style)           
        
    return

def simple_sizing(configs):

    base = configs.base
    base.pull_base()

    # zero fuel weight
    base.mass_properties.max_zero_fuel = 0.9 * base.mass_properties.max_takeoff 

    # Areas
    wetted_areas = get_vsp_measurements(base.tag)

    for wing in base.wings:
        wing.areas.wetted   = wetted_areas[wing.tag]
        wing.areas.exposed  = wetted_areas[wing.tag]
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
    id(su2_object.geometry._base)
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
    
    #configs, analyses = full_setup()
    #su2_object = analyses.configs.base.aerodynamics  # NILS: DIFFERENT FROM run_suave_avl_wrapper_nils.py!
    
    # print('analyses.configs.base =', analyses.configs.base.aerodynamics.keys())
    # sys.exit('Stop here.')
    
    # NOTE: `._base` required to avoid `AttributeError: 'NoneType' object has no attribute 'items'`
    # To understand why I have to use `aircraft = su2_object.geometry._base` instead of
    # `aircraft = SUAVE.Vehicle()`, compare `print(id(aircraft))` with `print(id(vehicle))`
    # in `full_setup()` above - they are not the same object instance!
    
    aircraft_is_b737 = False  # False
    keep_b737_defaults = False  # True
    wing_has_segments = True  # False
    htail_has_segments = True  # False
    ac_segment = "narrowbody"  # "narrowbody" or "regional"
    
    if not keep_b737_defaults:
        """
        NOTE: configs.base is a Config, not a Vehicle -> see configs_setup()
        """
        
        #assert id(configs.base) == id(su2_object.geometry)
        #id_old = id(configs.base)
        
        # Define clean-slate aircraft

        aircraft = SUAVE.Vehicle()   # new clean vehicle
        aircraft_base = SUAVE.Vehicle()   # base version for diffing
    
        # SUAVE requirement: `_base` holds the baseline configuration
        aircraft._base = aircraft_base
    
        # Attach this clean-slate "aircraft" to configs AND su2_object
        #configs.base = aircraft
        #su2_object.geometry = aircraft
        
        #id_new = id(configs.base)
        #assert id_new != id_old
        #assert id(configs.base) == id(su2_object.geometry)
        
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
                df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_narrowbody_kerosene_211225.csv')
                # df_geom = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_narrowbody_kerosene_{date_str}.csv')
                df_mass = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_narrowbody_kerosene_211225_6.csv')
                # df_mass = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_narrowbody_kerosene_{date_str}_{study_idx}.csv')
                aircraft.tag = 'Airbus_A220-100'
                t_tail_bool = False
            elif ac_segment == "regional":
                df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_kerosene_211225.csv')
                # df_geom = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_kerosene_{date_str}.csv')
                df_mass = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_regional_kerosene_211225_6.csv')
                # df_mass = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_regional_kerosene_{date_str}_{study_idx}.csv')
                aircraft.tag = 'ATR_72-600'
                t_tail_bool = True
            
        # sys.exit('Stop here.')
            
        # Aircraft mass properties
        
        for counter, row in df_mass.iterrows():
            print('counter =', counter)
            
            # =============================================================================
            #configs, analyses = full_setup()
            #su2_object = analyses.configs.base.aerodynamics  # NILS: DIFFERENT FROM run_suave_avl_wrapper_nils.py!
            
            # Define clean-slate aircraft

            aircraft = SUAVE.Vehicle()   # new clean vehicle
            aircraft_base = SUAVE.Vehicle()   # base version for diffing
        
            # SUAVE requirement: `_base` holds the baseline configuration
            aircraft._base = aircraft_base
        
            # Attach this clean-slate "aircraft" to configs AND su2_object
            #configs.base = aircraft
            #su2_object.geometry = aircraft
            
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
            
            ## for suave_body in aircraft.fuselages:  # write_geometry.py
            #fuselage.lengths.total = df_geom['body_lengths_total'][0]
            #fuselage.lengths.nose = df_geom['body_lengths_nose'][0]
            #fuselage.lengths.tail = df_geom['body_lengths_tail'][0]
            #fuselage.width = df_geom['body_widths_maximum'][0]
            #fuselage.heights.maximum = df_geom['body_heights_maximum'][0]
            
            ### @NILS: THEN RUN SINGLE OPERATING POINT WITH 2000 ITERATIONS, THEN IMPLEMENT
            # DISTRIBUTED PROPULSION BASED ON TASOPT GEOMETRY AND PROP LOADING!
            
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
            fuselage.lengths.fore_space    = 0.0
            fuselage.lengths.aft_space     = 0.0
            #fuselage.lengths.fore_space    = 6.    * Units.meter
            #fuselage.lengths.aft_space     = 5.    * Units.meter
            fuselage.width                 = df_geom['body_widths_maximum'][0]
            fuselage.heights.maximum       = df_geom['body_heights_maximum'][0]
            fuselage.effective_diameter    = df_geom['body_heights_maximum'][0]
            #fuselage.areas.side_projected  = 142.1948 * Units['meters**2'] 
            #fuselage.areas.wetted          = 446.718  * Units['meters**2'] 
            #fuselage.areas.front_projected = 12.57    * Units['meters**2'] 
            #fuselage.differential_pressure = 5.0e4 * Units.pascal # Maximum differential pressure
            
            fuselage.heights.at_quarter_length          = df_geom['body_heights_maximum'][0]
            fuselage.heights.at_three_quarters_length   = df_geom['body_heights_maximum'][0] * 3.65 / 3.74  # NILS: might have to be reduced in accordance with 
            fuselage.heights.at_wing_root_quarter_chord = df_geom['body_heights_maximum'][0]
            
            ###
            
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
                    
                    center_airfoil                          = SUAVE.Components.Airfoils.Airfoil()
                    center_airfoil.coordinate_file          = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737a.txt'
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
                    inboard_airfoil.coordinate_file       = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737b.txt'
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
                    outboard_airfoil.coordinate_file      = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737c.txt'
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
                    tip_airfoil.coordinate_file      = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737d.txt'
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
                '''
                # control surfaces -------------------------------------------
                slat                          = SUAVE.Components.Wings.Control_Surfaces.Slat()
                slat.tag                      = 'slat'
                slat.span_fraction_start      = df_geom['wing_slat_span_fraction_start'][0]
                slat.span_fraction_end        = df_geom['wing_slat_span_fraction_end'][0]
                slat.deflection               = df_geom['wing_slat_deflection'][0]
                slat.chord_fraction           = df_geom['wing_slat_chord_fraction'][0]
                wing.append_control_surface(slat)
    
                flap                          = SUAVE.Components.Wings.Control_Surfaces.Flap()
                flap.tag                      = 'flap'
                flap.span_fraction_start      = df_geom['wing_flap_span_fraction_start'][0]
                flap.span_fraction_end        = df_geom['wing_flap_span_fraction_end'][0]
                flap.deflection               = df_geom['wing_flap_deflection'][0]
                flap.configuration_type       = df_geom['wing_flap_configuration_type'][0]
                flap.chord_fraction           = df_geom['wing_flap_chord_fraction'][0]
                wing.append_control_surface(flap)
    
                aileron                       = SUAVE.Components.Wings.Control_Surfaces.Aileron()
                aileron.tag                   = 'aileron'
                aileron.span_fraction_start   = df_geom['wing_aileron_span_fraction_start'][0]
                aileron.span_fraction_end     = df_geom['wing_aileron_span_fraction_end'][0]
                aileron.deflection            = df_geom['wing_aileron_deflection'][0]
                aileron.chord_fraction        = df_geom['wing_aileron_chord_fraction'][0]
                wing.append_control_surface(aileron)
                '''
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
                '''
                # control surfaces -------------------------------------------
                elevator                       = SUAVE.Components.Wings.Control_Surfaces.Elevator()
                elevator.tag                   = 'elevator'
                elevator.span_fraction_start   = df_geom['htail_elevator_span_fraction_start'][0]
                elevator.span_fraction_end     = df_geom['htail_elevator_span_fraction_end'][0]
                elevator.deflection            = df_geom['htail_elevator_deflection'][0]
                elevator.chord_fraction        = df_geom['htail_elevator_chord_fraction'][0]
                wing.append_control_surface(elevator)
                '''
                
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
            nacelle.origin = np.expand_dims(df_geom['nacelle_origin'].to_numpy(), axis=0)  # NILS: added np.expand_dims() to avoid "IndexError: invalid index to scalar variable." in trunk\SUAVE\Input_Output\OpenVSP\vsp_nacelle.py
            nacelle.flow_through = df_geom['nacelle_flow_through'][0]
            nacelle_airfoil = SUAVE.Components.Airfoils.Airfoil() 
            nacelle_airfoil.naca_4_series_airfoil = '2410'
            nacelle.append_airfoil(nacelle_airfoil)
            
            # =============================================================================
            # n_segments = 18
            
            # # Wing Segments
            # for i_segs in range(n_segments):
            #     root_airfoil                          = SUAVE.Components.Airfoils.Airfoil()
            #     root_airfoil.coordinate_file          = r'/Documents/GitHub/SUAVE/regression/scripts/Vehicles/Airfoils/B737a.txt'
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
            
            nacelle_2 = deepcopy(nacelle)
            nacelle_2.tag = 'nacelle_2'
            nacelle_2_origin = deepcopy(nacelle.origin)
            nacelle_2_origin[0][1] *= -1  # NILS: added [0] relative to run_suave_avl_wrapper_nils.py
            nacelle_2.origin = nacelle_2_origin
            
            aircraft.append_component(nacelle)  
            aircraft.append_component(nacelle_2)
            # =============================================================================
            
            # print(aircraft.nacelles.nacelle_1.Airfoil)
            # sys.exit('Stop.')
            
            # =============================================================================
            # su2_object.sample_training(
            #     study_idx=study_idx, counter=counter,
            #     sigma_fcs = sigma_fcs,
            #     span_loc = span_loc,
            #     fcs_loc = fcs_loc,
            #     wing_frac = wing_frac,
            #     nacelle_frac = nacelle_frac,
            # )
            # su2_object.initialize()  # NILS: replaced above line from AVL analysis for SU2 analysis
            main(aircraft)  # NILS
            # =============================================================================
            
        # elif keep_b737_defaults:
            
        #     # aircraft = su2_object.geometry._base  # instance of SUAVE.Vehicle()
        #     # print('id(aircraft) =', id(aircraft))
        #     # # OR
        #     # aircraft = configs.base._base  # instance of SUAVE.Vehicle()
        #     # print('id(aircraft) =', id(aircraft))
            
        #     pass
            
        # sys.exit('Stop.')
        
        #%% Run sample_training() only
        
        # # su2_object.sample_training()
        # # =============================================================================
        # su2_object.sample_training(study_idx=study_idx, counter=counter)
        # # =============================================================================

