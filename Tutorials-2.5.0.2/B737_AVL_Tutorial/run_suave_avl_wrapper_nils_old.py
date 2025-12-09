# SUAVE imports
import SUAVE
from SUAVE.Core import Units, Data
from SUAVE.Core import redirect

from SUAVE.Analyses.Mission.Segments.Conditions.Aerodynamics import Aerodynamics
from SUAVE.Analyses.Mission.Segments.Conditions.Conditions   import Conditions

from SUAVE.Methods.Aerodynamics.AVL.write_geometry           import write_geometry
from SUAVE.Methods.Aerodynamics.AVL.write_mass_file          import write_mass_file
from SUAVE.Methods.Aerodynamics.AVL.write_run_cases          import write_run_cases
from SUAVE.Methods.Aerodynamics.AVL.write_input_deck         import write_input_deck
from SUAVE.Methods.Aerodynamics.AVL.run_analysis             import run_analysis
from SUAVE.Methods.Aerodynamics.AVL.translate_data           import translate_conditions_to_cases, translate_results_to_conditions
from SUAVE.Methods.Aerodynamics.AVL.purge_files              import purge_files
from SUAVE.Methods.Aerodynamics.AVL.Data.Settings            import Settings
from SUAVE.Methods.Aerodynamics.AVL.Data.Cases               import Run_Case
from SUAVE.Methods.Geometry.Two_Dimensional.Planform.populate_control_sections import populate_control_sections  
from SUAVE.Methods.Flight_Dynamics.Dynamic_Stability.compute_dynamic_flight_modes import  compute_dynamic_flight_modes
from SUAVE.Components.Wings.Control_Surfaces import Aileron , Elevator , Slat , Flap , Rudder 

# local imports 
from SUAVE.Analyses.Stability import Stability

# Package imports 
import os
import numpy as np
import sys  
from shutil import rmtree 
from scipy.interpolate import  RectBivariateSpline 


## @ingroup Analyses-Stability
class AVL(Stability):
    """This builds a surrogate and computes moment using AVL.

    Assumptions:
    None

    Source:
    None
    """  

    def __defaults__(self):
        """This sets the default values and methods for the analysis.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        N/A
        
        """
        self.tag                                    = 'avl' 
        
        self.current_status                         = Data()        
        self.current_status.batch_index             = 0
        self.current_status.batch_file              = None
        self.current_status.deck_file               = None
        self.current_status.cases                   = None      
        self.geometry                               = None   
                                                    
        self.settings                               = Settings()
        self.settings.filenames.log_filename        = sys.stdout
        self.settings.filenames.err_filename        = sys.stderr        
        self.settings.number_spanwise_vortices      = 20
        self.settings.number_chordwise_vortices     = 10
        self.settings.trim_aircraft                 = False 
                                                    
        # Conditions table, used for surrogate model training
        self.training                               = Data()   
        
        # Standard subsonic/transonic aircarft
        self.training.angle_of_attack               = np.array([-2.,0., 2.,5., 7., 10.])*Units.degrees
        self.training.Mach                          = np.array([0.05,0.15,0.25, 0.45,0.65,0.85]) 
        
        self.training.moment_coefficient            = None
        self.training.Cm_alpha_moment_coefficient   = None
        self.training.Cn_beta_moment_coefficient    = None
        self.training.neutral_point                 = None
        self.training_file                          = None
        
        # Initialize quantities
        self.configuration                          = Data()    
        self.geometry                               = Data()
                                                    
        # Regression Status      
        self.keep_files                             = False
        self.save_regression_results                = True  # False  # NILS: save these for visualisation
        self.regression_flag                        = False 


    def sample_training(self):
    
        # Unpack
        run_folder    = os.path.abspath(self.settings.filenames.run_folder)
        geometry      = self.geometry
        training      = self.training 
        trim_aircraft = self.settings.trim_aircraft  
        AoA           = training.angle_of_attack
        Mach          = training.Mach
        atmosphere    = SUAVE.Analyses.Atmospheric.US_Standard_1976()
        atmo_data     = atmosphere.compute_values(altitude = 0.0)         
        cg            = geometry.mass_properties.center_of_gravity[0][0]  # copied over from SUAVE 2.5.2
        MAC           = geometry.wings.main_wing.chords.mean_aerodynamic  # copied over from SUAVE 2.5.2
                      
        CM            = np.zeros((len(AoA),len(Mach)))
        Cm_alpha      = np.zeros_like(CM)
        Cn_beta       = np.zeros_like(CM)
        NP            = np.zeros_like(CM)
        static_margin = np.zeros_like(CM)  # copied over from SUAVE 2.5.2
       
        # remove old files in run directory  
        if os.path.exists('avl_files'):
            if not self.regression_flag:
                rmtree(run_folder)
                
        results_list = []  # NILS: added to store dynamic stability analysis results
        for i,_ in enumerate(Mach):
            # Set training conditions
            run_conditions = Aerodynamics()
            run_conditions.freestream.density           = atmo_data.density[0,0] 
            run_conditions.freestream.gravity           = 9.81               
            run_conditions.aerodynamics.angle_of_attack = AoA 
            run_conditions.freestream.speed_of_sound    = atmo_data.speed_of_sound[0,0] 
            run_conditions.aerodynamics.side_slip_angle = 0.0
            run_conditions.freestream.velocity          = Mach[i] * run_conditions.freestream.speed_of_sound
            run_conditions.freestream.mach_number       = Mach[i] 
            
            #Run Analysis at AoA[i] and Mach[i]
            results =  self.evaluate_conditions(run_conditions, trim_aircraft)
            results_list.append(results)  # NILS: added to store dynamic stability analysis results
    
            # Obtain CM Cm_alpha, Cn_beta and the Neutral Point 
            CM[:,i]       = results.aerodynamics.Cmtot[:,0]
            Cm_alpha[:,i] = results.stability.static.Cm_alpha[:,0]
            Cn_beta[:,i]  = results.stability.static.Cn_beta[:,0]
            NP[:,i]       = results.stability.static.neutral_point[:,0]
            static_margin[:,i] = (results.stability.static.neutral_point[:,0] - cg) / MAC  # copied over from SUAVE 2.5.2
        
        # Save the data for regression 
        # convert from 2D to 1D
        CM_1D       = CM.reshape([len(AoA)*len(Mach),1]) 
        Cm_alpha_1D = Cm_alpha.reshape([len(AoA)*len(Mach),1])  
        Cn_beta_1D  = Cn_beta.reshape([len(AoA)*len(Mach),1])         
        NP_1D       = NP.reshape([len(AoA)*len(Mach),1])  # NILS: TYPO - said `Cn_beta` like in line above
        static_margin_1D = static_margin.reshape([len(AoA)*len(Mach),1])  # copied over from SUAVE 2.5.2
        np.savetxt(geometry.tag+'_stability_data.txt',np.hstack([CM_1D,Cm_alpha_1D, Cn_beta_1D,NP_1D,static_margin_1D]),fmt='%10.8f',header='   CM       Cm_alpha       Cn_beta       NP       static_margin ')
            
        # NILS: save dynamic stability results too
            
        # Loop over dynamic stability analysis results across all Mach numbers
        final_out_list = []
        for results in results_list:
        
            # Extract SUAVE dynamic stability results
            dyn = results.dynamic_stability
            long = dyn.LongModes
            lat = dyn.LatModes
        
            # ---- Collect all fields (robust to scalars, 1D arrays, lists, eigenvalues) ---- #
            fields = [
                long.phugoidFreqHz,
                long.phugoidDamp,
                long.phugoidTimeDoubleHalf,
                long.shortPeriodFreqHz,
                long.shortPeriodDamp,
                long.shortPeriodTimeDoubleHalf,
        
                lat.dutchRollFreqHz,
                lat.dutchRollDamping,
                lat.dutchRollTimeDoubleHalf,
                lat.dutchRoll_mode_real,
                lat.rollSubsistenceFreqHz,
                lat.rollSubsistenceTimeConstant,
                lat.rollSubsistenceDamping,
                lat.spiralFreqHz,
                lat.spiralTimeDoubleHalf,
                lat.spiralDamping,
        
                dyn.pMax,
            ]
        
            # Convert each to a 1D numpy array
            cols = [np.atleast_1d(np.array(f)) for f in fields]
        
            # Determine max length
            L = max(len(c) for c in cols)
        
            # Pad columns to equal length with NaNs
            cols_padded = [
                np.pad(c, (0, L - len(c)), mode='constant', constant_values=np.nan)
                for c in cols
            ]
        
            # Build final 2D array (L rows × 17 columns)
            data_out = np.column_stack(cols_padded)
        
            # ----------------------- Write to text file ----------------------------- #
            header = (
                "phugoidFreqHz  phugoidDamp  phugoidTimeDoubleHalf  "
                "shortPeriodFreqHz  shortPeriodDamp  shortPeriodTimeDoubleHalf  "
                "dutchRollFreqHz  dutchRollDamping  dutchRollTimeDoubleHalf  dutchRoll_mode_real  "
                "rollSubsistenceFreqHz  rollSubsistenceTimeConstant  rollSubsistenceDamping  "
                "spiralFreqHz  spiralTimeDoubleHalf  spiralDamping  "
                "pMax  "
                "Long_Re1  Long_Im1  Long_Re2  Long_Im2  Long_Re3  Long_Im3  Long_Re4  Long_Im4  "
                "Lat_Re1   Lat_Im1   Lat_Re2   Lat_Im2   Lat_Re3   Lat_Im3   Lat_Re4   Lat_Im4"
            )
            # ---------- Save LongModes and LatModes in same output file ---------- #

            # Convert complex matrices to real+imag pairs
            def split_complex_matrix(M):
                M = np.asarray(M)
                real = M.real
                imag = M.imag
                # stack as [Re1 Im1 Re2 Im2 ...]
                cols = []
                for k in range(M.shape[1]):
                    cols.append(real[:, k])
                    cols.append(imag[:, k])
                return np.column_stack(cols)
            
            # Longitudinal modes (Nx4 complex → Nx8 real)
            long_modes_mat = split_complex_matrix(long.LongModes)
            
            # Lateral-directional modes (Nx4 complex → Nx8 real)
            lat_modes_mat = split_complex_matrix(lat.LatModes)
            
            # Pad dynamic scalar/vector data to same number of rows as modes
            n_long = long_modes_mat.shape[0]
            n_lat = lat_modes_mat.shape[0]
            nrows = max(n_long, n_lat, data_out.shape[0])
            
            def pad_rows(A, n):
                if A.shape[0] == n:
                    return A
                pad = np.full((n - A.shape[0], A.shape[1]), np.nan)
                return np.vstack((A, pad))
            
            data_out_padded = pad_rows(data_out, nrows)
            long_modes_padded = pad_rows(long_modes_mat, nrows)
            lat_modes_padded = pad_rows(lat_modes_mat, nrows)
            
            # Final array: [dynamic stability scalars | LongModes | LatModes]
            final_out = np.hstack([data_out_padded, long_modes_padded, lat_modes_padded])
            final_out_list.append(final_out)
            
        # Stack dynamic stability analysis results across all Mach numbers
        final_out_array = np.vstack(final_out_list)
            
        np.savetxt(
            geometry.tag + "_dynamic_stability_data.txt",
            final_out_array,  # data_out,
            fmt="%12.6f",
            header=header,
            comments="",  # prevents '#' from being added
        )
        
        return
    
    
    def evaluate_conditions(self,run_conditions, trim_aircraft  ):
        """Process vehicle to setup geometry, condititon, and configuration.
    
        Assumptions:
        None
    
        Source:
        N/A
    
        Inputs:
        run_conditions <SUAVE data type> aerodynamic conditions; until input
                method is finalized, will assume mass_properties are always as 
                defined in self.features
    
        Outputs:
        results        <SUAVE data type>
    
        Properties Used:
        self.settings.filenames.
          run_folder
          output_template
          batch_template
          deck_template
        self.current_status.
          batch_index
          batch_file
          deck_file
          cases
        """
        
        # unpack
        run_folder                       = os.path.abspath(self.settings.filenames.run_folder)
        run_script_path                  = run_folder.rstrip('avl_files').rstrip('/')
        aero_results_template_1          = self.settings.filenames.aero_output_template_1       # 'stability_axis_derivatives_{}.dat' 
        aero_results_template_2          = self.settings.filenames.aero_output_template_2       # 'surface_forces_{}.dat'
        aero_results_template_3          = self.settings.filenames.aero_output_template_3       # 'strip_forces_{}.dat'   
        aero_results_template_4          = self.settings.filenames.aero_output_template_4       # 'body_axis_derivatives_{}.dat'     
        dynamic_results_template_1       = self.settings.filenames.dynamic_output_template_1    # 'eigen_mode_{}.dat'
        dynamic_results_template_2       = self.settings.filenames.dynamic_output_template_2    # 'system_matrix_{}.dat'
        batch_template                   = self.settings.filenames.batch_template
        deck_template                    = self.settings.filenames.deck_template 
        
        # rename defaul avl aircraft tag
        self.tag                         = 'avl_analysis_of_{}'.format(self.geometry.tag) 
        print('self.geometry._base =', self.geometry._base)
        import sys
        sys.exit()
        self.settings.filenames.features = self.geometry._base.tag + '.avl'
        self.settings.filenames.mass_file= self.geometry._base.tag + '.mass'
        
        # update current status
        self.current_status.batch_index += 1
        batch_index                      = self.current_status.batch_index
        self.current_status.batch_file   = batch_template.format(batch_index)
        self.current_status.deck_file    = deck_template.format(batch_index)
               
        # control surfaces
        num_cs       = 0
        cs_names     = []
        cs_functions = []
        for wing in self.geometry.wings: # this parses through the wings to determine how many control surfaces does the vehicle have 
            if wing.control_surfaces:
                wing = populate_control_sections(wing)     
                num_cs_on_wing = len(wing.control_surfaces)
                num_cs +=  num_cs_on_wing
                for ctrl_surf in wing.control_surfaces:
                    cs_names.append(ctrl_surf.tag)  
                    if (type(ctrl_surf) ==  Slat):
                        ctrl_surf_function  = 'slat'
                    elif (type(ctrl_surf) ==  Flap):
                        ctrl_surf_function  = 'flap' 
                    elif (type(ctrl_surf) ==  Aileron):
                        ctrl_surf_function  = 'aileron'                          
                    elif (type(ctrl_surf) ==  Elevator):
                        ctrl_surf_function  = 'elevator' 
                    elif (type(ctrl_surf) ==  Rudder):
                        ctrl_surf_function = 'rudder'                      
                    cs_functions.append(ctrl_surf_function)   
        
        # translate conditions
        cases                            = translate_conditions_to_cases(self, run_conditions)    
        for case in cases:
            case.stability_and_control.number_control_surfaces = num_cs
            case.stability_and_control.control_surface_names   = cs_names
        self.current_status.cases        = cases  
        
       # write casefile names using the templates defined in MACE/Analyses/AVL/AVL_Data_Classes/Settings.py 
        for case in cases:  
            case.aero_result_filename_1     = aero_results_template_1.format(case.tag)      # 'stability_axis_derivatives_{}.dat'  
            case.aero_result_filename_2     = aero_results_template_2.format(case.tag)      # 'surface_forces_{}.dat'
            case.aero_result_filename_3     = aero_results_template_3.format(case.tag)      # 'strip_forces_{}.dat'  
            case.aero_result_filename_4     = aero_results_template_4.format(case.tag)      # 'body_axis_derivatives_{}.dat'
            case.eigen_result_filename_1    = dynamic_results_template_1.format(case.tag)   # 'eigen_mode_{}.dat'
            case.eigen_result_filename_2    = dynamic_results_template_2.format(case.tag)   # 'system_matrix_{}.dat'
        
        # write the input files
        with redirect.folder(run_folder,force=False):
            write_geometry(self,run_script_path)
            write_mass_file(self,run_conditions)
            write_run_cases(self,trim_aircraft)
            write_input_deck(self, trim_aircraft)
    
            # RUN AVL!
            results_avl = run_analysis(self)
    
        # translate results
        results = translate_results_to_conditions(cases,results_avl)
        
        # -----------------------------------------------------------------------------------------------------------------------                     
        # Dynamic Stability & System Matrix Computation
        # -----------------------------------------------------------------------------------------------------------------------      
        # Dynamic Stability
        if np.count_nonzero(self.geometry.mass_properties.moments_of_inertia.tensor) > 0:  
            results = compute_dynamic_flight_modes(results,self.geometry,run_conditions,cases)        
        
        if not self.keep_files:
            rmtree( run_folder )   
        
        # sys.exit('Stop here.')
        
        return results
    
#%%
    
if __name__ == '__main__':
    
    avl_object = AVL()
    
    import pandas as pd
    
    df = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl\suave_avl_wrapper_tasopt_inputs.csv')
    
    avl_object.geometry = SUAVE.Vehicle()  # tut_mission_B737_AVL.py
    aircraft = avl_object.geometry  # correspons to `vehicle`
    
    def configs_setup(vehicle):
        
        # ------------------------------------------------------------------
        #   Initialize Configurations
        # ------------------------------------------------------------------
        configs = SUAVE.Components.Configs.Config.Container()

        base_config = SUAVE.Components.Configs.Config(vehicle)
        base_config.tag = 'base'
        configs.append(base_config)

        return configs
    
    configs = configs_setup(aircraft)
    
    
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
    
    
    def analyses_setup(configs):

        analyses = SUAVE.Analyses.Analysis.Container()

        # build a base analysis for each config
        for tag,config in configs.items():
            analysis = base_analysis(config)
            analyses[tag] = analysis

        return analyses
    
    
    configs_analyses = analyses_setup(configs)
    analyses = SUAVE.Analyses.Analysis.Container()
    analyses.configs  = configs_analyses
    
    aircraft.mass_properties.center_of_gravity[0][0] = df['x_cg']
    aircraft.mass_properties.center_of_gravity[0][1] = df['y_cg']
    aircraft.mass_properties.center_of_gravity[0][2] = df['z_cg']
    aircraft.mass_properties.mass = df['mass']
    moments_of_inertia = aircraft.mass_properties.moments_of_inertia.tensor
    moments_of_inertia[0][0] = df['Ixx']
    moments_of_inertia[1][1] = df['Iyy']
    moments_of_inertia[2][2] = df['Izz']
    moments_of_inertia[0][1] = df['Ixy']
    moments_of_inertia[1][2] = df['Iyz']
    moments_of_inertia[2][0] = df['Izx']
    
    fuselage = SUAVE.Components.Fuselages.Fuselage()  # tut_mission_B737_AVL.py
    fuselage.tag = 'fuselage'
    aircraft.append_component(fuselage)
    
    wing = SUAVE.Components.Wings.Main_Wing()
    wing.tag = 'main_wing'
    aircraft.append_component(wing)
    
    for suave_body in aircraft.fuselages:  # write_geometry.py
        suave_body.lengths.total = df['body_lengths_total']
        suave_body.lengths.nose = df['body_lengths_nose']
        suave_body.lengths.tail = df['body_lengths_tail']
        suave_body.width = df['body_widths_maximum']
        suave_body.heights.maximum = df['body_heights_maximum']
    
    for suave_wing in aircraft.wings:  # write_geometry.py
        suave_wing.spans.projected = df['wing_spans_projected']
        suave_wing.origin = df['wing_origin']
        suave_wing.dihedral = df['wing_dihedral']
        
        for i in range(3):
            suave_wing.Segments.append(SUAVE.Components.Wings.Segment())
        center_segment = suave_wing.Segments[0]
        inboard_segment = suave_wing.Segments[1]
        outboard_segment = suave_wing.Segments[2]
                
        center_segment.sweeps.leading_edge = df['wing_segment_sweep_leading_edge_center']
        center_segment.root_chord_percent = df['wing_segment_root_chord_percent_center']
        center_segment.percent_span_location = df['wing_segment_percent_span_location_center']
        center_segment.sweeps.quarter_chord = df['wing_segment_sweep_quarter_chord_center']
        # suave_wing.Segment.twist = df['wing_segment_twist_center']
        center_segment.twist = df['wing_segment_twist_center']
        
        inboard_segment.sweeps.leading_edge = df['wing_segment_sweep_leading_edge_inboard']
        inboard_segment.root_chord_percent = df['wing_segment_root_chord_percent_inboard']
        inboard_segment.percent_span_location = df['wing_segment_percent_span_location_inboard']
        inboard_segment.sweeps.quarter_chord = df['wing_segment_sweep_quarter_chord_inboard']
        # suave_wing.Segment.twist = df['wing_segment_twist_inboard']
        inboard_segment.twist = df['wing_segment_twist_inboard']
        
        outboard_segment.sweeps.leading_edge = df['wing_segment_sweep_leading_edge_outboard']
        outboard_segment.root_chord_percent = df['wing_segment_root_chord_percent_outboard']
        outboard_segment.percent_span_location = df['wing_segment_percent_span_location_outboard']
        outboard_segment.sweeps.quarter_chord = df['wing_segment_sweep_quarter_chord_outboard']
        # suave_wing.Segment.twist = df['wing_segment_twist_outboard']
        outboard_segment.twist = df['wing_segment_twist_outboard']
    
    #%%    
    
    avl_object.sample_training()



    
    
    
    