## @ingroup Analyses-Stability
# AVL.py
#
# Created:  Apr 2017, M. Clarke 
# Modified: Apr 2020, M. Clarke

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

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
from .Stability import Stability

# Package imports 
import os
import numpy as np
import sys  
from shutil import rmtree 
from scipy.interpolate import  RectBivariateSpline 

# ----------------------------------------------------------------------
#  Class
# ----------------------------------------------------------------------

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
        self.settings.trim_aircraft                 = True  # NILS: toggle to trim or not
        self.settings.print_output                  = True  # NILS: added to match SUAVE 2.5.2 (toggle to print or not AVL console output)
        
        # Regression Status      
        self.settings.keep_files                    = True  # NILS: added to match SUAVE 2.5.2
        self.settings.save_regression_results       = False  # NILS: added to match SUAVE 2.5.2
        self.settings.regression_flag               = False  # NILS: added to match SUAVE 2.5.2
                                                    
        # Conditions table, used for surrogate model training
        self.training                               = Data()   
        
        # # Standard subsonic/transonic aircarft
        # # self.training.angle_of_attack               = np.array([-2.,0., 2.,5., 7., 10.])*Units.degrees  # NILS: default
        # # self.training.angle_of_attack = np.array([-2.0, 0.0, 5.0]) * Units.degrees  # NILS
        # self.training.angle_of_attack = np.linspace(-2.0, 10.0, 4) * Units.degrees  # NILS
        # # self.training.Mach                          = np.array([0.05,0.15,0.25, 0.45,0.65,0.85])   # NILS: default
        # # self.training.Mach = np.array([0.05, 0.45, 0.85])  # NILS
        # self.training.Mach = np.linspace(0.05, 0.85, 4)  # NILS
        
        # # NILS: added training parameters
        # # self.training.side_slip_angle = np.array([-10.0, -5.0, 0.0, 5.0, 10.0]) * Units.degrees
        # # self.training.side_slip_angle = np.array([-5.0, 0.0, 10.0]) * Units.degrees
        # self.training.side_slip_angle = np.linspace(-10.0, 10.0, 4) * Units.degrees
        # # self.training.altitude = np.linspace(0, 11e3, 5)
        # self.training.altitude = np.linspace(0, 11e3, 4)
        
        # NILS: longitudinal test cases from Table 6.2 in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris
        # NOTE: training inputs always have to be at least 1D, otherwise get
        # `TypeError: object of type 'float' has no len()` in
        # Documents\GitHub\SUAVE\trunk\SUAVE\Methods\Aerodynamics\AVL\translate_data.py
        self.training.Mach = np.array([0.2, 0.5, 0.7, 0.7, 0.5, 0.2])
        self.training.altitude = np.array([0.0, 10e3, 35e3, 35e3, 10e3, 0.0]) * 0.3048
        self.training.mass = None
        self.training.side_slip_angle = np.zeros_like(self.training.Mach) * Units.degrees
        # NILS: 6x faster if use `np.array([0])` instead of `np.zeros_like(self.training.Mach)` (6x duplication)
        self.training.angle_of_attack = np.array([0]) * Units.degrees  #  np.zeros_like(self.training.Mach) * Units.degrees  # to be trimmed
        self.training.load_factor = np.array([1.0, 2.5, 1.0, 1.0, 2.5, 1.0])
        
        self.settings.side_slip_angle               = 0.0  # NILS: added to match SUAVE 2.5.2 (can remain set to 0 as vary self.training.side_slip_angle in sample_training() below)
        self.settings.roll_rate_coefficient         = 0.0  # NILS: added to match SUAVE 2.5.2
        self.settings.pitch_rate_coefficient        = 0.0  # NILS: added to match SUAVE 2.5.2
        self.settings.lift_coefficient              = None  # NILS: added to match SUAVE 2.5.2 (computed in sample_training)
        self.settings.load_factor = None  # NILS: added for variability in Documents\GitHub\SUAVE\trunk\SUAVE\Methods\Aerodynamics\AVL\write_run_cases.py
        
        self.training.moment_coefficient            = None
        self.training.Cm_alpha_moment_coefficient   = None
        self.training.Cn_beta_moment_coefficient    = None
        self.training.neutral_point                 = None
        self.training_file                          = None
                                                    
        # Surrogate model
        self.surrogates                             = Data()
        self.surrogates.moment_coefficient          = None
        self.surrogates.Cm_alpha_moment_coefficient = None
        self.surrogates.Cn_beta_moment_coefficient  = None      
        self.surrogates.neutral_point               = None
    
        # Initialize quantities
        self.configuration                          = Data()    
        self.geometry                               = Data()
                                                    
        # Regression Status      
        self.keep_files                             = True  # NILS: toggle to save these for plotting geometry
        self.save_regression_results                = True  # NILS: toggle to save these for visualisation
        self.regression_flag                        = False 

    def finalize(self):
        """Drives functions to get training samples and build a surrogate.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        self.tag = 'avl_analysis_of_{}'.format( )

        Properties Used:
        self.geometry.tag
        """
        geometry                       = self.geometry
        self.tag                       = 'avl_analysis_of_{}'.format(geometry.tag) 
            
        # Sample training data
        self.sample_training()
        
        # Build surrogate
        self.build_surrogate()
    
        return

    def __call__(self,conditions):
        """Evaluates moment coefficient, stability and body axis deriviatives and neutral point using available surrogates.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        state.conditions.
          mach_number      [-]
          angle_of_attack  [radians]

        Outputs:
        results
            results.stability.static
            results.stability.dynamic
        

        Properties Used:
        self.surrogates.
           pitch_moment_coefficient [-] CM
           cm_alpha                 [-] Cm_alpha
           cn_beta                  [-] Cn_beta
           neutral_point            [-] NP

        """
        # Unpack
        surrogates          = self.surrogates       
        Mach                = conditions.freestream.mach_number
        AoA                 = conditions.aerodynamics.angle_of_attack 
        moment_model        = surrogates.moment_coefficient
        Cm_alpha_model      = surrogates.Cm_alpha_moment_coefficient
        Cn_beta_model       = surrogates.Cn_beta_moment_coefficient      
        neutral_point_model = surrogates.neutral_point
        cg                  = self.geometry.mass_properties.center_of_gravity[0][0]  # NILS: added to match SUAVE 2.5.2
        MAC                 = self.geometry.wings.main_wing.chords.mean_aerodynamic  # NILS: added to match SUAVE 2.5.2
        
        # set up data structures
        static_stability    = Data()
        dynamic_stability   = Data()    

        #Run Analysis
        data_len            = len(AoA)
        CM                  = np.zeros([data_len,1])
        Cm_alpha            = np.zeros([data_len,1])
        Cn_beta             = np.zeros([data_len,1])
        NP                  = np.zeros([data_len,1]) 

        for i,_ in enumerate(AoA):           
            CM[i]       = moment_model(AoA[i][0],Mach[i][0])[0]  
            Cm_alpha[i] = Cm_alpha_model(AoA[i][0],Mach[i][0])[0]  
            Cn_beta[i]  = Cn_beta_model(AoA[i][0],Mach[i][0])[0]  
            NP[i]       = neutral_point_model(AoA[i][0],Mach[i][0])[0]    
            
        static_stability.CM            = CM
        static_stability.Cm_alpha      = Cm_alpha 
        static_stability.Cn_beta       = Cn_beta   
        static_stability.neutral_point = NP
        static_stability.static_margin = (NP - cg)/MAC
 
        results         = Data()
        results.static  = static_stability
        results.dynamic = dynamic_stability
    
        return results   


    def sample_training(
        self,
        study_idx=None, counter=None,
        sigma_fcs=None,
        span_loc=None,
        fcs_loc=None,
        wing_frac=None,
        nacelle_frac=None,
    ):  # NILS: added 2nd and 3rd argument to allow running in parallel for different mass study cases and runs within
        """Call methods to run AVL for sample point evaluation.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        see properties used

        Outputs:
        self.training.
          coefficients     [-] CM, Cm_alpha, Cn_beta
          neutral point    [-] NP
          grid_points      [radians,-] angles of attack and Mach numbers 

        Properties Used:
        self.geometry.tag  <string>
        self.training.     
          angle_of_attack  [radians]
          Mach             [-]
        self.training_file (optional - file containing previous AVL data)
        """
        # =============================================================================
        # NILS: set different run folder for each study_idx to
        # prevent parallel processes from overwriting each other
        run_folder = 'avl_files_' + str(study_idx)
        os.makedirs(run_folder, exist_ok=True)
        self.settings.filenames.run_folder = run_folder
        # =============================================================================
        
        # Unpack
        run_folder    = os.path.abspath(self.settings.filenames.run_folder)
        geometry      = self.geometry
        training      = self.training
        trim_aircraft = self.settings.trim_aircraft
        AoA           = training.angle_of_attack
        Mach          = training.Mach
        
        # NILS: added training parameters
        Beta = training.side_slip_angle
        h = training.altitude
        n = training.load_factor
        training.mass = np.array([
            geometry.mass_properties.takeoff, geometry.mass_properties.takeoff, geometry.mass_properties.takeoff,
            geometry.mass_properties.max_zero_fuel, geometry.mass_properties.max_zero_fuel, geometry.mass_properties.max_zero_fuel
        ])
        W = training.mass * 9.81
        
        side_slip_angle        = self.settings.side_slip_angle  # NILS: added to match SUAVE 2.5.2
        roll_rate_coefficient  = self.settings.roll_rate_coefficient  # NILS: added to match SUAVE 2.5.2
        pitch_rate_coefficient = self.settings.pitch_rate_coefficient  # NILS: added to match SUAVE 2.5.2
        # lift_coefficient       = self.settings.lift_coefficient  # NILS: added to match SUAVE 2.5.2 (commented as included in loop below)
        atmosphere    = SUAVE.Analyses.Atmospheric.US_Standard_1976()
        # atmo_data     = atmosphere.compute_values(altitude = 0.0)  # NILS: commented as included in loop below
        cg            = geometry.mass_properties.center_of_gravity[0][0]  # NILS: added to match SUAVE 2.5.2
        MAC           = geometry.wings.main_wing.chords.mean_aerodynamic  # NILS: added to match SUAVE 2.5.2
                      
        # # CM            = np.zeros((len(AoA),len(Mach)))  # NILS: default
        # CM            = np.zeros((len(AoA),len(Mach), len(Beta), len(h)))  # NILS
        # Cm_alpha      = np.zeros_like(CM)
        # Cn_beta       = np.zeros_like(CM)
        # NP            = np.zeros_like(CM)
        # static_margin = np.zeros_like(CM)  # NILS: added to match SUAVE 2.5.2
        # Cl_beta = np.zeros_like(CM)  # NILS: added independently
        # Cn_r = np.zeros_like(CM)  # NILS: added independently
        # Cl_r = np.zeros_like(CM)  # NILS: added independently
        CM            = []  # NILS
        Cm_alpha      = []
        Cn_beta       = []
        NP            = []
        static_margin = []  # NILS: added to match SUAVE 2.5.2
        Cl_beta = []  # NILS: added independently
        Cn_r = []  # NILS: added independently
        Cl_r = []  # NILS: added independently
        
        # remove old files in run directory  
        # if os.path.exists('avl_files'):
        if os.path.exists(run_folder):  # NILS: use run_folder name to allow parallelisation
            if not self.regression_flag:
                rmtree(run_folder)
                
        results_list = []  # NILS: added to store dynamic stability analysis results
        AoA_list = []
        Mach_list = []
        Beta_list = []
        h_list = []
        n_list = []
        W_list = []
        # for i, _Mach in enumerate(Mach):
        #     for j, _Beta in enumerate(Beta):
        #         for k, _h in enumerate(h):
            
        for i, (_Mach, _Beta, _h, _n, _W) in enumerate(zip(Mach, Beta, h, n, W)):
            
            print('_Mach, _Beta, _h, _n, _W =', _Mach, _Beta, _h, _n, _W)
                    
            # atmo_data = atmosphere.compute_values(altitude = h[k])  # NILS: moved here to account for differences in altitude
            atmo_data = atmosphere.compute_values(altitude = _h)  # NILS: moved here to account for differences in altitude
            
            # Set training conditions
            run_conditions = Aerodynamics()
            run_conditions.freestream.density           = atmo_data.density[0,0] 
            run_conditions.freestream.gravity           = 9.81               
            run_conditions.freestream.speed_of_sound    = atmo_data.speed_of_sound[0,0]
            # run_conditions.freestream.velocity          = Mach[i] * run_conditions.freestream.speed_of_sound
            run_conditions.freestream.velocity          = _Mach * run_conditions.freestream.speed_of_sound
            # run_conditions.freestream.mach_number       = Mach[i]
            run_conditions.freestream.mach_number       = _Mach
            # run_conditions.aerodynamics.side_slip_angle = Beta[j]  # NILS: see `conditions.aerodynamics.angle_of_attack[i]/Units.deg` in trunk/SUAVE/Methods/Aerodynamics/AVL/translate_data.py (not applied to beta)
            run_conditions.aerodynamics.side_slip_angle = _Beta  # NILS: see `conditions.aerodynamics.angle_of_attack[i]/Units.deg` in trunk/SUAVE/Methods/Aerodynamics/AVL/translate_data.py (not applied to beta)
            run_conditions.aerodynamics.angle_of_attack = AoA
            run_conditions.aerodynamics.roll_rate_coefficient  = roll_rate_coefficient  # NILS: added to match SUAVE 2.5.2
            # run_conditions.aerodynamics.lift_coefficient       = lift_coefficient  # NILS: added to match SUAVE 2.5.2
            run_conditions.aerodynamics.pitch_rate_coefficient = pitch_rate_coefficient  # NILS: added to match SUAVE 2.5.2
            
            # Update aircraft mass
            geometry.mass_properties.mass = _W / 9.81
            # Update aircraft load factor
            run_conditions.aerodynamics.load_factor = _n
            # Calculate required lift coefficient
            CL = 2 * _n * _W / (run_conditions.freestream.density * run_conditions.freestream.velocity**2 * geometry.wings.main_wing.areas.reference)
            run_conditions.aerodynamics.lift_coefficient = CL  # NILS: added to match SUAVE 2.5.2
            
            #Run Analysis at AoA[i] and Mach[i]
            results =  self.evaluate_conditions(run_conditions, trim_aircraft)
            results_list.append(results)  # NILS: added to store dynamic stability analysis results

            # Obtain CM Cm_alpha, Cn_beta and the Neutral Point 
            # CM[:,i,j,k]       = results.aerodynamics.Cmtot[:,0]
            # Cm_alpha[:,i,j,k] = results.stability.static.Cm_alpha[:,0]
            # Cn_beta[:,i,j,k]  = results.stability.static.Cn_beta[:,0]
            # NP[:,i,j,k]       = results.stability.static.neutral_point[:,0]
            # static_margin[:,i,j,k] = (results.stability.static.neutral_point[:,0] - cg) / MAC  # NILS: added to match SUAVE 2.5.2
            # Cl_beta[:,i,j,k]  = results.stability.static.Cl_beta[:,0]  # NILS: added independently
            # Cn_r[:,i,j,k]  = results.stability.static.Cn_r[:,0]  # NILS: added independently
            # Cl_r[:,i,j,k]  = results.stability.static.Cl_r[:,0]  # NILS: added independently
            CM.append(results.aerodynamics.Cmtot[0,0])
            Cm_alpha.append(results.stability.static.Cm_alpha[0,0])
            Cn_beta.append(results.stability.static.Cn_beta[0,0])
            NP.append(results.stability.static.neutral_point[0,0])
            static_margin.append((results.stability.static.neutral_point[0,0] - cg) / MAC)  # NILS: added to match SUAVE 2.5.2
            Cl_beta.append(results.stability.static.Cl_beta[0,0])  # NILS: added independently
            Cn_r.append(results.stability.static.Cn_r[0,0])  # NILS: added independently
            Cl_r.append(results.stability.static.Cl_r[0,0])  # NILS: added independently
            
            # NILS: added independently
            # AoA_list.append(AoA)
            # Mach_list.append(_Mach * np.ones_like(AoA))
            # Beta_list.append(_Beta * np.ones_like(AoA))
            # h_list.append(_h * np.ones_like(AoA))
            AoA_list.append(results.aerodynamics.AoA[0,0])  # NILS: this is an OUTPUT of the trim analysis (deg)
            Mach_list.append(_Mach)
            Beta_list.append(_Beta)
            h_list.append(_h)
            n_list.append(_n)
            W_list.append(_W)
                    
            #         break
            #     break
            # break
        # sys.exit('Done.')
        
        # NILS: convert lists of row arrays into stacked column arrays
        # AoA_col = np.concatenate(AoA_list)[:, None]
        # Mach_col = np.concatenate(Mach_list)[:, None]
        # Beta_col = np.concatenate(Beta_list)[:, None]
        # h_col = np.concatenate(h_list)[:, None]
        # n_col = np.concatenate(n_list)[:, None]
        # W_col = np.concatenate(W_list)[:, None]
        # =============================================================================
        sigma_fcs_col = np.ones_like(np.array(AoA_list)[:, None]) * sigma_fcs
        span_loc_col = np.ones_like(np.array(AoA_list)[:, None]) * span_loc
        fcs_loc_col = np.ones_like(np.array(AoA_list)[:, None]) * fcs_loc
        wing_frac_col = np.ones_like(np.array(AoA_list)[:, None]) * wing_frac
        nacelle_frac_col = np.ones_like(np.array(AoA_list)[:, None]) * nacelle_frac
        # =============================================================================
        AoA_col = np.array(AoA_list)[:, None]
        Mach_col = np.array(Mach_list)[:, None]
        Beta_col = np.array(Beta_list)[:, None]
        h_col = np.array(h_list)[:, None]
        n_col = np.array(n_list)[:, None]
        W_col = np.array(W_list)[:, None]
        inputs_stack = np.column_stack([
            sigma_fcs_col, span_loc_col, fcs_loc_col, wing_frac_col, nacelle_frac_col,
            AoA_col, Mach_col, Beta_col, h_col, n_col, W_col,
        ])
            
        if self.training_file:
            # load data 
            data_array   = np.loadtxt(self.training_file)  
            CM_1D        = np.atleast_2d(data_array[:,0]) 
            Cm_alpha_1D  = np.atleast_2d(data_array[:,1])            
            Cn_beta_1D   = np.atleast_2d(data_array[:,2])
            NP_1D        = np.atleast_2d(data_array[:,3])
            
            # convert from 1D to 2D
            CM        = np.reshape(CM_1D, (len(AoA),-1))
            Cm_alpha  = np.reshape(Cm_alpha_1D, (len(AoA),-1))
            Cn_beta   = np.reshape(Cn_beta_1D , (len(AoA),-1))
            NP        = np.reshape(NP_1D , (len(AoA),-1))
        
        # Save the data for regression 
        if self.save_regression_results:
            # convert from 2D to 1D
            # CM_1D       = CM.reshape([len(AoA)*len(Mach)*len(Beta)*len(h),1]) 
            # Cm_alpha_1D = Cm_alpha.reshape([len(AoA)*len(Mach)*len(Beta)*len(h),1])  
            # Cn_beta_1D  = Cn_beta.reshape([len(AoA)*len(Mach)*len(Beta)*len(h),1])         
            # NP_1D       = NP.reshape([len(AoA)*len(Mach)*len(Beta)*len(h),1])  # NILS: TYPO - said `Cn_beta` like in line above
            # static_margin_1D = static_margin.reshape([len(AoA)*len(Mach)*len(Beta)*len(h),1])  # NILS: added to match SUAVE 2.5.2
            # Cl_beta_1D = Cl_beta.reshape([len(AoA)*len(Mach)*len(Beta)*len(h),1])  # NILS: added independently
            # Cn_r_1D = Cn_r.reshape([len(AoA)*len(Mach)*len(Beta)*len(h),1])  # NILS: added independently
            # Cl_r_1D = Cl_r.reshape([len(AoA)*len(Mach)*len(Beta)*len(h),1])  # NILS: added independently
            CM_1D       = np.array(CM)[:, None] 
            Cm_alpha_1D = np.array(Cm_alpha)[:, None]  
            Cn_beta_1D  = np.array(Cn_beta)[:, None]         
            NP_1D       = np.array(NP)[:, None]  # NILS: TYPO - said `Cn_beta` like in line above
            static_margin_1D = np.array(static_margin)[:, None]  # NILS: added to match SUAVE 2.5.2
            Cl_beta_1D = np.array(Cl_beta)[:, None]  # NILS: added independently
            Cn_r_1D = np.array(Cn_r)[:, None]  # NILS: added independently
            Cl_r_1D = np.array(Cl_r)[:, None]  # NILS: added independently
            
            static_stability_file_name = (
                'suave_static_stability_outputs_' + str(study_idx) + '_' + str(counter) + '.txt'
                if (study_idx != None and counter != None)
                else 'suave_static_stability_outputs.txt'
            )  # NILS: added for better logging control
            np.savetxt(
                # geometry.tag+'_stability_data.txt',
                static_stability_file_name,  # NILS: replaces line above
                np.hstack([
                    inputs_stack,
                    CM_1D,
                    Cm_alpha_1D,
                    Cn_beta_1D,
                    NP_1D,static_margin_1D,
                    Cl_beta_1D,
                    Cn_r_1D,
                    Cl_r_1D,
                ]),
                fmt='%10.8f',
                # header='   AoA       Mach       Beta       hCM       Cm_alpha       Cn_beta       NP       static_margin       Cl_beta       Cn_r       Cl_r'
                header=(
                    '   sigma_fcs       span_loc       fcs_loc       wing_frac       nacelle_frac       '
                    'AoA       Mach       Beta       h       n       W       '
                    'hCM       Cm_alpha       Cn_beta       NP       static_margin       Cl_beta       Cn_r       Cl_r'
                )
            )
        
        # <<< NILS: save dynamic stability results too
        if self.save_regression_results:
            
            # Loop over dynamic stability analysis results across all Mach numbers
            ev_list = []  # eigenvalues
            sm_list = []  # system matrix
            for results in results_list:
                
                ''' NILS: uncomment this to perform modal analysis with trunk\SUAVE\Methods\Flight_Dynamics\Dynamic_Stability\compute_dynamic_flight_modes.py
                # Extract SUAVE dynamic stability results
                dyn = results.dynamic_stability
                long = dyn.LongModes
                lat = dyn.LatModes
            
                # ---- Collect all fields (robust to scalars, 1D arrays, lists, eigenvalues) ---- #
                fields = [
                    long.phugoidInd,
                    long.phugoidFreqHz,
                    long.phugoidDamp,
                    long.phugoidTimeDoubleHalf,
                    long.shortPeriodInd,
                    long.shortPeriodFreqHz,
                    long.shortPeriodDamp,
                    long.shortPeriodTimeDoubleHalf,
            
                    lat.dutchRollInd,
                    lat.dutchRollFreqHz,
                    lat.dutchRollDamping,
                    lat.dutchRollTimeDoubleHalf,
                    lat.dutchRoll_mode_real,
                    lat.rollSubsistenceInd,
                    lat.rollSubsistenceFreqHz,
                    lat.rollSubsistenceTimeConstant,
                    lat.rollSubsistenceDamping,
                    lat.spiralInd,
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
            
                # Build final 2D array
                data_out = np.column_stack(cols_padded)
            
                # ----------------------- Write to text file ----------------------------- #
                header = (
                    # "AoA  Mach  Beta  h  "
                    "AoA  Mach  Beta  h  n  W  "
                    "phugoidInd  phugoidFreqHz  phugoidDamp  phugoidTimeDoubleHalf  "
                    "shortPeriodInd  shortPeriodFreqHz  shortPeriodDamp  shortPeriodTimeDoubleHalf  "
                    "dutchRollInd  dutchRollFreqHz  dutchRollDamping  dutchRollTimeDoubleHalf  dutchRoll_mode_real  "
                    "rollSubsistenceInd  rollSubsistenceFreqHz  rollSubsistenceTimeConstant  rollSubsistenceDamping  "
                    "spiralInd  spiralFreqHz  spiralTimeDoubleHalf  spiralDamping  "
                    "pMax  "
                    "Long_Re1  Long_Im1  Long_Re2  Long_Im2  Long_Re3  Long_Im3  Long_Re4  Long_Im4  "
                    "Lat_Re1   Lat_Im1   Lat_Re2   Lat_Im2   Lat_Re3   Lat_Im3   Lat_Re4   Lat_Im4  "
                    "Long_A  Long_B  Long_C  Long_D  Long_E  "
                    "Lat_A  Lat_B  Lat_C  Lat_D  Lat_E"
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
                long_poly_padded = np.array(long.polyLon)
                lat_poly_padded = np.array(lat.polyLat)
                
                # Final array: [dynamic stability scalars | LongModes | LatModes]
                final_out = np.hstack([
                    data_out_padded,
                    long_modes_padded,
                    lat_modes_padded,
                    long_poly_padded,
                    lat_poly_padded
                ])
                '''
                
                # NILS: use this to perform modal analysis with AVL
                final_out = np.hstack((
                    np.squeeze(np.array(results.stability.dynamic.eigenvalues_real)),
                    np.squeeze(np.array(results.stability.dynamic.eigenvalues_imag)),
                ))
                ev_list.append(final_out)
                sm_list.append(results.stability.dynamic.system_matrix)
                
            ''' NILS: uncomment this to perform modal analysis with trunk\SUAVE\Methods\Flight_Dynamics\Dynamic_Stability\compute_dynamic_flight_modes.py
            # Stack dynamic stability analysis results across all Mach numbers
            ev_list_stacked = np.hstack((inputs_stack, np.vstack(ev_list)))
            '''
            
            # NILS: use this to perform modal analysis with AVL
            header = (
                'sigma_fcs  span_loc  fcs_loc  wing_frac  nacelle_frac  '
                'AoA  Mach  Beta  h  n  W  '
                'Re1  Re2  Re3  Re4  Re5  Re6  Re7  Re8  Im1  Im2  Im3  Im4  Im5  Im6  Im7  Im8'
            )
            # ev_list_stacked = np.vstack(ev_list)
            print('ev_list =', ev_list)
            # ev_list_stacked = np.hstack((inputs_stack, np.vstack(ev_list)))
            
            # =============================================================================
            try:
                ev_stacked = np.vstack(ev_list)
            
            except ValueError:
                # --- Pad rows with NaNs so all arrays have equal length ---
                # Find the longest eigenvalue vector
                max_len = max(arr.size for arr in ev_list)
            
                ev_padded = []
                for arr in ev_list:
                    if arr.size < max_len:
                        pad = np.full(max_len - arr.size, np.nan)
                        arr = np.concatenate((arr, pad))
                    ev_padded.append(arr)
            
                # Now safe to vstack
                ev_stacked = np.vstack(ev_padded)
            
            # Finally:
            ev_list_stacked = np.hstack((inputs_stack, ev_stacked))
            # =============================================================================
                
            dynamic_stability_file_name = (
                'suave_dynamic_stability_outputs_' + str(study_idx) + '_' + str(counter) + '.txt'
                if (study_idx != None and counter != None) else 'suave_dynamic_stability_outputs.txt'
            )  # added for better logging control
            np.savetxt(
                # geometry.tag + "_dynamic_stability_data.txt",
                dynamic_stability_file_name,  # replaces line above
                ev_list_stacked,  # data_out,
                fmt="%12.6f",
                header=header,
                comments="",  # prevents '#' from being added
            )
            
            # =============================================================================
            # NILS: use this to perform modal analysis with AVL
            header = (
                'u  w  q  the  v  p  r  phi  x  y  z  psi  |  slat  flap  aileron  elevator'
            )
            sm_list_stacked = np.vstack(sm_list)
            # print()
            # print('sm_list_stacked =', sm_list_stacked)
            # print()
            # print('np.shape(sm_list_stacked) =', np.shape(sm_list_stacked))
            # sys.exit()
                
            dynamic_stability_file_name = (
                'suave_dynamic_stability_matrix_' + str(study_idx) + '_' + str(counter) + '.txt'
                if (study_idx != None and counter != None) else 'suave_dynamic_stability_outputs.txt'
            )  # added for better logging control
            # np.savetxt(
            #     # geometry.tag + "_dynamic_stability_data.txt",
            #     dynamic_stability_file_name,  # replaces line above
            #     sm_list_stacked,  # data_out,
            #     fmt="%12.6f",
            #     header=header,
            #     comments="",  # prevents '#' from being added
            # )
            with open(dynamic_stability_file_name, "w") as fh:
                fh.write(header + "\n")
                for k in range(sm_list_stacked.shape[0]):
                    if k > 0:
                        fh.write("\n")  # empty line between 2D blocks
                    np.savetxt(
                        fh,
                        sm_list_stacked[k],
                        fmt="%12.6f"
                    )

            # =============================================================================
            
        # >>>
        
        # Store training data
        # Save the data for regression
        training_data = np.zeros((4,len(AoA),len(Mach)))
        training_data[0,:,:] = CM       
        training_data[1,:,:] = Cm_alpha 
        training_data[2,:,:] = Cn_beta  
        training_data[3,:,:] = NP      
            
        # Store training data
        training.coefficients = training_data
 
        return        

    def build_surrogate(self):
        """Builds a surrogate based on sample evalations using a Guassian process.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        self.training.
          coefficients     [-] CM, Cm_alpha, Cn_beta 
          neutral point    [meters] NP
          grid_points      [radians,-] angles of attack and Mach numbers 

        Outputs:
        self.surrogates.
          moment_coefficient           
          Cm_alpha_moment_coefficient  
          Cn_beta_moment_coefficient   
          neutral_point                      

        Properties Used:
        No others
        """
        # Unpack data
        training                                    = self.training
        AoA_data                                    = training.angle_of_attack
        mach_data                                   = training.Mach
        CM_data                                     = training.coefficients[0,:,:]
        Cm_alpha_data                               = training.coefficients[1,:,:]
        Cn_beta_data                                = training.coefficients[2,:,:]
        NP_data                                     = training.coefficients[3,:,:]
        
        self.surrogates.moment_coefficient          = RectBivariateSpline(AoA_data, mach_data, CM_data      ) 
        self.surrogates.Cm_alpha_moment_coefficient = RectBivariateSpline(AoA_data, mach_data, Cm_alpha_data) 
        self.surrogates.Cn_beta_moment_coefficient  = RectBivariateSpline(AoA_data, mach_data, Cn_beta_data ) 
        self.surrogates.neutral_point               = RectBivariateSpline(AoA_data, mach_data, NP_data      )  
        
        return

    
# ----------------------------------------------------------------------
#  Helper Functions
# ----------------------------------------------------------------------
        
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
        # run_script_path                  = run_folder.rstrip('avl_files').rstrip('/')
        run_script_path = run_folder.rstrip(self.settings.filenames.run_folder).rstrip('/')
        aero_results_template_1          = self.settings.filenames.aero_output_template_1       # 'stability_axis_derivatives_{}.dat' 
        aero_results_template_2          = self.settings.filenames.aero_output_template_2       # 'surface_forces_{}.dat'
        aero_results_template_3          = self.settings.filenames.aero_output_template_3       # 'strip_forces_{}.dat'   
        aero_results_template_4          = self.settings.filenames.aero_output_template_4       # 'body_axis_derivatives_{}.dat'     
        dynamic_results_template_1       = self.settings.filenames.dynamic_output_template_1    # 'eigen_mode_{}.dat'
        dynamic_results_template_2       = self.settings.filenames.dynamic_output_template_2    # 'system_matrix_{}.dat'
        batch_template                   = self.settings.filenames.batch_template
        deck_template                    = self.settings.filenames.deck_template
        print_output                     = self.settings.print_output  # NILS: added to match SUAVE 2.5.2
        
        # rename defaul avl aircraft tag
        self.tag                         = 'avl_analysis_of_{}'.format(self.geometry.tag) 
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
        control_surfaces = False  # NILS: added to match SUAVE 2.5.2
        for wing in self.geometry.wings: # this parses through the wings to determine how many control surfaces does the vehicle have 
            if wing.control_surfaces:
                control_surfaces = True
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
        
        # NILS: raise Exception if control surfaces have not been defined
        # I am only considering trimmed cases, for which these are always needed
        if not control_surfaces:
            raise Exception("No control surfaces have been defined.")
        
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
            write_input_deck(self, trim_aircraft, control_surfaces, run_modal=True)  # NILS: added last argument to match SUAVE 2.5.2

            # RUN AVL!
            results_avl = run_analysis(self, print_output)  # NILS: added last argument to match SUAVE 2.5.2
    
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
        
        return results
