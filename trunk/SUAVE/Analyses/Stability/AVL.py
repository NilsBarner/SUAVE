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
        self.keep_files                             = True  # False  # NILS: save these for plotting geometry
        self.save_regression_results                = True  # False  # NILS: save these for visualisation
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
        # print('__call__')
        
        # Unpack
        surrogates          = self.surrogates       
        Mach                = conditions.freestream.mach_number
        AoA                 = conditions.aerodynamics.angle_of_attack 
        moment_model        = surrogates.moment_coefficient
        Cm_alpha_model      = surrogates.Cm_alpha_moment_coefficient
        Cn_beta_model       = surrogates.Cn_beta_moment_coefficient      
        neutral_point_model = surrogates.neutral_point
        cg                  = self.geometry.mass_properties.center_of_gravity[0][0]  # copied over from SUAVE 2.5.2
        MAC                 = self.geometry.wings.main_wing.chords.mean_aerodynamic  # copied over from SUAVE 2.5.2
        
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


    def sample_training(self):
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
        # print('sample_training')
        # Unpack
        run_folder    = os.path.abspath(self.settings.filenames.run_folder)
        geometry      = self.geometry
        # print("geometry.wings['main_wing'].Segments.keys() =", geometry.wings['main_wing'].Segments.keys())
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
        Cl_beta = np.zeros_like(CM)  # added by NILS
        Cn_r = np.zeros_like(CM)  # added by NILS
        Cl_r = np.zeros_like(CM)  # added by NILS
       
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
            Cl_beta[:,i]  = results.stability.static.Cl_beta[:,0]
            Cn_r[:,i]  = results.stability.static.Cn_r[:,0]
            Cl_r[:,i]  = results.stability.static.Cl_r[:,0]
            
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
            CM_1D       = CM.reshape([len(AoA)*len(Mach),1]) 
            Cm_alpha_1D = Cm_alpha.reshape([len(AoA)*len(Mach),1])  
            Cn_beta_1D  = Cn_beta.reshape([len(AoA)*len(Mach),1])         
            NP_1D       = NP.reshape([len(AoA)*len(Mach),1])  # NILS: TYPO - said `Cn_beta` like in line above
            static_margin_1D = static_margin.reshape([len(AoA)*len(Mach),1])  # copied over from SUAVE 2.5.2
            Cl_beta_1D = Cl_beta.reshape([len(AoA)*len(Mach),1])  # added by NILS
            Cn_r_1D = Cn_r.reshape([len(AoA)*len(Mach),1])  # added by NILS
            Cl_r_1D = Cl_r.reshape([len(AoA)*len(Mach),1])  # added by NILS
            # print("geometry.tag+'_stability_data.txt' =", geometry.tag+'_stability_data.txt')
            np.savetxt(
                geometry.tag+'_stability_data.txt',
                np.hstack([
                    CM_1D,Cm_alpha_1D, Cn_beta_1D,NP_1D,static_margin_1D, Cl_beta_1D, Cn_r_1D, Cl_r_1D,
                ]),fmt='%10.8f',header='   CM       Cm_alpha       Cn_beta       NP       static_margin       Cl_beta       Cn_r       Cl_r')
            
        # =============================================================================
        # NILS: save dynamic stability results too
        if self.save_regression_results:
            
            # Loop over dynamic stability analysis results across all Mach numbers
            final_out_list = []
            for results in results_list:
            
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
                
                # =============================================================================
                # @ NILS: next, add long.polyLon and lat.polyLat to header below!
                # =============================================================================
            
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
                final_out = np.hstack([data_out_padded, long_modes_padded, lat_modes_padded, long_poly_padded, lat_poly_padded])
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
        
            sys.exit("Saved SUAVE dynamic stability results to dynamic_stability_data.txt")
        # =============================================================================
        
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
        # print('build_surrogate')
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
                    
        # raise Exception('Follow traceback from here.')                                   
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
        # print('evaluate_conditions')
        
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
        # print('self.__class__ =', self.__class__)
        # print('self.geometry._base =', self.geometry._base)
        # raise Exception
        # import sys
        # sys.exit()
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
        # print('self.geometry.mass_properties.moments_of_inertia.tensor =', self.geometry.mass_properties.moments_of_inertia.tensor)
        if np.count_nonzero(self.geometry.mass_properties.moments_of_inertia.tensor) > 0:
            results = compute_dynamic_flight_modes(results,self.geometry,run_conditions,cases)        
            
        # print('TEST results.dynamic_stability.LongModes =', results.dynamic_stability.LongModes)
        # sys.exit()
             
        if not self.keep_files:
            rmtree( run_folder )   
        
        # sys.exit('Stop here.')
        
        return results
