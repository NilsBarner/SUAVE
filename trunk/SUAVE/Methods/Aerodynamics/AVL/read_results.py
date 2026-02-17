## @ingroup Methods-Aerodynamics-AVL
#read_results.py
# 
# Created:  Mar 2015, T. Momose
# Modified: Jan 2016, E. Botero
#           Dec 2017, M. Clarke
#           Aug 2019, M. Clarke
# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

"""
NILS: INDICES IN THIS FILE HAVE BEEN UPDATED SUBSTANTIALLY
BY NILS TO COMPLY WITH AVL 3.52 (2025)! EVEN
trunk/SUAVE/Methods/Aerodynamics/AVL OF SUAVE 2.5.2 IS
OUT-OF-DATE!
"""

from SUAVE.Core import Data
from SUAVE.Methods.Aerodynamics.AVL.Data.Wing import Control_Surface_Data ,  Control_Surface_Results 
import numpy as np 

## @ingroup Methods-Aerodynamics-AVL
def read_results(avl_object, backend='AVL'):  # NILS: skip read_results() when JVL backend is used (else encounter read error in surface forces result file)
    """ This functions reads the results from the results text file created 
    at the end of an AVL function call

    Assumptions:
        None
        
    Source:
        Drela, M. and Youngren, H., AVL, http://web.mit.edu/drela/Public/web/avl

    Inputs:
        None

    Outputs:
        results     

    Properties Used:
        N/A
    """    
    # unpack
    aircraft = avl_object.geometry
    results  = Data()
    case_idx = 0  
    for case in avl_object.current_status.cases:
        num_ctrl =  case.stability_and_control.number_control_surfaces
        # open newly written result files and read in aerodynamic properties 
        with open(case.aero_result_filename_1,'r') as stab_der_vile:
            # Extract results from stability axis derivatives file                                                                
            case_res                                                        = Data()  
            case_res.aerodynamics                                           = Data()
            case_res.stability                                              = Data()
            case_res.stability.control_surfaces                             = Control_Surface_Data()   
            case_res.stability.alpha_derivatives                            = Data()
            case_res.stability.beta_derivatives                             = Data()   
                                                                            
            case_res.tag                                                    = case.tag 
            lines                                                           = stab_der_vile.readlines()

            if backend == 'AVL':  # NILS

                case_res.S_ref                                                  = float(lines[8][8:16].strip())
                case_res.c_ref                                                  = float(lines[8][29:37].strip())
                case_res.b_ref                                                  = float(lines[8][50:58].strip())
                case_res.X_ref                                                  = float(lines[9][8:16].strip())
                case_res.Y_ref                                                  = float(lines[9][29:37].strip())
                case_res.Z_ref                                                  = float(lines[9][50:58].strip())

                case_res.aerodynamics.AoA                                       = float(lines[15][9:19].strip())
                case_res.aerodynamics.CX                                        = float(lines[19][9:19].strip())
                case_res.aerodynamics.CY                                        = float(lines[20][9:19].strip()) 
                case_res.aerodynamics.CZ                                        = float(lines[21][9:19].strip())

                case_res.aerodynamics.Cltot                                     = float(lines[19][31:41].strip())
                case_res.aerodynamics.Cmtot                                     = float(lines[20][31:41].strip()) 
                case_res.aerodynamics.Cntot                                     = float(lines[21][31:41].strip())

                case_res.aerodynamics.roll_moment_coefficient                   = float(lines[19][31:42].strip())  # NILS: repeat of above
                case_res.aerodynamics.pitch_moment_coefficient                  = float(lines[20][31:42].strip())  # NILS: repeat of above
                case_res.aerodynamics.yaw_moment_coefficient                    = float(lines[21][31:42].strip())  # NILS: repeat of above
                case_res.aerodynamics.total_lift_coefficient                    = float(lines[23][9:19].strip())
                case_res.aerodynamics.total_drag_coefficient                    = float(lines[24][9:19].strip())
                case_res.aerodynamics.viscous_drag_coefficient                  = float(lines[25][9:19].strip())
                case_res.aerodynamics.induced_drag_coefficient                  = float(lines[25][31:41].strip())
                case_res.aerodynamics.oswald_efficiency                         = float(lines[27][31:41].strip())

                case_res.stability.alpha_derivatives.lift_curve_slope           = float(lines[36+num_ctrl][23:34].strip()) # CL_a
                case_res.stability.alpha_derivatives.side_force_derivative      = float(lines[37+num_ctrl][23:34].strip()) # CY_a
                # NILS: note that CDa is not used by SUAVE (probably not available in the version of AVL used at the time)
                case_res.stability.alpha_derivatives.roll_moment_derivative     = float(lines[39+num_ctrl][23:34].strip()) # Cl_a
                case_res.stability.alpha_derivatives.pitch_moment_derivative    = float(lines[40+num_ctrl][23:34].strip()) # Cm_a
                case_res.stability.alpha_derivatives.yaw_moment_derivative      = float(lines[41+num_ctrl][23:34].strip()) # Cn_a
                case_res.stability.beta_derivatives.lift_coefficient_derivative = float(lines[36+num_ctrl][43:54].strip()) # CL_b
                case_res.stability.beta_derivatives.side_force_derivative       = float(lines[37+num_ctrl][43:54].strip()) # CY_b
                # NILS: note that CDb is not used by SUAVE (probably not available in the version of AVL used at the time)
                case_res.stability.beta_derivatives.roll_moment_derivative      = float(lines[39+num_ctrl][43:54].strip()) # Cl_b
                case_res.stability.beta_derivatives.pitch_moment_derivative     = float(lines[40+num_ctrl][43:54].strip()) # Cm_b
                case_res.stability.beta_derivatives.yaw_moment_derivative       = float(lines[41+num_ctrl][43:54].strip()) # Cn_b

                case_res.stability.CL_p                                         = float(lines[45+num_ctrl][23:34].strip())
                case_res.stability.CL_q                                         = float(lines[45+num_ctrl][43:54].strip())
                case_res.stability.CL_r                                         = float(lines[45+num_ctrl][63:74].strip())
                case_res.stability.CY_p                                         = float(lines[46+num_ctrl][23:34].strip())  # NILS: gets overwritten further down!
                case_res.stability.CY_q                                         = float(lines[46+num_ctrl][43:54].strip())
                case_res.stability.CY_r                                         = float(lines[46+num_ctrl][63:74].strip())
                # NILS: note that CDp, CDq, and CDr are not used by SUAVE (probably not available in the version of AVL used at the time)
                case_res.stability.Cl_p                                         = float(lines[48+num_ctrl][23:34].strip())
                case_res.stability.Cl_q                                         = float(lines[48+num_ctrl][43:54].strip())
                case_res.stability.Cl_r                                         = float(lines[48+num_ctrl][63:74].strip())
                case_res.stability.Cm_p                                         = float(lines[49+num_ctrl][23:34].strip())
                case_res.stability.Cm_q                                         = float(lines[49+num_ctrl][43:54].strip())
                case_res.stability.Cm_r                                         = float(lines[49+num_ctrl][63:74].strip())
                case_res.stability.Cn_p                                         = float(lines[50+num_ctrl][23:34].strip())
                case_res.stability.Cn_q                                         = float(lines[50+num_ctrl][43:54].strip())
                case_res.stability.Cn_r                                         = float(lines[50+num_ctrl][63:74].strip())

            elif backend == 'JVL':  # NILS

                case_res.S_ref = float(lines[8][8:17].strip())
                case_res.c_ref = float(lines[8][30:39].strip())
                case_res.b_ref = float(lines[8][52:61].strip())
                case_res.X_ref = float(lines[9][8:17].strip())
                case_res.Y_ref = float(lines[9][30:39].strip())
                case_res.Z_ref = float(lines[9][52:61].strip())
                case_res.X_mom = float(lines[10][8:17].strip())
                case_res.Y_mom = float(lines[10][30:39].strip())
                case_res.Z_mom = float(lines[10][52:61].strip())
                # print('case_res.S_ref, case_res.c_ref, case_res.b_ref, case_res.X_ref, case_res.Y_ref, case_res.Z_ref, case_res.X_mom, case_res.Y_mom, case_res.Z_mom =', case_res.S_ref, case_res.c_ref, case_res.b_ref, case_res.X_ref, case_res.Y_ref, case_res.Z_ref, case_res.X_mom, case_res.Y_mom, case_res.Z_mom)
            
                case_res.aerodynamics.AoA = float(lines[16][9:19].strip())
                case_res.aerodynamics.CX = float(lines[20][9:19].strip())
                case_res.aerodynamics.CY = float(lines[21][9:19].strip()) 
                case_res.aerodynamics.CZ = float(lines[22][9:19].strip())
                # print('case_res.aerodynamics.AoA, case_res.aerodynamics.CX, case_res.aerodynamics.CY, case_res.aerodynamics.CZ =', case_res.aerodynamics.AoA, case_res.aerodynamics.CX, case_res.aerodynamics.CY, case_res.aerodynamics.CZ)

                case_res.aerodynamics.Cltot = float(lines[20][31:41].strip())
                case_res.aerodynamics.Cmtot = float(lines[21][31:41].strip()) 
                case_res.aerodynamics.Cntot = float(lines[22][31:41].strip())
                # print('case_res.aerodynamics.Cltot, case_res.aerodynamics.Cmtot, case_res.aerodynamics.Cntot =', case_res.aerodynamics.Cltot, case_res.aerodynamics.Cmtot, case_res.aerodynamics.Cntot)
                
                case_res.aerodynamics.total_lift_coefficient = float(lines[24][9:19].strip())
                case_res.aerodynamics.total_drag_coefficient = float(lines[26][9:19].strip())
                case_res.aerodynamics.viscous_drag_coefficient = float(lines[33][9:19].strip())
                case_res.aerodynamics.induced_drag_coefficient = float(lines[31][9:19].strip())
                case_res.aerodynamics.oswald_efficiency = float(lines[26][52:62].strip())
                # print('case_res.aerodynamics.total_lift_coefficient, case_res.aerodynamics.total_drag_coefficient, case_res.aerodynamics.viscous_drag_coefficient, case_res.aerodynamics.induced_drag_coefficient, case_res.aerodynamics.oswald_efficiency =', case_res.aerodynamics.total_lift_coefficient, case_res.aerodynamics.total_drag_coefficient, case_res.aerodynamics.viscous_drag_coefficient, case_res.aerodynamics.induced_drag_coefficient, case_res.aerodynamics.oswald_efficiency)

                # JVL-specific output
                case_res.aerodynamics.CLjet = float(lines[29][9:19].strip())
                case_res.aerodynamics.CDjet = float(lines[32][9:19].strip())
                case_res.aerodynamics.DCQtot = float(lines[35][9:21].strip())
                case_res.aerodynamics.DCJtot = float(lines[35][32:44].strip()) 
                case_res.aerodynamics.DCEtot = float(lines[35][55:67].strip()) 
                case_res.aerodynamics.CTtot = float(lines[36][32:44].strip()) 
                case_res.aerodynamics.CPtot = float(lines[36][55:67].strip()) 
                # print('case_res.aerodynamics.CLjet, case_res.aerodynamics.CDjet, case_res.aerodynamics.DCQtot, case_res.aerodynamics.DCJtot, case_res.aerodynamics.DCEtot, case_res.aerodynamics.CTtot, case_res.aerodynamics.CPtot =', case_res.aerodynamics.CLjet, case_res.aerodynamics.CDjet, case_res.aerodynamics.DCQtot, case_res.aerodynamics.DCJtot, case_res.aerodynamics.DCEtot, case_res.aerodynamics.CTtot, case_res.aerodynamics.CPtot)
            
            # this block of text reads in aerodynamic results related to the defined control surfaces
            if num_ctrl != 0: 
                for ctrl_idx in range(num_ctrl):
                    ctrl_surf = Control_Surface_Results()
                    ctrl_surf.tag                 = str(lines[29+ctrl_idx][2:11].strip())
                    ctrl_surf.deflection          = float(lines[29+ctrl_idx][20:30].strip())
                    ctrl_surf.CL                  = float(lines[54+num_ctrl][(21*ctrl_idx + 24):(21*ctrl_idx + 35)].strip())
                    ctrl_surf.CY                  = float(lines[55+num_ctrl][(21*ctrl_idx + 24):(21*ctrl_idx + 35)].strip())
                    # NILS: note that CDd01, CDd02, CDd03, and CDd04 are not used by SUAVE (probably not available in the version of AVL used at the time)
                    ctrl_surf.Cl                  = float(lines[57+num_ctrl][(21*ctrl_idx + 24):(21*ctrl_idx + 35)].strip())
                    ctrl_surf.Cm                  = float(lines[58+num_ctrl][(21*ctrl_idx + 24):(21*ctrl_idx + 35)].strip())
                    ctrl_surf.Cn                  = float(lines[59+num_ctrl][(21*ctrl_idx + 24):(21*ctrl_idx + 35)].strip())
                    ctrl_surf.CDff                = float(lines[60+num_ctrl][(21*ctrl_idx + 24):(21*ctrl_idx + 35)].strip())
                    ctrl_surf.e                   = float(lines[61+num_ctrl][(21*ctrl_idx + 24):(21*ctrl_idx + 35)].strip())
                    case_res.stability.control_surfaces.append_control_surface_result(ctrl_surf)

            if backend == 'AVL':  # NILS
                case_res.stability.neutral_point      = float(lines[54+11*(num_ctrl>0)+num_ctrl][22:33].strip())    
                case_res.stability.spiral_criteria    = float(lines[56+11*(num_ctrl>0)+num_ctrl][22:33].strip())
            elif backend == 'JVL':  # NILS
                case_res.stability.neutral_point      = float(lines[75+11*(num_ctrl>0)+num_ctrl][21:31].strip())    
                case_res.stability.spiral_criteria    = float(lines[77+11*(num_ctrl>0)+num_ctrl][21:32].strip())
                # print('case_res.stability.neutral_point, case_res.stability.spiral_criteria =', case_res.stability.neutral_point, case_res.stability.spiral_criteria)
        
        # get number of wings, spanwise discretization for surface and strip force result extraction
        n_sw    = avl_object.settings.number_spanwise_vortices
        n_wings = 0 
        for wing in aircraft.wings:
            n_wings += 1
            if wing.symmetric:
                n_wings += 1
        n_fus_sec = 0
        for fuselage in aircraft.fuselages:
            n_fus_sec += 2
            
        # NILS: add nacelles (symmetric - same number on each wing)
        n_nacelles = 0
        for nacelle in aircraft.nacelles:
            if nacelle.flow_through == True:
                n_nacelles += 2  # total of two surfaces (one per wing)
            elif nacelle.flow_through == False:
                n_nacelles += 4  # total of four surfaces (two per wing: one horizontal, one vertical)
        
        wing_area            = np.zeros(n_wings)
        wing_CL              = np.zeros(n_wings)
        wing_CD              = np.zeros(n_wings)  
        wing_local_span      = np.zeros((n_wings,n_sw))
        wing_sectional_chord = np.zeros((n_wings,n_sw))
        wing_cl              = np.zeros((n_wings,n_sw))
        alpha_i              = np.zeros((n_wings,n_sw))
        wing_cd              = np.zeros((n_wings,n_sw))   
        
        if backend == 'AVL':  # NILS

            # Extract resulst from surface forces result file
            with open(case.aero_result_filename_2,'r') as aero_res_file:
                aero_lines   = aero_res_file.readlines()
                line_idx     = 0
                header       = 12 + n_wings + n_fus_sec + n_nacelles  # NILS: added final term
                for i in range(n_wings):  # NILS: note that this EXCLUDES the fuselage (see surface_forces_cas_01_01.txt and below quote from https://suave.stanford.edu/tutorials/avl.html)
                    """Despite AVL having the capability of modelling bodies, a decision was made to model the fuselage
                    as a wake-producing, lifting surface. The entire body is defined by a series of vertical and horizontal
                    chords that create a cross when viewed from the front."""
                    wing_area[i] = float(aero_lines[header + line_idx][4:14].strip())
                    wing_CL[i]   = float(aero_lines[header + line_idx][24:33].strip())
                    wing_CD[i]   = float(aero_lines[header + line_idx][33:42].strip())
                    line_idx += 1                   
                case_res.aerodynamics.wing_areas = wing_area 
                case_res.aerodynamics.wing_CLs   = wing_CL 
                case_res.aerodynamics.wing_CDs   = wing_CD
            
        # Extract resulst from  strip forces result file
        with open(case.aero_result_filename_3,'r') as aero_res_file_2:

            if backend == 'AVL':  # NILS

                aero_lines_2     = aero_res_file_2.readlines()
                line_idx         = 0
                header           = 22
                divider_header   = 17
                
                for i in range(n_wings): 
                    for j in range(n_sw):
                        wing_local_span[i,j]      = float(aero_lines_2[header + j + line_idx][15:24].strip())
                        wing_sectional_chord[i,j] = float(aero_lines_2[header + j + line_idx][33:42].strip()) 
                        wing_cl[i,j]              = float(aero_lines_2[header + j + line_idx][78:87].strip())  
                        # At high angle of attacks, AVL does not give an answer 
                        try:
                            alpha_i[i,j]              = float(aero_lines_2[header + j + line_idx][60:69].strip())
                            wing_cd[i,j]              = float(aero_lines_2[header + j + line_idx][87:96].strip())
                        except:
                            alpha_i[i,j]              = 0.
                            wing_cd[i,j]              = 0.
                    line_idx = divider_header +  n_sw + line_idx            
                case_res.aerodynamics.wing_local_spans         = wing_local_span
                case_res.aerodynamics.wing_section_chords      = wing_sectional_chord 
                case_res.aerodynamics.wing_section_cls         = wing_cl 
                case_res.aerodynamics.wing_section_aoa_i       = alpha_i 
                case_res.aerodynamics.wing_section_cds         = wing_cd 

            elif backend == 'JVL':  # NILS

                aero_lines_2 = aero_res_file_2.readlines()
                line_idx = 0
                header = 19
                divider_header = 17

                for i in range(n_wings): 
                    
                    surface_info_line = aero_lines_2[header + line_idx - 13]  # NILS: e.g. "  Surface # 1     main_wing_1                             " or "  Surface #10     horizontal_stabilizer (YDUP)            "
                    if 'main' in surface_info_line:
                        n_sw = avl_object.settings.Nspanwise_main_wing  # NILS: changed w.r.t. AVL to allow varying Nspanwise_main_wing of MAIN WING only; done to avoid
                        # *** Cannot adjust spanwise spacing at SECTION  2, on SURFACE main_wing_1
                        # *** Insufficient number of spanwise vortices to work with
                    else:
                        n_sw = avl_object.settings.number_spanwise_vortices  # NILS: keep original value for other surfaces

                    # NILS: copied here from outside elif-statement to use local n_sw
                    wing_local_span = np.zeros((n_wings,n_sw))
                    wing_sectional_chord = np.zeros((n_wings,n_sw))
                    wing_cl = np.zeros((n_wings,n_sw))
                    alpha_i = np.zeros((n_wings,n_sw))
                    wing_cd = np.zeros((n_wings,n_sw))

                    for j in range(n_sw):
                        wing_local_span[i,j] = float(aero_lines_2[header + j + line_idx][14:30].strip())
                        # print('wing_local_span[i,j] =', wing_local_span[i,j])
                        wing_sectional_chord[i,j] = float(aero_lines_2[header + j + line_idx][45:50].strip()) 
                        # print('wing_sectional_chord[i,j] =', wing_sectional_chord[i,j])
                        wing_cl[i,j] = float(aero_lines_2[header + j + line_idx][102:110].strip())  
                        # print('wing_cl[i,j] =', wing_cl[i,j])
                        # At high angle of attacks, AVL does not give an answer 
                        try:
                            alpha_i[i,j] = float(aero_lines_2[header + j + line_idx][86:102].strip())
                            wing_cd[i,j] = float(aero_lines_2[header + j + line_idx][126:139].strip())
                        except:
                            alpha_i[i,j] = 0.
                            wing_cd[i,j] = 0.
                    line_idx = divider_header +  n_sw + line_idx            
                case_res.aerodynamics.wing_local_spans = wing_local_span
                case_res.aerodynamics.wing_section_chords = wing_sectional_chord 
                case_res.aerodynamics.wing_section_cls = wing_cl 
                case_res.aerodynamics.wing_section_aoa_i = alpha_i 
                case_res.aerodynamics.wing_section_cds = wing_cd
                # print('wing_local_span, wing_sectional_chord, wing_cl, alpha_i, wing_cd =', wing_local_span, wing_sectional_chord, wing_cl, alpha_i, wing_cd)

        if backend == 'AVL':  # NILS

            with open(case.aero_result_filename_4,'r') as bod_der_vile:
                # Extract results from body axis derivatives file                         
                                                            
                # NILS: every column has length 11
                lines_2                  = bod_der_vile.readlines() 
                case_res.stability.CX_u  = float(lines_2[36+num_ctrl][23:34].strip())
                case_res.stability.CX_v  = float(lines_2[36+num_ctrl][43:54].strip())
                case_res.stability.CX_w  = float(lines_2[36+num_ctrl][63:74].strip())
                case_res.stability.CY_u  = float(lines_2[37+num_ctrl][23:34].strip())
                case_res.stability.CY_v  = float(lines_2[37+num_ctrl][43:54].strip())
                case_res.stability.CY_w  = float(lines_2[37+num_ctrl][63:74].strip())
                case_res.stability.CZ_u  = float(lines_2[38+num_ctrl][23:34].strip())
                case_res.stability.CZ_v  = float(lines_2[38+num_ctrl][43:54].strip())
                case_res.stability.CZ_w  = float(lines_2[38+num_ctrl][63:74].strip())
                case_res.stability.Cl_u  = float(lines_2[39+num_ctrl][23:34].strip())
                case_res.stability.Cl_v  = float(lines_2[39+num_ctrl][43:54].strip())
                case_res.stability.Cl_w  = float(lines_2[39+num_ctrl][63:74].strip())
                case_res.stability.Cm_u  = float(lines_2[40+num_ctrl][23:34].strip())
                case_res.stability.Cm_v  = float(lines_2[40+num_ctrl][43:54].strip())
                case_res.stability.Cm_w  = float(lines_2[40+num_ctrl][63:74].strip())
                case_res.stability.Cn_u  = float(lines_2[41+num_ctrl][23:34].strip())
                case_res.stability.Cn_v  = float(lines_2[41+num_ctrl][43:54].strip())
                case_res.stability.Cn_w  = float(lines_2[41+num_ctrl][63:74].strip())
                
                case_res.stability.CX_p  = float(lines_2[45+num_ctrl][23:34].strip())
                case_res.stability.CX_q  = float(lines_2[45+num_ctrl][43:54].strip())
                case_res.stability.CX_r  = float(lines_2[45+num_ctrl][63:74].strip())
                case_res.stability.CY_p  = float(lines_2[46+num_ctrl][23:34].strip())
                case_res.stability.CY_q  = float(lines_2[46+num_ctrl][43:54].strip())
                case_res.stability.CY_r  = float(lines_2[46+num_ctrl][63:74].strip())
                case_res.stability.CZ_p  = float(lines_2[47+num_ctrl][23:34].strip())
                case_res.stability.CZ_q  = float(lines_2[47+num_ctrl][43:54].strip())
                case_res.stability.CZ_r  = float(lines_2[47+num_ctrl][63:74].strip())
                case_res.stability.Cl_p  = float(lines_2[48+num_ctrl][23:34].strip())
                case_res.stability.Cl_q  = float(lines_2[48+num_ctrl][43:54].strip())
                case_res.stability.Cl_r  = float(lines_2[48+num_ctrl][63:74].strip())
                case_res.stability.Cm_p  = float(lines_2[49+num_ctrl][23:34].strip())
                case_res.stability.Cm_q  = float(lines_2[49+num_ctrl][43:54].strip())
                case_res.stability.Cm_r  = float(lines_2[49+num_ctrl][63:74].strip())
                case_res.stability.Cn_p  = float(lines_2[50+num_ctrl][23:34].strip())
                case_res.stability.Cn_q  = float(lines_2[50+num_ctrl][43:54].strip())
                case_res.stability.Cn_r  = float(lines_2[50+num_ctrl][63:74].strip())
                
            # NILS: added to extract content of case.eigen_result_filename_1 (not previously implemented)
            with open(case.eigen_result_filename_1,'r') as eigval_file:
                # Extract results from eigenvalues file
            
                lines_eig = eigval_file.readlines()
            
                eig_real = []
                eig_imag = []
            
                for line in lines_eig:
                    # skip comments / empty lines
                    if line.strip() == '' or line.strip().startswith('#'):
                        continue
            
                    eig_real.append(float(line[8:26].strip()))
                    eig_imag.append(float(line[26:40].strip()))
            
                case_res.stability.eigenvalues_real = np.array(eig_real)
                case_res.stability.eigenvalues_imag = np.array(eig_imag)
                
            # NILS: added to extract content of case.eigen_result_filename_2 (not previously implemented)
            with open(case.eigen_result_filename_2,'r') as systmat_file:
                # Extract results from system matrix file
                
                lines_sys = systmat_file.readlines()
                
                # first non-empty line is header with variable names
                header_idx = 0
                while lines_sys[header_idx].strip() == '':
                    header_idx += 1
                
                header = lines_sys[header_idx].split()
                numeric_lines = lines_sys[header_idx+1:]
                
                system_matrix = []
                
                for line in numeric_lines:
                    if line.strip() == '':
                        continue
                    system_matrix.append([float(val) for val in line.split()])
                
                case_res.stability.system_matrix_labels = header
                case_res.stability.system_matrix = np.array(system_matrix)
            
        results.append(case_res)

    return results