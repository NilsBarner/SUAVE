## @ingroup Methods-Aerodynamics-AVL
#create_avl_datastructure.py
# 
# Created:  Oct 2014, T. Momose
# Modified: Jan 2016, E. Botero
#           Apr 2017, M. Clarke
#           Jul 2017, T. MacDonald
#           Aug 2019, M. Clarke
#           Mar 2020, M. Clarke

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import scipy
import numpy as np

from copy import deepcopy

# SUAVE Imports
from SUAVE.Core import Data , Units

# SUAVE-AVL Imports
from .Data.Inputs                                                  import Inputs
from .Data.Wing                                                    import Wing, Section, Control_Surface
from .Data.Body                                                    import Body
from .Data.Aircraft                                                import Aircraft
from .Data.Cases                                                   import Run_Case
from .Data.Configuration                                           import Configuration
from SUAVE.Components.Wings.Control_Surfaces                       import Aileron , Elevator , Slat , Flap , Rudder 
from SUAVE.Methods.Aerodynamics.AVL.write_avl_airfoil_file         import write_avl_airfoil_file  
from SUAVE.Methods.Geometry.Two_Dimensional.Planform.wing_planform import wing_planform

## @ingroup Methods-Aerodynamics-AVL
def translate_avl_wing(suave_wing):
    """ Translates wing geometry from the vehicle setup to AVL format

    Assumptions:
        None

    Source:
        None

    Inputs:
        suave_wing.tag                                                          [-]
        suave_wing.symmetric                                                    [boolean]
        suave_wing.verical                                                      [boolean]
        suave_wing - passed into the populate_wing_sections function            [data stucture]

    Outputs:
        w - aircraft wing in AVL format                                         [data stucture] 

    Properties Used:
        N/A
    """         
    w                 = Wing()
    w.tag             = suave_wing.tag
    w.symmetric       = suave_wing.symmetric
    w.vertical        = suave_wing.vertical
    w                 = populate_wing_sections(w,suave_wing)

    return w

def translate_avl_body(suave_body):
    """ Translates body geometry from the vehicle setup to AVL format

    Assumptions:
        None

    Source:
        None

    Inputs:
        body.tag                                                       [-]
        suave_wing.lengths.total                                       [meters]    
        suave_body.lengths.nose                                        [meters]
        suave_body.lengths.tail                                        [meters]
        suave_wing.verical                                             [meters]
        suave_body.width                                               [meters]
        suave_body.heights.maximum                                     [meters]
        suave_wing - passed into the populate_body_sections function   [data stucture]

    Outputs:
        b - aircraft body in AVL format                                [data stucture] 

    Properties Used:
        N/A
    """  
    b                 = Body()
    b.tag             = suave_body.tag
    b.symmetric       = True
    b.lengths.total   = suave_body.lengths.total
    b.lengths.nose    = suave_body.lengths.nose
    b.lengths.tail    = suave_body.lengths.tail
    b.widths.maximum  = suave_body.width
    b.heights.maximum = suave_body.heights.maximum
    # b                 = populate_body_sections(b,suave_body)
    b                 = populate_body_sections_nils(b,suave_body)  # NILS: refined fuselage contour shape

    return b

# NILS: added this function from scratch to support modelling nacelles with
# (turbofan/ducted fan) and without (turboprop/e-prop) flow_through = True
def translate_avl_nacelle_nils(suave_nacelle):
    """ Translates wing geometry from the vehicle setup to AVL format

    Assumptions:
        None

    Source:
        None

    Inputs:
        suave_wing.tag                                                          [-]
        suave_wing.symmetric                                                    [boolean]
        suave_wing.verical                                                      [boolean]
        suave_wing - passed into the populate_wing_sections function            [data stucture]

    Outputs:
        w - aircraft wing in AVL format                                         [data stucture] 

    Properties Used:
        N/A
    """         
    if suave_nacelle.flow_through == True:
        n                 = Wing()
    elif suave_nacelle.flow_through == False:
        n                 = Body()
    n.tag             = suave_nacelle.tag
    n.flow_through    = suave_nacelle.flow_through
    n.symmetric       = True  # NILS: must be True to model both nacelles in AVL
    if n.flow_through == True:
        n = populate_turbofan_nacelle_sections_nils(n, suave_nacelle)
    elif n.flow_through == False:
        n = populate_turboprop_nacelle_sections_nils(n, suave_nacelle)

    return n

def populate_wing_sections(avl_wing,suave_wing): 
    """ Creates sections of wing geometry and populates the AVL wing data structure

    Assumptions:
        None

    Source:
        None

    Inputs:
        avl_wing.symmetric                         [boolean]
        suave_wing.spans.projected                 [meters]
        suave_wing.origin                          [meters]
        suave_wing.dihedral                        [radians]
        suave_wing.Segments.sweeps.leading_edge    [radians]
        suave_wing.Segments.root_chord_percent     [-]
        suave_wing.Segments.percent_span_location  [-]
        suave_wing.Segments.sweeps.quarter_chord   [radians]
        suave_wing.Segment.twist                   [radians]

    Outputs:
        avl_wing - aircraft wing in AVL format     [data stucture] 

    Properties Used:
        N/A
    """         

    if len(suave_wing.Segments.keys())>0:
        # obtain the geometry for each segment in a loop                                            
        symm                 = avl_wing.symmetric
        semispan             = suave_wing.spans.projected*0.5 * (2 - symm)
        avl_wing.semispan    = semispan   
        root_chord           = suave_wing.chords.root
        segment_percent_span = 0    
        segments             = suave_wing.Segments
        n_segments           = len(segments.keys())
        segment_sweeps       = []
        origin               = []
        
        origin.append(suave_wing.origin)

        for i_segs in range(n_segments):
            if (i_segs == n_segments-1):
                segment_sweeps.append(0)                                  
            else: # this converts all sweeps defined by the quarter chord to leading edge sweep since AVL needs the start of each wing section
                #from the leading edge coordinate and not the quarter chord coordinate
                if segments[i_segs].sweeps.leading_edge is not None: 
                    # if leading edge sweep is defined 
                    segment_sweep       = segments[i_segs].sweeps.leading_edge  
                else:   
                    # if quarter chord sweep is defined, convert it to leading edge sweep
                    sweep_quarter_chord = segments[i_segs].sweeps.quarter_chord 
                    chord_fraction      = 0.25                          
                    segment_root_chord  = root_chord*segments[i_segs].root_chord_percent
                    segment_tip_chord   = root_chord*segments[i_segs+1].root_chord_percent
                    segment_span        = semispan*(segments[i_segs+1].percent_span_location - segments[i_segs].percent_span_location )
                    segment_sweep       = np.arctan(((segment_root_chord*chord_fraction) + (np.tan(sweep_quarter_chord )*segment_span - chord_fraction*segment_tip_chord)) /segment_span)
                segment_sweeps.append(segment_sweep)
            dihedral       = segments[i_segs].dihedral_outboard  
            ctrl_surf_at_seg = False 

            # condition for the presence of control surfaces in segment 
            if getattr(segments[i_segs],'control_surfaces',False):    
                dihedral_ob   = segments[i_segs-1].dihedral_outboard 
                section_spans = []
                for cs in segments[i_segs].control_surfaces:     
                    # create a vector if all the section breaks in a segment. sections include beginning and end of control surfaces and end of segment      
                    control_surface_start = semispan*cs.span_fraction_start
                    control_surface_end   = semispan*cs.span_fraction_end
                    section_spans.append(control_surface_start)
                    section_spans.append(control_surface_end)                                
                ordered_section_spans = sorted(list(set(section_spans)))     # sort the section_spans in order to create sections in spanwise order
                num_sections = len(ordered_section_spans)                    # count the number of sections breaks that the segment will contain    \

                for section_count in range(num_sections):        
                    # create and append sections onto avl wing structure  
                    if ordered_section_spans[section_count] == semispan*segments[i_segs-1].percent_span_location:  
                        # if control surface begins at beginning of segment, redundant section is removed
                        section_tags = list(avl_wing.sections.keys())
                        del avl_wing.sections[section_tags[-1]]

                    # create section for each break in the wing        
                    section                   = Section()              
                    section.tag               = segments[i_segs].tag + '_section_'+ str(ordered_section_spans[section_count]) + 'm'
                    root_section_chord        = root_chord*segments[i_segs-1].root_chord_percent
                    tip_section_chord         = root_chord*segments[i_segs].root_chord_percent
                    semispan_section_fraction = (ordered_section_spans[section_count] - semispan*segments[i_segs-1].percent_span_location)/(semispan*(segments[i_segs].percent_span_location - segments[i_segs-1].percent_span_location ))   
                    section.chord             = np.interp(semispan_section_fraction,[0.,1.],[root_section_chord,tip_section_chord])
                    root_section_twist        = segments[i_segs-1].twist/Units.degrees 
                    tip_section_twist         = root_chord*segments[i_segs].twist/Units.degrees  
                    section.twist             = np.interp(semispan_section_fraction,[0.,1.],[root_section_twist,tip_section_twist]) 

                    # if wing is a vertical wing, the y and z coordinates are swapped 
                    if avl_wing.vertical:
                        dz = ordered_section_spans[section_count] -  semispan*segments[i_segs-1].percent_span_location 
                        dy = dz*np.tan(dihedral_ob)
                        l  = dz/np.cos(dihedral_ob)
                        dx = l*np.tan(segment_sweeps[i_segs-1])                                                            
                    else:
                        dy = ordered_section_spans[section_count] - semispan*segments[i_segs-1].percent_span_location 
                        dz = dy*np.tan(dihedral_ob)
                        l  = dy/np.cos(dihedral_ob)
                        dx = l*np.tan(segment_sweeps[i_segs-1])
                    section.origin = [[origin[i_segs-1][0][0] + dx , origin[i_segs-1][0][1] + dy, origin[i_segs-1][0][2] + dz]]              

                    # this loop appends all the control surfaces within a particular wing section
                    for index  , ctrl_surf in enumerate(segments[i_segs].control_surfaces):
                        if  (semispan*ctrl_surf.span_fraction_start == ordered_section_spans[section_count]) or \
                                                    (ordered_section_spans[section_count] == semispan*ctrl_surf.span_fraction_end):
                            c                     = Control_Surface()
                            c.tag                 = ctrl_surf.tag                # name of control surface   
                            c.sign_duplicate      = '+1'                         # this float indicates control surface deflection symmetry
                            c.x_hinge             = 1 - ctrl_surf.chord_fraction # this float is the % location of the control surface hinge on the wing 
                            c.deflection          = ctrl_surf.deflection / Units.degrees 
                            c.order               = index

                            # if control surface is an aileron, the deflection is asymmetric. This is standard convention from AVL
                            if (type(ctrl_surf) ==  Aileron):
                                c.sign_duplicate = '-1'
                                c.function       = 'aileron'
                                c.gain           = -1.0
                            # if control surface is a slat, the hinge is taken from the leading edge        
                            elif (type(ctrl_surf) ==  Slat):
                                c.x_hinge   =  -ctrl_surf.chord_fraction
                                c.function  = 'slat'
                                c.gain      = -1.0
                            elif (type(ctrl_surf) ==  Flap):
                                c.function  = 'flap'    
                                c.gain      = 1.0
                            elif (type(ctrl_surf) ==  Elevator):
                                c.function  = 'elevator'
                                c.gain      = 1.0
                            elif (type(ctrl_surf) ==  Rudder):
                                c.function  = 'rudder'
                                c.gain      = 1.0
                            else:
                                raise AttributeError("Define control surface function as 'slat', 'flap', 'elevator' , 'aileron' or 'rudder'")
                            section.append_control_surface(c) 

                        elif  (semispan*ctrl_surf.span_fraction_start < ordered_section_spans[section_count]) and \
                                                      (ordered_section_spans[section_count] < semispan*ctrl_surf.span_fraction_end):
                            c                     = Control_Surface()
                            c.tag                 = ctrl_surf.tag                # name of control surface   
                            c.sign_duplicate      = '+1'                         # this float indicates control surface deflection symmetry
                            c.x_hinge             = 1 - ctrl_surf.chord_fraction # this float is the % location of the control surface hinge on the wing 
                            c.deflection          = ctrl_surf.deflection / Units.degrees 
                            c.order               = index

                            # if control surface is an aileron, the deflection is asymmetric. This is standard convention from AVL
                            if (type(ctrl_surf) ==  Aileron):
                                c.sign_duplicate = '-1'
                                c.function       = 'aileron'
                                c.gain           = -1.0
                            # if control surface is a slat, the hinge is taken from the leading edge        
                            elif (type(ctrl_surf) ==  Slat):
                                c.x_hinge   =  -ctrl_surf.chord_fraction
                                c.function  = 'slat'
                                c.gain      = -1.0
                            elif (type(ctrl_surf) ==  Flap):
                                c.function  = 'flap'    
                                c.gain      = 1.0
                            elif (type(ctrl_surf) ==  Elevator):
                                c.function  = 'elevator'
                                c.gain      = 1.0
                            elif (type(ctrl_surf) ==  Rudder):
                                c.function  = 'rudder'
                                c.gain      = 1.0
                            else:
                                raise AttributeError("Define control surface function as 'slat', 'flap', 'elevator' , 'aileron' or 'rudder'")
                            section.append_control_surface(c)                                                  

                    if segments[i_segs].Airfoil:
                        if segments[i_segs].Airfoil.airfoil.coordinate_file is not None:
                            section.airfoil_coord_file   = write_avl_airfoil_file(segments[i_segs].Airfoil.airfoil.coordinate_file)
                        elif segments[i_segs].Airfoil.airfoil.naca_airfoil is not None:
                            section.naca_airfoil         = segments[i_segs].Airfoil.airfoil.naca_airfoil 

                    avl_wing.append_section(section)   

                # check if control surface ends at end of segment         
                if ordered_section_spans[section_count] == semispan*segments[i_segs].percent_span_location:  
                    ctrl_surf_at_seg = True

            if ctrl_surf_at_seg:  # if a control surface ends at the end of the segment, there is not need to append another segment
                pass
            else: # if there is no control surface break at the end of the segment, this block appends a segment
                section        = Section() 
                section.tag    = segments[i_segs].tag
                section.chord  = root_chord*segments[i_segs].root_chord_percent 
                section.twist  = segments[i_segs].twist/Units.degrees    
                section.origin = origin[i_segs]
                if segments[i_segs].Airfoil:
                    if segments[i_segs].Airfoil.airfoil.coordinate_file is not None:
                        section.airfoil_coord_file   = write_avl_airfoil_file(segments[i_segs].Airfoil.airfoil.coordinate_file)
                    elif segments[i_segs].Airfoil.airfoil.naca_airfoil is not None:
                        section.naca_airfoil         = segments[i_segs].Airfoil.airfoil.naca_airfoil     
                # append section to wing
                avl_wing.append_section(section)                               

            # update origin for next segment
            if (i_segs == n_segments-1):                                          
                return avl_wing

            segment_percent_span =    segments[i_segs+1].percent_span_location - segments[i_segs].percent_span_location     
            if avl_wing.vertical:
                dz = semispan*segment_percent_span
                dy = dz*np.tan(dihedral)
                l  = dz/np.cos(dihedral)
                dx = l*np.tan(segment_sweep)
            else:
                dy = semispan*segment_percent_span
                dz = dy*np.tan(dihedral)
                l  = dy/np.cos(dihedral)
                dx = l*np.tan(segment_sweep)
            origin.append( [[origin[i_segs][0][0] + dx , origin[i_segs][0][1] + dy, origin[i_segs][0][2] + dz]])               

    else:    
        symm                  = avl_wing.symmetric  
        dihedral              = suave_wing.dihedral
        span                  = suave_wing.spans.projected
        semispan              = suave_wing.spans.projected * 0.5 * (2 - symm) 
        if suave_wing.sweeps.leading_edge  is not None: 
            sweep      = suave_wing.sweeps.leading_edge
        else: 
            suave_wing = wing_planform(suave_wing)
            sweep      = suave_wing.sweeps.leading_edge
        avl_wing.semispan     = semispan
        origin                = suave_wing.origin[0]  
        
        # define root section 
        root_section          = Section()
        root_section.tag      = 'root_section'
        root_section.origin   = [origin]
        root_section.chord    = suave_wing.chords.root 
        root_section.twist    = suave_wing.twists.root/Units.degrees 
        root_section.semispan  = semispan

        # define tip section
        tip_section           = Section()
        tip_section.tag       = 'tip_section'
        tip_section.chord     = suave_wing.chords.tip 
        tip_section.twist     = suave_wing.twists.tip/Units.degrees 
        tip_section.semispan  = 0

        # assign location of wing tip         
        if avl_wing.vertical:
            tip_section.origin    = [[origin[0]+semispan*np.tan(sweep),origin[1]+semispan*np.tan(dihedral),origin[2]+semispan]]
        else: 
            tip_section.origin    = [[origin[0]+semispan*np.tan(sweep),origin[1]+semispan,origin[2]+semispan*np.tan(dihedral)]]

        # assign wing airfoil
        if suave_wing.Airfoil:
            root_section.airfoil_coord_file  = suave_wing.Airfoil.airfoil.coordinate_file          
            tip_section.airfoil_coord_file   = suave_wing.Airfoil.airfoil.coordinate_file    


        avl_wing.append_section(root_section)
        avl_wing.append_section(tip_section)

    return avl_wing

# NILS: added this function from scratch to support modelling nacelles with
# (turbofan/ducted fan) and without (turboprop/e-prop) flow_through = True
def populate_turbofan_nacelle_sections_nils(avl_nacelle, suave_nacelle):
    """ Creates sections of wing geometry and populates the AVL wing data structure

    Assumptions:
        None

    Source:
        None

    Inputs:
        avl_nacelle.symmetric                         [boolean]
        suave_nacelle.spans.projected                 [meters]
        suave_nacelle.origin                          [meters]
        suave_nacelle.dihedral                        [radians]
        suave_nacelle.Segments.sweeps.leading_edge    [radians]
        suave_nacelle.Segments.root_chord_percent     [-]
        suave_nacelle.Segments.percent_span_location  [-]
        suave_nacelle.Segments.sweeps.quarter_chord   [radians]
        suave_nacelle.Segment.twist                   [radians]

    Outputs:
        avl_nacelle - aircraft wing in AVL format     [data stucture] 

    Properties Used:
        N/A
    """
    
    # obtain the geometry for each segment in a loop                                            
    semispan             = np.pi * suave_nacelle.diameter / 2  # suave_nacelle.spans.projected*0.5 * (2 - symm)
    avl_nacelle.semispan = semispan   
    root_chord           = suave_nacelle.length  # suave_nacelle.chords.root
    origin               = []  # NILS: kept for now as need for future studies with neng > 2
    origin.append(suave_nacelle.origin)
    
    diameter = suave_nacelle.diameter
    radius = diameter * 0.5
    
    # Compute section origin by wrapping around circle
    # X is fixed: X_le is origin.x
    X0 = origin[0][0]
    Yc = origin[0][1]
    Zc = origin[0][2]
    
    avl_nacelle.x_scale = 1.0
    avl_nacelle.y_scale = radius
    avl_nacelle.z_scale = radius
    avl_nacelle.x_transl = X0
    avl_nacelle.y_transl = Yc
    avl_nacelle.z_transl = Zc
    
    # Number of unique points around the circle
    n_segments = 12      # not 13
    
    # Start at top (π/2), go clockwise (negative direction), avoid duplicate wrap
    thetas = np.linspace(0, 2 * np.pi, n_segments, endpoint=False)
    
    # Rotate so 0 → π/2
    thetas = (np.pi/2 - thetas) % (2*np.pi)
    
    for i_segs in range(n_segments):
        
        # Circular coordinates
        theta = thetas[i_segs]
        dy = np.cos(theta)
        dz = np.sin(theta)
        
        # if there is no control surface break at the end of the segment, this block appends a segment
        section        = Section() 
        section.tag    = f'segment_{i_segs}'
        section.chord  = root_chord
        section.twist  = 0.0
        section.origin = [[0.0, dy, dz]]
        
        # NILS: use NACA airfoil and append to section (same for all)
        section.naca_airfoil = suave_nacelle.Airfoil.airfoil.naca_4_series_airfoil
        avl_nacelle.append_section(section)            

        # update origin for next segment
        if (i_segs == n_segments-1):                                          
            return avl_nacelle

    return avl_nacelle

# =============================================================================
# NILS: added this function from scratch to support modelling nacelles with
# (turbofan/ducted fan) and without (turboprop/e-prop) flow_through = True
def populate_turboprop_nacelle_sections_nils(avl_nacelle, suave_nacelle):
    """ Creates sections of wing geometry and populates the AVL wing data structure

    Assumptions:
        None

    Source:
        None

    Inputs:
        avl_nacelle.symmetric                         [boolean]
        suave_nacelle.spans.projected                 [meters]
        suave_nacelle.origin                          [meters]
        suave_nacelle.dihedral                        [radians]
        suave_nacelle.Segments.sweeps.leading_edge    [radians]
        suave_nacelle.Segments.root_chord_percent     [-]
        suave_nacelle.Segments.percent_span_location  [-]
        suave_nacelle.Segments.sweeps.quarter_chord   [radians]
        suave_nacelle.Segment.twist                   [radians]

    Outputs:
        avl_nacelle - aircraft wing in AVL format     [data stucture] 

    Properties Used:
        N/A
    """
    
    symm = avl_nacelle.symmetric   
    semispan_h = suave_nacelle.diameter * 0.5 * (2 - symm)
    semispan_v = suave_nacelle.diameter * 0.5
    origin = suave_nacelle.origin
    
    # from nils.parametric_geometry import generate_streamlined_body_geometry
    
    import sys
    import numpy as np
    from scipy.integrate import simpson
    from typing import Annotated, Tuple
    from scipy.optimize import minimize
    from scipy.interpolate import Akima1DInterpolator
    
    def generate_streamlined_body_geometry(
        R_le_over_c: Annotated[float, "[-]"],
        beta_tail: Annotated[float, "[deg]"],
        Psi_zeta_max: Annotated[float, "[-]"], # this value does not fully work as expected - to be investigated once less pressing stuff has been dealt with
        zeta_max: Annotated[float, "[-]"],
        zeta_te: Annotated[float, "[-]"],
        dimensional_known_dict: Annotated[dict, "[m] or [m^3]"] # length or volume
    ) -> Tuple[float, float, float, float, float, float]:
        """
        Generate nacelle geometry from airfoil-shaped body of revolution based on low_paper_2008_univparamgeomreprmeth_kulfan
        by either specifying the volume and fineness ratio or the length and fineness ratio (diameter).
        """
        
        def _define_quadratic(x0, y0, x1, y1):
            """
            Generate second-order polynomial tangent to point (x0, y0)
            and secant through point (x1, y1).
            """
            a = (y1 - y0) / (x1 - x0)**2
            b = -2 * a * x0
            c = y0 + a * x0**2
            return np.poly1d([a, b, c])
        
        # Fixed parameters
        N_1, N_2 = 0.5, 1 # for round-nose airfoil
        Psi_le = 0
        Psi_te = 1
        
        l_over_d = 1/zeta_max/2 # concerns HALF of airfoil (see definition of zeta_max on page 2 and in Fig. 1 of low_paper_2008_univparamgeomreprmeth_kulfan)
        
        Vol = None
        x_list_sorted = None
        y_list_sorted = None
        Psi_list = None
        zeta_list = None
        zeta_list_sorted = None
        
        def _calculate_volume_from_length(l):
            nonlocal Vol, y_list_sorted, x_list_sorted, Psi_list, zeta_list, zeta_list_sorted
            
            d_max = l[0]/zeta_max
            
            S_le = np.sqrt(2 * R_le_over_c) # (4)
            S_te = np.tan(beta_tail * np.pi/180) + zeta_te # (5)
            S_zeta_max = (zeta_max - Psi_zeta_max * zeta_te) / (np.sqrt(Psi_zeta_max) * (1 - Psi_zeta_max)) # (2)

            quadratic_front = _define_quadratic(Psi_zeta_max, S_zeta_max, Psi_le, S_le)
            quadratic_aft = _define_quadratic(Psi_zeta_max, S_zeta_max, Psi_te, S_te)
            
            zeta_list = []
            Psi_list = np.linspace(1e-4, 1, 100)
            
            for Psi in Psi_list:
                if Psi >= Psi_zeta_max:
                    S = quadratic_aft(Psi)
                elif Psi < Psi_zeta_max:
                    S = quadratic_front(Psi)
                C = Psi**N_1 * (1 - Psi)**N_2 # (6)
                zeta = C * S + Psi * zeta_te # (7)
                zeta_list.append(zeta)
            
            # Sort points in ascending x-order for integration
            sorted_indices = np.argsort(Psi_list)
            Psi_list_sorted = Psi_list[sorted_indices]
            zeta_list_sorted = np.array(zeta_list)[sorted_indices]
            
            # Go from nondimensional coordinates to dimensional lengths
            x_list_sorted = Psi_list_sorted * l[0]
            y_list_sorted = zeta_list_sorted * l[0]  # zeta_max refers to max. ndim. radius, hence x2
            
            # Calculate the volume of revolution around the x-axis using the disk method: V = pi * integral( y^2 dx )
            Vol = np.pi * simpson(y_list_sorted**2, x=x_list_sorted)
            
            f = abs(Vol - _Vol)
            
            return f
        
        if list(dimensional_known_dict.keys())[0] == 'length':
            _Vol = np.inf # dummy value
            l = list(dimensional_known_dict.values())[0]
            _calculate_volume_from_length([l]) # input must be of type list
        elif list(dimensional_known_dict.keys())[0] == 'volume':
            _Vol = list(dimensional_known_dict.values())[0]
            bnds = [(0, np.inf)]
            x0 = (_Vol * 4 * l_over_d**2 / np.pi)**(1/3) # use cylinder as initial guess
            outputs = minimize(
                _calculate_volume_from_length, x0=x0, method='SLSQP', bounds=bnds
            )
            l = outputs.x[0]
        
        # Calculate the surface area of revolution around the x-axis: S = 2 * pi * integral(f * sqrt(1 + fprime**2))
        dy_dx_list = np.gradient(y_list_sorted, x_list_sorted)
        integrand = 2 * np.pi * y_list_sorted * np.sqrt(1 + dy_dx_list**2)
        surface_area = simpson(integrand, x=x_list_sorted)
        
        # Calculate the surface area up to the point of maximum thickness
        idx_zeta_max = np.argmin(abs(zeta_list_sorted - zeta_max))
        y_list_masked = y_list_sorted[:idx_zeta_max + 1]
        x_list_masked = x_list_sorted[:idx_zeta_max + 1]
        dy_dx_list_masked = np.gradient(y_list_masked, x_list_masked)
        integrand_masked = 2 * np.pi * y_list_masked * np.sqrt(1 + dy_dx_list_masked**2)
        surface_area_le_to_zeta_max = simpson(integrand_masked, x=x_list_masked)
        
        # Calculate the frontal area
        frontal_area = np.pi * (zeta_max * l)**2
        
        # Calculate the circumference (Euclidean distance between consecutive points) in 2D
        circumference = np.sum(np.sqrt(np.diff(x_list_sorted)**2 + np.diff(y_list_sorted)**2)) * 2 # x2 since there is an upper and a lower half of the airfoil
        
        # Create closed area from half-airfoil contour
        Psi_list_closed = np.append(Psi_list, np.flip(Psi_list))
        zeta_list_closed = np.append(np.array(zeta_list), -1 * np.flip(np.array(zeta_list)))
        
        # Check: plot side-view of nacelle contour
        # fig, ax = plt.subplots()
        # ax.plot(Psi_list_closed, zeta_list_closed, color = 'black')
        # ax.set_aspect('equal')
        # ax.set_xlabel(r'$\Psi=x/c$')
        # ax.set_ylabel(r'$\zeta=z/c$')
        # ax.set_xlim(0, 1)
        # plt.show()
        
        return Vol, l, circumference, frontal_area, surface_area, Psi_list_closed, zeta_list_closed, surface_area_le_to_zeta_max
        # return Vol, l, circumference, frontal_area, surface_area, Psi_list, zeta_list, surface_area_le_to_zeta_max
    
    zeta_max = suave_nacelle.diameter / suave_nacelle.length / 2
    
    Vol, l, _, _, surface_area, Psi_list_closed, zeta_list_closed, _ = \
        generate_streamlined_body_geometry(
            R_le_over_c=0.06, beta_tail=30, Psi_zeta_max=0.35, zeta_max=zeta_max, zeta_te=0,
            dimensional_known_dict={'length': suave_nacelle.length}
        )
    
    x_fore = np.hstack((
        Psi_list_closed[np.nanargmin(zeta_list_closed):],
        Psi_list_closed[:np.nanargmax(zeta_list_closed) + 1],  # NILS: add `+ 1` to ensure that Akima1DInterpolator does not return nan when queried at boundary of data
    )) * suave_nacelle.length
    y_fore = np.hstack((
        zeta_list_closed[np.nanargmin(zeta_list_closed):],
        zeta_list_closed[:np.nanargmax(zeta_list_closed) + 1],  # NILS: add `+ 1` to ensure that Akima1DInterpolator does not return nan when queried at boundary of data
    )) * suave_nacelle.length
    
    x_aft = Psi_list_closed[np.nanargmax(zeta_list_closed) : np.nanargmin(zeta_list_closed) + 1][::-1] * suave_nacelle.length  # NILS: add `+ 1` to ensure that Akima1DInterpolator does not return nan when queried at boundary of data
    y_aft = zeta_list_closed[np.nanargmax(zeta_list_closed) : np.nanargmin(zeta_list_closed) + 1][::-1] * suave_nacelle.length  # NILS: add `+ 1` to ensure that Akima1DInterpolator does not return nan when queried at boundary of data
    
    # Remove duplicate y-values while preserving order
    _, idx = np.unique(y_aft, return_index=True)
    idx = np.sort(idx)
    x_aft = x_aft[idx]
    y_aft = y_aft[idx]
    
    origin = []  # NILS: kept for now as need for future studies with neng > 2
    origin.append(suave_nacelle.origin)
    
    diameter = suave_nacelle.diameter
    radius = diameter * 0.5
    
    # Compute section origin by wrapping around circle
    # X is fixed: X_le is origin.x
    X0 = origin[0][0]
    Yc = origin[0][1]
    Zc = origin[0][2]
    
    avl_nacelle.x_scale = 1.0
    avl_nacelle.y_scale = radius
    avl_nacelle.z_scale = radius
    avl_nacelle.x_transl = X0
    avl_nacelle.y_transl = Yc
    avl_nacelle.z_transl = Zc
    
    x_makima_fore_interp = Akima1DInterpolator(
        y_fore,
        x_fore,
        method="makima",
        # extrapolate=True,  # NILS: not available in my version of SciPy (added in version 1.13.0 - see https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.Akima1DInterpolator.html)
    )
    x_makima_aft_interp = Akima1DInterpolator(
        y_aft,
        x_aft,
        method="makima",
        # extrapolate=True,  # NILS: not available in my version of SciPy (added in version 1.13.0 - see https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.Akima1DInterpolator.html)
    )
    
    # Horizontal Sections of Fuselage
    if semispan_h != 0.0:                
        width_array = np.linspace(-semispan_h, semispan_h, num=11,endpoint=True)
        for section_width in width_array:
            nacelle_h_section               = Section()
            
            nacelle_h_section_fore_length   = np.nanmax(x_fore) - x_makima_fore_interp(section_width)
            nacelle_h_section_aft_length    = x_makima_aft_interp(section_width) - np.nanmax(x_fore)
            nacelle_h_section_nose_origin   = np.nanmax(x_fore) - nacelle_h_section_fore_length
            
            nacelle_h_section.tag           =  'nacelle_horizontal_section_at_' +  str(section_width) + '_m'
            nacelle_h_section.origin        = [ origin[0][0] + nacelle_h_section_nose_origin , origin[0][1] + section_width, origin[0][2]]
            nacelle_h_section.chord         = nacelle_h_section_fore_length + nacelle_h_section_aft_length
            
            avl_nacelle.append_section(nacelle_h_section,'horizontal')
            
    # Vertical Sections of Fuselage 
    if semispan_v != 0:               
        height_array = np.linspace(-semispan_v, semispan_v, num=11,endpoint=True)
        for section_height in height_array :
            nacelle_v_section               = Section()
            
            nacelle_v_section_fore_length   = np.nanmax(x_fore) - x_makima_fore_interp(section_height)
            nacelle_v_section_aft_length    = x_makima_aft_interp(section_height) - np.nanmax(x_fore)
            nacelle_v_section_nose_origin   = np.nanmax(x_fore) - nacelle_v_section_fore_length
            
            nacelle_v_section.tag           = 'nacelle_vertical_top_section_at_' +  str(section_height) + '_m'        
            nacelle_v_section.origin        = [ origin[0][0] + nacelle_v_section_nose_origin,  origin[0][1],  origin[0][2] + section_height ]
            nacelle_v_section.chord         = nacelle_v_section_fore_length + nacelle_v_section_aft_length
            
            avl_nacelle.append_section(nacelle_v_section,'vertical')
            
    return avl_nacelle
# =============================================================================

def populate_body_sections(avl_body,suave_body):
    """ Creates sections of body geometry and populates the AVL body data structure

    Assumptions:
        None

    Source:
        None

    Inputs:
        avl_wing.symmetric                       [boolean]
        avl_body.widths.maximum                  [meters]
        avl_body.heights.maximum                 [meters]
        suave_body.fineness.nose                 [meters]
        suave_body.fineness.tail                 [meters]
        avl_body.lengths.total                   [meters]
        avl_body.lengths.nose                    [meters] 
        avl_body.lengths.tail                    [meters]  

    Outputs:
        avl_body - aircraft body in AVL format   [data stucture] 

    Properties Used:
        N/A
    """  

    symm = avl_body.symmetric   
    semispan_h = avl_body.widths.maximum * 0.5 * (2 - symm)
    semispan_v = avl_body.heights.maximum * 0.5
    origin = suave_body.origin[0]

    # Compute the curvature of the nose/tail given fineness ratio. Curvature is derived from general quadratic equation
    # This method relates the fineness ratio to the quadratic curve formula via a spline fit interpolation
    vec1 = [2 , 1.5, 1.2 , 1]
    vec2 = [1  ,1.57 , 3.2,  8]
    x = np.linspace(0,1,4)
    fuselage_nose_curvature =  np.interp(np.interp(suave_body.fineness.nose,vec2,x), x , vec1)
    fuselage_tail_curvature =  np.interp(np.interp(suave_body.fineness.tail,vec2,x), x , vec1) 


    # Horizontal Sections of Fuselage
    if semispan_h != 0.0:                
        width_array = np.linspace(-semispan_h, semispan_h, num=11,endpoint=True)
        for section_width in width_array:
            fuselage_h_section               = Section()
            fuselage_h_section_cabin_length  = avl_body.lengths.total - (avl_body.lengths.nose + avl_body.lengths.tail)
            fuselage_h_section_nose_length   = ((1 - ((abs(section_width/semispan_h))**fuselage_nose_curvature ))**(1/fuselage_nose_curvature))*avl_body.lengths.nose
            fuselage_h_section_tail_length   = ((1 - ((abs(section_width/semispan_h))**fuselage_tail_curvature ))**(1/fuselage_tail_curvature))*avl_body.lengths.tail
            fuselage_h_section_nose_origin   = avl_body.lengths.nose - fuselage_h_section_nose_length
            fuselage_h_section.tag           =  'fuselage_horizontal_section_at_' +  str(section_width) + '_m'
            fuselage_h_section.origin        = [ origin[0] + fuselage_h_section_nose_origin , origin[1] + section_width, origin[2]]
            fuselage_h_section.chord         = fuselage_h_section_cabin_length + fuselage_h_section_nose_length + fuselage_h_section_tail_length
            avl_body.append_section(fuselage_h_section,'horizontal')

    # Vertical Sections of Fuselage 
    if semispan_v != 0:               
        height_array = np.linspace(-semispan_v, semispan_v, num=11,endpoint=True)
        for section_height in height_array :
            fuselage_v_section               = Section()
            fuselage_v_section_cabin_length  = avl_body.lengths.total - (avl_body.lengths.nose + avl_body.lengths.tail)
            fuselage_v_section_nose_length   = ((1 - ((abs(section_height/semispan_v))**fuselage_nose_curvature ))**(1/fuselage_nose_curvature))*avl_body.lengths.nose
            fuselage_v_section_tail_length   = ((1 - ((abs(section_height/semispan_v))**fuselage_tail_curvature ))**(1/fuselage_tail_curvature))*avl_body.lengths.tail
            fuselage_v_section_nose_origin   = avl_body.lengths.nose - fuselage_v_section_nose_length
            fuselage_v_section.tag           = 'fuselage_vertical_top_section_at_' +  str(section_height) + '_m'        
            fuselage_v_section.origin        = [ origin[0] + fuselage_v_section_nose_origin,  origin[1],  origin[2] + section_height ]
            fuselage_v_section.chord         = fuselage_v_section_cabin_length + fuselage_v_section_nose_length + fuselage_v_section_tail_length
            avl_body.append_section(fuselage_v_section,'vertical')

    return avl_body

# ====================================== NILS =======================================
def populate_body_sections_nils(avl_body,suave_body):
    """ Creates sections of body geometry and populates the AVL body data structure

    Assumptions:
        None

    Source:
        None

    Inputs:
        avl_wing.symmetric                       [boolean]
        avl_body.widths.maximum                  [meters]
        avl_body.heights.maximum                 [meters]
        suave_body.fineness.nose                 [meters]
        suave_body.fineness.tail                 [meters]
        avl_body.lengths.total                   [meters]
        avl_body.lengths.nose                    [meters] 
        avl_body.lengths.tail                    [meters]  

    Outputs:
        avl_body - aircraft body in AVL format   [data stucture] 

    Properties Used:
        N/A
    """  

    symm = avl_body.symmetric   
    semispan_h = avl_body.widths.maximum * 0.5 * (2 - symm)
    semispan_v = avl_body.heights.maximum * 0.5
    origin = suave_body.origin[0]
    
    # =============================================================================
    import copy
    from scipy.interpolate import Akima1DInterpolator
    
    N_points_long = 500
    N_points_vert = 50
    N_points_horz = 50
    r_fuse = semispan_h
    x = np.array([
        0,
        avl_body.lengths.nose,
        (avl_body.lengths.nose + (avl_body.lengths.total - avl_body.lengths.tail)) / 2,
        avl_body.lengths.total - avl_body.lengths.tail,
        avl_body.lengths.total,
    ])
    z = np.array([-r_fuse / 2, 0, 0, 0, r_fuse / 2])
    x_nose_start = x[0]
    x_nose_end = x[1]
    x_tail_start = x[-2]
    x_tail_end = x[-1]
    z_centreline = z[2]
    xs = np.linspace(min(x), max(x), num=N_points_long)
    xs_nose = copy.deepcopy(xs)[xs < x_nose_end]
    xs_tail = copy.deepcopy(xs)[xs > x_tail_start]
    nnose = len(xs_nose)
    ntail = len(xs_tail)
    
    # Horizontal fuselage cross-section

    # Nose left

    a_nose_left = 1.5

    x_left_nose = []
    y_left_nose = []
    for i, _x in enumerate(xs_nose, start=1):
        fraci = (i - 1) / (nnose - 1)
        fracx = np.cos(0.5 * np.pi * fraci)
        x_left_nose.append(x_nose_end + (x_nose_start - x_nose_end) * fracx)
        y_left_nose.append(r_fuse * (1.0 - fracx**a_nose_left)**(1.0 / a_nose_left))

    y_makima_nose = np.linspace(-r_fuse, r_fuse, num=N_points_vert,endpoint=True)
    x_makima_horz_nose_interp = Akima1DInterpolator(
        np.concatenate((-np.array(y_left_nose[::-1][:-1]), y_left_nose)),
        np.concatenate((x_left_nose[::-1][:-1], x_left_nose)),
        method="makima",
    )#(y_makima_nose)

    # Tail left

    btail = 2

    x_left_tail = []
    y_left_tail = []
    for i, _x in enumerate(xs_tail, start=1):
        fraci = (i - 1) / (ntail - 1)
        fracx = np.cos(0.5 * np.pi * fraci)
        x_left_tail.append(x_tail_start + (x_tail_end - x_tail_start) * fracx)
        y_left_tail.append(r_fuse + (-0.7 * r_fuse) * fracx**btail)
        
    y_makima_tail = np.linspace(-r_fuse, r_fuse, num=N_points_vert,endpoint=True)
    x_makima_horz_tail_interp = Akima1DInterpolator(
        np.concatenate((
            -np.array(y_left_tail[::-1][:-1]),
            [-y_left_tail[0], 0, y_left_tail[0]],
            y_left_tail[1:],
        )),
        np.concatenate((
            x_left_tail[::-1][:-1],
            [x_left_tail[0], x_left_tail[0], x_left_tail[0]],
            x_left_tail[1:],
        )),
        method="makima",
    )#(y_makima_tail)
    
    # Vertical fuselage cross-section
    
    # z_akima = Akima1DInterpolator(x, z, method="akima")(xs)
    z_makima = Akima1DInterpolator(x, z, method="makima")(xs)
    z_offset_nose = z_makima[0]
    z_offset_tail = z_makima[-1]

    # Nose top

    a_nose_top = 1.5
    r_nose_top = r_fuse + z[1] - z[0]

    x_top_nose = []
    z_top_nose = []
    for i, _x in enumerate(xs_nose, start=1):
        fraci = (i - 1) / (nnose - 1)
        fracx = np.cos(0.5 * np.pi * fraci)
        x_top_nose.append(x_nose_end + (x_nose_start - x_nose_end) * fracx)
        z_top_nose.append(z_offset_nose + (r_nose_top * (1.0 - fracx**a_nose_top)**(1.0 / a_nose_top)))

    # Nose bottom

    a_nose_bottom = 3
    r_nose_bottom = r_fuse - (z[1] - z[0])

    x_bottom_nose = []
    z_bottom_nose = []
    for i, _x in enumerate(xs_nose, start=1):
        fraci = (i - 1) / (nnose - 1)
        fracx = np.cos(0.5 * np.pi * fraci)
        x_bottom_nose.append(x_nose_end + (x_nose_start - x_nose_end) * fracx)
        z_bottom_nose.append(z_offset_nose - (r_nose_bottom * (1.0 - fracx**a_nose_bottom)**(1.0 / a_nose_bottom)))

    z_makima_nose = np.linspace(-r_fuse, r_fuse, num=N_points_vert,endpoint=True)
    x_makima_vert_nose_interp = Akima1DInterpolator(
        np.concatenate((z_bottom_nose[::-1][:-1], z_top_nose)),
        np.concatenate((x_bottom_nose[::-1][:-1], x_top_nose)),
        method="makima",
    )#(z_makima_nose)
    # sys.exit()

    # Tail top

    b_tail_top = 3
    r_tail_top = r_fuse - (z[-1] - z[-2])

    x_top_tail = []
    z_top_tail = []
    for i, _x in enumerate(xs_tail, start=1):
        fraci = (i - 1) / (ntail - 1)
        fracx = np.cos(0.5 * np.pi * fraci)
        x_top_tail.append(x_tail_start + (x_tail_end - x_tail_start) * fracx)
        z_top_tail.append(z_offset_tail + (r_tail_top + (-0.4*r_tail_top) * fracx**b_tail_top))

    # Tail bottom

    b_tail_bottom = 2
    r_tail_bottom = r_fuse + (z[-1] - z[-2])

    x_bottom_tail = []
    z_bottom_tail = []
    for i, _x in enumerate(xs_tail, start=1):
        fraci = (i - 1) / (ntail - 1)
        fracx = np.cos(0.5 * np.pi * fraci)
        x_bottom_tail.append(x_tail_start + (x_tail_end - x_tail_start) * fracx)
        z_bottom_tail.append(z_offset_tail - (r_tail_bottom + (-0.9 * r_tail_bottom) * fracx**b_tail_bottom))
        
    z_makima_tail = np.linspace(-r_fuse, r_fuse, num=N_points_horz,endpoint=True)
    x_makima_vert_tail_interp = Akima1DInterpolator(
        np.concatenate((
            z_bottom_tail[::-1][:-1],
            [z_bottom_tail[0], (z_bottom_tail[0] + z_top_tail[0]) / 2, z_top_tail[0]],
            z_top_tail[1:],
        )),
        np.concatenate((
            x_bottom_tail[::-1][:-1],
            [x_bottom_tail[0], (x_bottom_tail[0] + x_top_tail[0]) / 2, x_top_tail[0]],
            x_top_tail[1:],
        )),
        method="makima",
    )#(z_makima_tail)
    
    # =============================================================================

    # Horizontal Sections of Fuselage
    if semispan_h != 0.0:                
        width_array = np.linspace(-semispan_h, semispan_h, num=11,endpoint=True)
        for section_width in width_array:
            fuselage_h_section               = Section()
            
            fuselage_h_section_cabin_length  = avl_body.lengths.total - (avl_body.lengths.nose + avl_body.lengths.tail)
            fuselage_h_section_nose_length   = avl_body.lengths.nose - x_makima_horz_nose_interp(section_width)
            fuselage_h_section_tail_length   = x_makima_horz_tail_interp(section_width) - (avl_body.lengths.total - avl_body.lengths.tail)
            fuselage_h_section_nose_origin   = avl_body.lengths.nose - fuselage_h_section_nose_length
            
            fuselage_h_section.tag           =  'fuselage_horizontal_section_at_' +  str(section_width) + '_m'
            fuselage_h_section.origin        = [ origin[0] + fuselage_h_section_nose_origin , origin[1] + section_width, origin[2]]
            fuselage_h_section.chord         = fuselage_h_section_cabin_length + fuselage_h_section_nose_length + fuselage_h_section_tail_length
            
            avl_body.append_section(fuselage_h_section,'horizontal')

    # Vertical Sections of Fuselage 
    if semispan_v != 0:               
        height_array = np.linspace(-semispan_v, semispan_v, num=11,endpoint=True)
        for section_height in height_array :
            fuselage_v_section               = Section()
            
            fuselage_v_section_cabin_length  = avl_body.lengths.total - (avl_body.lengths.nose + avl_body.lengths.tail)
            fuselage_v_section_nose_length   = avl_body.lengths.nose - x_makima_vert_nose_interp(section_height)
            # print('fuselage_v_section_nose_length =', fuselage_v_section_nose_length)
            fuselage_v_section_tail_length   = x_makima_vert_tail_interp(section_height) - (avl_body.lengths.total - avl_body.lengths.tail)
            # print('fuselage_v_section_tail_length =', fuselage_v_section_tail_length)
            fuselage_v_section_nose_origin   = avl_body.lengths.nose - fuselage_v_section_nose_length
            
            fuselage_v_section.tag           = 'fuselage_vertical_top_section_at_' +  str(section_height) + '_m'        
            fuselage_v_section.origin        = [ origin[0] + fuselage_v_section_nose_origin,  origin[1],  origin[2] + section_height ]
            fuselage_v_section.chord         = fuselage_v_section_cabin_length + fuselage_v_section_nose_length + fuselage_v_section_tail_length
            
            avl_body.append_section(fuselage_v_section,'vertical')
            
    return avl_body
# ====================================== NILS =======================================

