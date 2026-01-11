## @ingroup Methods-Aerodynamics-AVL
#write_geometry.py
# 
# Created:  Oct 2015, T. Momose
# Modified: Jan 2016, E. Botero
#           Oct 2018, M. Clarke
#           Aug 2019, M. Clarke
#           Apr 2020, M. Clarke

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import re  # NILS: added
from .purge_files import purge_files
from SUAVE.Methods.Aerodynamics.AVL.Data.Settings    import Settings
import numpy as np
import shutil
from .create_avl_datastructure import translate_avl_wing, translate_avl_body, translate_avl_nacelle_nils  # NILS: added 'translate_avl_nacelle

## @ingroup Methods-Aerodynamics-AVL
def write_geometry(avl_object,run_script_path):
    """This function writes the translated aircraft geometry into text file read 
    by AVL when it is called

    Assumptions:
        None
        
    Source:
        Drela, M. and Youngren, H., AVL, http://web.mit.edu/drela/Public/web/avl

    Inputs:
        avl_object

    Outputs:
        None

    Properties Used:
        N/A
    """    
    
    # unpack inputs
    aircraft                   = avl_object.geometry
    geometry_file              = avl_object.settings.filenames.features
    number_spanwise_vortices   = avl_object.settings.number_spanwise_vortices
    number_chordwise_vortices  = avl_object.settings.number_chordwise_vortices
    # Open the geometry file after purging if it already exists
    purge_files([geometry_file]) 
    geometry             = open(geometry_file,'w')

    with open(geometry_file,'w') as geometry:
        header_text       = make_header_text(avl_object)
        geometry.write(header_text)
        
        for w in aircraft.wings:
            avl_wing      = translate_avl_wing(w)
            wing_text     = make_surface_text(avl_wing,number_spanwise_vortices,number_chordwise_vortices)
            ### NILS
            Nspanwise_main_wing = avl_object.settings.Nspanwise_main_wing  # NILS: reduce number of spanwise vortices to avoid SPUPL error
            if avl_wing.tag == 'main_wing':
                wing_text = adjust_wing_text_for_jvl(avl_wing, wing_text, Nspanwise_main_wing)
            ###
            geometry.write(wing_text)
            
        for b in aircraft.fuselages:
            avl_body  = translate_avl_body(b)
            body_text = make_body_text(avl_body,number_spanwise_vortices,number_chordwise_vortices)  # NILS: added second argument
            geometry.write(body_text)
            
        # NILS: TODO implement separate chord- and spanwise number of vortices for nacelles
        for n in aircraft.nacelles:
            avl_nacelle = translate_avl_nacelle_nils(n)
            
            if n.flow_through == True:
                nacelle_text = make_surface_text(avl_nacelle, 12, 6)
            elif n.flow_through == False:
                nacelle_text = make_body_text(avl_nacelle, 11, 6)  # NILS: added second argument
            geometry.write(nacelle_text)
        
    return


def make_header_text(avl_object):  
    """This function writes the header using the template required for the AVL executable to read

    Assumptions:
        None
        
    Source:
        None

    Inputs:
        avl_object.settings.flow_symmetry.xz_plane                      [-]
        avl_object.settings.flow_symmetry.xy_parallel                   [-]
        avl_object.settings.flow_symmetry.z_symmetry_plane              [-]
        avl_object.geometry.wings['main_wing'].areas.reference          [meters**2]
        avl_object.geometry.wings['main_wing'].chords.mean_aerodynamic  [meters]
        avl_object.geometry.wings['main_wing'].spans.projected          [meters]
        avl_object.geometry.mass_properties.center_of_gravity           [meters]
        avl_object.geometry.tag                                         [-]
    
    Outputs:
        header_text                                                     [-]

    Properties Used:
        N/A
    """      
    header_base = \
'''{0}

#Mach
 {1}
 
#Iysym   IZsym   Zsym
  {2}      {3}     {4}
  
#Sref    Cref    Bref 	<meters>
{5}      {6}     {7}

#Xref    Yref    Zref   <meters>
{8}      {9}     {10}

'''

    # Unpack inputs
    Iysym = avl_object.settings.flow_symmetry.xz_plane
    Izsym = avl_object.settings.flow_symmetry.xy_parallel
    Zsym  = avl_object.settings.flow_symmetry.z_symmetry_plane
    Sref  = avl_object.geometry.wings['main_wing'].areas.reference
    Cref  = avl_object.geometry.wings['main_wing'].chords.mean_aerodynamic
    Bref  = avl_object.geometry.wings['main_wing'].spans.projected
    Xref  = avl_object.geometry.mass_properties.center_of_gravity[0][0]
    Yref  = avl_object.geometry.mass_properties.center_of_gravity[0][1]
    Zref  = avl_object.geometry.mass_properties.center_of_gravity[0][2]
    name  = avl_object.geometry.tag

    mach = 0.0

    # Insert inputs into the template
    header_text = header_base.format(name,mach,Iysym,Izsym,Zsym,Sref,Cref,Bref,Xref,Yref,Zref)

    return header_text


def make_surface_text(avl_wing,number_spanwise_vortices,number_chordwise_vortices):
    """This function writes the surface text using the template required for the AVL executable to read

    Assumptions:
        None
        
    Source:
        None

    Inputs:
       avl_wing.symmetric
       avl_wing.tag
        
    Outputs:
        surface_text                                                 

    Properties Used:
        N/A
    """       
    ordered_tags = []         
    surface_base = \
        '''

#---------------------------------------------------------
SURFACE
{0}
#Nchordwise  Cspace   Nspanwise  Sspace
{1}         {2}         {3}      {4}{5}{6}{7}{8}
'''  # NILS: added {6}{7}{8} and reassigned {5}
    # Unpack inputs
    symm = avl_wing.symmetric
    name = avl_wing.tag
    
    # NILS
    if avl_wing.tag == 'nacelle':
        ycomp = '\n\nCOMPONENT\n1\n'
        yscale = f'\n\nSCALE\n{avl_wing.x_scale}  {avl_wing.y_scale}  {avl_wing.z_scale}\n'
        ytransl = f'\n\nTRANSLATE\n{avl_wing.x_transl}  {avl_wing.y_transl}  {avl_wing.z_transl}\n'
    else:
        ycomp     = ' ' 
        yscale     = ' ' 
        ytransl     = ' ' 

    if symm:
        ydup = '\n\nYDUPLICATE\n0.0\n' # Duplication of wing about xz plane
    else:
        ydup     = ' ' 
    
    # NILS: Circular nacelles
    if avl_wing.tag == 'nacelle':
        # Define precision of analysis. See AVL documentation for reference
        chordwise_vortex_spacing = 1.0        
        spanwise_vortex_spacing  = 0.0
        ordered_tags = sorted(
            avl_wing.sections,
            key=lambda s: np.arctan2(s.origin[0][2], s.origin[0][1])
        )  # sort by angle rather than spanwise coordinate
    
        # Write text  
        # surface_text = surface_base.format(name,number_chordwise_vortices,chordwise_vortex_spacing,number_spanwise_vortices ,spanwise_vortex_spacing,ydup)
        surface_text = surface_base.format(name,number_chordwise_vortices,chordwise_vortex_spacing,number_spanwise_vortices ,spanwise_vortex_spacing,ycomp,ydup,yscale,ytransl)
        for i in range(len(ordered_tags)):
            section_text    = make_turbofan_nacelle_section_text_nils(ordered_tags[i])
            surface_text    = surface_text + section_text
            
    else:
    
        # Vertical Wings
        if avl_wing.vertical:
            # Define precision of analysis. See AVL documentation for reference 
            chordwise_vortex_spacing = 1.0
            spanwise_vortex_spacing  = -1.1                              # cosine distribution i.e. || |   |    |    |  | ||
            ordered_tags = sorted(avl_wing.sections, key = lambda x: x.origin[0][2])
            
            # Write text 
            # surface_text = surface_base.format(name,number_chordwise_vortices,chordwise_vortex_spacing,number_spanwise_vortices ,spanwise_vortex_spacing,ydup)
            surface_text = surface_base.format(name,number_chordwise_vortices,chordwise_vortex_spacing,number_spanwise_vortices ,spanwise_vortex_spacing,ycomp,ydup,yscale,ytransl)
            for i in range(len(ordered_tags)):
                section_text    = make_wing_section_text(ordered_tags[i])
                surface_text    = surface_text + section_text
        
        # Horizontal Wings        
        else:        
            # Define precision of analysis. See AVL documentation for reference
            chordwise_vortex_spacing = 1.0        
            spanwise_vortex_spacing  = 1.0                              # cosine distribution i.e. || |   |    |    |  | ||
            ordered_tags = sorted(avl_wing.sections, key = lambda x: x.origin[0][1])
        
            # Write text  
            # surface_text = surface_base.format(name,number_chordwise_vortices,chordwise_vortex_spacing,number_spanwise_vortices ,spanwise_vortex_spacing,ydup)
            surface_text = surface_base.format(name,number_chordwise_vortices,chordwise_vortex_spacing,number_spanwise_vortices ,spanwise_vortex_spacing,ycomp,ydup,yscale,ytransl)
            for i in range(len(ordered_tags)):
                section_text    = make_wing_section_text(ordered_tags[i])
                surface_text    = surface_text + section_text
    
    return surface_text


def make_body_text(avl_body,number_spanwise_vortices,number_chordwise_vortices):   
    """This function writes the body text using the template required for the AVL executable to read

    Assumptions:
        None
        
    Source:
        None

    Inputs:
        avl_body.sections.horizontal
        avl_body.sections.vertical
    
    Outputs:
        body_text                                                 

    Properties Used:
        N/A
    """      
    surface_base = \
'''

#---------------------------------------------------------
SURFACE
{0}
#Nchordwise  Cspace   Nspanwise  Sspace
{1}           {2}        {3}      {4}{5}
'''  # NILS: added {3}{4}{5}
    # Unpack inputs
    name = avl_body.tag
    
    # Define precision of analysis. See AVL documentation for reference 
    chordwise_vortex_spacing = 1.0 
    
    if name == 'nacelle':
        ydup = '\n\nYDUPLICATE\n0.0\n' # Duplication of wing about xz plane
    else:
        ydup = ' ' 
    
    # Form the horizontal part of the + shaped fuselage    
    hname           = name + '_horizontal'
    spanwise_vortex_spacing  = 0.0  # NILS: impose spanwise symmetry
    horizontal_text = surface_base.format(hname,number_chordwise_vortices,chordwise_vortex_spacing,number_spanwise_vortices,spanwise_vortex_spacing,ydup)  # NILS: added last three arguments
       
    ordered_tags = []
    ordered_tags = sorted(avl_body.sections.horizontal, key = lambda x: x.origin[1])
    for i in range(len(ordered_tags)):
        section_text    = make_body_section_text(ordered_tags[i])
        horizontal_text = horizontal_text + section_text
        
    # Form the vertical part of the + shaped fuselage
    vname         = name + '_vertical'
    spanwise_vortex_spacing  = 0.0  # NILS: impose spanwise symmetry
    vertical_text = surface_base.format(vname,number_chordwise_vortices,chordwise_vortex_spacing,number_spanwise_vortices,spanwise_vortex_spacing,ydup)  # NILS: added last three arguments
    ordered_tags = []
    ordered_tags = sorted(avl_body.sections.vertical, key = lambda x: x.origin[2])
    for i in range(len(ordered_tags)):
        section_text    = make_body_section_text(ordered_tags[i])
        vertical_text = vertical_text + section_text
        
    body_text = horizontal_text + vertical_text
    return body_text  


def make_wing_section_text(avl_section):
    """This function writes the wing text using the template required for the AVL executable to read

    Assumptions:
        None
        
    Source:
        None

    Inputs:
       avl_section.origin             [meters]
       avl_section.chord              [meters]
       avl_section.twist              [radians]
       avl_section.airfoil_coord_file [-] 
        
    Outputs:
        wing_section_text                                                 

    Properties Used:
        N/A
    """      
#     section_base = \
# '''
# SECTION
# #Xle    Yle      Zle      Chord     Ainc  Nspanwise  Sspace
# {0}  {1}    {2}    {3}    {4}     
# '''
    section_base = \
'''
SECTION
# {0}
#Xle    Yle      Zle      Chord     Ainc  Nspanwise  Sspace
{1}  {2}    {3}    {4}    {5}     
'''  # NILS
    airfoil_base = \
'''AFILE
{}
'''
    naca_airfoil_base = \
'''NACA
{}
'''
    # Unpack inputs
    section_tag = avl_section.tag  # NILS
    x_le          = avl_section.origin[0][0]
    y_le          = avl_section.origin[0][1]
    z_le          = avl_section.origin[0][2]
    chord         = avl_section.chord
    ainc          = avl_section.twist
    airfoil_coord = avl_section.airfoil_coord_file
    naca_airfoil  = avl_section.naca_airfoil
     
    # wing_section_text = section_base.format(round(x_le,4),round(y_le,4), round(z_le,4),round(chord,4),round(ainc,4))
    wing_section_text = section_base.format(section_tag,round(x_le,4),round(y_le,4), round(z_le,4),round(chord,4),round(ainc,4))  # NILS
    if airfoil_coord:
        wing_section_text = wing_section_text + airfoil_base.format(airfoil_coord)
    if naca_airfoil:
        wing_section_text = wing_section_text + naca_airfoil_base.format(naca_airfoil)        
    
    ordered_cs = []
    ordered_cs = sorted(avl_section.control_surfaces, key = lambda x: x.order)
    for i in range(len(ordered_cs)):
        control_text = make_controls_text(ordered_cs[i])
        wing_section_text = wing_section_text + control_text

    return wing_section_text

# NILS
def make_turbofan_nacelle_section_text_nils(avl_section):
    """This function writes the wing text using the template required for the AVL executable to read

    Assumptions:
        None
        
    Source:
        None

    Inputs:
       avl_section.origin             [meters]
       avl_section.chord              [meters]
       avl_section.twist              [radians]
       avl_section.airfoil_coord_file [-] 
        
    Outputs:
        wing_section_text                                                 

    Properties Used:
        N/A
    """      
    section_base = \
'''
SECTION
#Xle     Yle      Zle      Chord     Ainc  Nspanwise  Sspace
{0}    {1}     {2}     {3}     {4}      1        0.
'''
    airfoil_base = \
'''AFILE
{}
'''
    naca_airfoil_base = \
'''NACA
{}
'''
    # Unpack inputs
    x_le          = avl_section.origin[0][0]
    y_le          = avl_section.origin[0][1]
    z_le          = avl_section.origin[0][2]
    chord         = avl_section.chord
    ainc          = avl_section.twist
    airfoil_coord = avl_section.airfoil_coord_file
    naca_airfoil  = avl_section.naca_airfoil
     
    wing_section_text = section_base.format(round(x_le,4),round(y_le,4), round(z_le,4),round(chord,4),round(ainc,4))
    if airfoil_coord:
        wing_section_text = wing_section_text + airfoil_base.format(airfoil_coord)
    if naca_airfoil:
        wing_section_text = wing_section_text + naca_airfoil_base.format(naca_airfoil)        
    
    return wing_section_text
    
def make_body_section_text(avl_body_section):
    """This function writes the body text using the template required for the AVL executable to read

    Assumptions:
        None
        
    Source:
        None

    Inputs:
       avl_section.origin             [meters]
       avl_section.chord              [meters]
       avl_section.twist              [radians]
       avl_section.airfoil_coord_file [-] 
                  
    Outputs:
        body_section_text                                                 

    Properties Used:
        N/A
    """    
    section_base = \
'''
SECTION
#Xle     Yle      Zle      Chord     Ainc  Nspanwise  Sspace
{0}    {1}     {2}     {3}     {4}      1        0
'''
    airfoil_base = \
'''AFILE
{}
'''

    # Unpack inputs
    x_le    = avl_body_section.origin[0]
    y_le    = avl_body_section.origin[1]
    z_le    = avl_body_section.origin[2]
    chord   = avl_body_section.chord
    ainc    = avl_body_section.twist
    airfoil = avl_body_section.airfoil_coord_file

    body_section_text = section_base.format(round(x_le,4),round(y_le,4), round(z_le,4),round(chord,4),round(ainc,4))
    if airfoil:
        body_section_text = body_section_text + airfoil_base.format(airfoil)
    
    return body_section_text

    
def make_controls_text(avl_control_surface):
    """This function writes the control surface text using the template required 
    for the AVL executable to read

    Assumptions:
        None
        
    Source:
        None

    Inputs:
        avl_control_surface.tag             [-]
        avl_control_surface.gain            [-]
        avl_control_surface.x_hinge         [-]
        avl_control_surface.hinge_vector    [-]
        avl_control_surface.sign_duplicate  [-]
                  
    Outputs:
        control_text                                                 

    Properties Used:
        N/A
    """    
    control_base = \
'''CONTROL
{0}    {1}   {2}   {3}  {4}
'''

    # Unpack inputs
    name     = avl_control_surface.tag
    gain     = avl_control_surface.gain
    xhinge   = avl_control_surface.x_hinge
    hv       = avl_control_surface.hinge_vector
    sign_dup = avl_control_surface.sign_duplicate

    control_text = control_base.format(name,gain,xhinge,hv,sign_dup)

    return control_text


### NILS
def adjust_wing_text_for_jvl(avl_wing, wing_text, Nspanwise_main_wing):
    
    import re
    from textwrap import dedent
    
    def parse_sections(text):
        blocks = re.split(r"\n(?=SECTION\n)", text)
        out = []
        for b in blocks:
            if b.startswith("SECTION"):
                name = re.search(r"#\s*(\S+)", b).group(1)
                out.append((name, b.strip()))
        return out


    def is_prop_in(name):
        return re.fullmatch(r"prop_\d+_in", name)


    def is_prop_out(name):
        return re.fullmatch(r"prop_\d+_out", name)


    def surface_header(idx, first):
        return dedent(f"""
        #---------------------------------------------------------
        SURFACE
        main_wing_{idx}
        #Nchordwise  Cspace   Nspanwise  Sspace
        10         1.0         {Nspanwise_main_wing}      1.0 

        YDUPLICATE
        {"0.0" if first else "-0.0"}
        """).strip()


    def transform(source):
        sections = parse_sections(source)
        surfaces = []
        current = []

        for name, block in sections:
            current.append(block)

            if is_prop_in(name):
                surfaces.append(current)
                current = [block]

            elif is_prop_out(name):
                surfaces.append(current)
                current = [block]

        if current:
            surfaces.append(current)

        out = []
        for i, surf in enumerate(surfaces, 1):
            out.append(surface_header(i, first=(i == 1)))
            out.append("")
            out.append("\n\n".join(surf))
            out.append("")

        return "\n".join(out).strip()
    
    wing_text = transform(wing_text)
    
    # print('wing_text =', wing_text)
    # import sys
    # sys.exit('Stop here.')
    
    jet_base = \
'''
JETCONTROL
#Jname   Jgain    SgnDup
DVjet    {0}      1.0
JETPARAM
#hdisk   fh      djet0   djet1   djet3
{1}    {2}     {3}    {4}   {5}
'''
    
    # Unpack inputs
    Njet = avl_wing.Njet
    # Jgain = avl_wing.Jgain
    # hdisk = avl_wing.hdisk
    fh = avl_wing.fh
    djet0 = avl_wing.djet0
    djet1 = avl_wing.djet1
    djet3 = avl_wing.djet3
    
    surfaces = wing_text.split("#---------------------------------------------------------")
    new_surfaces = [surfaces[0]]  # keep header part untouched
    
    for surf in surfaces[1:]:
        # Check whether this is a prop-SURFACE:
        # exactly two SECTIONs, both with 'prop' in the tag line
        section_tags = re.findall(r'SECTION\s*\n#\s*(.*)', surf)
        
        # Check whether current surf is a "prop-surf"
        
        section_tags = re.findall(r'SECTION\s*\n#\s*(.*)', surf)

        # Require that `id` be the same in subsequent prop_{id}_in and prop_{id}_out section
        prop_ids = []
        for tag in section_tags:
            m = re.match(r'prop_(\d+)_', tag)
            if m:
                prop_ids.append(m.group(1))
        
        is_prop_surface = (
            len(section_tags) == 2 and
            len(prop_ids) == 2 and
            prop_ids[0] == prop_ids[1]
        )
    
        if is_prop_surface:
            # 1) Add Njet to header and value line
            surf = surf.replace(
                "#Nchordwise  Cspace   Nspanwise  Sspace",
                "#Nchordwise  Cspace   Nspanwise  Sspace  Njet"
            )
    
            surf = re.sub(
                r'(\n\d+[^\n]*)',
                r'\1      {}'.format(Njet),
                surf,
                count=1
            )
    
            # 2) Insert jet_text after EACH *.dat line
            Jgain = avl_wing.Jgain
            hdisk = avl_wing.hdisk
            jet_text = jet_base.format(
                round(Jgain,4),round(hdisk,4), round(fh,4),round(djet0,4),round(djet1,4),round(djet3,4)
            )
            surf = re.sub(
                r'(AFILE\s*\n.*\.dat)',
                r'\1\n' + jet_text,
                surf
            )
    
        new_surfaces.append("#---------------------------------------------------------" + surf)
    
    wing_text = "".join(new_surfaces)
    
    # print('wing_text =', wing_text)
    # import sys
    # sys.exit('Stop here.')
    
    return wing_text
###


