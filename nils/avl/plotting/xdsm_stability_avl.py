from pyxdsm.XDSM import XDSM, SOLVER, FUNC, LEFT, DOE

""" Brainstorming

# TASOPT.jl
# Perform retrofit analysis
# Impose FCS loads
# Calculate c.g. location and MMOIs

# Range, FCS weight -> Aircraft sizing -> tasopt_geometry_file
# FCS location -> Aircraft loading -> tasopt_mass_file (c.g., MMOIs)

# SUAVE
# Create SUAVE Aircraft() object from TASOPT.jl geometry
# Create AVL vortex lattice model from SUAVE Aircraft() object

# tasopt_geometry_file, tasopt_mass_file -> SUAVE-AVL wrapper -> avl_geometry_file, avl_mass_file

# AVL
# Trim aircraft in set of TASOPT mission points
# Perform flow analysis
# Perform modal analysis

# avl_mass_file -> Trim analysis -> alpha, delta_elev, delta_flap
# avl_mass_file, avl_geometry_file -> Flow analysis -> aero & stability derivatives
# avl_mass_file, avl_geometry_file -> Modal analysis -> eigenvalues and system matrix

# Python
# Evaluate static stability
# Evaluate dynamic stability
# Create surrogate model

# aero & stability derivatives -> Static stability analysis -> yes/no
# eigenvalues and system matrix -> Dynamic stability analysis -> yes/no
# aero & stability derivatives -> Surrogate model -> stability margin, pitching moment derivatives

# Top-level inputs
# MIL-STD-1797, Range, sigma_FCS, FCS locations
# Top-level outputs
# Space of stability-feasible FCS locations

# Post-process eigenmodes

# "retrofit_tasopt"  -> geometry
# "load_spec" -> magnitude and location
# "cg_mmoi_calcs" -> 
# "suave_ac_obj"
# "avl_vlm_model"
# "trim_avl"
# "oper_avl"
# "mode_avl"
# "stat_stab"
# "dyn_stab"
# "surr_model"
"""

#%% Retrofit under on-design conditions

x = XDSM(use_sfmath=False)

x.add_system(
    "doe",
    DOE,
    "\mathrm{DOE}",
)

# TASOPT
x.add_system(
    "tasopt_sizing",
    FUNC,
    r'\mathrm{Aircraft\,sizing}',
)
x.add_system(
    "mass_analysis",
    FUNC,
    r'\mathrm{Mass\,analysis}',
)

# SUAVE
x.add_system(
    "suave_wrapper",
    FUNC,
    r'\mathrm{VLM\,wrapper}',
)

# AVL
x.add_system(
    "avl_trimming",
    FUNC,
    r'\mathrm{Trim\,analysis}',
)
x.add_system(
    "avl_modal",
    FUNC,
    r'\mathrm{Modal\,analysis}',
)

# Python
x.add_system(
    "python_static_stability",
    FUNC,
    r'\mathrm{Static\,stability\,test}',
)
x.add_system(
    "python_dynamic_stability",
    FUNC,
    r'\mathrm{Dynamic\,stability\,test}',
)
x.add_system(
    "python_surrogate",
    FUNC,
    r'\mathrm{Surrogate\,model}'
)

x.connect(
    "doe",
    "mass_analysis",
    r"\mathrm{FCS\,location}",
)
x.connect(
    "tasopt_sizing",
    "suave_wrapper",
    r"\mathrm{Solid\,geometry}",
)
x.connect(
    "tasopt_sizing",
    "mass_analysis",
    r"W_\mathrm{PL,des},\,\mathrm{c.g.,\,no\,FCS}",
)
x.connect(
    "mass_analysis",
    "suave_wrapper",
    r"\mathrm{c.g.,\,MMOIs}",
)
x.connect(
    "tasopt_sizing",
    "avl_trimming",
    r"M,\,h,\,n,\,W",
)
x.connect(
    "tasopt_sizing",
    "avl_modal",
    r"M,\,h,\,n,\,W",
)
x.connect(
    "tasopt_sizing",
    "python_surrogate",
    r"M,\,h,\,n,\,W",
)
x.connect(
    "suave_wrapper",
    "avl_trimming",
    r'\mathrm{VLM\,geometry}',
)
x.connect(
    "suave_wrapper",
    "avl_modal",
    r'\mathrm{VLM\,geometry}',
)
x.connect(
    "avl_trimming",
    "avl_modal",
    r"\alpha,\,\delta_\mathrm{flap},\,\delta_\mathrm{elevator}",
)
x.connect(
    "avl_trimming",
    "python_static_stability",
    r"\begin{array}{c}"
    r"\mathrm{Stability\,and}\\"
    r"\mathrm{aero\,coefficients}"
    r"\end{array}"
)
x.connect(
    "avl_trimming",
    "python_surrogate",
    r"\begin{array}{c}"
    r"\mathrm{Stability\,and}\\"
    r"\mathrm{aero\,coefficients}"
    r"\end{array}"
)
x.connect(
    "avl_modal",
    "python_dynamic_stability",
    r"\begin{array}{c}"
    r"\mathrm{Eigenvalues\,and}\\"
    r"\mathrm{system\,matrix}"
    r"\end{array}"
)
x.connect(
    "python_static_stability",
    "doe",
    r"\mathrm{Yes/no}",
)
x.connect(
    "python_dynamic_stability",
    "doe",
    r"\mathrm{Yes/no}",
)

x.add_input(
    "tasopt_sizing",
    "R_\mathrm{des},\,\sigma_\mathrm{FCS}",
)
x.add_input(
    "python_static_stability",
    "\mathrm{MIL-STD-1797}",
)
x.add_input(
    "python_dynamic_stability",
    "\mathrm{MIL-STD-1797}",
)
x.add_input(
    "doe",
    r"\begin{array}{c}"
    r"\mathrm{Conceivable}\\"
    r"\mathrm{FCS\,locations}"
    r"\end{array}"
)

x.add_output(
    "doe",
    r"\begin{array}{c}"
    r"\mathrm{Feasible}\\"
    r"\mathrm{FCS\,locations}"
    r"\end{array}",
    side=LEFT,
)
x.add_output(
    "python_surrogate",
    r"\mathrm{Stability\,surrogate}",
    side=LEFT,
)

x.write("xdsm_stability_avl_on_design")  # BE SURE TO DELETE EXISTING FILE WITH THAT NAME AS LATEX-ERROR IN THAT ONE WILL PREVENT NEW FILE FROM SAVING

#%% Retrofit under off-design conditions

x = XDSM(use_sfmath=False)

x.add_system(
    "doe",
    DOE,
    "\mathrm{DOE}",
)

# TASOPT
x.add_system(
    "tasopt_sizing",
    FUNC,
    r'\mathrm{Aircraft\,sizing}',
)
x.add_system(
    "tasopt_rating",
    FUNC,
    r'\mathrm{Aircraft\,rating}',
)
x.add_system(
    "mass_analysis",
    FUNC,
    r'\mathrm{Mass\,analysis}',
)

# SUAVE
x.add_system(
    "suave_wrapper",
    FUNC,
    r'\mathrm{VLM\,wrapper}',
)

# AVL
x.add_system(
    "avl_trimming",
    FUNC,
    r'\mathrm{Trim\,analysis}',
)
x.add_system(
    "avl_modal",
    FUNC,
    r'\mathrm{Modal\,analysis}',
)

# Python
x.add_system(
    "python_static_stability",
    FUNC,
    r'\mathrm{Static\,stability\,test}',
)
x.add_system(
    "python_dynamic_stability",
    FUNC,
    r'\mathrm{Dynamic\,stability\,test}',
)
x.add_system(
    "python_surrogate",
    FUNC,
    r'\mathrm{Surrogate\,model}'
)

x.connect(
    "doe",
    "mass_analysis",
    r"\mathrm{FCS\,location}",
)
x.connect(
    "tasopt_sizing",
    "suave_wrapper",
    r"\mathrm{Solid\,geometry}",
)
x.connect(
    "tasopt_sizing",
    "tasopt_rating",
    r"W_\mathrm{PL,des}",
)
x.connect(
    "tasopt_rating",
    "mass_analysis",
    r"\mathrm{c.g.,\,no\,FCS}",
)
x.connect(
    "mass_analysis",
    "suave_wrapper",
    r"\mathrm{c.g.,\,MMOIs}",
)
x.connect(
    "tasopt_sizing",
    "avl_trimming",
    r"M,\,h,\,n,\,W",
)
x.connect(
    "tasopt_sizing",
    "avl_modal",
    r"M,\,h,\,n,\,W",
)
x.connect(
    "tasopt_sizing",
    "python_surrogate",
    r"M,\,h,\,n,\,W",
)
x.connect(
    "suave_wrapper",
    "avl_trimming",
    r'\mathrm{VLM\,geometry}',
)
x.connect(
    "suave_wrapper",
    "avl_modal",
    r'\mathrm{VLM\,geometry}',
)
x.connect(
    "avl_trimming",
    "avl_modal",
    r"\alpha,\,\delta_\mathrm{flap},\,\delta_\mathrm{elevator}",
)
x.connect(
    "avl_trimming",
    "python_static_stability",
    r"\begin{array}{c}"
    r"\mathrm{Stability\,and}\\"
    r"\mathrm{aero\,coefficients}"
    r"\end{array}"
)
x.connect(
    "avl_trimming",
    "python_surrogate",
    r"\begin{array}{c}"
    r"\mathrm{Stability\,and}\\"
    r"\mathrm{aero\,coefficients}"
    r"\end{array}"
)
x.connect(
    "avl_modal",
    "python_dynamic_stability",
    r"\begin{array}{c}"
    r"\mathrm{Eigenvalues\,and}\\"
    r"\mathrm{system\,matrix}"
    r"\end{array}"
)
x.connect(
    "python_static_stability",
    "doe",
    r"\mathrm{Yes/no}",
)
x.connect(
    "python_dynamic_stability",
    "doe",
    r"\mathrm{Yes/no}",
)

x.add_input(
    "tasopt_sizing",
    "R_\mathrm{des},\,\sigma_\mathrm{FCS}",
)
x.add_input(
    "tasopt_rating",
    "R_\mathrm{off},\,W_\mathrm{PL,off}",
)
x.add_input(
    "python_static_stability",
    "\mathrm{MIL-STD-1797}",
)
x.add_input(
    "python_dynamic_stability",
    "\mathrm{MIL-STD-1797}",
)
x.add_input(
    "doe",
    r"\begin{array}{c}"
    r"\mathrm{Conceivable}\\"
    r"\mathrm{FCS\,locations}"
    r"\end{array}"
)

x.add_output(
    "doe",
    r"\begin{array}{c}"
    r"\mathrm{Feasible}\\"
    r"\mathrm{FCS\,locations}"
    r"\end{array}",
    side=LEFT,
)
x.add_output(
    "python_surrogate",
    r"\mathrm{Stability\,surrogate}",
    side=LEFT,
)

x.write("xdsm_stability_avl_off_design")  # BE SURE TO DELETE EXISTING FILE WITH THAT NAME AS LATEX-ERROR IN THAT ONE WILL PREVENT NEW FILE FROM SAVING

