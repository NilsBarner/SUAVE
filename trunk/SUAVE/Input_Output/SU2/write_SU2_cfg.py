## @ingroup Input_Output-SU2
def write_SU2_cfg(tag, SU2_settings):
    """Creates an SU2 .cfg file compatible with SU2 v8.3.0"""

    ref_area = SU2_settings.reference_area
    mach     = SU2_settings.mach_number
    AOA      = SU2_settings.angle_of_attack
    iters    = SU2_settings.maximum_iterations
    p0 = SU2_settings.freestream_pressure  # NILS: added to support analysis at different altitudes
    T0 = SU2_settings.freestream_temperature  # NILS: added to support analysis at different altitudes

    filename = tag + '.cfg'
    f = open(filename, mode='w')

    # ------------------------------------------------------------------
    # Problem definition
    # ------------------------------------------------------------------
    f.write('SOLVER = EULER\n\n')
    f.write('KIND_TURB_MODEL = NONE\n\n')
    f.write('MATH_PROBLEM = DIRECT\n\n')
    f.write('AXISYMMETRIC = NO\n\n')
    f.write('RESTART_SOL = NO\n\n')
    f.write('DISCARD_INFILES = NO\n\n')
    f.write('SYSTEM_MEASUREMENTS = SI\n\n')

    # ------------------------------------------------------------------
    # Freestream definition
    # ------------------------------------------------------------------
    f.write(f'MACH_NUMBER = {float(mach)}\n\n')
    f.write(f'AOA = {float(AOA)}\n\n')
    f.write('SIDESLIP_ANGLE = 0.0\n\n')
    f.write(f'FREESTREAM_PRESSURE = {float(p0)}\n\n')  # NILS: previously fixed at 101325.0
    f.write(f'FREESTREAM_TEMPERATURE = {float(T0)}\n\n')  # NILS: previously fixed at 288.15

    # ------------------------------------------------------------------
    # Reference definition
    # ------------------------------------------------------------------
    f.write('REF_ORIGIN_MOMENT_X = 0.25\n\n')
    f.write('REF_ORIGIN_MOMENT_Y = 0.0\n\n')
    f.write('REF_ORIGIN_MOMENT_Z = 0.0\n\n')
    f.write('REF_LENGTH = 1.0\n\n')
    f.write(f'REF_AREA = {float(ref_area)}\n\n')
    f.write('REF_DIMENSIONALIZATION = FREESTREAM_VEL_EQ_ONE\n\n')
    
    # ------------------------------------------------------------------
    # Boundary conditions
    # ------------------------------------------------------------------
    f.write('MARKER_EULER = ( VEHICLE )\n\n')
    f.write('MARKER_FAR = ( FARFIELD )\n\n')
    f.write('MARKER_SYM = ( SYMPLANE )\n\n')

    f.write('MARKER_PLOTTING = ( VEHICLE )\n\n')
    f.write('MARKER_MONITORING = ( VEHICLE )\n\n')
    f.write('MARKER_DESIGNING = ( VEHICLE )\n\n')  # NILS

    # ------------------------------------------------------------------
    # Numerical methods
    # ------------------------------------------------------------------
    f.write('NUM_METHOD_GRAD = WEIGHTED_LEAST_SQUARES\n\n')
    f.write('OBJECTIVE_FUNCTION = DRAG\n\n')

    f.write('CFL_NUMBER = 5.0\n\n')
    #f.write('CFL_ADAPT = NO\n\n')
    f.write('CFL_ADAPT = YES\n\n')
    f.write('CFL_ADAPT_PARAM = ( 0.5, 1.5, 1.0, 100.0 )\n\n')  # NILS: string format has changed from `f.write('CFL_ADAPT_PARAM = ( 1.5, 0.5, 1.0, 100.0 )\n\n')`
    f.write('RK_ALPHA_COEFF = ( 0.66667, 0.66667, 1.000000 )\n\n')  # NILS

    f.write(f'INNER_ITER = {int(iters)}\n\n')

    f.write('LINEAR_SOLVER = FGMRES\n\n')
    f.write('LINEAR_SOLVER_ERROR = 1E-6\n\n')
    f.write('LINEAR_SOLVER_ITER = 2\n\n')

    # ------------------------------------------------------------------
    # Multigrid
    # ------------------------------------------------------------------
    f.write('MGLEVEL = 3\n\n')
    f.write('MGCYCLE = W_CYCLE\n\n')
    f.write('MG_PRE_SMOOTH = ( 1, 2, 3, 3 )\n\n')
    f.write('MG_POST_SMOOTH = ( 0, 0, 0, 0 )\n\n')
    f.write('MG_DAMP_RESTRICTION = 0.9\n\n')
    f.write('MG_DAMP_PROLONGATION = 0.9\n\n')

    # ------------------------------------------------------------------
    # Flow discretization
    # ------------------------------------------------------------------
    f.write('CONV_NUM_METHOD_FLOW = JST\n\n')
    f.write('MUSCL_FLOW = NO\n\n')
    # f.write('MUSCL_FLOW = YES\n\n')  # NILS: cannot use MUSCL with JST (Error Exit: "Centered schemes do not use MUSCL reconstruction (use MUSCL_FLOW= NO).")
    f.write('SLOPE_LIMITER_FLOW = VENKATAKRISHNAN\n\n')
    f.write('JST_SENSOR_COEFF = ( 0.5, 0.02 )\n\n')
    f.write('TIME_DISCRE_FLOW = EULER_IMPLICIT\n\n')

    # ------------------------------------------------------------------
    # Convergence (SU2 v8+ compliant)
    # ------------------------------------------------------------------
    f.write('CONV_FIELD = ( RMS_DENSITY, RMS_ENERGY )\n\n')
    f.write('CONV_RESIDUAL_MINVAL = -12\n\n')

    # ------------------------------------------------------------------
    # Output control (modern)
    # ------------------------------------------------------------------
    f.write('SCREEN_OUTPUT = ( INNER_ITER, RMS_DENSITY, RMS_ENERGY, LIFT, DRAG )\n\n')
    f.write('SCREEN_WRT_FREQ_INNER = 10\n\n')

    f.write('HISTORY_OUTPUT = ( ITER, RMS_RESIDUAL, LIFT, DRAG )\n\n')
    f.write('HISTORY_WRT_FREQ_INNER = 10\n\n')

    f.write(f'MESH_FILENAME = {tag}.su2\n\n')
    f.write('MESH_FORMAT = SU2\n\n')

    f.write(f'CONV_FILENAME = {tag}_history\n\n')

    f.write('OUTPUT_FILES = ( RESTART, PARAVIEW, SURFACE_PARAVIEW )\n\n')

    f.write(f'SOLUTION_FILENAME = solution_flow.dat\n\n')
    f.write(f'RESTART_FILENAME = {tag}_restart_flow.dat\n\n')
    f.write(f'VOLUME_FILENAME = {tag}_flow\n\n')
    f.write(f'SURFACE_FILENAME = {tag}_surface_flow\n\n')

    f.close()
