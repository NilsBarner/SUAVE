"""
This script is based on suave_ac_object_converter.py.
"""

import sys
import ast
import numpy as np
import pandas as pd
import pylab as plt
from copy import deepcopy
from datetime import date

import SUAVE
assert SUAVE.__version__=='2.5.0', 'These tutorials only work with the SUAVE 2.5.0 release'
from SUAVE.Core import Units
from SUAVE.Methods.Geometry.Two_Dimensional.Planform import segment_properties

# NILS: added below lines
from SUAVE.Methods.Propulsion import propeller_design
from SUAVE.Components.Energy.Networks.Battery_Propeller import Battery_Propeller

# NILS: added below lines
import os, re, glob
import pandas as pd
import sys
os.chdir(r'C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\B737_AVL_Tutorial\jvl')

from nils.suave_base_ac_object import full_setup
from nils.jvl.jvl_wrapper import wrap_jvl
from nils.jvl.shared_methods import build_jvl_grids

configs, analyses = full_setup()
jvl_object = analyses.configs.base.stability

# NOTE: `._base` required to avoid `AttributeError: 'NoneType' object has no attribute 'items'`
# To understand why I have to use `aircraft = jvl_object.geometry._base` instead of
# `aircraft = SUAVE.Vehicle()`, compare `print(id(aircraft))` with `print(id(vehicle))`
# in `full_setup()` above - they are not the same object instance!

ac_segment = "regional"  # "narrowbody" or "regional"
study_idx = 6

"""
NOTE: configs.base is a Config, not a Vehicle -> see configs_setup()
"""

assert id(configs.base) == id(jvl_object.geometry)
id_old = id(configs.base)

# Define clean-slate aircraft

aircraft = SUAVE.Vehicle()   # new clean vehicle
aircraft_base = SUAVE.Vehicle()   # base version for diffing

# SUAVE requirement: `_base` holds the baseline configuration
aircraft._base = aircraft_base

# Attach this clean-slate "aircraft" to configs AND jvl_object
configs.base = aircraft
jvl_object.geometry = aircraft

id_new = id(configs.base)
assert id_new != id_old
assert id(configs.base) == id(jvl_object.geometry)

# =============================================================================
folder = os.path.join(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\jvl\data')  # or any other folder path

# Find all .csv files in above folder containing geometry and mass distribution information

# geometry_files = glob.glob(os.path.join(folder, 'geometry_*.csv'))
# mass_distr_files = glob.glob(os.path.join(folder, 'mass_distr_*.csv'))

geom_df, mass_df = build_jvl_grids(folder)
N_eng_range = list(geom_df.head().index)  # list of ints
# Prop_PR_des_range = [float(element) for element in list(geom_df.columns)]
Prop_PR_des_range = list(geom_df.columns)  # list of strings
N_eng_grid, Prop_PR_des_grid = np.meshgrid(
    N_eng_range, Prop_PR_des_range, indexing='ij',
)

# Aircraft mass properties

shape = N_eng_grid.shape
CTtot_grid = np.zeros(shape, dtype=float)
CLtot_grid = np.zeros(shape, dtype=float)
dY_grid = np.zeros(shape, dtype=float)
dy_grid = np.zeros(shape, dtype=float)

Nspanwise_main_wing_list = [100, 50, 50, 20, 20]

for _i, N_eng in enumerate(N_eng_range):
    Nspanwise_main_wing = Nspanwise_main_wing_list[_i]
    
    for _j, Prop_PR_des in enumerate(Prop_PR_des_range):
        print('N_eng, Prop_PR_des =', N_eng, Prop_PR_des)
        
        # Read geometry data from .csv file
        date_str = date.today().strftime("%d%m%y")
        try:
            df_geom = pd.read_csv(geom_df.loc[N_eng, Prop_PR_des])
            df_mass = pd.read_csv(mass_df.loc[N_eng, Prop_PR_des])
        except ValueError:
            print("This input pair did not result in a valid TASOPT.jl design. Proceed to next input pair.")
            continue
        aircraft.tag = 'ATR_72-600'
        t_tail_bool = True
        
        for counter, row in df_mass.iterrows():
        
            # NILS: to avoid
            # *** Cannot adjust spanwise spacing at SECTION  2, on SURFACE main_wing_1
            # *** Insufficient number of spanwise vortices to work with
            # jvl_object.settings.Nspanwise_main_wing = Nspanwise_main_wing
            
            try:
                results_list = wrap_jvl(
                    df_geom, counter, row, t_tail_bool,
                )
            except FileNotFoundError:
                print("This input pair crashed JVL. Proceed to next input pair.")
                continue
            
            # Extract overall lift coefficient and overall drag coefficient
            
            CTtot = results_list[0]['case_01_01'].aerodynamics.CTtot
            CLtot = results_list[0]['case_01_01'].aerodynamics.total_lift_coefficient
            dY = df_geom['dY'][0]
            dy = df_geom['dy'][0]
            
            CTtot_grid[_i, _j] = CTtot
            CLtot_grid[_i, _j] = CLtot
            dY_grid[_i, _j] = dY
            dy_grid[_i, _j] = dy
            
            print('CTtot, CLtot =', CTtot, CLtot)
            
            
