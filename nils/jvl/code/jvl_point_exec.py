r"""
This file is adapted from tut_mission_B737_AVL.py
to contain the minimum set of SUAVE methods needed to
perform a blown-wing analysis using JVL on
aircraft geometry data imported as a .csv from TASOPT.jl,
intended for verification of the surrogate model called
in TASOPT-priv.jl.

To be run as module from Anaconda Prompt via:
(suave) C:\Users\nmb48\Documents\GitHub\SUAVE>python -m nils.jvl.run_suave_jvl_wrapper_nils_new

For further comments, see bottom of file.
"""

__all__ = []

import re
import sys
import ast
import numpy as np
import pandas as pd
from copy import deepcopy
from datetime import date

import SUAVE
assert SUAVE.__version__=='2.5.0', 'These tutorials only work with the SUAVE 2.5.0 release'
from SUAVE.Core import Units
from SUAVE.Methods.Geometry.Two_Dimensional.Planform import segment_properties
from SUAVE.Methods.Propulsion import propeller_design
from SUAVE.Methods.Aerodynamics.AVL.read_results import read_results
from SUAVE.Components.Energy.Networks.Battery_Propeller import Battery_Propeller
from nils.suave_base_ac_object import full_setup
from nils.jvl.jvl_wrapper import wrap_jvl

#%% Command-line parsing
# `python run_suave_avl_wrapper_nils.py study_idx`, where study_idx is an integer 1..7

if len(sys.argv) > 1:
    study_idx = int(sys.argv[1])
else:
    study_idx = 6  # default (choose sensible default) — change if you want no default

#%% Aircraft object

configs, analyses = full_setup()
avl_object = analyses.configs.base.stability

assert id(configs.base) == id(avl_object.geometry)
id_old = id(configs.base)

# Define clean-slate aircraft

aircraft = SUAVE.Vehicle()   # new clean vehicle
aircraft_base = SUAVE.Vehicle()   # base version for diffing

# SUAVE requirement: `_base` holds the baseline configuration
aircraft._base = aircraft_base

# Attach this clean-slate "aircraft" to configs AND avl_object
configs.base = aircraft
avl_object.geometry = aircraft

id_new = id(configs.base)
assert id_new != id_old
assert id(configs.base) == id(avl_object.geometry)

# Read geometry data from .csv file
date_str = date.today().strftime("%d%m%y")
df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\mass\data\mass_distr_geometry_regional_kerosene_230126.csv')
# df_geom = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_kerosene_090126.csv')
# df_geom = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\mass\data\mass_distr_geometry_regional_kerosene_{date_str}.csv')
df_mass = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\mass\data\mass_distr_results_regional_kerosene_230126_6.csv')
# df_mass = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_regional_kerosene_090126_6.csv')
# df_mass = pd.read_csv(rf'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\mass\data\mass_distr_results_regional_kerosene_{date_str}_{study_idx}.csv')
aircraft.tag = 'ATR_72-600'
t_tail_bool = True

#%% Loop over mass configurations of selected study case
# NOTE that unlike in suave_avl_wrapper_exec.py I am only considering one
# mass configuration here, but to use the same format of .csv files I shall
# keep this dummy iteration

for counter, row in df_mass.iterrows():

    results_list = wrap_jvl(
        df_geom, counter, row, t_tail_bool,
    )
    
    
"""
NOTE: the definition of a wing segment always requires at
least two instances of SUAVE.Components.Wings.Segment() -
one at segment.percent_span_location = 0.0 and
one at segment.percent_span_location = 1.0, so that the
wing itself can be interpolated in between.


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

configs.base._base.keys()
Out[19]: dict_keys(['tag', 'fuselages', 'wings', 'networks', 'nacelles', 'systems', 'mass_properties', 'payload', 'costs', 'envelope', 'landing_gear', 'reference_area', 'passengers', 'performance'])

configs.base.keys()
Out[20]: dict_keys(['tag', 'fuselages', 'wings', 'networks', 'nacelles', 'systems', 'mass_properties', 'payload', 'costs', 'envelope', 'landing_gear', 'reference_area', 'passengers', 'performance', '_base', '_diff'])

id(avl_object.geometry._base)
Out[29]: 1532937625664

id(configs.base._base)
Out[30]: 1532937625664


# NOTE: `._base` required to avoid `AttributeError: 'NoneType' object has no attribute 'items'`
# To understand why I have to use `aircraft = avl_object.geometry._base` instead of
# `aircraft = SUAVE.Vehicle()`, compare `print(id(aircraft))` with `print(id(vehicle))`
# in `full_setup()` above - they are not the same object instance!

NOTE: configs.base is a Config, not a Vehicle -> see configs_setup()
"""