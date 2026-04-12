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

from suave_base_ac_object import full_setup

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

def build_jvl_grids(folder, load=False):
    pat = re.compile(r'_(\d+)_([0-9]+\.[0-9]+)\.csv$')

    def parse(files):
        rows = []
        for f in files:
            m = pat.search(os.path.basename(f))
            if m:
                rows.append((int(m.group(1)), m.group(2), f))
        return rows

    geom = parse(glob.glob(os.path.join(folder, 'geometry_*.csv')))
    mass = parse(glob.glob(os.path.join(folder, 'mass_distr_*.csv')))

    N_eng = sorted({n for n, _, _ in geom})
    PRs   = sorted({p for _, p, _ in geom}, key=float)

    def make_df(entries):
        df = pd.DataFrame(index=N_eng, columns=PRs, dtype=object)
        for n, p, f in entries:
            df.at[n, p] = pd.read_csv(f) if load else f
        return df

    return make_df(geom), make_df(mass)


geom_df, mass_df = build_jvl_grids(folder)

N_eng_range = list(geom_df.head().index)  # list of ints
# Prop_PR_des_range = [float(element) for element in list(geom_df.columns)]
Prop_PR_des_range = list(geom_df.columns)  # list of strings

N_eng_grid, Prop_PR_des_grid = np.meshgrid(N_eng_range, Prop_PR_des_range, indexing='ij')

# print('geom_df =', geom_df)
# print()
# print('mass_df =', mass_df)

# geom_df.loc[12, '1.028']
# mass_df.loc[12, '1.028']

# sys.exit()

# =============================================================================

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
            
            configs, analyses = full_setup()
            jvl_object = analyses.configs.base.stability
            
            # Define clean-slate aircraft
        
            aircraft = SUAVE.Vehicle()   # new clean vehicle
            aircraft_base = SUAVE.Vehicle()   # base version for diffing
        
            # SUAVE requirement: `_base` holds the baseline configuration
            aircraft._base = aircraft_base
        
            # Attach this clean-slate "aircraft" to configs AND jvl_object
            configs.base = aircraft
            jvl_object.geometry = aircraft
            aircraft.tag = 'ATR_72-600'
            
            sigma_fcs = row['sigma_fcs']
            span_loc = row['span_loc']
            fcs_loc = row['fcs_loc']
            wing_frac = row['wing_frac']
            nacelle_frac = row['nacelle_frac']
            
            aircraft.mass_properties.center_of_gravity[0][0] = row['x_cg']
            aircraft.mass_properties.center_of_gravity[0][1] = row['y_cg']
            aircraft.mass_properties.center_of_gravity[0][2] = row['z_cg']
            aircraft.mass_properties.mass = row['mass']
            aircraft.mass_properties.max_takeoff = row['max_takeoff']
            aircraft.mass_properties.takeoff = row['takeoff']
            aircraft.mass_properties.max_zero_fuel = row['max_zero_fuel']
            moments_of_inertia = aircraft.mass_properties.moments_of_inertia.tensor
            moments_of_inertia[0][0] = row['Ixx']
            moments_of_inertia[1][1] = row['Iyy']
            moments_of_inertia[2][2] = row['Izz']
            moments_of_inertia[0][1] = row['Ixy']
            moments_of_inertia[1][2] = row['Iyz']
            moments_of_inertia[2][0] = row['Izx']
        
            # Fuselage
            
            fuselage = SUAVE.Components.Fuselages.Fuselage()
            fuselage.tag = 'fuselage'
            
            # for suave_body in aircraft.fuselages:  # write_geometry.py
            fuselage.lengths.total = df_geom['body_lengths_total'][0]
            fuselage.lengths.nose = df_geom['body_lengths_nose'][0]
            fuselage.lengths.tail = df_geom['body_lengths_tail'][0]
            fuselage.width = df_geom['body_widths_maximum'][0]
            fuselage.heights.maximum = df_geom['body_heights_maximum'][0]
            
            aircraft.append_component(fuselage)
                
            # Main wing
            
            wing = SUAVE.Components.Wings.Main_Wing()
            wing.tag = 'main_wing'
            
            wing.aspect_ratio = df_geom['wing_aspect_ratio'][0]
            wing.sweeps.quarter_chord = df_geom['wing_sweeps_quarter_chord'][0]
            # wing.thickness_to_chord = df_geom['wing_thickness_to_chord'][0]  # defined per segment below
            # wing.taper = df_geom['wing_taper'][0]  # defined per segment below
            wing.spans.projected = df_geom['wing_spans_projected'][0]
            wing.chords.root = df_geom['wing_chords_root'][0]
            wing.chords.tip = df_geom['wing_chords_tip'][0]
            wing.chords.mean_aerodynamic = df_geom['wing_chords_mean_aerodynamic'][0]
            wing.areas.reference = df_geom['wing_areas_reference'][0]
            wing.twists.root = df_geom['wing_twists_root'][0]
            wing.twists.tip = df_geom['wing_twists_tip'][0]
            wing.origin = [[
                df_geom['wing_origin_x'][0],
                df_geom['wing_origin_y'][0],
                df_geom['wing_origin_z'][0],
            ]]
            wing.vertical = df_geom['wing_vertical'][0]
            wing.symmetric = df_geom['wing_symmetric'][0]
            wing.high_lift = df_geom['wing_high_lift'][0]
            # wing.dihedral = df_geom['wing_dihedral'][0]  # defined per segment below
            
            ###
            import re
        
            section_indices = sorted(
                {int(re.search(r"sect_(\d+)_", c).group(1))
                 for c in df_geom.columns
                 if c.startswith("sect_")}
            )
            
            for i in section_indices:
                segment_airfoil                          = SUAVE.Components.Airfoils.Airfoil()
                segment_airfoil.coordinate_file          = r'C:\Users\nmb48\Documents\GitHub\SUAVE\regression\scripts\Vehicles\Airfoils\B737a.txt'
                segment                               = SUAVE.Components.Wings.Segment()
                segment.tag                           = df_geom[f'sect_{i}_label'][0]
                # print('segment.tag =', segment.tag)
                segment.percent_span_location         = df_geom[f'sect_{i}_percent_span_location'][0]
                segment.twist                         = 0.0  # not modelled in TASOPT.jl
                segment.root_chord_percent            = df_geom[f'sect_{i}_root_chord_percent'][0]
                segment.thickness_to_chord            = df_geom[f'sect_{i}_thickness_to_chord'][0]
                segment.dihedral_outboard             = df_geom[f'sect_{i}_dihedral_outboard'][0]
                segment.sweeps.quarter_chord          = df_geom[f'sect_{i}_sweeps_quarter_chord'][0]
                segment.append_airfoil(segment_airfoil)
                wing.append_segment(segment)
            
            aircraft.append_component(wing)
            
            # Horizontal stabiliser
            
            wing = SUAVE.Components.Wings.Horizontal_Tail()
            wing.tag = 'horizontal_stabilizer'
            
            wing.aspect_ratio = df_geom['htail_aspect_ratio'][0]
            wing.sweeps.quarter_chord = df_geom['htail_sweeps_quarter_chord'][0]
            # wing.thickness_to_chord = df_geom['htail_thickness_to_chord'][0]  # defined per segment below
            # wing.taper = df_geom['htail_taper'][0]  # defined per segment below
            wing.spans.projected = df_geom['htail_spans_projected'][0]
            wing.chords.root = df_geom['htail_chords_root'][0]
            wing.chords.tip = df_geom['htail_chords_tip'][0]
            # wing.chords.mean_aerodynamic = df_geom['htail_chords_mean_aerodynamic'][0]
            wing.areas.reference = df_geom['htail_areas_reference'][0]
            wing.twists.root = df_geom['htail_twists_root'][0]
            wing.twists.tip = df_geom['htail_twists_tip'][0]
            wing.origin = [[
                df_geom['htail_origin_x'][0],
                df_geom['htail_origin_y'][0],
                df_geom['htail_origin_z'][0],
            ]]
            wing.vertical = df_geom['htail_vertical'][0]
            wing.symmetric = df_geom['htail_symmetric'][0]
            wing.high_lift = df_geom['htail_high_lift'][0]
            # wing.dihedral = df_geom['htail_dihedral'][0]  # defined per segment below
            
            segment                               = SUAVE.Components.Wings.Segment()
            segment.tag                           = 'root_segment'
            segment.percent_span_location         = df_geom['htail_root_percent_span_location'][0]
            segment.twist                         = df_geom['htail_root_twist'][0]
            segment.root_chord_percent            = df_geom['htail_root_root_chord_percent'][0]
            segment.thickness_to_chord            = df_geom['htail_root_thickness_to_chord'][0]
            segment.dihedral_outboard             = df_geom['htail_root_dihedral_outboard'][0]
            segment.sweeps.quarter_chord          = df_geom['htail_root_sweeps_quarter_chord'][0]
            wing.append_segment(segment)
            
            segment                               = SUAVE.Components.Wings.Segment()
            segment.tag                           = 'tip_segment'
            segment.percent_span_location         = df_geom['htail_tip_percent_span_location'][0]
            segment.twist                         = df_geom['htail_tip_twist'][0]
            segment.root_chord_percent            = df_geom['htail_tip_root_chord_percent'][0]
            segment.thickness_to_chord            = df_geom['htail_tip_thickness_to_chord'][0]
            segment.dihedral_outboard             = df_geom['htail_tip_dihedral_outboard'][0]
            segment.sweeps.quarter_chord          = df_geom['htail_tip_sweeps_quarter_chord'][0]
            wing.append_segment(segment)
            
            """ Commented on 23.01.2026 as do not want to trim aircraft for blown-wing design space exploration
            # control surfaces -------------------------------------------
            elevator                       = SUAVE.Components.Wings.Control_Surfaces.Elevator()
            elevator.tag                   = 'elevator'
            elevator.span_fraction_start   = df_geom['htail_elevator_span_fraction_start'][0]
            elevator.span_fraction_end     = df_geom['htail_elevator_span_fraction_end'][0]
            elevator.deflection            = df_geom['htail_elevator_deflection'][0]
            elevator.chord_fraction        = df_geom['htail_elevator_chord_fraction'][0]
            wing.append_control_surface(elevator)
            """
                
            aircraft.append_component(wing)
            
            # Vertical stabiliser
            wing = SUAVE.Components.Wings.Vertical_Tail()
            wing.tag = 'vertical_stabilizer'
            wing.aspect_ratio = df_geom['vtail_aspect_ratio'][0]
            wing.sweeps.quarter_chord = df_geom['vtail_sweeps_quarter_chord'][0]
            wing.thickness_to_chord = df_geom['vtail_thickness_to_chord'][0]
            wing.taper = df_geom['vtail_taper'][0]
            wing.spans.projected = df_geom['vtail_spans_projected'][0]
            wing.chords.root = df_geom['vtail_chords_root'][0]
            wing.chords.tip = df_geom['vtail_chords_tip'][0]
            # wing.chords.mean_aerodynamic = df_geom['vtail_chords_mean_aerodynamic'][0]
            wing.areas.reference = df_geom['vtail_areas_reference'][0]
            wing.twists.root = df_geom['vtail_twists_root'][0]
            wing.twists.tip = df_geom['vtail_twists_tip'][0]
            wing.origin = [[
                df_geom['vtail_origin_x'][0],
                df_geom['vtail_origin_y'][0],
                df_geom['vtail_origin_z'][0],
            ]]
            wing.vertical = df_geom['vtail_vertical'][0]
            wing.symmetric = df_geom['vtail_symmetric'][0]
            wing.high_lift = df_geom['vtail_high_lift'][0]
            wing.dihedral = df_geom['vtail_dihedral'][0]
            wing.t_tail = t_tail_bool
            aircraft.append_component(wing)
            
            # Nacelles
            
            #####
            # Number of nacelles (must be even)
            N_nacelles = len(ast.literal_eval(df_geom['nacelle_origin'].iloc[0])) * 2  # .csv file only contains nacelle positions for one wing half
            assert N_nacelles % 2 == 0, "Number of nacelles must be even"
            N_half = N_nacelles // 2
            
            # --- Base nacelle ---
            base_nacelle = SUAVE.Components.Nacelles.Nacelle()
            base_nacelle.tag = 'nacelle'
            base_nacelle.length = df_geom['nacelle_length'][0]
            base_nacelle.inlet_diameter = df_geom['nacelle_inlet_diameter'][0]
            base_nacelle.diameter = df_geom['nacelle_diameter'][0]
            base_nacelle.areas.wetted = df_geom['nacelle_areas_wetted'][0]
            base_nacelle.flow_through = df_geom['nacelle_flow_through'][0]
            
            nacelle_airfoil = SUAVE.Components.Airfoils.Airfoil()
            nacelle_airfoil.naca_4_series_airfoil = '2410'
            base_nacelle.append_airfoil(nacelle_airfoil)
            
            # --- Create nacelles ---
            nacelles = []
            
            # Positive-Y side
            for i in range(N_half):
                nac = deepcopy(base_nacelle)
                nac.tag = 'nacelle' if i == 0 else f'nacelle_{i+1}'
                nac.origin = np.array([
                    ast.literal_eval(row)[i]
                    for row in df_geom['nacelle_origin']
                ])
                nacelles.append(nac)
            
            # Mirrored negative-Y side
            for i in range(N_half):
                nac = deepcopy(nacelles[i])
                nac.tag = f'nacelle_{i + 1 + N_half}'
                nac.origin[1] *= -1
                nacelles.append(nac)
            
            # --- Append to aircraft ---
            for nac in nacelles:
                aircraft.append_component(nac)
            #####
            
            ##### NILS: from X57_Maxwell_Mod2.py
            #---------------------------------------------------------------------------------------------
            # DEFINE PROPELLER
            #---------------------------------------------------------------------------------------------
            # build network
            net = Battery_Propeller()
            net.number_of_propeller_engines  = N_nacelles
            net.identical_propellers         = True
            
            for nacelle in nacelles:
            
                # Component 2 the Propeller 
                prop = SUAVE.Components.Energy.Converters.Propeller()
                prop.tag = 'propeller_1'
                prop.tip_radius = df_geom['Dprop'][0] / 2
                prop.Jgain = df_geom['Jgain'][0]
                net.propellers.append(prop)
            
            # add the network to the vehicle
            aircraft.append_component(net)
            #####
            
            # NILS: to avoid
            # *** Cannot adjust spanwise spacing at SECTION  2, on SURFACE main_wing_1
            # *** Insufficient number of spanwise vortices to work with
            jvl_object.settings.Nspanwise_main_wing = Nspanwise_main_wing
            
            try:
                results_list = jvl_object.sample_training(
                    study_idx=study_idx, counter=counter,
                    sigma_fcs=sigma_fcs,
                    span_loc=span_loc,
                    fcs_loc=fcs_loc,
                    wing_frac=wing_frac,
                    nacelle_frac=nacelle_frac,
                    N_eng=N_eng,
                    Prop_PR_des=Prop_PR_des,
                )
            except FileNotFoundError:
                print("This input pair crashed JVL. Proceed to next input pair.")
                continue
            
            # print('Hello world')
            # sys.exit()
            
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
            
            # sys.exit('Stop after first run.')
            
#%%

np.savetxt("CTtot_grid.csv", CTtot_grid, delimiter=",", fmt="%.6g")
np.savetxt("CLtot_grid.csv", CLtot_grid, delimiter=",", fmt="%.6g")
np.savetxt("dY_upper_grid.csv", dY_grid, delimiter=",", fmt="%.6g")  # otherwise Windows will overwrite files as does not distinguish between lower and upper case
np.savetxt("dy_lower_grid.csv", dy_grid, delimiter=",", fmt="%.6g")  # otherwise Windows will overwrite files as does not distinguish between lower and upper case

print('CTtot_grid =', CTtot_grid)
print('CLtot_grid =', CLtot_grid)
print('dY_grid =', dY_grid)
print('dy_grid =', dy_grid)

sys.exit('Stop here.')

fig, ax = plt.subplots(figsize=(8, 5.5))

left_ends_list = []
for i, N_eng in enumerate(N_eng_vals):
    
    # For S_wing
    # ax.plot(Prop_PR_des_vals, S_wing[i, :], color=colors[0], label=f"N_eng = {N_eng}", linewidth=0.8)
    # left_end_idx = np.nanargmax(S_wing[i, :])
    # left_end = [Prop_PR_des_vals[left_end_idx], S_wing[i, left_end_idx]]
    
    # For P_prop_max
    ax.plot(Prop_PR_des_vals, P_prop_max[i, :], color=colors[0], label=f"N_eng = {N_eng}", linewidth=0.8)
    finite_mask = np.isfinite(P_prop_max[i, :])
    left_end_idx = np.where(P_prop_max[i, :] == P_prop_max[i, :][finite_mask][0])[0]
    left_end = [Prop_PR_des_vals[left_end_idx], P_prop_max[i, left_end_idx]]
    
    left_ends_list.append(left_end)
    
left_ends_array = np.array(left_ends_list)
ax.plot(left_ends_array[:, 0], left_ends_array[:, 1], marker='.', color='k')

# levels = np.linspace(np.nanmin(blown_wing_fraction), np.nanmax(blown_wing_fraction), 5)
# levels[-1] -= 1e-6
levels = [0.2, 0.3, 0.4, 0.5]

# For S_wing
cs = ax.contour(
    Prop_PR_des_grid, S_wing, blown_wing_fraction, levels=levels, colors=colors[1], zorder=10, extend='max', linewidths=0.8,
)

# For P_prop_max
cs = ax.contour(
    Prop_PR_des_grid, P_prop_max, blown_wing_fraction, levels=levels, colors=colors[1], zorder=10, extend='max', linewidths=0.8,
)

ax.clabel(cs, fmt='%0.1f', fontsize=12, inline=True, inline_spacing=50)

ax.set_xlabel("Propeller stagnation pressure ratio, $\\Pi_\mathrm{o,des}$", labelpad=15)
# ax.set_ylabel("Wing area, $S$ ($\mathrm{m}^2$)", labelpad=15)
ax.set_ylabel("Overall installed power, $P_\mathrm{des}$ (MW)", labelpad=15)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)  #, pad=40)

custom_lines = [
    Line2D([0], [0], color=colors[0], linewidth=0.8),
    Line2D([0], [0], color=colors[1], linewidth=0.8),
    Line2D([0], [0], color='k'),
]
ax.legend(
    custom_lines, [
        r'$N_\mathrm{prop}$',
        r'Blown wing fraction',
        r'Zero prop spacing'
    ],
    frameon=False,
    # loc='lower left',  # for S_wing
    loc='upper left',  # for P_prop_max
    ncols=1,
)

add_margin(ax, m=0.05)

plt.tight_layout()
plt.show()
            
        
    