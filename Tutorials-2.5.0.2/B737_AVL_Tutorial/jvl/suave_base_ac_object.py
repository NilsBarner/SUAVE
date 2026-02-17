"""

"""

__all__ = []

import SUAVE
assert SUAVE.__version__=='2.5.0', 'These tutorials only work with the SUAVE 2.5.0 release'

# ----------------------------------------------------------------------
#   Analysis Setup
# ----------------------------------------------------------------------

def full_setup(vehicle=None):

    # vehicle data
    if vehicle == None:
        vehicle  = vehicle_setup()
        # print(id(vehicle))
    configs  = configs_setup(vehicle)

    # vehicle analyses
    configs_analyses = analyses_setup(configs)

    analyses = SUAVE.Analyses.Analysis.Container()
    analyses.configs  = configs_analyses
    
    return configs, analyses

# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ----------------------------------------------------------------------

def analyses_setup(configs):

    analyses = SUAVE.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

def base_analysis(vehicle):

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = SUAVE.Analyses.Vehicle()

    # ------------------------------------------------------------------
    #  Stability Analysis
    stability = SUAVE.Analyses.Stability.AVL()
    # # stability.settings.filenames.avl_bin_name = r"C:\Users\nmb48\Documents\GitHub\SUAVE\avl3.52\avl.exe"  # NILS: set AVL executable name
    # stability.settings.filenames.avl_bin_name = r"C:\Users\nmb48\Documents\GitHub\SUAVE\jvl2.16\jvl.exe"
    # =============================================================================
    stability.settings.filenames.avl_bin_name = r"C:\Users\nmb48\Downloads\JVL2.16\bin\jvl.exe"
    # =============================================================================
    #stability.settings.spanwise_vortex_density                  = 3
    stability.geometry = vehicle
    analyses.append(stability)
    
    return analyses

def vehicle_setup():
    
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------    
    
    vehicle = SUAVE.Vehicle()
    
    return vehicle

# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def configs_setup(vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------
    configs = SUAVE.Components.Configs.Config.Container()

    base_config = SUAVE.Components.Configs.Config(vehicle)
    base_config.tag = 'base'
    configs.append(base_config)

    return configs

    
    