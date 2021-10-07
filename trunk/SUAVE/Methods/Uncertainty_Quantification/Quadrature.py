## @ingroup Methods-Uncertainty_Quantification
# Quadrature.py
#
# Created:  Oct 2021, R. Erhard
# Modified: 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Improve the quality of the figures
mpl.rcParams['figure.dpi'] = 200
import seaborn as sns
sns.set_style('white')
sns.set_context('talk')

def Quadrature(f , k=5):
    # approximate f by a family of orthogonal polynomials of degree k
    get_orthogonal_polynomials(k)
    
    return


def get_orthogonal_polynomials(k):
    PI = np.ndarray((k,))
    PI[0]  = 1
    PI[-1] = 0
    
    for i in range(k-1)+1:
        PI[i] = (1/s)*(beta_k*PI[i-1] + alpha_k*PI[i] + beta_kplus*PI[i+1] )
    return