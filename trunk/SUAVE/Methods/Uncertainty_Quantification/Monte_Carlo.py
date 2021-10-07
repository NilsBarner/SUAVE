## @ingroup Methods-Uncertainty_Quantification
# Monte_Carlo.py
#
# Created:  Sep 2021, R. Erhard
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


def mean_and_variance_plot(f):
    
    # Number of independent MC runs
    num_mc = 1
    
    # Maximum number of samples to take per MC run
    max_n  = 100
    
    # Array holding # of samples
    all_ns = np.arange(1, max_n+1)
    fig, ax = plt.subplots()
    
    # Plot the true expected value for reference
    true_mean, = ax.plot(all_ns, [0.09] * max_n, color='black', linewidth=3.0)
    for k in range(num_mc):
        S_n_vals, V_n_vals = get_mc_estimate(f, max_n = max_n)
        ax.plot(all_ns, S_n_vals, color = sns.color_palette()[3], lw=1)
        
        # Give an initial error bar just for the first MC run
        if k == 0:
            # The lower bound
            l = S_n_vals - 2. / np.sqrt(all_ns) * np.sqrt(V_n_vals)
            
            # The upper bound
            u = S_n_vals + 2. / np.sqrt(all_ns) * np.sqrt(V_n_vals)
            var_rng = ax.fill_between(all_ns, l, u, color = sns.color_palette()[3], alpha = 0.25)
    
    ax.set_xlabel('$n$')
    ax.set_ylabel('Sample Mean, $\\overline{C_T}$')# usetex = True);
    ax.set_title('Uncertainty in Inflow Angle')
    err_bar_lab = 'Confidence Region'
    ax.legend([true_mean, var_rng],['Experimental Result',err_bar_lab], loc='best')
    
    plt.savefig('UQ_PSI_CT.png',bbox_inches='tight', dpi=300)
    
    plt.show()
    return



def get_mc_estimate(f, max_n =1000,  sampler = np.random.rand):
    """
    Return the MC estimate of the 1D integral of a user-defined f function.
    
    :param func:    The user-defined function to integrate
    :param max_n:   Maximum number of sample
    :param sampler: A function that samples the X_i;
                    the default is np.random.rand, which is the uniform distribution on [0,1)
    """
    S_n = np.ndarray((max_n,))  # Sample mean; an array with the running estimates of the expecation
    V_n = np.ndarray((max_n,))  # Sample variance;  an array with the running estimates of Var(f(X))
    s  = 0.                     # Variable to keep track of the sum (for expected val)
    s2 = 0.                     # Variable to keep track of the sum (for var)
    for i in range(max_n):
        # Sample X_i
        X_i = sampler(1)
        
        # Update the sum
        s += f(X_i)
        
        # Current approx of the original integral (expected value)
        S_n[i] = s / (i + 1)
        s2 += (f(X_i) - S_n[i]) ** 2
        
        # The current approximation of the epistemic variance
        V_n[i] = s2 / (i + 1)
        
    return(S_n, V_n)


