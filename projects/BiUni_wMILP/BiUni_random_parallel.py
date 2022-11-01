
from open.optim.LP_MILP_random import *
import time
import pandas as pd
import numpy as np
from sys import argv
from open.utils.postprocess import convert_results_to_df, remove_log_calc_overallgamma,plot_convergence_graph

"This function was when splitting ratio is fixed"
"random bi uni mechanism is studied in test split_ratio script" \
"this script is not included in the manuscript"
_, gamma_overall,alpha = argv

gamma_overall = float(gamma_overall)
alpha = float(alpha)

t=time.time()
q_equilibrium = 2.0
limit=5.0

A_conc = np.logspace(np.log10(0.001), np.log10(5), 10)
# A_conc=np.array([5.0])

# to sample equally spaced points in isolines and depending on the range reduce points
limit = 5.0
P_concentration = np.linspace(0.001,5,10)
#P_concentration = np.linspace(1,1,1)

q_equilibrium = 2.0

for P_conc in P_concentration:
    B_min = P_conc / q_equilibrium / gamma_overall / limit

    if B_min < 2.5:
        A_conc = np.logspace(np.log10(B_min), np.log10(limit), 20)

    elif ((B_min > 4.0) & (B_min < 4.8)):
        A_conc = np.logspace(np.log10(B_min), np.log10(limit), 5)

    elif ((B_min >= 4.8) & (B_min <= 5.0)):
        A_conc = np.logspace(np.log10(B_min), np.log10(limit), 3)

    elif ((B_min > 5.0)):
        A_conc = np.array([])

    else:
        A_conc = np.logspace(np.log10(B_min), np.log10(limit), 10)


    B_conc = P_conc / q_equilibrium / gamma_overall / A_conc

    for A_concentration in A_conc[(B_conc >= 0.001) & (B_conc <= 5.0)]:

        B_concentration = P_conc/(A_concentration*q_equilibrium*gamma_overall)

        obj_primal_milp_4step, variables_primal_milp_4step,var_analysis = milp_problem_4step_biuni_random(gamma_overall, q=q_equilibrium,
                                                                                  S=A_concentration, P=P_conc,alpha=alpha,
                                                                                  variability_analysis=False)

        df_milp_variables = pd.DataFrame(variables_primal_milp_4step,index=[0])
    #print(df_milp_variables.T)
        df_milp_variables.to_hdf(
    './output_BiUni_parallel_random_05_090120/maximize_flux_milp_{}_A_{}_B_{}_P_{}_q_{}_alpha_{}.h5'.format(
        gamma_overall, A_concentration,B_concentration,P_conc,
        q_equilibrium,alpha), key='s')


elapsed = time.time() - t
print('time for optim', elapsed)