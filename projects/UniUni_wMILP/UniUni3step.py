
from open.optim.LP_MILP_wpiecewise import *
import time
import pandas as pd
import numpy as np
from sys import argv
import os
from open.utils.postprocess import convert_results_to_df, remove_log_calc_overallgamma,plot_convergence_graph


_, gamma_overall = argv

gamma_overall = float(gamma_overall)

t=time.time()
q_equilibrium =2.0
limit=5.0
S_conc = np.linspace((0.001), (limit), 50)
P_conc = S_conc* q_equilibrium *gamma_overall

'put your output file name'
output_file='./output_3step_fwd_170321_q_{}'.format(q_equilibrium)
if not os.path.exists(output_file):
    os.makedirs(output_file)

for S_concentration in S_conc[(P_conc >= 0.001) & (P_conc <= 6.0)]:
    obj_primal_milp, variables_primal_milp= milp_problem_3step(gamma_overall, q=q_equilibrium, S=S_concentration, variability_analysis=False)

    df_milp_variables = pd.DataFrame(variables_primal_milp,index=[0])
    #print(df_milp_variables.T)
    df_milp_variables.to_hdf(
    output_file+'/maximize_flux_milp_{}_S_{}_q_{}.h5'.format(
        gamma_overall, S_concentration,
        q_equilibrium), key='s')


elapsed = time.time() - t
print('time for optim', elapsed)