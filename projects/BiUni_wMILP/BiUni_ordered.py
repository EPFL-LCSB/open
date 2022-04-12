
from open.optim.LP_MILP_wpiecewise import *
import time
import pandas as pd
import numpy as np
from sys import argv
from open.utils.postprocess import convert_results_to_df, remove_log_calc_overallgamma,plot_convergence_graph


_, gamma_overall = argv

gamma_overall = float(gamma_overall)

t=time.time()
q_equilibrium = 2.0
limit=5.0
#S_conc = np.linspace((0.001), (limit), 20)
#P_conc = S_conc* q_equilibrium *gamma_overall


A_conc = np.logspace(np.log10(0.001), np.log10(5), 10)
# A_conc=np.array([5.0])

# to sample equally spaced points in isolines and depending on the range reduce points
limit = 5.0
P_concentration = np.linspace(0.001,5,15)
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

        obj_primal_milp_4step, variables_primal_milp_4step = milp_problem_4step_biuni(gamma_overall, q=q_equilibrium,
                                                                                  S=A_concentration, P=P_conc,
                                                                                  variability_analysis=False)

        df_milp_variables = pd.DataFrame(variables_primal_milp_4step,index=[0])
    #print(df_milp_variables.T)
        df_milp_variables.to_hdf(
    './output_BiUni_parallel_161120/maximize_flux_milp_{}_A_{}_B_{}_P_{}_q_{}.h5'.format(
        gamma_overall, A_concentration,B_concentration,P_conc,
        q_equilibrium), key='s')


elapsed = time.time() - t
print('time for optim', elapsed)



#try here bi-bi
from open.optim.MILP_BiBi_ordered import *
import time
import pandas as pd
import numpy as np
from sys import argv

S1_conc=1.0
S2_conc=1,0
P1_conc=1.0
q_equilibrium=1.1
gamma_ov=1/1/q_equilibrium
t=time.time()

obj_primal_milp_5step, variables_primal_milp_5step,var_analysis = milp_problem_5step_bibi(gamma_ov, q=2.0, S1=S1_conc,S2=S2_conc,P1=P1_conc, variability_analysis=False)

elapsed = time.time() - t
print('time for optim', elapsed)
df_milp= pd.DataFrame(variables_primal_milp_5step,index=[0])

df_milp.k1f*df_milp.k2f*df_milp.k3f*df_milp.k4f*df_milp.k5f/(df_milp.k1b*df_milp.k2b*df_milp.k3b*df_milp.k4b*df_milp.k5b)

(df_milp.z6*df_milp.z7*df_milp.z8*df_milp.z9*df_milp.z10)/(df_milp.z1*df_milp.z2*df_milp.z3*df_milp.z4*df_milp.z5)

df_milp.E+df_milp.EA+df_milp.EAB+df_milp.EPQ+df_milp.EP
