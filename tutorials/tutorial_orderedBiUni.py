
from open.optim.LP_MILP_wpiecewise import *
import time
import pandas as pd
import numpy as np
from sys import argv

"Tutorial for a ordered BiUni mechanism for given inputs"


t=time.time()
P_conc=1.0
q_equilibrium=2.0
A_concentration=3.0
gamma_overall=0.5

B_concentration = P_conc/(A_concentration*q_equilibrium*gamma_overall)

obj_primal_milp_4step, variables_primal_milp_4step,var_analysis = milp_problem_4step_biuni(gamma_overall, q=q_equilibrium,
                                                                                  S=A_concentration, P=P_conc,
                                                                                  variability_analysis=False)

#convert to dataframe
df_milp_variables = pd.DataFrame(variables_primal_milp_4step,index=[0])

elapsed = time.time() - t
print('time for optim', elapsed)
