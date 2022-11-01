
from open.optim.LP_MILP_wpiecewise import *
import time
import pandas as pd
import numpy as np
from sys import argv
import os


"""This is a tutorial to show how to estimate optimal operating conditions 
for a given set up inputs
In this case inputs include Concentrations for both the substrate and the product and equilibrium constant
as a result overall thermodynamic displacement (gamma) is also defined"""

"Characteristic concentration is  1e-4 M  1e4s-1/1e8s-1M-1"
"Hence S=3.0 is equivalent to 3e-4M 0.3 mM"

t=time.time()
S_concentration=3.0
P_concentration=2.0
q_equilibrium=2.0
gamma_overall=P_concentration/S_concentration/q_equilibrium

if gamma_overall>1:
    print('Reaction proceeds in the reverse direction please reaarrange inputs')
obj_primal_milp, variables_primal_milp, var_analysis = milp_problem_3step(gamma_overall, q=q_equilibrium,
                                                                          S=S_concentration, variability_analysis=False)
df_milp_variables = pd.DataFrame(variables_primal_milp, index=[0])

elapsed = time.time() - t
print('time for optim', elapsed)


