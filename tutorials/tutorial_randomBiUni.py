from open.optim.LP_MILP_random import *
import time
import pandas as pd


"random ordered bi-uni mechanism"
P_concentration=1.0
A_concentration=3.0
B_concentration=3.0
q_equilibrium=2.0
df_st = pd.DataFrame(columns=['A', 'B', 'P', 'q', 'alpha_max', 'alpha_min', 'v_net', 'gamma_ov'])
gamma_overall = P_concentration / A_concentration / q_equilibrium / B_concentration

t = time.time()
obj_primal_milp_4step_random_split, variables_primal_milp_4step_random_split, var_analysis_feasibility = milp_problem_4step_biuni_random_split_ratio(
    gamma_overall, q=q_equilibrium,
    S=A_concentration, P=P_concentration,
    variability_analysis=True)
elapsed = time.time() - t
print('time for optimization', elapsed)
split_max = var_analysis_feasibility.loc['v_upper', 'maximum'] / obj_primal_milp_4step_random_split
split_min = var_analysis_feasibility.loc['v_upper', 'minimum'] / obj_primal_milp_4step_random_split
df_milp_variables_random_split = pd.DataFrame(variables_primal_milp_4step_random_split, index=[0])

values_to_add = {'A': A_concentration, 'B': B_concentration, 'P': P_concentration, 'q': q_equilibrium,
                 'alpha_max': split_max, \
                 'alpha_min': split_min, 'v_net': obj_primal_milp_4step_random_split, 'gamma_ov': gamma_overall}
row_to_add = pd.Series(values_to_add, name=str(0))
df_st = df_st.append(row_to_add)
