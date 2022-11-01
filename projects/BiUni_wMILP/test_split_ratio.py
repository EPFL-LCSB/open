from open.optim.LP_MILP_wpiecewise import *
from open.optim.LP_MILP_random import *
import time
from sys import argv
import os


''' to run for general bi uni mechanism
with variability analysis on the splitting ratio
To generate data for Figure 3
'''

t=time.time()
#
"Inputs as arguments and run the bash script"
_, p_conc = argv

P_concentration=float(p_conc)

"or else you can give it in the script"
#P_concentration=1.0


"Here maximal values given directly as normalized values" \
"if not check normalize_parameters function in utils/preprocessing and manuscript"
xmax=5.0
ymax=5.0
xx, yy = np.meshgrid(np.linspace((0.001), (xmax), 25),
                     np.linspace((0.001), (ymax), 25))


xx, yy = np.meshgrid(np.linspace((0.001), (xmax), 3),
                     np.linspace((0.001), (ymax), 3))
q_equilibrium=2.0
df=pd.DataFrame()
df['A']=xx.ravel()
df['B']=yy.ravel()
df['AB']=df.A*df.B
#P_concentration = 1.0
df_test=df[df.AB>=P_concentration/q_equilibrium/0.95]
n=df_test.shape[0]
df_st = pd.DataFrame(columns=['A', 'B', 'P', 'q', 'alpha_max', 'alpha_min', 'v_net', 'gamma_ov'])


output_file='./output_130322_P_{}_0.001_res'.format(P_concentration)
#output_file='./milp_3step_011120'
output_var=output_file+'/var_analysis'
output_df=output_file+'/df'

if not os.path.exists(output_file):
    os.makedirs(output_file)
if not os.path.exists(output_var):
    os.makedirs(output_var)
if not os.path.exists(output_df):
    os.makedirs(output_df)


for i in range(n):
    print(i)
    S_concentration = df_test['A'].iloc[i]
    B_concentration = df_test['B'].iloc[i]
    #
    #P_concentration = 4.0
    gamma_overall = P_concentration / S_concentration / q_equilibrium / B_concentration
    t3 = time.time()
    obj_primal_milp_4step_random_split, variables_primal_milp_4step_random_split,var_analysis_feasibility = milp_problem_4step_biuni_random_split_ratio(
        gamma_overall, q=q_equilibrium,
        S=S_concentration, P=P_concentration,
        variability_analysis=True)
    elapsed = time.time() - t3
    print('time for optim', elapsed)



    split_max = var_analysis_feasibility.loc['v_upper', 'maximum'] / obj_primal_milp_4step_random_split
    split_min = var_analysis_feasibility.loc['v_upper', 'minimum'] /obj_primal_milp_4step_random_split
    print(split_max)


    df_milp_variables_random_split = pd.DataFrame(variables_primal_milp_4step_random_split, index=[0])
    df_milp_variables_random_split.to_hdf(output_df+'/df_milp_A_{}_B_{}_P_{}_gamma_{}.h5'.format(S_concentration,B_concentration,P_concentration,gamma_overall), key='s')
    var_analysis_feasibility.to_hdf(output_var+'/df_milp_var_analysis_A_{}_B_{}_P_{}_gamma_{}.h5'.format(S_concentration,B_concentration,P_concentration,gamma_overall), key='s')
    values_to_add = {'A': S_concentration, 'B': B_concentration, 'P':P_concentration,'q':q_equilibrium, 'alpha_max':split_max, \
                     'alpha_min':split_min,'v_net':obj_primal_milp_4step_random_split,'gamma_ov':gamma_overall}
    row_to_add = pd.Series(values_to_add, name=str(i))

    df_st = df_st.append(row_to_add)

df_st.to_hdf(output_file+'/df_split_ratio.h5', key='s')

