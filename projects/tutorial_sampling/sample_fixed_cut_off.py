"workflow kind of to sample from close to optimal solutions"

from open.optim.LP_MILP_wpiecewise import *
from open.optim.feasibility_problem import *
from open.sampling.manipulation import *
from open.sampling.sampling import *
import time
import os
from open.utils.postprocess import calculate_kinetic_params_uniuni, calculate_rate_constants_uniuni
output_file='./samples_for_all'
if not os.path.exists(output_file):
    os.makedirs(output_file)
"set your parameters"
q_equilibrium = 2.0
S_concentration=1.0
P_concentration=1.0
gamma_overall=P_concentration/S_concentration/q_equilibrium
t=time.time()
'run milp problem and find the optimal solution'
obj_primal_milp_3step, variables_primal_milp_3step,var_analysis= milp_problem_3step(gamma_overall, q=q_equilibrium, S=S_concentration, variability_analysis=False)

'for each close to optimal'
'run variability analysis to have the ranges mainly for gammas when close to optimal  e.g. 0.90% of the optimal '

close_to_optimal=0.90
obj_feasibility, variables_primal_milp_3step_feasibility,model,\
var_analysis_milp= milp_problem_3step_feasibility(obj_primal_milp_3step*close_to_optimal,obj_primal_milp_3step*(close_to_optimal),
                                                gamma_overall, q=q_equilibrium, S=S_concentration, variability_analysis=True)


df=pd.DataFrame(variables_primal_milp_3step, index=[0])
#now run with gammas fixed and get the variability
obj, vars, model_lin, var_analysis = linear_problem_for_3_step_linearized_feasibility(obj_primal_milp_3step*close_to_optimal, obj_primal_milp_3step*(close_to_optimal),float(df.gamma_1), float(df.gamma_2),
                                                                        float(df.gamma_3), gamma_overall, q=q_equilibrium, S=S_concentration,
                                                                        variability_analysis=True)




'sample gammas from the variability analysis of the milp'
n_sample_gammas=10
gamma_1=np.random.uniform(var_analysis_milp.loc['gamma_1'].minimum  , var_analysis_milp.loc['gamma_1'].maximum  , n_sample_gammas)
a=np.random.uniform(var_analysis_milp.loc['a'].minimum  , var_analysis_milp.loc['a'].maximum  , n_sample_gammas)
'here meshing takes place'
x, y = np.meshgrid(gamma_1, a)
gamma_1=x.reshape(n_sample_gammas**2,1)
a=y.reshape(n_sample_gammas**2,1)
list_=[]
gamma_3=a/gamma_1
gamma_2=gamma_overall/a

for i in range(n_sample_gammas**2):
    print(i)
    #i=15
    gamma1=gamma_1[i]
    gamma2=gamma_2[i]
    gamma3=gamma_3[i]
    #how many samples from the feasible set
    obj_2, vars_2, model_lin_2, var_analysis_2 = linear_problem_for_3_step_linearized_feasibility(
        obj_primal_milp_3step * close_to_optimal, obj_primal_milp_3step * (close_to_optimal), gamma1[0],
        gamma2[0],
        gamma3[0], gamma_overall, q=q_equilibrium, S=S_concentration,
        variability_analysis=True)
    if model_lin_2.status!='infeasible':
        'sampling does not work with integer variables'
        'fix the upper and lower bounds'
        tight_model = apply_generic_variability(model_lin_2, var_analysis_2)
        continuous_model = strip_from_integer_variables(tight_model)
        n_sample_hr=10
        sampling = sample(model_lin_2, n_sample_hr, processes = 1,method='achr')
        sampling['gamma_1']=gamma1*np.ones(n_sample_hr)
        sampling['gamma_2']=gamma2*np.ones(n_sample_hr)
        sampling['gamma_3']=gamma3*np.ones(n_sample_hr)

    #linear_problem_for_3_step_linearized_feasibility(obj_primal_milp_3step * 0.90, gamma1, gamma2,

     #                                                gamma3, gamma_overall, q=2.0, S=S_concentration,
      #                                               variability_analysis=True)
        list_.append(sampling)

    else:
        continue

'df containing samples'
df_samples = pd.concat(list_, ignore_index=True)

elapsed = time.time() - t
print('time for optim', elapsed)

df_samples['overall_gamma']=df_samples.gamma_1*df_samples.gamma_2*df_samples.gamma_3
df_samples['P']=P_concentration*np.ones(df_samples.shape[0])
df_samples['S']=S_concentration*np.ones(df_samples.shape[0])
calculate_rate_constants_uniuni(df_samples)
calculate_kinetic_params_uniuni(df_samples)
df_samples.to_hdf(
    output_file+'/sample_close_to_optimal_cutoff_{}_S_{}_P_{}_q_{}_gamma_{}_size_{}.h5'.format(
        close_to_optimal, S_concentration, P_concentration, q_equilibrium,df_samples.shape[0],\
        gamma_overall), key='s')
