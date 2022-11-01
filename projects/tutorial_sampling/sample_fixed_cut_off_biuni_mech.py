"""workflow to sample the random ordered mechanism
at equal concentrations of both substrates"""

from open.optim.LP_MILP_random import *
from open.optim.feasibility_problem import *
from open.sampling.manipulation import *
from open.sampling.sampling import *
import time
import os
from open.utils.postprocess import calculate_kinetic_params_uniuni, calculate_rate_constants_uniuni

"set your parameters concentrations"
"gamma overall should be <1"
S_concentration =  3.0004
B_concentration =  3.0004
q_equilibrium = 2.0
P_concentration = 1.0
gamma_overall=P_concentration/S_concentration/q_equilibrium/B_concentration

'make the required directory'
#output_file='./samples_for_variability_split_ratio_cut_off_A_{}_B_{}_P_{}_q_{}_190421'.format(S_concentration,B_concentration,P_concentration,q_equilibrium)
output_file='./samples_for_variability_split_ratio_cut_off_A_{}_B_{}_P_{}_q_{}_211221'.format(S_concentration,B_concentration,P_concentration,q_equilibrium)

if not os.path.exists(output_file):
    os.makedirs(output_file)


t=time.time()
'run milp problem and find the optimal solution'

obj_primal_milp_4step_random_split, variables_primal_milp_4step_random_split,df_var_analysis = milp_problem_4step_biuni_random_split_ratio(
    gamma_overall, q=q_equilibrium,
    S=S_concentration, P=P_concentration,
    variability_analysis=False)

df_optimal = pd.DataFrame(variables_primal_milp_4step_random_split, index=[0])

'for each close to optimal'
'run variability analysis to have the ranges mainly for gammas when close to optimal  e.g. 0.90% of the optimal '

close_to_optimal_list=[1.0]

for close_to_optimal in close_to_optimal_list:
    flexibility=0.0
    obj_feasibility, variables_primal_milp_biuni_random,model,\
    var_analysis_milp= milp_problem_4step_biuni_random_split_ratio_feasibility(
        obj_primal_milp_4step_random_split * close_to_optimal,obj_primal_milp_4step_random_split * (close_to_optimal+flexibility), gamma_overall, q=q_equilibrium,
        S=S_concentration, P=P_concentration,
        variability_analysis=True)
    df=pd.DataFrame(variables_primal_milp_biuni_random, index=[0])

#here first get out of the loop if split ratio is not variablee
    if (var_analysis_milp['maximum'].v_upper-var_analysis_milp['minimum'].v_upper)<=1e-9:
        print('The solution is unique no need to sample')
        break

    else:
        var_analysis_milp=var_analysis_milp.drop(['v_upper','v_lower'])
    n_sample_gammas=10
    if (var_analysis_milp['maximum']-var_analysis_milp['minimum']).all()<=1e-9:
        n_sample_gammas=1

    print('chosen n gammas is')
    print(n_sample_gammas)
    gamma_1=np.random.uniform(var_analysis_milp.loc['gamma_1'].minimum  , var_analysis_milp.loc['gamma_1'].maximum  , n_sample_gammas)
    a=np.random.uniform(var_analysis_milp.loc['a'].minimum  , var_analysis_milp.loc['a'].maximum  , n_sample_gammas)
    gamma_3=np.random.uniform(var_analysis_milp.loc['gamma_3'].minimum  , var_analysis_milp.loc['gamma_3'].maximum  , n_sample_gammas)
    gamma_5=np.random.uniform(var_analysis_milp.loc['gamma_5'].minimum  , var_analysis_milp.loc['gamma_5'].maximum  , n_sample_gammas)

    'here meshing takes place'
    'to have mesh of each possible gammas'
    x,y,z,t = np.meshgrid(gamma_1, a, gamma_3,gamma_5)
    gamma_1=x.reshape(n_sample_gammas**4,1)
    a=y.reshape(n_sample_gammas**4,1)
    gamma_3=z.reshape(n_sample_gammas**4,1)
    gamma_5=t.reshape(n_sample_gammas**4,1)
    gamma_2=a/gamma_1
    gamma_6=a/gamma_5
    gamma_4=gamma_overall/a/gamma_3

    list_=[]


    for i in range(n_sample_gammas**4):
        print(i)
        #i=15
        gamma1=float(gamma_1[i])
        gamma2=float(gamma_2[i])
        gamma3=float(gamma_3[i])
        gamma4=float(gamma_4[i])
        gamma5=float(gamma_5[i])
        gamma6=float(gamma_6[i])
        #how many samples from the feasible set
        obj_2, vars_2, model_lin_2, var_analysis_2 = milp_problem_4step_biuni_random_split_ratio_feasibility_with_gammas(obj_primal_milp_4step_random_split * close_to_optimal,obj_primal_milp_4step_random_split *( close_to_optimal+flexibility),gamma1,gamma2, \
                                                                        gamma3,gamma4,gamma5, gamma6,gamma_overall, q=q_equilibrium, S=S_concentration, P=P_concentration,variability_analysis=True)

        if model_lin_2.status!='infeasible':
            'sampling does not work with integer variables'
            'fix the upper and lower bounds'
            tight_model = apply_generic_variability(model_lin_2, var_analysis_2)
            continuous_model = strip_from_integer_variables(tight_model)
            n_sample_hr=50
            sampling = sample(model_lin_2, n_sample_hr, processes = 1,method='achr')
            sampling['gamma_1']=gamma1*np.ones(n_sample_hr)
            sampling['gamma_2']=gamma2*np.ones(n_sample_hr)
            sampling['gamma_3']=gamma3*np.ones(n_sample_hr)
            sampling['gamma_4']=gamma4*np.ones(n_sample_hr)
            sampling['gamma_5']=gamma5*np.ones(n_sample_hr)
            sampling['gamma_6']=gamma6*np.ones(n_sample_hr)


            list_.append(sampling)

        else:
            continue

    'df containing samples'
    df_samples = pd.concat(list_, ignore_index=True)

    elapsed = time.time() - t
    print('time for optim', elapsed)

    df_samples['overall_gamma']=df_samples.gamma_1*df_samples.gamma_2*df_samples.gamma_3*df_samples.gamma_4
    df_samples['P']=P_concentration*np.ones(df_samples.shape[0])
    df_samples['B']=B_concentration*np.ones(df_samples.shape[0])
    df_samples['S']=S_concentration*np.ones(df_samples.shape[0])
    df_samples['split_ratio']=df_samples.v_upper/df_samples.v
#to save the data
    df_samples.to_hdf(
        output_file+'/sample_close_to_optimal_cutoff_{}_A_{}_B_{}_P_{}_q_{}_gamma_{}_size_{}.h5'.format(
            close_to_optimal, S_concentration,B_concentration, P_concentration, q_equilibrium,df_samples.shape[0],\
            gamma_overall), key='s')


