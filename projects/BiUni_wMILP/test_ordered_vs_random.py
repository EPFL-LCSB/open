from open.optim.LP_MILP_wpiecewise import *
from open.optim.LP_MILP_random import *
import time

'script to compare ordered and general bi-uni mechanism' \
'not included in the manuscript'
t=time.time()
S_concs =  np.linspace(300,1000,10)
B_concentration =  5.0
q_equilibrium = 2.0
P_concentration = 3.0
list_=[]
df_st = pd.DataFrame(columns=['A', 'B', 'P', 'q', 'random', 'ordered_A', 'ordered_B','alpha', 'gamma_ov'])

for i in S_concs:
    df=pd.DataFrame()
    S_concentration=i
    gamma_overall = P_concentration / S_concentration / q_equilibrium / B_concentration
    t3 = time.time()
    obj_primal_milp_4step_random_split, variables_primal_milp_4step_random_split,var_analysis = milp_problem_4step_biuni_random_split_ratio(
    gamma_overall, q=q_equilibrium,
    S=S_concentration, P=P_concentration,
    variability_analysis=False)
    df_milp_variables_random_split = pd.DataFrame(variables_primal_milp_4step_random_split,index=[0])
    elapsed = time.time() - t3
    print('time for optim', elapsed)

    obj_primal_milp_4step_ordered_A, variables_primal_milp_4step_ordered_A = milp_problem_4step_biuni(gamma_overall, q=q_equilibrium,
                                                                              S=S_concentration, P=P_concentration,
                                                                               variability_analysis=False)

    obj_primal_milp_4step_ordered_B, variables_primal_milp_4step_ordered_B = milp_problem_4step_biuni(gamma_overall, q=q_equilibrium,
                                                                              S=B_concentration, P=P_concentration,
                                                                               variability_analysis=False)

    alpha=round(variables_primal_milp_4step_random_split['alpha'],2)
    values_to_add = {'A': S_concentration, 'B': B_concentration, 'P':P_concentration,'q':q_equilibrium, 'random':obj_primal_milp_4step_random_split, \
                     'ordered_A':obj_primal_milp_4step_ordered_A, 'ordered_B':obj_primal_milp_4step_ordered_B,'alpha':alpha , 'gamma_ov':gamma_overall}
    row_to_add = pd.Series(values_to_add, name='x')

    df_st = df_st.append(row_to_add)

df_st.to_hdf('./compare_ordered_vs_random/df_comparison_more_P3_A_highest_range_300_1000.h5', key='s')


df_st_1=pd.read_hdf('/open/projects/efficiency_tests_and_sampling/compare_ordered_vs_random/df_comparison_more_P3.h5')
df_st_2=pd.read_hdf('./compare_ordered_vs_random/df_comparison_more_P3_A_bigger.h5')
df_st_3=pd.read_hdf('./compare_ordered_vs_random/df_comparison_more_P3_A_high_range.h5')
df_st_4=pd.read_hdf('./compare_ordered_vs_random/df_comparison_more_P3_A_highest_range_300_1000.h5')



df_st=df_st_1.append(df_st_2)
df_st=df_st.append(df_st_3)
df_st=df_st.append(df_st_4)
df_st=df_st_1
plt.rcParams.update({'font.size':14,
                     'font.family':'Calibri',
                     'legend.fontsize':10,
                    'xtick.labelsize' : 14,
                    'ytick.labelsize' : 14
                     })

import matplotlib.pyplot as plt
plt.figure()
plt.plot(df_st.A,df_st.ordered_A,ls='--',marker='o',label='ordered A',alpha=0.5)
plt.plot(df_st.A,df_st.ordered_B,ls='--',marker='o',label='ordered B',alpha=0.5)
#plt.plot(df_st.A,df_st.random,ls='--',marker='o',label='random',alpha=0.5)
plt.xlabel(r'[A]')
plt.ylabel(r'$\nu_{net}$')
plt.legend()
n=df_st.shape[0]
# for i in range(n):
#     plt.annotate(str(df_st.alpha.iloc[i]),  # this is the text
#              (df_st['A'].iloc[i], df_st['random'].iloc[i]),  # this is the point to label
#              textcoords="offset points",  # how to position the text
#              size=7,
#              xytext=(0, 5),  # distance from text to points (x,y)
#              ha='center')  # horizontal alignment can be left, right or center
#plt.xscale('log')
#plt.yscale('log')
plt.xlim([0,10.1])
plt.ylim([0,0.27])
plt.xticks(np.linspace(0,10,11))
plt.tight_layout()
plt.savefig('./compare_ordered_vs_random/comparison_plot_more_P3_low_range_ordered.svg')


plt.figure()
plt.plot(df_st.A,df_st.alpha,marker='o',ls='dashed')
plt.xlabel('[A]')
#plt.yscale('log')
plt.xscale('log')
plt.ylim([0,1])
#plt.yticks(np.linspace(0.1,1,10))
plt.ylabel(r'$\alpha$')
plt.tight_layout()
plt.savefig('./compare_ordered_vs_random/comparison_plot_P_05_split_ratio_vs_conc.svg')
