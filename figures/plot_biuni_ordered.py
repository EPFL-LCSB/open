import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from figures.main_figures import *
import os
import numpy as np
from figures.main_figures_for_random import *
import os
plt.rcParams.update({'font.size':16,
                     'font.family':'Calibri',
                     'legend.fontsize':16,
                    'xtick.labelsize' : 18,
                    'ytick.labelsize' : 18
                     })



map_dict={1:27,
2:23,
3:11,
4:12,
5:24,
6:13,
7:14,
8:31,
9:19,
10:5,
11:21,
12:7,
13:22,
14:8,
15:2,
16:10,
17:3,
18:4}


df = pd.read_hdf('/open/data/BiUni/df_biuni_milp_q2_040321.h5')
df=pd.read_hdf('/open/data/BiUni/df_biuni_milp_q2_P1_100321.h5')
df=pd.read_hdf('/open/data/BiUni/df_biuni_milp_q2_P1_251221.h5')
df=pd.read_hdf('/open/data/BiUni/df_biuni_milp_q2_P1_261221.h5')

#to screen the one with most results
interested_A=[df.A.unique()[i] for i in np.linspace(1,df.A.nunique()-1,20).astype(int)]
df=df[df.A.isin(interested_A)]

q_equilibrium=2.0
df['kcat']=df.k3f*df.k4f/(df.k3f+df.k4f+df.k3b)
output_file='./bi_uni_figures/q_{}'.format(q_equilibrium)
output_file='./biuni_ordered/biuni_040321'

output_file='./biuni_ordered/biuni_100321_P1'

output_file='./biuni_ordered/biuni_paper_verif_100422'

if not os.path.exists(output_file):
    os.makedirs(output_file)
"DATA FROM OPTIMIZATION"
P_conc=1.0
df_all_3step_unique_rounded=df[df.P==P_conc]
limit_cut=0.99
"Assign ranges based on the limit you define and classify accordingly"
#gr_general=['k1f_range', 'k2f_range', 'k3f_range', 'k4f_range','k5f_range', 'k6f_range', 'k1b_range', 'k2b_range', 'k3b_range', 'k4b_range','k5b_range', 'k6b_range']
#gr_general_according=['k1f_range', 'k3f_range', 'k5f_range', 'k6f_range', 'k1b_range', 'k3b_range','k5b_range', 'k6b_range']
gr_general=['k1f_range', 'k2f_range', 'k3f_range', 'k4f_range', 'k1b_range', 'k2b_range', 'k3b_range', 'k4b_range']


assign_k_ranges_random(df_all_3step_unique_rounded,lim=(limit_cut))
add_label(df_all_3step_unique_rounded,gr_general)
print(df_all_3step_unique_rounded.label.nunique())

"Classification from data points ONLY"
plot_colored_data_random(df_all_3step_unique_rounded, label='label', qeq=q_equilibrium)
plt.xlim([0,5.1])
plt.ylim([0, 5.1])
plt.tight_layout()
plt.savefig(output_file+'/261221_k_ranges_datapoints_q{}_lim{}_P_{}.svg'.format(q_equilibrium,limit_cut,P_conc  ))
plt.close()

#
print(df_all_3step_unique_rounded.label.nunique())
knn_k=11
df_trained=return_knn_trained_model_k_ranges_random(df_all_3step_unique_rounded,q_equilibrium,k=knn_k,group=gr_general)
'to normalize the added numbers '
df_trained['label_old']=df_trained.label
add_label(df_trained,'label_old')
print(df_trained.label.nunique())
plot_k_ranges_with_knn_random(df_trained,q_equilibrium,mapdict=map_dict,map=True)
plt.scatter(df_all_3step_unique_rounded.A,df_all_3step_unique_rounded.B,alpha=0.2)
plt.xlim([0,5.1])
plt.ylim([0,5.1])
plt.tight_layout()
plt.savefig(output_file+'/030122_NEW_NEW_k_ranges_knn_plot_q{}_lim{}_P_{}_k_{}.svg'.format(q_equilibrium,limit_cut,P_conc,knn_k))
plt.close()


plot_colored_data_random(df_trained, label='label', qeq=q_equilibrium)
plt.xlim([0,5.1])
plt.ylim([0, 5.1])
plt.tight_layout()
plt.savefig(output_file+'/030122_trained_k_ranges_datapoints_q{}_lim{}_P_{}.svg'.format(q_equilibrium,limit_cut,P_conc  ))
plt.close()




level_number=15
P_concentration=df_all_3step_unique_rounded.P.iloc[0]
base=np.linspace(0,0.27,level_number)
#df_all_3step_unique_rounded['sat']=1-df_all_3step_unique_rounded['E']
main_contour_plot_fluxes_random(df_all_3step_unique_rounded,base,q_equilibrium,P_conc=P_concentration,label='v',level_num=level_number,sampled_points=True)
# plt.xlim([-0.5,5.1])
# plt.ylim([-0.5,5.1])
#plt.ylim([0,10])
plt.tight_layout()
plt.savefig(output_file+'/main_contour_plot_FLUX_q{}_ordered_P_{}.svg'.format(q_equilibrium,P_concentration))
plt.close()


main_contour_plot_fluxes_random(df_all_3step_unique_rounded,df_all_3step_unique_rounded.kcat,q_equilibrium,label='kcat',level_num=level_number,sampled_points=True)
#plt.xlim([0,10])
#plt.ylim([0,10])
plt.savefig(output_file+'/main_contour_plot_kcat_q{}_ordered_P_{}.svg'.format(q_equilibrium,limit_cut,P_conc))
plt.close()

df_all = pd.read_hdf('/open/data/BiUni/df_biuni_milp_q2_040321.h5')
P_concentration=5.0
df_all_3step_unique_rounded=df_all[df_all.P==P_concentration]
level_number=21
base=np.linspace(0.5,1,15)
df_all_3step_unique_rounded['sat']=1-df_all_3step_unique_rounded['E']
main_contour_plot_fluxes_random(df_all_3step_unique_rounded,base,q_equilibrium,P_conc=P_concentration,label='sat',level_num=level_number,sampled_points=True)
# plt.xlim([-0.5,5.1])
# plt.ylim([-0.5,5.1])
#plt.ylim([0,10])
plt.tight_layout()
plt.savefig(output_file+'/main_contour_plot_OLD_sat_q{}_ordered_P_{}.svg'.format(q_equilibrium,P_concentration))
plt.close()

#
# fig, axes = plt.subplots()
#     #
#     # for i in a:
#     #     l = df[df[x] == i][y].values
#     #     list_.append(l)
#     #     x_axis.append(str(i))
# axes.violinplot(df.kcat, showmedians=True)
# axes.violinplot(df.kcat, showmedians=True)
# axes.yaxis.grid(True)
# axes.set_xlabel('Mechanisms')
# #axes.set_xticks(np.linspace(1, df[x].nunique(), df[x].nunique()))
# #axes.set_xticklabels(x_axis)
# plt.savefig(output_file+'/violin_q{}_ordered.png'.format(q_equilibrium))
# #   axes.set_ylabel('c_2nd')
from MCA.mca_plots import plot_violin_MCA_with_gamma_matplotlib

fig,axes=plot_violin_MCA_with_gamma_matplotlib(df,x='P',y='kcat')
#axes.set_ylim([0,1])
axes.set_ylabel(r'$k_{cat}$')
axes.set_xlabel(r'[P]')
axes.set_ylim([-1e-1,1.01])
plt.savefig(output_file+'/kcat_effect_violin_qeq_{}.svg')


## from here plot random
'plot random rate constant distributions'
'neglect this for now'
df_01=pd.read_hdf('/open/projects/efficiency_tests_and_sampling/output_190321_P_0.1_0.001_res/df_combined.h5')
df_05=pd.read_hdf('/open/projects/efficiency_tests_and_sampling/output_190321_P_0.5_0.001_res/df_combined.h5')
df_1=pd.read_hdf('/open/projects/efficiency_tests_and_sampling/output_190321_P_1.0_0.001_res/df_combined.h5')
df_2=pd.read_hdf('/open/projects/efficiency_tests_and_sampling/output_190321_P_2.0_0.001_res/df_combined.h5')
df_5=pd.read_hdf('/open/projects/efficiency_tests_and_sampling/output_190321_P_5.0_0.001_res/df_combined.h5')

df_st_01=pd.read_hdf('/open/projects/efficiency_tests_and_sampling/output_190321_P_0.1_0.001_res/df_split_ratio.h5')
df_st_05=pd.read_hdf('/open/projects/efficiency_tests_and_sampling/output_190321_P_0.5_0.001_res/df_split_ratio.h5')
df_st_1=pd.read_hdf('/open/projects/efficiency_tests_and_sampling/output_190321_P_1.0_0.001_res/df_split_ratio.h5')
df_st_2=pd.read_hdf('/open/projects/efficiency_tests_and_sampling/output_190321_P_2.0_0.001_res/df_split_ratio.h5')
df_st_5=pd.read_hdf('/open/projects/efficiency_tests_and_sampling/output_190321_P_5.0_0.001_res/df_split_ratio.h5')


df=df_05.append(df_1)
df=df.append(df_2)
df=df.append(df_5)
df=df.append(df_01)
df=df.reset_index()


df=df_01.append(df_1)
df=df.append(df_5)
df=df.reset_index()

df_st=df_st_05.append(df_st_1)
df_st=df_st.append(df_st_2)
df_st=df_st.append(df_st_5)
df_st=df_st.append(df_st_01)
df_st=df_st.reset_index()

from open.utils.postprocess import calculate_rate_constants_random_biuni
from figures.main_figures_for_random import *
calculate_rate_constants_random_biuni(df)

gr_general=['k1f_range', 'k3f_range', 'k5f_range', 'k6f_range', 'k1b_range', 'k3b_range','k5b_range', 'k6b_range']
gr_general=['k1f_range', 'k5f_range', 'k4b_range']

df['k1b']=df['k1b'].fillna(0)
df['k5b']=df['k5b'].fillna(0)
df['k2f']=df['k2f'].fillna(0)
df['k6f']=df['k6f'].fillna(0)

limit_cut=0.95

assign_k_ranges_random(df,lim=(limit_cut))
add_label(df,gr_general)
df.label.nunique()

output_file='./biuni_random_k_classes'
if not os.path.exists(output_file):
    os.makedirs(output_file)

P_conc=5.0
plot_colored_data_random(df[df.P==P_conc], label='label', qeq=q_equilibrium,mapdict=dict(),P_conc=P_conc)
plt.xlim([0,5.1])
plt.ylim([0, 5.1])
plt.tight_layout()
plt.savefig(output_file+'/k_ranges_binding_all_datapoints_q{}_lim{}_P_{}.svg'.format(q_equilibrium,limit_cut,P_conc))
plt.close()

assign_k_ranges_random(df_01,lim=(limit_cut))
add_label(df_01,gr_general)
df_01.label.nunique()
plot_colored_data_random(df_01, label='label', qeq=q_equilibrium,P_conc=0.1)
plt.xlim([0,5.1])
plt.ylim([0, 5.1])
plt.tight_layout()
plt.savefig(output_file+'/k_ranges_datapoints_q{}_lim{}_P_{}.svg'.format(q_equilibrium,limit_cut,0.1  ))
plt.close()


plt.figure()
plt.scatter(df[df.A!=df.B].gamma_ov,df[df.A!=df.B].A/df[df.A!=df.B].B,c=df[df.A!=df.B].alpha,cmap=cm.seismic,s=(df.P)*10)
plt.xlabel(r'$\Gamma$')
plt.ylabel(r'$\frac{[A]}{[B]}$')
plt.yscale('log')
plt.xscale('log')
plt.xlim([1e-3,1.1])
plt.colorbar()
plt.tight_layout()
plt.savefig(output_file+'/TRY_q{}_lim{}_P.svg'.format(q_equilibrium,limit_cut ))
plt.close()


plt.figure()
plt.scatter(df[df.P>=1.0].gamma_ov,df[df.P>=1.0].A/df[df.P>=1.0].B,c=df[df.P>=1.0].alpha,cmap=cm.seismic)
plt.xlabel(r'$\Gamma$')
plt.ylabel(r'$\frac{[A]}{[B]}$')
plt.yscale('log')
plt.xscale('log')
plt.xlim([1e-3,1.1])
plt.colorbar()
plt.tight_layout()
plt.savefig(output_file+'/TRY_high_conc_q{}_lim{}_P.svg'.format(q_equilibrium,limit_cut ))
plt.close()

plt.figure()
plt.scatter(df[df.P<=1.0].gamma_ov,df[df.P<=1.0].A/df[df.P<=1.0].B,c=df[df.P<=1.0].alpha,cmap=cm.seismic)
plt.xlabel(r'$\Gamma$')
plt.ylabel(r'$\frac{[A]}{[B]}$')
plt.yscale('log')
plt.xscale('log')
plt.xlim([1e-3,1.1])
plt.colorbar()
plt.tight_layout()
plt.savefig(output_file+'/TRY_low_conc_q{}_lim{}_P.svg'.format(q_equilibrium,limit_cut ))
plt.close()


plt.figure()
plt.scatter(df[df.P==1.0].gamma_ov,df[df.P==1.0].A/df[df.P==1.0].B,c=df[df.P==1.0].alpha,cmap=cm.seismic)
plt.xlabel(r'$\Gamma$')
plt.ylabel(r'$\frac{[A]}{[B]}$')
plt.yscale('log')
plt.xscale('log')
plt.xlim([1e-3,1.1])
plt.colorbar()
plt.tight_layout()
plt.savefig(output_file+'/TRY_p1_conc_q{}_lim{}_P.svg'.format(q_equilibrium,limit_cut ))
plt.close()

plt.figure()
df.loc[df.A==df.B,'alpha']=0.5*np.ones(df[df.A==df.B].shape[0])
df=df[df.P!=5]
plt.scatter(df.A/df.P,df.B/df.P,c=df.alpha,cmap=cm.seismic,alpha=0.5,s=10+(df.P)*10)
plt.xlabel(r'$\frac{[A]}{[P]}$')
plt.ylabel(r'$\frac{[B]}{[P]}$')
plt.yscale('log')
plt.xscale('log')
plt.xlim([1e-1,70])
plt.colorbar()
plt.tight_layout()
plt.savefig(output_file+'/TRY_pb_pa_conc_limited_less_q{}_lim{}_P.svg'.format(q_equilibrium,limit_cut ))
plt.close()

plt.figure()
plt.scatter(df.A/df.P,df.B/df.P,c=df.label,cmap=cm.jet,alpha=0.3,s=10+(df.P)*10)
plt.xlabel(r'$\frac{[A]}{[P]}$')
plt.ylabel(r'$\frac{[B]}{[P]}$')
plt.yscale('log')
plt.xscale('log')
#plt.xlim([1e-3,1.1])
plt.colorbar()
plt.tight_layout()
plt.savefig(output_file+'/TRY_clored_data_pb_pa_conc_llimited_less_q{}_lim{}_P.svg'.format(q_equilibrium,limit_cut ))
plt.close()


plt.figure()
df['sat']=1-df.E
plt.scatter(df.A/df.P,df.B/df.P,c=df.sat,cmap=cm.jet,alpha=0.3,s=10+(df.P)*10)
plt.xlabel(r'$\frac{[A]}{[P]}$')
plt.ylabel(r'$\frac{[B]}{[P]}$')
plt.yscale('log')
plt.xscale('log')
#plt.xlim([1e-3,1.1])
plt.colorbar()
plt.tight_layout()
plt.savefig(output_file+'/TRY_clored_data_saturation_pb_pa_conc_llimited_less_q{}_lim{}_P.svg'.format(q_equilibrium,limit_cut ))
plt.close()