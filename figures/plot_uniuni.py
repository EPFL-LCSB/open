import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from figures.main_figures import *
import os
import numpy as np

'''read data frame'''
"k_eq=1"
df_q2 = pd.read_hdf('/open/data/UniUni/df_uninuni_milp_q1_161120_finer.h5')
"k_eq=2"
df_q2 = pd.read_hdf('/open/data/UniUni/df_uninuni_milp_011120.h5')
"k_eq=3"
df_q2 = pd.read_hdf('/open/data/UniUni/df_uninuni_milp_q3_161120_finer.h5')
"k_eq=5"
df_q2 = pd.read_hdf('/open/data/UniUni/df_uninuni_milp_q5_161120_finer.h5')


q_equilibrium=(df_q2.k1f*df_q2.k2f*df_q2.k3f/(df_q2.k1b*df_q2.k2b*df_q2.k3b)).iloc[0].round(3)

df_q2_base = pd.read_hdf('/open/data/UniUni/df_uninuni_milp_011120.h5')


output_file='./uni_uni_figures_paper'

if not os.path.exists(output_file):
    os.makedirs(output_file)
"DATA FROM OPTIMIZATION"
"put the path for the data"

df_q2['kcat'] = df_q2.k2f * df_q2.k3f / (
            df_q2.k2f + df_q2.k3f + df_q2.k2b)
df_q2['nu_calculated'] = df_q2.v / df_q2['kcat']
df_q2['residual'] = (1 - df_q2.overall_gamma) * (1 - df_q2.E)
df_q2['KM_S']=(df_q2.k2f*df_q2.k3f+df_q2.k1b*df_q2.k3f+df_q2.k1b*df_q2.k2b)/ \
(df_q2.k1f*(df_q2.k2f + df_q2.k3f + df_q2.k2b))
df_q2['KM_P']=(df_q2.k2f*df_q2.k3f+df_q2.k1b*df_q2.k3f+df_q2.k1b*df_q2.k2b)/ \
(df_q2.k3b*(df_q2.k2f + df_q2.k1b + df_q2.k2b))
df_q2['kcat_kms']=df_q2.kcat/df_q2.KM_S
df_q2['kcat_back']= df_q2.k1b*df_q2.k2b/(df_q2.k2f+df_q2.k1b+df_q2.k2b)

limit_cut=0.99
gr_general_3step=['k1f_range','k2f_range','k3f_range','k1b_range','k2b_range','k3b_range']

assign_k_ranges(df_q2,lim=(limit_cut))
add_label(df_q2,gr_general_3step)
assign_k_ranges(df_q2,lim=(limit_cut))
add_label(df_q2,gr_general_3step)
df_q2.label.nunique()

level_number=11
df_q2['sat']=1-df_q2.E
df_q2_base['sat']=1-df_q2_base.E
base=np.linspace(0,1,15)
plt.figure()
main_contour_plot_fluxes_3step(df_q2,base,q_equilibrium,label='sat',level_num=level_number,sampled_points=True)
plt.xlim([0,5.1])
plt.ylim([0,5.1])
plt.savefig(output_file+'/main_contour_plot_sat_q{}_limit_notscaled.svg'.format(q_equilibrium))
plt.close()

level_number=11
df_q2['sat']=1-df_q2.E
#df_q2_base['sat']=1-df_q2_base.E
base=np.linspace(df_q2.v.min(),df_q2.v.max(),15)
plt.figure()
main_contour_plot_fluxes_3step(df_q2,base,q_equilibrium,label='v',level_num=level_number,sampled_points=True)
plt.xlim([0,5.1])
plt.ylim([0,5.1])
plt.savefig(output_file+'/main_contour_plot_flux_q{}.svg'.format(q_equilibrium))
plt.close()
