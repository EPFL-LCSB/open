import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from figures.main_figures import *
import os
import numpy as np
from open.utils.postprocess import calculate_rate_constants_random_biuni
from figures.main_figures_for_random import *
import os
plt.rcParams.update({'font.size':18,
                     'font.family':'Calibri',
                     'legend.fontsize':16,
                    'xtick.labelsize' : 18,
                    'ytick.labelsize' : 18
                     })
'split ratio dataframes'
plt.rcParams['svg.fonttype'] = 'none'



df_1=pd.read_hdf('/open/data/BiUni_general/output_010321_more_P1_2/df_split_ratio_010321.h5')
df_1_comb=pd.read_hdf('/open/data/BiUni_general/output_010321_more_P1_2/df_combined.h5')

df_01=pd.read_hdf('/open/data/BiUni_general/output_130322_P_0.1_0.001_res/df_split_ratio.h5')
df_01_comb=pd.read_hdf('/open/data/BiUni_general/output_130322_P_0.1_0.001_res/df_combined.h5')


df_5=pd.read_hdf('/open/data/BiUni_general/output_190321_P_5.0_0.001_res/df_split_ratio.h5')
df_5_comb=pd.read_hdf('/open/data/BiUni_general/output_190321_P_5.0_0.001_res/df_combined.h5')


"more resolved"
df_1_new_comb=pd.read_hdf('/open/data/BiUni_general/output_130322_P_1.0_0.001_res/df_combined.h5')
df_1_new=pd.read_hdf('/open/data/BiUni_general/output_130322_P_1.0_0.001_res/df_split_ratio.h5')


df_1_comb['tf1']=np.log(df_1_comb.gamma_1)/np.log(df_1_comb.gamma_ov)
df_1_comb['tf2']=np.log(df_1_comb.gamma_2)/np.log(df_1_comb.gamma_ov)
df_1_comb['tf3']=np.log(df_1_comb.gamma_3)/np.log(df_1_comb.gamma_ov)
df_1_comb['tf4']=np.log(df_1_comb.gamma_4)/np.log(df_1_comb.gamma_ov)
df_1_comb['tf5']=np.log(df_1_comb.gamma_5)/np.log(df_1_comb.gamma_ov)
df_1_comb['tf6']=np.log(df_1_comb.gamma_6)/np.log(df_1_comb.gamma_ov)

df_1_comb['tfcycle']=df_1_comb['tf1']+df_1_comb['tf2']

df_1_comb['v1f']=1/(1-df_1_comb['gamma_1'])*df_1_comb.alpha
df_1_comb['v2f']=1/(1-df_1_comb['gamma_2'])*df_1_comb.alpha
df_1_comb['v3f']=1/(1-df_1_comb['gamma_3'])#*df_q2.v
df_1_comb['v4f']=1/(1-df_1_comb['gamma_4'])#*df_q2.v
df_1_comb['v5f']=1/(1-df_1_comb['gamma_5'])*(1-df_1_comb.alpha)
df_1_comb['v6f']=1/(1-df_1_comb['gamma_6'])*(1-df_1_comb.alpha)


df_1_comb['v1b']=df_1_comb['gamma_1']/(1-df_1_comb['gamma_1'])*df_1_comb.alpha
df_1_comb['v2b']=df_1_comb['gamma_2']/(1-df_1_comb['gamma_2'])*df_1_comb.alpha
df_1_comb['v3b']=df_1_comb['gamma_3']/(1-df_1_comb['gamma_3'])#*df_q2.v
df_1_comb['v4b']=df_1_comb['gamma_4']/(1-df_1_comb['gamma_4'])#*df_q2.v
df_1_comb['v5b']=df_1_comb['gamma_5']/(1-df_1_comb['gamma_5'])*(1-df_1_comb.alpha)
df_1_comb['v6b']=df_1_comb['gamma_6']/(1-df_1_comb['gamma_6'])*(1-df_1_comb.alpha)




q_equilibrium = 2.0
#gamma_overall=P_concentration/S_concentration/q_equilibrium/B_concentration
'start from here'

fig, ax = plt.subplots(figsize=(10, 6))
#df_st=df_1
#indicate which data to plot
df_st=df_1
P_concentration = df_st.P.iloc[0].round(2)
scaled=True
output_file='./biuni_general/paper_review'


#output_file='./milp_3step_011120'
if not os.path.exists(output_file):
    os.makedirs(output_file)

xs = df_st['A']
ys = df_st['B']
list_color=[]
n=df_st.shape[0]
list_facecol=[]
list_edgecol=[]
list_linest=[]
list_linewidth=[]
difference=df_st.alpha_max-df_st.alpha_min

'second for loop for diverging color maps'
for i in range (n):
    if difference.iloc[i]<=1e-2 :
        #c='red'
        facecol = 'black'
        edgecol = 'black'
        c = (df_st.iloc[i].alpha_min+df_st.iloc[i].alpha_max)/2
        #c=0.5
        linest='-'
        linewidth=0.01
    else :
        #c='green'
        c=0.5
        facecol = 'none'
        edgecol = 'black'
        linest='--'
        linewidth=0.1


    # if (df_st.iloc[i].A>df_st.iloc[i].B):
    #     c=df_st.iloc[i].alpha_min
    # elif (df_st.iloc[i].A<df_st.iloc[i].B):
    #     c=1-df_st.iloc[i].alpha_min
    # else:
    #     c=0
    list_linest.append(linest)
    list_color.append(c)
    list_facecol.append(facecol)
    list_edgecol.append(edgecol)
    list_linewidth.append(linewidth)

#plt.figure()
# fig, ax = plt.subplots(figsize=(5.5, 4))
# fig, ax = plt.subplots(figsize=(8, 8))

# plot the points
x=np.linspace(0.0001,7,1000)
plt.plot(x,P_concentration/q_equilibrium/x,c='black',linewidth=1)
#overall_gamma_list=[0.02,0.03,0.05,0.07,0.12,0.2,0.3,0.7]
#p1
if P_concentration==1:
    overall_gamma_list=[0.02,0.03,0.05,0.12,0.3,0.7]
if P_concentration==5:
    overall_gamma_list=[0.12,0.18,0.3,0.7]
if P_concentration==0.01:
    overall_gamma_list=[0.002,0.003,0.005,0.01,0.02,0.07]

#overall_gamma_list=[0.002,0.005,0.008,0.02,0.03,0.05,0.12,0.3,0.9]

#for k in df_st[df_st.A==df_st.B].gamma_ov:
i=0
colors=['midnightblue','mediumblue','blue','royalblue','slategrey','lightsteelblue']
for k in overall_gamma_list:
    plt.plot(x,P_concentration/q_equilibrium/x/k,ls='dashed',label=r'$\Gamma=$'+str(round(k,2)),alpha=0.5,c='grey')
    i=i+1

plt.xlabel(r'$\widetilde{A}$')
plt.ylabel(r'$\widetilde{B}$')
# plt.ylabel('B')
plt.xlim([-1e-1,5.5])
plt.ylim([-1e-1,5.5])
plt.xticks(np.linspace(0,5,6))
plt.yticks(np.linspace(0,5,6))
n=df_st.shape[0]
if scaled:
    #before seismic
    cax=ax.scatter(xs,ys,label=None,facecolors=list_facecol,edgecolors=list_edgecol,s=100,c=list_color,cmap=cm.seismic,linestyle=list_linest, vmin=0,vmax=1.0)#cmap=cm.seismic
    cbar=fig.colorbar(cax,ticks=[0,0.5, 1])
    cbar.ax.set_yticklabels([str(0), '0.5', str(1)])
else:
    cax=ax.scatter(xs,ys,label=None,facecolors=list_facecol,edgecolors=list_edgecol,s=100,c=list_color,cmap=cm.seismic,linestyle=list_linest)#cmap=cm.seismic
    cbar=fig.colorbar(cax,ticks=[min(list_color),0.5, max(list_color)])
    cbar.ax.set_yticklabels([str(min(list_color).round(2)), '0.5', str(max(list_color).round(2))])

'if not scaled'
#cbar=fig.colorbar(cax,ticks=[min(list_color),0.5, max(list_color)])
#cbar.ax.set_yticklabels([str(min(list_color).round(2)), '0.5', str(max(list_color).round(2))])

'if scaled'
#cbar=fig.colorbar(cax,ticks=[0,0.5, 1])
#cbar.ax.set_yticklabels([str(0), '0.5', str(1)])

'to annotate with splitting ratios'
for i in range(n):

    if difference.iloc[i]<=1e-2 :
        label = str(df_st.alpha_max.iloc[i].round(2))
        # label=' '
    else:
        label = str(df_st.alpha_min.iloc[i].round(2)) + '-'+ str(df_st.alpha_max.iloc[i].round(2))
        # label = ' '

    #if (df_st['A'].iloc[i]==df_st['B'].iloc[i]):
    plt.annotate(label, # this is the text
                 (df_st['A'].iloc[i],df_st['B'].iloc[i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 size=8,
                 xytext=(0,7), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center


#plt.legend(bbox_to_anchor=(1.02, 1.0), loc='upper left')
plt.tight_layout()
plt.savefig(output_file+'/LABELED_split_ratio_isolines_plot_not_annotated_P{}_diverging_scaled_{}.svg'.format(P_concentration,scaled))
# plt.savefig(output_file+'/split_ratio_isolines_plot_not_annotated_P{}_diverging_scaled_{}.svg'.format(P_concentration,scaled))
plt.close()


# #A= 2.857571 b= 2.857571 P=1.0 q=2.0
# #plt.subplots(figsize=(8, 6))
# fig, ax = plt.subplots(figsize=(8, 6))
# a=df_st.A/df_st.P
# plt.hlines(1.0,a.min(),a.max()+5,linestyles='--',alpha=0.3)
# plt.vlines(1.0,a.min(),a.max()+5,linestyles='--',alpha=0.3)
# #plt.figure()
# plt.scatter(df_st.A/df_st.P,df_st.B/df_st.P,c=list_color,cmap=cm.coolwarm,alpha=0.8,s=30+(df_st.P)*15,edgecolors=list_edgecol,linestyle=list_linest,facecolors=list_facecol,linewidth=list_linewidth)
# plt.xlabel(r'$\frac{[A]}{[P]}$')
# plt.ylabel(r'$\frac{[B]}{[P]}$')
# plt.yscale('log')
# plt.xscale('log')
# #plt.xlim([1e-1,70])
# #fig.colorbar()
# cbar=fig.colorbar(cax,ticks=[0,0.5,1.0])
# plt.tight_layout()
# plt.savefig(output_file+'/TRY_pb_pa_conc_limited_less_q{}_lim{}_P.svg'.format(q_equilibrium,limit_cut ))
# plt.close()


df_A_eq_B=pd.read_csv('/open/data/BiUni_general/data_figure3_var_analysis_A_3.0004_B_3.0004_P_1.0_keq_2.csv',index_col=0)
calculate_rate_constants_random_biuni(df_A_eq_B)
df_A_eq_B['gamma_ov']=df_A_eq_B.gamma_1*df_A_eq_B.gamma_2*df_A_eq_B.gamma_3*df_A_eq_B.gamma_4
#
# df_A_eq_B['tf1']=np.log(df_A_eq_B.gamma_1)/np.log(df_A_eq_B.gamma_ov)
# df_A_eq_B['tf2']=np.log(df_A_eq_B.gamma_2)/np.log(df_A_eq_B.gamma_ov)
# df_A_eq_B['tf3']=np.log(df_A_eq_B.gamma_3)/np.log(df_A_eq_B.gamma_ov)
# df_A_eq_B['tf4']=np.log(df_A_eq_B.gamma_4)/np.log(df_A_eq_B.gamma_ov)
# df_A_eq_B['tf5']=np.log(df_A_eq_B.gamma_5)/np.log(df_A_eq_B.gamma_ov)
# df_A_eq_B['tf6']=np.log(df_A_eq_B.gamma_6)/np.log(df_A_eq_B.gamma_ov)
#
# df_1_comb=df_A_eq_B
# df_1_comb['alpha']=df_1_comb.split_ratio
# df_1_comb['v1f']=1/(1-df_1_comb['gamma_1'])*df_1_comb.alpha
# df_1_comb['v2f']=1/(1-df_1_comb['gamma_2'])*df_1_comb.alpha
# df_1_comb['v3f']=1/(1-df_1_comb['gamma_3'])#*df_q2.v
# df_1_comb['v4f']=1/(1-df_1_comb['gamma_4'])#*df_q2.v
# df_1_comb['v5f']=1/(1-df_1_comb['gamma_5'])*(1-df_1_comb.alpha)
# df_1_comb['v6f']=1/(1-df_1_comb['gamma_6'])*(1-df_1_comb.alpha)
#
#
# df_1_comb['v1b']=df_1_comb['gamma_1']/(1-df_1_comb['gamma_1'])*df_1_comb.alpha
# df_1_comb['v2b']=df_1_comb['gamma_2']/(1-df_1_comb['gamma_2'])*df_1_comb.alpha
# df_1_comb['v3b']=df_1_comb['gamma_3']/(1-df_1_comb['gamma_3'])#*df_q2.v
# df_1_comb['v4b']=df_1_comb['gamma_4']/(1-df_1_comb['gamma_4'])#*df_q2.v
# df_1_comb['v5b']=df_1_comb['gamma_5']/(1-df_1_comb['gamma_5'])*(1-df_1_comb.alpha)
# df_1_comb['v6b']=df_1_comb['gamma_6']/(1-df_1_comb['gamma_6'])*(1-df_1_comb.alpha)






level_number_v=10
df=df_st
q_equilibrium=2.0
P_concentration=df_st.P.iloc[0]
base_v=np.linspace(0.005,0.3,15)
fig, ax = plt.subplots(figsize=(5.5, 4))
main_contour_plot_fluxes_random(df,base_v,q_equilibrium,label='v_net',level_num=level_number_v,sampled_points=True,colormap=plt.cm.Reds)
#plt.xlim([0,10])
#plt.ylim([0,10])
plt.xlim([-1e-1,5.3])
plt.ylim([-1e-1,5.3])
plt.xticks(np.linspace(0,5,6))
plt.yticks(np.linspace(0,5,6))
plt.savefig(output_file+'/main_contour_plot_flux_q{}_random_P_{}.svg'.format(q_equilibrium,P_concentration))
plt.close()

#def main_contour_plot_fluxes_random_double(df,Z2,Z3,qeq,P_conc=1.0,label='v',label_z3='v',level_num=7,level_num_Z3=7,sampled_points=False,colormap=plt.cm.coolwarm):

#'to plot saturation graphs'
df_1_comb=df_1_comb
level_number=15
base=np.linspace(0.0,1,15)
base=np.linspace(0.70,1,3)
P_concentration=df_1_comb.P.iloc[0]
df_1_comb['sat']=1-df_1_comb['E']
q_equilibrium=2.0
main_contour_plot_fluxes_random_double(df_1_comb,base,base_v,q_equilibrium,P_conc=P_concentration,label='sat',label_z3='v',level_num=level_number,level_num_Z3=level_number_v,sampled_points=False,colormap=plt.cm.Reds)
#main_contour_plot_fluxes_random(df_1_comb,base,q_equilibrium,P_conc=P_concentration,label='sat',level_num=level_number,sampled_points=True,colormap=plt.cm.Reds)
# plt.xlim(-0.1,5.1)
# plt.ylim(-0.1,5.1)
plt.tight_layout()
plt.savefig(output_file+'/main_contour_plot_sat_q{}_random_P_{}.svg'.format(q_equilibrium,P_concentration),format='svg')
plt.close()





# "illustration figure in Figure1 piecewise constant approximation"
# fig, ax = plt.subplots(figsize=(6, 5))
# gamma_ov=0.2
# bin_size=0.1
# lb=np.array([0.2,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95])
# ub=np.array([0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.0])
# y_values=np.linspace(0.2,1,9)
# y_values_2=np.linspace(0.15,1.3,9)
# plt.hlines(y_values,lb,ub,linestyles='-',alpha=1.0,label=r'$\hat{\gamma}_{i}$',linewidth=3.5)
# plt.plot(y_values_2,y_values_2,c='grey',linestyle='dashed',alpha=0.3,label='y=x',linewidth=3.5)
# #plt.vlines(1.0,a.min(),a.max()+5,linestyles='--',alpha=0.3)
# #plt.figure()
# #plt.scatter(df_st.A/df_st.P,df_st.B/df_st.P,c=list_color,cmap=cm.seismic,alpha=0.8,s=30+(df_st.P)*15,edgecolors=list_edgecol,linestyle=list_linest,facecolors=list_facecol,linewidth=list_linewidth)
# plt.xlabel(r'$\gamma_{i}$',labelpad=20)
# plt.ylabel(r'$\hat{\gamma}_{i}$',rotation=0,labelpad=20)
# ax.set_xticks(np.linspace(0.2,1,9))
# ax.set_yticks(np.linspace(0.2,1,9))
# ax.tick_params(which='both',labelbottom=False,labelleft=False, length=4, width=1.5,)
# plt.xlim([0.15,1.02])
# plt.ylim([0.15,1.02])
# #fig.colorbar()
# #cbar=fig.colorbar(cax,ticks=[0,0.5,1.0])
# plt.tight_layout()
# plt.legend()
# from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
# ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# #ax.xaxis.set_major_locator(MultipleLocator(1))
# #ax.yaxis.set_major_locator(MultipleLocator(1))
# ax.grid(which='major',alpha=0.1,c='gray')
# plt.savefig(output_file+'/piecewise_linear_figre.svg'.format(q_equilibrium))
# plt.close()