from figures.main_figures import *
import os

"all figures " \
"generated for Michaelis-Menten mechanism for K_eq=2.0"


plt.rcParams.update({'font.size':16,
                     'font.family':'Calibri',
                     'legend.fontsize':16,
                    'xtick.labelsize' : 18,
                    'ytick.labelsize' : 18
                     })
plt.rcParams['svg.fonttype'] = 'none'

map_dict={1:7,
2:8,
3:9,
4:10,
5:5,
6:4,
7:1,
8:6,
9:2,
10:3
}


"Function to train for labeling k ranges"
def return_knn_trained_model_k_ranges_3step_2(df,qeq,k=6,group=['']):
    X=df.S
    Y=df.P
    #xmax=X.max().round()
    #ymax=Y.max().round()
    xmax=math.ceil(X.max())
    ymax=math.ceil(Y.max())
    Z = df.label
    k = k
    X_train = pd.DataFrame([X, Y]).T

    knn_model = KNeighborsClassifier(n_neighbors=k)
    knn_model.fit(X_train, Z)

    xx, yy = np.meshgrid(np.linspace((0.001), (xmax), 500),
                         np.linspace((0.001), (ymax), 500))

    Z_test = knn_model.predict(np.c_[xx.ravel(), yy.ravel()])
    X_test = xx.ravel()
    Y_test = yy.ravel()
    df_trained = pd.DataFrame(np.c_[X_test, Y_test, Z_test])
    df_trained.columns = ['S', 'P', 'label']
    #df_trained['SP'] = df_trained['P'] / df_trained['S']
    #df_trained = df_trained[df_trained['SP'] <= qeq]

    return df_trained

'where to store the figures'

output_file='./milp_3step_q2_all/paper_latest'
#output_file='./milp_3step_011120'
if not os.path.exists(output_file):
    os.makedirs(output_file)
"DATA FROM OPTIMIZATION"
"put the path for the data"


df_all=pd.read_csv('/open/data/UniUni/df_all_uniuni_q_eq_2.0_latest.csv',index_col=0)
df_all=df_all.reset_index()


gr_general_3step=['k1f_range','k2f_range','k3f_range','k1b_range','k2b_range','k3b_range']
q_eq=df_all.k1f*df_all.k2f*df_all.k3f/(df_all.k1b*df_all.k2b*df_all.k3b)
q_equilibrium=q_eq.iloc[0].round(1)
limit_cut=0.99
assign_k_ranges(df_all,lim=(limit_cut))
add_label(df_all,gr_general_3step)

fig,ax=plot_colored_data(df_all, label='label', qeq=2,mapdict=map_dict,map=False)
plt.xlim([-0.05,5.1])
plt.ylim([-0.05, 5.1])
plt.savefig(output_file+'/all_k_ranges_datapoints_q{}_limit_{}_231221.svg'.format(q_equilibrium.round(),limit_cut))
plt.close()

df_trained=pd.DataFrame()
neighbors=13
df_trained=return_knn_trained_model_k_ranges_3step_2(df_all,q_equilibrium,k=neighbors,group=gr_general_3step)
print(df_trained.label.nunique())
plot_k_ranges_with_knn_3step(df_trained,q_equilibrium,mapdict=map_dict,map=True)
plt.scatter(df_all.S,df_all.P,s=15,alpha=0.05,edgecolor='black',facecolor='black')
plt.xlim([0,5.1])
plt.ylim([0,5.1])
plt.savefig(output_file+'/030122_all_NEW_k_ranges_knn_accent_plot_q{}_limit_{}_251221_k_{}.svg'.format(q_equilibrium.round(),limit_cut,neighbors),format='svg')
plt.close()

df_q2=df_all
df_q2['kcat'] = df_q2.k2f * df_q2.k3f / (
            df_q2.k2f + df_q2.k3f + df_q2.k2b)

df_q2['nu_calculated'] = df_q2.v / df_q2['kcat']

#df_q2['residual'] = (1 - df_q2.overall_gamma) * (1 - df_q2.E)
df_q2['KM_S']=(df_q2.k2f*df_q2.k3f+df_q2.k1b*df_q2.k3f+df_q2.k1b*df_q2.k2b)/ \
(df_q2.k1f*(df_q2.k2f + df_q2.k3f + df_q2.k2b))
df_q2['KM_P']=(df_q2.k2f*df_q2.k3f+df_q2.k1b*df_q2.k3f+df_q2.k1b*df_q2.k2b)/ \
(df_q2.k3b*(df_q2.k2f + df_q2.k1b + df_q2.k2b))
df_q2['kcat_back']= df_q2.k1b*df_q2.k2b/(df_q2.k2f+df_q2.k1b+df_q2.k2b)




import matplotlib.pyplot as plt
import matplotlib.cm as cm


df_q2=df_all
level_number=13
df_q2['overall_gamma']=df_q2['gamma_overall']
rang=np.linspace(-0.35,0.35,13)
main_contour_plot_fluxes_3step(df_q2,rang,q_equilibrium,label='v',level_num=level_number,sampled_points=True)
plt.xlim([-1e-1,5.1])
plt.ylim([-1e-1,5.1])
plt.xlim([0,5.1])
plt.ylim([0,5.1])
plt.tight_layout()
#cbar.set_label(r'$\sigma$', rotation=0, labelpad=-40, y=1.12)
plt.savefig(output_file+'/030122_contour_plot_NEW_fluxes_all_direction_new_cmap_q_{}_more_Dat.svg'.format(q_equilibrium))
plt.close()


#new
df_q2=df_all
level_number=15
df_q2['overall_gamma']=df_q2['gamma_overall']
rang=np.linspace(-0.35,0.35,15)
main_contour_plot_fluxes_3step(df_q2,rang,q_equilibrium,label='v',level_num=level_number,sampled_points=True)
plt.xlim([-1e-1,5.1])
plt.ylim([-1e-1,5.1])
plt.xlim([0,5.1])
plt.ylim([0,5.1])
plt.tight_layout()
#cbar.set_label(r'$\sigma$', rotation=0, labelpad=-40, y=1.12)
plt.savefig(output_file+'/030122_contour_plot_NEW_fluxes_all_direction_new_cmap_q_{}_more_Dat.svg'.format(q_equilibrium),format='svg')
plt.close()


df_q2=df_all
df_q2['sat']=1-df_q2.E
level_number=6
df_q2['overall_gamma']=df_q2['gamma_overall']
rang=np.linspace(-0.35,0.35,15)
main_contour_plot_fluxes_double_3step(df_q2,df_q2['sat'],rang,q_equilibrium,label='sat',label3='v',level_num=10,level_num_Z3=level_number,sampled_points=False)
plt.xlim([-1e-1,5.1])
plt.ylim([-1e-1,5.1])
plt.xlim([0,5.1])
plt.ylim([0,5.1])
plt.tight_layout()
#cbar.set_label(r'$\sigma$', rotation=0, labelpad=-40, y=1.12)
plt.savefig(output_file+'/130322_contour_double_plot_NEW_fluxes_all_direction_new_cmap_q_{}_more_Dat.svg'.format(q_equilibrium),format='svg')
plt.close()


level_number=10
df_q2['sat']=1-df_q2.E
df_q2['overall_gamma']=df_q2['gamma_overall']
main_contour_plot_fluxes_3step(df_q2,df_q2['sat'],q_equilibrium,label='sat',level_num=level_number,sampled_points=True)
plt.xlim([-1e-1,5.1])
plt.ylim([-1e-1,5.1])
plt.tight_layout()
#cbar.set_label(r'$\sigma$', rotation=0, labelpad=-40, y=1.12)
plt.savefig(output_file+'/contour_plot_saturation_all_direction_new_cmap_q_{}_more_Dat.svg'.format(q_equilibrium))
plt.close()

df_q2['v1f']=1/(1-df_q2['gamma_1'])#*df_q2.v
df_q2['v2f']=1/(1-df_q2['gamma_2'])#*df_q2.v
df_q2['v3f']=1/(1-df_q2['gamma_3'])#*df_q2.v


df_q2['v1b']=df_q2['gamma_1']/(1-df_q2['gamma_1'])#*df_q2.v
df_q2['v2b']=df_q2['gamma_2']/(1-df_q2['gamma_2'])#*df_q2.v
df_q2['v3b']=df_q2['gamma_3']/(1-df_q2['gamma_3'])#*df_q2.v


df_q2['overall_gamma']=df_q2['gamma_overall']
Z2=np.array([0.001,0.1,0.28,0.35,0.41])

main_contour_plot_kcat_3step(df_q2,Z2,q_equilibrium,label='kcat',level_num=level_number,sampled_points=True)
plt.xlim([-1e-1,5.1])
plt.ylim([-1e-1,5.1])
plt.tight_layout()
#cbar.set_label(r'$\sigma$', rotation=0, labelpad=-40, y=1.12)
plt.savefig(output_file+'/contour_plot_kcat_fwd_new_all_direction_q_{}.svg'.format(q_equilibrium))
plt.close()


main_contour_plot_kcat_3step(df_q2,Z2,q_equilibrium,label='kcat_back',level_num=level_number,sampled_points=True)
plt.xlim([-1e-1,5.1])
plt.ylim([-1e-1,5.1])
plt.tight_layout()
#cbar.set_label(r'$\sigma$', rotation=0, labelpad=-40, y=1.12)
plt.savefig(output_file+'/contour_plot_kcat_back_new_all_direction_q_{}.svg'.format(q_equilibrium))
plt.close()

fig,ax=plot_colored_data_KMS_KP(df_q2, label='label', qeq=2,mapdict=map_dict,map=True)
plt.xlim([-0.05,2.5])
plt.ylim([-0.05, 2.5])
plt.savefig(output_file+'/all_KMS_KMP_q_{}_231221.svg'.format(q_equilibrium.round()))
plt.close()





