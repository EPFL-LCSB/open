import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KNeighborsClassifier



def assign_k_ranges(df,lim=0.95,rounded=9):
    bins = pd.IntervalIndex.from_tuples([(0, lim),(lim, 1.0)])
    for col in df.columns:
        if ('k' in col) and ('max' not in col) and ('range' not in col):
            df[col]=df[col].round(rounded)
            df[col+'_range']=pd.cut(df[col], bins)

    return df

"Main contour plot"
def main_contour_plot_fluxes_3step(df,Z2,qeq,label='v',level_num=7,sampled_points=False):
    #2D
    fig=plt.figure()
    ax=plt.gca()
    fig, ax = plt.subplots(figsize=(8, 6))
    X=df.S
    Y=df.P
    xmax=X.max().round()
    ymax=Y.max().round()
    Z=df[label]
    #Different ways for levels
    level_number=np.linspace(Z2.min(),Z2.max(),level_num)
    cs=ax.tricontourf(X,Y,Z,level_number, cmap=plt.cm.coolwarm) #
    c=ax.tricontour(X,Y,Z,level_number,linewidths=0.2, colors='black')#  levels=level_number,
    ax.plot(np.linspace(0.001,xmax,100),qeq*np.linspace(0.001,xmax,100),ls='--',color='black', linewidth=2.5,alpha=0.5)#crimson #or black
    #plot also your data points
    if sampled_points==True:
        plt.scatter(X,Y,s=12,alpha=0.05,edgecolor='black',facecolor='black') #'midnightblue' #s=15c='black'
        #
        # for i in range(df.overall_gamma.nunique()):
        #     a=np.linspace(0.1,xmax,100)
        #     ax.plot(a,qeq*a*df.overall_gamma.unique()[i],ls='--',color='black', linewidth=1.0,alpha=0.1)#crimson #or black


    xmax=X.max().round()
    ymax=Y.max().round()
    plt.xlim([0,xmax])
    plt.ylim([0,ymax])
    ax.set_xlabel(r'$\mathrm{\widetilde{S}}$')
    ax.set_ylabel(r'$\mathrm{\widetilde{P}}$')
   # ax.set_ylabel(r'$P$',rotation=0)
    #ax.yaxis.set_label_coords(-0.1,0.95)
    ax.xaxis.set_label_coords(0.5,-0.1)
    plt.xticks(np.linspace(0, int(xmax), int(xmax)+1))
    plt.yticks(np.linspace(0, int(ymax), int(ymax)+1))
    plt.tight_layout()
    cbar=fig.colorbar(cs,ax=ax,format='%.2f')
    #cbar.set_ticks(np.linspace(0,1,level_num))
    #new to add labels
    #for saturation plots uncomment this
    #cbar.set_ticks(np.linspace(0,1,level_num))
    #cbar.set_ticklabels(level_number.round(2))
    #cbar.set_ticklabels(np.linspace(-0.35, 0.35, 15))
    #cbar.set_label(r'$\mathrm{\widetilde{v}_{net}}$', rotation=0, labelpad=-40, y=1.12)
    #cbar.set_label(r'$\sigma$', rotation=0, labelpad=-40, y=1.12)



"Main contour plot for double contour"
def main_contour_plot_fluxes_double_3step(df,Z2,Z3,qeq,label='sat',label3='v',level_num=7,level_num_Z3=7,sampled_points=False):
    #2D
    fig=plt.figure()
    ax=plt.gca()
    fig, ax = plt.subplots(figsize=(8, 6))
    X=df.S
    Y=df.P
    xmax=X.max().round()
    ymax=Y.max().round()
    Z=df[label]
    Z3_=df[label3]
    #Different ways for levels
    level_number=np.linspace(Z2.min(),Z2.max(),level_num)
    #flux
    level_number_Z3 = np.linspace(Z3.min(), Z3.max(), level_num_Z3)
    cs=ax.tricontourf(X,Y,Z,level_number, cmap=plt.cm.coolwarm) #
    cs_=ax.tricontour(X,Y,Z,level_number,linewidths=0.2, colors='black')#  levels=level_number,
    loc = matplotlib.ticker.MaxNLocator(7)
    lvls = loc.tick_values(Z3.min(), Z3.max())
    c=ax.tricontour(X,Y,Z3_,level_number_Z3,linewidths=1.0, colors='black',alpha=0.5,linestyles=np.where(lvls >= 0, "-", "-"))#  levels=level_number,
    ax.plot(np.linspace(0.001,xmax,100),qeq*np.linspace(0.001,xmax,100),ls='--',color='black', linewidth=2.5,alpha=0.5)#crimson #or black
    #plot also your data points
    if sampled_points==True:
        plt.scatter(X,Y,s=12,alpha=0.05,edgecolor='black',facecolor='black') #'midnightblue' #s=15c='black'

        # for i in range(df.overall_gamma.nunique()):
        #     a=np.linspace(0.1,xmax,100)
        #     ax.plot(a,qeq*a*df.overall_gamma.unique()[i],ls='--',color='black', linewidth=1.0,alpha=0.1)#crimson #or black


    xmax=X.max().round()
    ymax=Y.max().round()
    plt.xlim([0,xmax])
    plt.ylim([0,ymax])
    ax.set_xlabel(r'$\mathrm{\widetilde{S}}$')
    ax.set_ylabel(r'$\mathrm{\widetilde{P}}$')
   # ax.set_ylabel(r'$P$',rotation=0)
    #ax.yaxis.set_label_coords(-0.1,0.95)
    ax.xaxis.set_label_coords(0.5,-0.1)
    plt.xticks(np.linspace(0, int(xmax), int(xmax)+1))
    plt.yticks(np.linspace(0, int(ymax), int(ymax)+1))
    plt.tight_layout()
    cbar=fig.colorbar(cs,ax=ax,format='%.2f')
    #cbar2=fig.colorbar(c,ax=ax,format='%.2f')
    ax.clabel(c, inline=1, fontsize=18)
    #new to add labels
    #for saturation plots uncomment this
    #cbar.set_ticks(np.linspace(0,1,level_num))
    #cbar.set_ticklabels(level_number.round(2))
    #cbar.set_label(r'$\mathrm{\widetilde{v}_{net}}$', rotation=0, labelpad=-40, y=1.12)
    #cbar.set_label(r'$\sigma$', rotation=0, labelpad=-40, y=1.12)



def main_contour_plot_kcat_3step(df,Z2,qeq,label='v',level_num=7,sampled_points=False):
    #2D
    fig=plt.figure()
    ax=plt.gca()
    X=df.S
    Y=df.P
    xmax=X.max().round()
    ymax=Y.max().round()
    Z=df[label]
    #Different ways for levels
    level_number=Z2
    cs=ax.tricontourf(X,Y,Z,level_number, cmap=plt.cm.GnBu,alpha=1.0) #
    c=ax.tricontour(X,Y,Z,level_number,linewidths=0.2, colors='black')#  levels=level_number,
    ax.plot(np.linspace(0.1,xmax,100),qeq*np.linspace(0.1,xmax,100),ls='--',color='black', linewidth=2.5,alpha=0.5)#crimson #or black
    #plot also your data points
    if sampled_points==True:
        plt.scatter(X,Y,s=15,c='midnightblue',alpha=0.03)

        for i in range(df.overall_gamma.nunique()):
            a=np.linspace(0.1,xmax,100)
            ax.plot(a,qeq*a*df.overall_gamma.unique()[i],ls='--',color='black', linewidth=1.0,alpha=0.03)#crimson #or black


    xmax=X.max().round()
    ymax=Y.max().round()
    plt.xlim([0,xmax])
    plt.ylim([0,ymax])
    ax.set_xlabel(r'$\mathrm{\widetilde{S}}$')
    ax.set_ylabel(r'$\mathrm{\widetilde{P}}$')
   # ax.set_ylabel(r'$P$',rotation=0)
    #ax.yaxis.set_label_coords(-0.1,0.95)
    ax.xaxis.set_label_coords(0.5,-0.1)
    plt.xticks(np.linspace(0, int(xmax), int(xmax)+1))
    plt.yticks(np.linspace(0, int(ymax), int(ymax)+1))
    cbar=fig.colorbar(cs,ax=ax,format='%.2f')
    #new to add labels
    #cbar.set_ticks(np.linspace(0,1,level_num))
    #cbar.set_ticklabels(level_number.round(2))
    #cbar.set_label(r'$\mathrm{\widetilde{v_{net}}}$', rotation=0, labelpad=-40, y=1.12)
    cbar.set_label(r'$\widetilde{k}_{cat,f}$', rotation=0, labelpad=-40, y=1.12)



def add_label(df,group=[''],name='label'):
    df[name] = df.groupby(by=group).grouper.group_info[0] + 1

import math
"Function to train for labeling k ranges"
def return_knn_trained_model_k_ranges_3step(df,qeq,k=6,group=['']):
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

    xx, yy = np.meshgrid(np.linspace((0.001), (xmax), 100),
                         np.linspace((0.001), (ymax), 100))

    Z_test = knn_model.predict(np.c_[xx.ravel(), yy.ravel()])
    X_test = xx.ravel()
    Y_test = yy.ravel()
    df_trained = pd.DataFrame(np.c_[X_test, Y_test, Z_test])
    df_trained.columns = ['S', 'P', 'label']
    df_trained['SP'] = df_trained['P'] / df_trained['S']
    df_trained = df_trained[df_trained['SP'] <= qeq]

    return df_trained


import matplotlib
"Function to plot k ranges as in Heinrich Klipp et al."
def plot_k_ranges_with_knn_3step(df,qeq,mapdict=dict(),map=False):
        X=df.S
        Y=df.P
        xmax=X.max().round()
        ymax=Y.max().round()
        Z = df.label
        fig = plt.figure()
        ax = plt.gca()
        fig, ax = plt.subplots(figsize=(8, 6))

        level_number = np.linspace(1, Z.max() + 1, Z.nunique() + 1)
        "this is to put your own defined colours"
        colors =  ['mediumblue','blue','moccasin']
        cmap = matplotlib.colors.ListedColormap(colors, "", len(colors))
        "this one is when you want to put specific colors"

        #cs = ax.tricontourf(X, Y, Z, level_number, cmap=cmap)#cmap=plt.cm.jet)  #  ['mediumblue','blue','brown']
        cs = ax.tricontourf(X, Y, Z, level_number, cmap=plt.cm.ocean,alpha=0.85)  #  ['mediumblue','blue','brown']

        c = ax.tricontour(X, Y, Z, level_number, linewidths=0.5, colors='black')
        #equilibrium line
        ax.plot(np.linspace(0.1,xmax,100),qeq*np.linspace(0.1,xmax,100),ls='--',color='black', linewidth=3)#crimson #or black
        ax.set_xlabel(r'$\mathrm{\widetilde{S}}$')
        ax.set_ylabel(r'$\mathrm{\widetilde{P}}$')
        #ax.yaxis.set_label_coords(-0.2, 0.95)
        ax.xaxis.set_label_coords(0.5, -0.1)

        plt.xlim([0, xmax])
        plt.ylim([0, ymax])
        plt.xticks(np.linspace(0, int(xmax), int(xmax)+1))
        plt.yticks(np.linspace(0, int(ymax), int(ymax)+1))
        plt.tight_layout()
        cbar = fig.colorbar(cs, ax=ax,format='%d.')
        #cbar = fig.colorbar(c, ax=ax,format='%d.')

        lab = np.arange(1, Z.nunique() + 1, 1)
        loc = lab + .5
        #real_label
        cbar.set_ticks(loc)

        if map == True:
            true_label = [mapdict[value] for value in lab]
            cbar.set_ticklabels(true_label)

        else:
            cbar.set_ticklabels(lab.round(0))
        #cbar.set_label(r'$labels$', rotation=0, labelpad=-40, y=1.12)




import matplotlib
"Function to plot k ranges as in Heinrich Klipp et al."
def plot_k_ranges_with_knn_3step_2(df,qeq,mapdict=dict(),map=False):
    X=df.S
    Y=df.P
    xmax=X.max().round()
    ymax=Y.max().round()
    Z = df.label
    fig = plt.figure()
    ax = plt.gca()
    level_number = np.linspace(1, Z.max() + 1, Z.nunique() + 1)
    "this is to put your own defined colours"
    colors =  ['mediumblue','blue','moccasin']
    cmap = matplotlib.colors.ListedColormap(colors, "", len(colors))
    "this one is when you want to put specific colors"
    #cs = ax.tricontourf(X, Y, Z, level_number, cmap=cmap)#cmap=plt.cm.jet)  #  ['mediumblue','blue','brown']
    #cs = ax.tricontourf(X, Y, Z, level_number, cmap=plt.cm.jet)  #  ['mediumblue','blue','brown']

    c = ax.tricontour(X, Y, Z, level_number, linewidths=1.0, colors='black')
    #equilibrium line
    ax.plot(np.linspace(0.1,xmax,100),qeq*np.linspace(0.1,xmax,100),ls='--',color='crimson', linewidth=4)#crimson #or black
    ax.set_xlabel(r'$\mathrm{\widetilde{S}}$')
    ax.set_ylabel(r'$\mathrm{\widetilde{P}}$')
    #ax.yaxis.set_label_coords(-0.2, 0.95)
    ax.xaxis.set_label_coords(0.5, -0.1)

    plt.xlim([0, xmax])
    plt.ylim([0, ymax])
    plt.xticks(np.linspace(0, int(xmax), int(xmax)+1))
    plt.yticks(np.linspace(0, int(ymax), int(ymax)+1))

    #cbar = fig.colorbar(cs, ax=ax,format='%d.')
    cbar = fig.colorbar(c, ax=ax,format='%d.')

    lab = np.arange(1, Z.nunique() + 1, 1)
    loc = lab + .5
    #real_label
    cbar.set_ticks(loc)

    if map == True:
        true_label = [mapdict[value] for value in lab]
        cbar.set_ticklabels(true_label)

    else:
        cbar.set_ticklabels(lab.round(0))
    #cbar.set_label(r'$labels$', rotation=0, labelpad=-40, y=1.12)




def plot_colored_data(df,label,qeq,mapdict=dict(),map=False,cmap='jet'):
    tick = np.linspace(1, df[label].nunique(), df[label].nunique())
    fig = plt.figure()
    ax = plt.gca()
    lab = df[label]
    cmap = plt.cm.jet
    norm = matplotlib.colors.BoundaryNorm(np.arange(0.5, df[label].max() + 1, 1), cmap.N)
    X=df.S
    Y=df.P
    im = ax.scatter(X, Y, c=lab, cmap=cmap, norm=norm)
    ax.set_ylim([-0.05, 5.1])
    ax.set_xlim([-0.05, 5.1])
    ax.plot(np.linspace(0.1, df.S.max(), 100), qeq*np.linspace(0.1, df.S.max(), 100), ls='--', color='crimson',
            linewidth=2)
    cbar = fig.colorbar(im, ax=ax, ticks=tick)
    cbar.set_ticklabels(tick)
    plt.xlabel(r'$\mathrm{\widetilde{S}}$')
    plt.ylabel(r'$\mathrm{\widetilde{P}}$')
    plt.xticks(np.linspace(0, int(X.max().round()), int(X.max().round()) + 1))
    plt.yticks(np.linspace(0, int(Y.max().round()), int(Y.max().round()) + 1))
    if map == True:
        true_label = [mapdict[value] for value in tick]
        cbar.set_ticklabels(true_label)
    plt.xlim([0, X.max()*1.01])
    plt.ylim([0, Y.max()*1.01])
    return fig,ax


def plot_colored_data_KMS_KP(df,label,qeq,mapdict=dict(),map=False,cmap='jet'):
    tick = np.linspace(1, df[label].nunique(), df[label].nunique())
    fig = plt.figure()
    ax = plt.gca()
    fig, ax = plt.subplots(figsize=(8, 6))
    lab = df[label]
    import matplotlib
    cmap = plt.cm.ocean
    norm = matplotlib.colors.BoundaryNorm(np.arange(0.5, df[label].max() + 1, 1), cmap.N)
    X=df.KM_S
    Y=df.KM_P
    im = ax.scatter(X, Y, c=lab,edgecolor='black',linewidth=0.1, cmap=cmap, norm=norm,alpha=0.85)
    ax.set_ylim([-0.05, 5.1])
    ax.set_xlim([-0.05, 5.1])
    # ax.plot(np.linspace(0.1, df.S.max(), 100), qeq*np.linspace(0.1, df.S.max(), 100), ls='--', color='crimson',
    #         linewidth=2)
    cbar = fig.colorbar(im, ax=ax, ticks=tick)
    cbar.set_ticklabels(tick)

    plt.xlabel(r'$\mathrm{\widetilde{K}_{M,S}}$')
    plt.ylabel(r'$\mathrm{\widetilde{K}_{M,P}}$')
    # plt.xticks(np.linspace(0, int(X.max().round()), int(X.max().round()) + 1))
    # plt.yticks(np.linspace(0, int(Y.max().round()), int(Y.max().round()) + 1))
    plt.xticks(np.linspace(0, 2, 3))
    plt.yticks(np.linspace(0, 2, 3))
    if map == True:
        true_label = [mapdict[value] for value in tick]
        cbar.set_ticklabels(true_label)
    plt.xlim([0, X.max()*1.01])
    plt.ylim([0, Y.max()*1.01])
    plt.tight_layout()

    return fig,ax
#
#
def plot_colored_data_2(df,label,qeq,mapdict=dict(),map=False,cmap='jet'):
    tick = np.linspace(1, df[label].nunique(), df[label].nunique())
    fig = plt.figure()
    ax = plt.gca()
    lab = df[label]
    true_label=df['true_label']
    import matplotlib
    colors =  ['gold','green','red','darkred','yellow','aquamarine','orange','darkblue','navy', 'coral','dodgerblue']
    #cmap = matplotlib.colors.ListedColormap(colors, "", len(colors))
    #norm = matplotlib.colors.BoundaryNorm(np.arange(0.5, 10 + 1, 1), cmap.N)
    X=df.S
    Y=df.P
    im = ax.scatter(X, Y, c= [colors[int(i)] for i in true_label])

    ax.set_ylim([-0.05, 5.1])
    ax.set_xlim([-0.05, 5.1])
    ax.plot(np.linspace(0.1, df.S.max(), 100), qeq*np.linspace(0.1, df.S.max(), 100), ls='--', color='crimson',
            linewidth=2)
    #cbar = fig.colorbar(im, ax=ax, ticks=tick)
    #cbar.set_ticklabels(tick)
    plt.xlabel(r'$\mathrm{\widetilde{S}}$')
    plt.ylabel(r'$\mathrm{\widetilde{P}}$')
    plt.xticks(np.linspace(0, int(X.max().round()), int(X.max().round()) + 1))
    plt.yticks(np.linspace(0, int(Y.max().round()), int(Y.max().round()) + 1))
    if map == True:
        true_label = [mapdict[value] for value in tick]
     #   cbar.set_ticklabels(true_label)
    plt.xlim([0, X.max()*1.01])
    plt.ylim([0, Y.max()*1.01])
    return fig,ax


import matplotlib
"Function to plot k ranges as in Heinrich Klipp et al."
def plot_k_ranges_with_knn_3step_2(df,qeq,mapdict=dict(),map=False):
    X=df.S
    Y=df.P
    xmax=X.max().round()
    ymax=Y.max().round()
    Z = df.label
    fig = plt.figure()
    ax = plt.gca()
    level_number = np.linspace(1, Z.max() + 1, Z.nunique() + 1)
    "this is to put your own defined colours"
    colors =  ['blue','green','orange','red','brown','darkviolet','steelblue','darkturquoise','navy', 'coral','medium blue']
    cmap = matplotlib.colors.ListedColormap(colors, "", len(colors))
    "this one is when you want to put specific colors"
    #cs = ax.tricontourf(X, Y, Z, level_number, cmap=cmap)#cmap=plt.cm.jet)  #  ['mediumblue','blue','brown']
    #cs = ax.tricontourf(X, Y, Z, level_number, c=colors)  #  ['mediumblue','blue','brown']

    c = ax.tricontour(X, Y, Z, level_number, linewidths=1.0, colors='black')
    #equilibrium line
    ax.plot(np.linspace(0.1,xmax,100),qeq*np.linspace(0.1,xmax,100),ls='--',color='crimson', linewidth=4)#crimson #or black
    ax.set_xlabel(r'$\mathrm{\widetilde{S}}$')
    ax.set_ylabel(r'$\mathrm{\widetilde{P}}$')
    #ax.yaxis.set_label_coords(-0.2, 0.95)
    ax.xaxis.set_label_coords(0.5, -0.1)

    plt.xlim([0, xmax])
    plt.ylim([0, ymax])
    plt.xticks(np.linspace(0, int(xmax), int(xmax)+1))
    plt.yticks(np.linspace(0, int(ymax), int(ymax)+1))

    #cbar = fig.colorbar(cs, ax=ax,format='%d.')
    cbar = fig.colorbar(c, ax=ax,format='%d.')

    lab = np.arange(1, Z.nunique() + 1, 1)
    loc = lab + .5
    #real_label
    cbar.set_ticks(loc)

    if map == True:
        true_label = [mapdict[value] for value in lab]
        cbar.set_ticklabels(true_label)

    else:
        cbar.set_ticklabels(lab.round(0))
    #cbar.set_label(r'$labels$', rotation=0, labelpad=-40, y=1.12)

def plot_colored_data_variability(df,label,qeq,cmap='jet'):
    tick = np.linspace(1, df[label].nunique(), df[label].nunique())
    fig = plt.figure()
    ax = plt.gca()
    lab = df[label]
    cmap = plt.cm.jet
    import matplotlib
   # norm = matplotlib.colors.BoundaryNorm(np.arange(0.5, df[label].max() + 1, 1), cmap.N)
    X=df.S
    Y=df.P
    im = ax.scatter(X, Y, c=lab, cmap=cmap)#, norm=norm)
    ax.plot(np.linspace(0.1, df.S.max(), 100), qeq*np.linspace(0.1, df.S.max(), 100), ls='--', color='crimson',
            linewidth=2)
    cbar = fig.colorbar(im, ax=ax, ticks=tick)
    cbar.set_ticklabels(tick)
    plt.xlabel(r'$S$')
    plt.ylabel(r'$P$')
    plt.xticks(np.linspace(0, int(X.max().round()), int(X.max().round()) + 1))
    plt.yticks(np.linspace(0, int(Y.max().round()), int(Y.max().round()) + 1))
    plt.xlim([0, X.max()*1.01])
    plt.ylim([0, Y.max()*1.01])

from mpl_toolkits.mplot3d import Axes3D
def plot_colored_data_gammas(df,label,qeq,mapdict=dict(),map=False,cmap='jet'):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    df=df[df.overall_gamma==0.99]
    X=df.gamma_1
    Y=df.gamma_2
    Z=df.gamma_3
    cmap = plt.cm.jet
    lab = df[label]
    norm = matplotlib.colors.BoundaryNorm(np.arange(0.5, df[label].max() + 1, 1), cmap.N)
    im = ax.scatter(X, Y,Z, c=lab, cmap=cmap, norm=norm)
    ax.set_xlabel('gamma_1')
    ax.set_ylabel('gamma_2')
    ax.set_zlabel('gamma_3')

    ax.axes.set_xlim3d(left=0.99, right=1)
    ax.axes.set_ylim3d(bottom=0.99, top=1)
    ax.axes.set_zlim3d(bottom=0.99, top=1)

    ax.view_init(30, 180)

    for i in (df.overall_gamma.unique()):
        ln_gammas = np.random.dirichlet((1, 1, 1), 5000)
        gammas = np.exp((ln_gammas * np.log(i)))

        ax.scatter(gammas[:,0],gammas[:,1],gammas[:,2], c='grey', marker='x',
            alpha=0.1)  # crimson #or black

def plot_colored_data_gammas_2d(df,label,qeq,mapdict=dict(),map=False,cmap='jet'):
    tick = np.linspace(1, df[label].nunique(), df[label].nunique())
    fig = plt.figure()
    ax = plt.gca()
    lab = df[label]
    cmap = plt.cm.jet
    import matplotlib
    norm = matplotlib.colors.BoundaryNorm(np.arange(0.5, df[label].max() + 1, 1), cmap.N)
    Y=df.gamma_2
    X=1/df.gamma_1

    #X=1/(df.gamma_2**2*df.gamma_3)
    #X=1*2*df.overall_gamma/df.gamma_3


    #cbar = fig.colorbar(im, ax=ax, ticks=tick)
   # cbar.set_ticklabels(tick)
    #plt.xlabel(r'$\mathrm{\gamma_{1}\gamma_{2}}$')
    #plt.ylabel(r'$\mathrm{\gamma_{3}}$')
    plt.xticks(np.linspace(0, int(X.max().round()), int(X.max().round()) + 1))
    plt.yticks(np.linspace(0, int(Y.max().round()), int(Y.max().round()) + 1))

    #
    # for i in (df.overall_gamma.unique()):
    #     ln_gammas = np.random.dirichlet((1, 1, 1), 3000)
    #     gammas = np.exp((ln_gammas * np.log(i)))
    #
    #     #plt.plot(gammas[:,0]*gammas[:,1],gammas[:,2], color='crimson',
    #      #   alpha=0.1)  # crimson #or black
    #
    #     plt.scatter(gammas[:,0]*gammas[:,1],gammas[:,2],marker='o', color='grey',s=3,
    #         alpha=0.1)  # crimson #or black

    im = ax.scatter(X, Y, c=lab, cmap=cmap, norm=norm)

    cbar = fig.colorbar(im, ax=ax, ticks=tick)
    cbar.set_ticklabels(tick)

    if map == True:
        true_label = [mapdict[value] for value in tick]
        cbar.set_ticklabels(true_label)
    plt.xlim([0, X.max()*1.01])
    plt.ylim([0, Y.max()*1.01])

#
#
def plot_colored_data_grouped(df,label,qeq,mapdict=dict(),map=False,cmap='jet'):
    tick = np.linspace(1, df[label].nunique(), df[label].nunique())
    fig = plt.figure()
    ax = plt.gca()
    lab = df[label]
    true_label=df['true_label']
    import matplotlib
    colors =  ['gold','green','red','darkred','yellow','aquamarine','orange','darkblue','navy', 'coral','dodgerblue']
    #cmap = matplotlib.colors.ListedColormap(colors, "", len(colors))
    #norm = matplotlib.colors.BoundaryNorm(np.arange(0.5, 10 + 1, 1), cmap.N)
    X=df.S
    Y=df.overall_gamma
    im = ax.scatter(X, Y, c= [colors[int(i)] for i in true_label])

    #ax.set_ylim([-0.05, 5.1])
    #ax.set_xlim([-0.05, 5.1])
    #ax.plot(np.linspace(0.1, df.S.max(), 100), qeq*np.linspace(0.1, df.S.max(), 100), ls='--', color='crimson',
     #       linewidth=2)
    #cbar = fig.colorbar(im, ax=ax, ticks=tick)
    #cbar.set_ticklabels(tick)
    plt.ylabel(r'$\mathrm{\Gamma}$')
    plt.xlabel(r'$\mathrm{\widetilde{S}}$')
    plt.xticks(np.linspace(0, int(X.max().round()), int(X.max().round()) + 1))
    plt.yticks(np.linspace(0, int(Y.max().round()), int(Y.max().round()) + 1))
    if map == True:
        true_label = [mapdict[value] for value in tick]
     #   cbar.set_ticklabels(true_label)
    plt.xlim([-0.5, X.max()*1.1])
    plt.ylim([-0.01, Y.max()*1.1])
    return fig,ax

def return_knn_trained_model_k_ranges_3step_grouped(df,qeq,k=6,group=['']):
    X=df.S
    Y=df.overall_gamma
    #xmax=X.max().round()
    #ymax=Y.max().round()
    xmax=math.ceil(X.max())
    ymax=math.ceil(Y.max())
    Z = df.label
    k = k
    X_train = pd.DataFrame([X, Y]).T

    knn_model = KNeighborsClassifier(n_neighbors=k)
    knn_model.fit(X_train, Z)

    xx, yy = np.meshgrid(np.linspace((0.001), (xmax), 100),
                         np.linspace((0.0001), (ymax), 100))

    Z_test = knn_model.predict(np.c_[xx.ravel(), yy.ravel()])
    X_test = xx.ravel()
    Y_test = yy.ravel()
    df_trained = pd.DataFrame(np.c_[X_test, Y_test, Z_test])
    df_trained.columns = ['S', 'overall_gamma', 'label']

    #df_trained['SP'] = df_trained['P'] / df_trained['S']
    #df_trained = df_trained[df_trained['SP'] <= qeq]

    return df_trained