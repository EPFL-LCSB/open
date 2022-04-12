import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KNeighborsClassifier



def assign_k_ranges_random(df,lim=0.95,rounded=9):
    bins = pd.IntervalIndex.from_tuples([(-1e-5, lim),(lim, 1.0)])
    for col in df.columns:
        if ('k' in col) and ('max' not in col) and ('range' not in col) and ('_' not in col):
            df[col]=df[col].round(rounded)
            df[col+'_range']=pd.cut(df[col], bins)

    return df

"Main contour plot"
def main_contour_plot_fluxes_random(df,Z2,qeq,P_conc=1.0,label='v',level_num=7,sampled_points=False,colormap=plt.cm.coolwarm):
    #2D
    # fig=plt.figure()
    # ax=plt.gca()
    fig, ax = plt.subplots(figsize=(8, 6))
    X=df.A
    Y=df.B
    xmax=X.max().round()
    ymax=Y.max().round()
    Z=df[label]
    #Different ways for levels
    level_number=np.linspace(Z2.min(),Z2.max(),level_num)
    cs=ax.tricontourf(X,Y,Z,level_number, cmap=colormap,alpha=0.95) #coolwarm
    c=ax.tricontour(X,Y,Z,level_number,linewidths=0.2, colors='black',alpha=0.8)#  levels=level_number,
    ax.plot(np.linspace(0.1, 5, 100), P_conc/(qeq* np.linspace(0.1, 5, 100)), ls='--', color='black',
            linewidth=2.5,alpha=0.5)
    #plot also your data points
    if sampled_points==True:
        plt.scatter(X,Y,c='black',alpha=0.05)
        gamma_ov_range=[0.02,0.03,0.05,0.12,0.3,0.7]
       # for i in gamma_ov_range: #range(df.gamma_ov.nunique()):
        #    a=np.linspace(0.1,xmax,100)
         #   ax.plot(a,1/(qeq*a*i),ls='--',color='black', linewidth=1.0,alpha=0.3)#crimson #or black


    xmax=X.max().round()
    ymax=Y.max().round()
    plt.xlim([0,xmax])
    plt.ylim([0,ymax])
    ax.set_xlabel(r'$\mathrm{\widetilde{A}}$')
    ax.set_ylabel(r'$\mathrm{\widetilde{B}}$')
    #ax.yaxis.set_label_coords(-0.1,0.95)
    ax.xaxis.set_label_coords(0.5,-0.1)
    plt.xticks(np.linspace(0, int(xmax), int(xmax)+1))
    plt.yticks(np.linspace(0, int(ymax), int(ymax)+1))
    cbar=fig.colorbar(cs,ax=ax,format='%.2f')
    #to be aranged based on desired scaling
    level_num_colorbar=11
    cbar.set_ticks(np.linspace(Z2.min(),Z2.max(),level_num_colorbar))
    #cbar.set_label(r'$\mathrm{\widetilde{v_{net}}}$', rotation=0, labelpad=-40, y=1.12)

"Main contour plot for double contours"
def main_contour_plot_fluxes_random_double(df,Z2,Z3,qeq,P_conc=1.0,label='v',label_z3='v',level_num=7,level_num_Z3=7,sampled_points=False,colormap=plt.cm.coolwarm):
    #2D
    # fig=plt.figure()
    # ax=plt.gca()
    fig, ax = plt.subplots(figsize=(8, 6))
    X=df.A
    Y=df.B
    xmax=X.max().round()
    ymax=Y.max().round()
    Z=df[label]
    Z3_=df[label_z3]
    #Different ways for levels
    level_number=np.linspace(Z2.min(),Z2.max(),level_num)
    #level_number_z3=np.linspace(Z2.min(),Z2.max(),level_num)
    level_number_Z3 = np.linspace(Z3.min(), Z3.max(), level_num_Z3)
    cs=ax.tricontourf(X,Y,Z,level_number, cmap=colormap,alpha=0.95) #coolwarm
    c=ax.tricontour(X,Y,Z3_,level_number_Z3,linewidths=0.5, colors='black',alpha=0.8)#  levels=level_number,

    cs_=ax.tricontour(X,Y,Z,level_number,linewidths=0.02, colors='black',alpha=0.8)#  levels=level_number,
    ax.plot(np.linspace(0.1, 5, 100), P_conc/(qeq* np.linspace(0.1, 5, 100)), ls='--', color='black',
            linewidth=2.5,alpha=0.5)
    #plot also your data points
    if sampled_points==True:
        plt.scatter(X,Y,c='black',alpha=0.05)
        gamma_ov_range=[0.02,0.03,0.05,0.12,0.3,0.7]
       # for i in gamma_ov_range: #range(df.gamma_ov.nunique()):
        #    a=np.linspace(0.1,xmax,100)
         #   ax.plot(a,1/(qeq*a*i),ls='--',color='black', linewidth=1.0,alpha=0.3)#crimson #or black


    xmax=X.max().round()
    ymax=Y.max().round()
    plt.xlim([0,xmax])
    plt.ylim([0,ymax])
    ax.set_xlabel(r'$\mathrm{\widetilde{A}}$')
    ax.set_ylabel(r'$\mathrm{\widetilde{B}}$')
    #ax.yaxis.set_label_coords(-0.1,0.95)
    ax.xaxis.set_label_coords(0.5,-0.1)
    plt.xticks(np.linspace(0, int(xmax), int(xmax)+1))
    plt.yticks(np.linspace(0, int(ymax), int(ymax)+1))
    cbar=fig.colorbar(cs,ax=ax,format='%.2f')
    ax.clabel(c, inline=1, fontsize=12)
    #to be aranged based on desired scaling
    level_num_colorbar=11
    cbar.set_ticks(np.linspace(Z2.min(),Z2.max(),level_num_colorbar))
    #cbar.set_label(r'$\mathrm{\widetilde{v_{net}}}$', rotation=0, labelpad=-40, y=1.12)


def add_label(df,group=['']):
    df['label'] = df.groupby(by=group).grouper.group_info[0] + 1

"Function to train for labeling k ranges"
def return_knn_trained_model_k_ranges_random(df,qeq,k=6,group=['']):
    X=df.A
    Y=df.B
    xmax=X.max().round()
    ymax=Y.max().round()
    Z = df.label
    k = k
    X_train = pd.DataFrame([X, Y]).T

    knn_model = KNeighborsClassifier(n_neighbors=k)
    knn_model.fit(X_train, Z)

    n_sample_size=800
    xx, yy = np.meshgrid(np.linspace((0.001), (xmax), n_sample_size),
                         np.linspace((0.001), (ymax), n_sample_size))

    Z_test = knn_model.predict(np.c_[xx.ravel(), yy.ravel()])
    X_test = xx.ravel()
    Y_test = yy.ravel()
    df_trained = pd.DataFrame(np.c_[X_test, Y_test, Z_test])
    df_trained.columns = ['A', 'B', 'label']
    df_trained['AB'] = df_trained['A'] * df_trained['B']
    df_trained = df_trained[df_trained['AB'] >= 1/qeq]
    #df_trained['SP'] = df_trained['P'] / df_trained['S']
    #df_trained = df_trained[df_trained['SP'] <= qeq]

    return df_trained


import matplotlib
"Function to plot k ranges as in Heinrich Klipp et al."
def plot_k_ranges_with_knn_random(df,qeq,mapdict=dict(),map=False):
    X=df.A
    Y=df.B
    xmax=X.max().round()
    ymax=Y.max().round()
    Z = df.label
    fig = plt.figure()
    ax = plt.gca()
    fig, ax = plt.subplots(figsize=(8, 6))
    level_number = np.linspace(1, Z.max() + 1, Z.nunique() + 1)
    #level_number = np.linspace(1, Z.nunique(), Z.nunique())
    "this is to put your own defined colours"
    colors =  ['mediumblue','blue','moccasin']
    cmap = matplotlib.colors.ListedColormap(colors, "", len(colors))
    "this one is when you want to put specific colors"
    #cs = ax.tricontourf(X, Y, Z, level_number, cmap=cmap)#cmap=plt.cm.jet)  #  ['mediumblue','blue','brown']
    cs = ax.tricontourf(X, Y, Z, level_number, cmap=plt.cm.jet)  #  ['mediumblue','blue','brown']

    c = ax.tricontour(X, Y, Z, level_number, linewidths=0.2, colors='black')
    #equilibrium line
    ax.plot(np.linspace(0.1, 5, 100), 1/(qeq* np.linspace(0.1, 5, 100)), ls='--', color='crimson',
            linewidth=4)
    ax.set_xlabel(r'$\mathrm{\widetilde{A}}$')
    ax.set_ylabel(r'$\mathrm{\widetilde{B}}$')
    #ax.yaxis.set_label_coords(-0.2, 0.95)
    ax.xaxis.set_label_coords(0.5, -0.1)

    plt.xlim([0, xmax])
    plt.ylim([0, ymax])
    plt.xticks(np.linspace(0, int(xmax), int(xmax)+1))
    plt.yticks(np.linspace(0, int(ymax), int(ymax)+1))

    cbar = fig.colorbar(cs, ax=ax,format='%d.')
    #cbar=fig.colorbar(cs,ax=ax,format='%.0f')
    lab = np.arange(1, Z.nunique() + 1, 1)
    lab=np.linspace(1,Z.nunique(),Z.nunique())
    loc = lab + .5
    #real_label
    cbar.set_ticks(loc)


    if map == True:
        true_label = [mapdict[value] for value in lab]
        cbar.set_ticklabels(true_label)

    #else:
        #cbar.set_ticklabels(lab.round(0))
    #cbar.set_label(r'$labels$', rotation=0, labelpad=-40, y=1.12)



def plot_colored_data_random(df,label,qeq,mapdict=dict(),P_conc=1.0,map=False,cmap='jet'):
    max=df[label].nunique()
    #tick = np.linspace(1, df[label].nunique(), df[label].nunique())
    tick = np.linspace(1, max, max)

    fig = plt.figure()
    ax = plt.gca()
    lab = df[label]
    cmap = plt.cm.jet
    import matplotlib
    #norm = matplotlib.colors.BoundaryNorm(np.arange(0.5, df[label].max() + 1, 1), cmap.N)
    norm = matplotlib.colors.BoundaryNorm(np.arange(0.5, max + 1, 1), cmap.N)

    X=df.A
    Y=df.B
    im = ax.scatter(X, Y, c=lab, cmap=cmap, norm=norm)
    ax.plot(np.linspace(0.001, 5, 100), P_conc/(qeq* np.linspace(0.001, 5, 100)), ls='--', color='crimson',
            linewidth=4)
    cbar = fig.colorbar(im, ax=ax, ticks=tick)
    cbar.set_ticklabels(tick)
    plt.xlabel(r'$\mathrm{\widetilde{A}}$')
    plt.ylabel(r'$\mathrm{\widetilde{B}}$')
    plt.xticks(np.linspace(0, int(X.max().round()), int(X.max().round()) + 1))
    plt.yticks(np.linspace(0, int(Y.max().round()), int(Y.max().round()) + 1))
    plt.xlim([0, X.max()*1.01])
    plt.ylim([0, Y.max()*1.01])
    if map == True:
        true_label = [mapdict[value] for value in tick]
        cbar.set_ticklabels(true_label)


def plot_colored_data_variability_random(df,label,qeq,cmap='jet'):
    tick = np.linspace(1, df[label].nunique(), df[label].nunique())
    fig = plt.figure()
    ax = plt.gca()
    lab = df[label]
    cmap = plt.cm.jet
    import matplotlib
   # norm = matplotlib.colors.BoundaryNorm(np.arange(0.5, df[label].max() + 1, 1), cmap.N)
    X=df.A
    Y=df.B
    im = ax.scatter(X, Y, c=lab, cmap=cmap)#, norm=norm)
    ax.plot(np.linspace(0.1, 5, 100), 1/(qeq* np.linspace(0.1, 5, 100)), ls='--', color='crimson',
            linewidth=4)
    cbar = fig.colorbar(im, ax=ax, ticks=tick)
    cbar.set_ticklabels(tick)
    plt.xlabel(r'$\mathrm{\widetilde{A}}$')
    plt.ylabel(r'$\mathrm{\widetilde{B}}$')
    plt.xticks(np.linspace(0, int(X.max().round()), int(X.max().round()) + 1))
    plt.yticks(np.linspace(0, int(Y.max().round()), int(Y.max().round()) + 1))
    plt.xlim([0, X.max()*1.01])
    plt.ylim([0, Y.max()*1.01])



from mpl_toolkits.mplot3d import Axes3D
def plot_colored_data_random_3D(df,label,qeq,mapdict=dict(),map=False,cmap='jet'):
    tick = np.linspace(1, df[label].nunique(), df[label].nunique())
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    lab = df[label]
    cmap = plt.cm.jet
    import matplotlib
    norm = matplotlib.colors.BoundaryNorm(np.arange(0.5, df[label].max() + 1, 1), cmap.N)
    X=df.A
    Y=df.B
    Z=df.P
    im = ax.scatter(X, Y, Z,c=lab, cmap=cmap, norm=norm,alpha=0.5)
    ax.plot(np.linspace(0.1, 5, 100), 1/(qeq* np.linspace(0.1, 5, 100)), ls='--', color='crimson',
            linewidth=4)
    cbar = fig.colorbar(im, ax=ax, ticks=tick)
    cbar.set_ticklabels(tick)

    ax.view_init(30, 150)
    ax.set_xlabel(r'$\mathrm{\widetilde{A}}$')
    ax.set_ylabel(r'$\mathrm{\widetilde{B}}$')
    ax.set_zlabel(r'$\mathrm{\widetilde{P}}$')
    #plt.xticks(np.linspace(0, int(X.max().round()), int(X.max().round()) + 1))
    #plt.yticks(np.linspace(0, int(Y.max().round()), int(Y.max().round()) + 1))
    #ax.axes.set_xlim3d(left=0.0, right=1)
    #ax.axes.set_ylim3d(bottom=0.0, top=1)
    #ax.axes.set_zlim3d(bottom=0.0, top=1)
    if map == True:
        true_label = [mapdict[value] for value in tick]
        cbar.set_ticklabels(true_label)
    #
    # for i in (df.overall_gamma.unique()):
    #     sample=100
    #     A_conc = np.linspace(0.001,5,sample)
    #     P_conc = np.linspace(0.001,5,sample)
    #     B_conc=P_conc/A_conc/qeq/i
    #
    #
    #         #plt.plot(gammas[:,0]*gammas[:,1],gammas[:,2], color='crimson',
    #          #   alpha=0.1)  # crimson #or black
    #
    #     ax.scatter(A_conc, B_conc, P_conc, c='grey')


    ax.view_init(15, 180)

    return fig,ax




def return_knn_trained_model_k_ranges_random_3D(df,qeq,k=6,group=['']):
    X=df.A
    Y=df.B
    P=df.P
    xmax=X.max().round()
    ymax=Y.max().round()
    pmax=P.max().round()
    Z = df.label
    k = k
    X_train = pd.DataFrame([X, Y, P]).T

    knn_model = KNeighborsClassifier(n_neighbors=k)
    knn_model.fit(X_train, Z)

    xx, yy, pp = np.meshgrid(np.linspace((0.001), (xmax), 100),
                         np.linspace((0.001), (ymax), 100),
                            np.linspace((0.001), (ymax), 100))

    Z_test = knn_model.predict(np.c_[xx.ravel(), yy.ravel(),pp.ravel()])
    X_test = xx.ravel()
    Y_test = yy.ravel()
    P_test=pp.ravel()
    df_trained = pd.DataFrame(np.c_[X_test, Y_test,P_test, Z_test])
    df_trained.columns = ['A', 'B', 'P','label']
    df_trained['AB'] = df_trained['A'] * df_trained['B']
    df_trained = df_trained[df_trained['AB'] >= df_trained['P']/qeq]
    #df_trained['SP'] = df_trained['P'] / df_trained['S']
    #df_trained = df_trained[df_trained['SP'] <= qeq]

    return df_trained


import matplotlib
"Function to plot k ranges as in Heinrich Klipp et al."
def plot_k_ranges_with_knn_random_3D(df,qeq,mapdict=dict(),map=False):
    tick = np.linspace(1, df[label].nunique(), df[label].nunique())
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X=df.A
    Y=df.B
    P=df.P
    xmax=X.max().round()
    ymax=Y.max().round()
    pmax=P.max().round()
    Z = df.label
    fig = plt.figure()
    ax = plt.gca()
    level_number = np.linspace(1, Z.max() + 1, Z.nunique() + 1)
    "this is to put your own defined colours"
    colors =  ['mediumblue','blue','moccasin']
    cmap = matplotlib.colors.ListedColormap(colors, "", len(colors))
    "this one is when you want to put specific colors"
    #cs = ax.tricontourf(X, Y, Z, level_number, cmap=cmap)#cmap=plt.cm.jet)  #  ['mediumblue','blue','brown']
    cs = ax.tricontourf(X, Y, Z, level_number, cmap=plt.cm.jet)  #  ['mediumblue','blue','brown']

    c = ax.tricontour(X, Y, Z, level_number, linewidths=0.2, colors='black')
    #equilibrium line
    ax.plot(np.linspace(0.1, 5, 100), 1/(qeq* np.linspace(0.1, 5, 100)), ls='--', color='crimson',
            linewidth=4)
    ax.set_xlabel(r'$\mathrm{\widetilde{A}}$')
    ax.set_ylabel(r'$\mathrm{\widetilde{B}}$')
    #ax.yaxis.set_label_coords(-0.2, 0.95)
    ax.xaxis.set_label_coords(0.5, -0.1)

    plt.xlim([0, xmax])
    plt.ylim([0, ymax])
    plt.xticks(np.linspace(0, int(xmax), int(xmax)+1))
    plt.yticks(np.linspace(0, int(ymax), int(ymax)+1))

    cbar = fig.colorbar(cs, ax=ax,format='%d.')
    lab = np.arange(1, Z.nunique() + 1, 1)
    lab=np.linspace(1,Z.nunique(),Z.nunique())
    loc = lab + .5
    #real_label
    cbar.set_ticks(loc)

    if map == True:
        true_label = [mapdict[value] for value in lab]
        cbar.set_ticklabels(true_label)

    else:
        cbar.set_ticklabels(lab.round(0))
    #cbar.set_label(r'$labels$', rotation=0, labelpad=-40, y=1.12)

'this is a trial to be seen if this works '

"First idea was to P/B vs A so each line --> isoline but not unique always" \
"so now do overall gamma vs A/B"
def plot_colored_data_random_grouped(df, label, qeq, mapdict=dict(), map=False, cmap='jet'):
    tick = np.linspace(1, df[label].nunique(), df[label].nunique())
    fig = plt.figure()
    ax = plt.gca()
    lab = df[label]
    cmap = plt.cm.jet
    import matplotlib
    norm = matplotlib.colors.BoundaryNorm(np.arange(0.5, df[label].max() + 1, 1), cmap.N)
    X = (df.A/df.B)
    Y = df.overall_gamma
    #X = (df.A*df.B)
    #Y = (df.P)
    #X = (df.A**2*df.B**2)
    #Y = (df.A*df.P*df.B)
    im = ax.scatter(X, Y, c=lab, cmap=cmap, norm=norm)
    # ax.plot(np.linspace(0.1, 5, 100), 1 / (qeq * np.linspace(0.1, 5, 100)), ls='--', color='crimson',
    #       linewidth=4)
    cbar = fig.colorbar(im, ax=ax, ticks=tick)
    cbar.set_ticklabels(tick)
    ax.set_xscale('log')
    plt.xlabel(r'$log[\frac{\mathrm{\widetilde{A}}}{\mathrm{\widetilde{B}}}]$')
    plt.ylabel(r'$\mathrm{\Gamma}$')
    #plt.xticks(np.linspace(0, int(X.max().round()), int(X.max().round()) + 1))
    #plt.yticks(np.linspace(0, int(Y.max().round()), int(Y.max().round()) + 1))
    plt.xlim([X.min(), X.max() * 1.01])
    #plt.xlim([-10, 10])
    #plt.ylim([0, Y.max() * 1.01])
    #plt.xlim([0, X.max() * 1.01])
    plt.ylim([0, Y.max() * 1.01])

    if map == True:
        true_label = [mapdict[value] for value in tick]
        cbar.set_ticklabels(true_label)

    return fig, ax
"First idea was to P/B vs A so each line --> isoline but not unique always" \
"also partition into A-B regions"
def plot_colored_data_random_grouped_concentrationwise(df, label, qeq, mapdict=dict(), map=False, cmap='jet'):
    tick = np.linspace(1, df[label].nunique(), df[label].nunique())
    fig = plt.figure()
    ax = plt.gca()
    lab = df[label]
    cmap = plt.cm.jet
    import matplotlib
    norm = matplotlib.colors.BoundaryNorm(np.arange(0.5, df[label].max() + 1, 1), cmap.N)
    #X = (df.A/df.B)
    #Y = df.overall_gamma
    X = (df.A*df.B)
    Y = (df.P)
    #X = (df.A**2*df.B**2)
    #Y = (df.A*df.P*df.B)
    im = ax.scatter(X, Y, c=lab, cmap=cmap, norm=norm)
    # ax.plot(np.linspace(0.1, 5, 100), 1 / (qeq * np.linspace(0.1, 5, 100)), ls='--', color='crimson',
    #       linewidth=4)
    cbar = fig.colorbar(im, ax=ax, ticks=tick)
    cbar.set_ticklabels(tick)
    #ax.set_xscale('log')
    plt.xlabel(r'$\mathrm{\widetilde{A}}\mathrm{\widetilde{B}}$')
    plt.ylabel(r'$\mathrm{\widetilde{P}}$')
    #plt.xticks(np.linspace(0, int(X.max().round()), int(X.max().round()) + 1))
    #plt.yticks(np.linspace(0, int(Y.max().round()), int(Y.max().round()) + 1))
   # plt.xlim([0, X.max() * 1.01])
    #plt.xlim([-10, 10])
    #plt.ylim([0, Y.max() * 1.01])
    #plt.xlim([0, X.max() * 1.01])
    plt.ylim([0, Y.max() * 1.01])

    if map == True:
        true_label = [mapdict[value] for value in tick]
        cbar.set_ticklabels(true_label)

    return fig, ax

from sympy.solvers.solvers import check_assumptions