import numpy as np

'''in tutorials optimization is performed directly on normalized values of the parameters 
this function can be used to normalize metabolite concentrations and equilibrium constants''' \

def normalize_parameters(second_order_limit=1e8,first_order_limit=1e4,param='concentration',num_s=1,num_p=1,param_value=0.01):
    '''uppper_limit_2:second order rate_constant upper limit default 1e8 M^-1s^-1
    upper_limit_1:1st order rate constant upper limit default 1e4 s^-2
    param=string of either "concentration","eq_constant"
    num_s=number of substrates in the mechanism
    num_p=number of products in the mechanism
    param_value=value of the parameter before scaling put in SI units [M] for concentrations f([M,s]) 0.01M-->10 mM
    for equilibrium constaant
    '''
    charact_conc=first_order_limit/second_order_limit #[M]
    diff_s_p=int(num_s-num_p)
    if param=='concentration':
        scaled=param_value/charact_conc
    if param=='eq_constant':
        scaled = param_value*charact_conc**(diff_s_p)

    return scaled


def calculate_feasible_gamma_ranges(df,c_min,c_max):
    df['gamma_min'] = np.zeros(df.shape[0])
    df['gamma_max'] = np.zeros(df.shape[0])

    'take gamma_min and gamma_max values '
    for k in df.index:

        q_equilibrium = float(df.loc[k].k_eq_correct)

        sample_size = 50
        subs_conc = np.linspace(c_min, c_max, sample_size)
        prods_conc = np.linspace(c_min, c_max, sample_size)

        'here meshing takes place'
        'to have mesh of each possible gammas'

        df['gamma_max'].loc[k] = c_max/c_min/q_equilibrium
        df['gamma_min'].loc[k] = c_min/c_max/q_equilibrium

        if df['gamma_max'].loc[k]>1.0:
            df['gamma_max'].loc[k] = 0.99

        if df['gamma_min'].loc[k]<=1e-6:
            df['gamma_min'].loc[k] = 1e-6

        if df['gamma_min'].loc[k]>1.0:
            df['gamma_max'].loc[k] = 'not_feasible'
            df['gamma_min'].loc[k] = 'not_feasible'



    drop_rows = df[df['gamma_min'] == 'not_feasible'].index
    df = df.drop(drop_rows, axis=0)

    return df