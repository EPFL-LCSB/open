
from optlang.cplex_interface import Model, Variable, Constraint, Objective
import numpy as np
from open.optim.variability_analysis import *
import math
import pandas as pd


"for  3 step MILP version"

def milp_problem_3step_feasibility(v_net_lb,v_net_ub,gamma_ov, q=2.0, S=1.0, variability_analysis=False):
    "ideally this big M is e_TOT but here as it is scaled it is 1 E_tot=1"
    big_M=1.0
    'irreversible number of steps'
    n_step=6
    model_max = Model(name='Max model')
    "E+S-->EP-->E+P"
    #N=math.ceil((1-gamma_ov)/0.005)
    #0.001
    bin_size=1e-4
    N=math.ceil((1-gamma_ov)/bin_size)
    binary_variables=math.ceil(math.log2(N))
    print(binary_variables)
    P = S * q * gamma_ov
    epsilon=0.0
    model_max.add(Variable('E' ,lb=epsilon, ub=1))
    model_max.add(Variable('EP' ,lb=epsilon, ub=1))
    model_max.add(Variable('ES' ,lb=epsilon, ub=1))


    "binary variables  ceil(log_2_N)"
    #although it is to be seen if you should add the same number of binary variables for different approximated variables
    for i in range(binary_variables + 1):
        model_max.add(Variable('d' + str(i), type='binary'))
        model_max.add(Variable('d_a' + str(i), type='binary'))

    # for z1 to z6
    for i in range(n_step):
        model_max.add(Variable('z' + str(i + 1), lb=0))
        for j in range(binary_variables + 1):
            # todo you check if it can be 0 or not
            #Lets say k is fir gamma
            model_max.add(Variable('k_' + str(i + 1) + '_' + str(j), lb=0))

            #lets m is for the multiplocation variable gamma1*gamma2
            model_max.add(Variable('m_' + str(i + 1) + '_' + str(j), lb=0))

            # add the three constraints for peterson linearization
            model_max.add(Constraint(
                model_max.variables['k_' + str(i + 1) + '_' + str(j)] - big_M * model_max.variables['d' + str(j)],
                ub=0.0))
            model_max.add(Constraint(
                model_max.variables['k_' + str(i + 1) + '_' + str(j)] - model_max.variables['z' + str(i + 1)], ub=0.0))
            model_max.add(Constraint(model_max.variables['k_' + str(i + 1) + '_' + str(j)] - (
                        model_max.variables['z' + str(i + 1)] + big_M * model_max.variables['d' + str(j)]), lb=-1*big_M))

            #and three for multiplication variable
            model_max.add(Constraint(
                model_max.variables['m_' + str(i + 1) + '_' + str(j)] - big_M * model_max.variables['d_a' + str(j)],
                ub=0.0))
            model_max.add(Constraint(
                model_max.variables['m_' + str(i + 1) + '_' + str(j)] - model_max.variables['z' + str(i + 1)], ub=0.0))
            model_max.add(Constraint(model_max.variables['m_' + str(i + 1) + '_' + str(j)] - (
                        model_max.variables['z' + str(i + 1)] + big_M * model_max.variables['d_a' + str(j)]), lb=-1*big_M))


    #to fix it
    v=Variable('v', lb=v_net_lb, ub=v_net_ub)
    model_max.add(v)


    # for v related constraints
    for j in range(binary_variables + 1):
        model_max.add(Variable('k_' + str(n_step+1) + '_' + str(j), lb=0))
        model_max.add(
            Constraint(
                model_max.variables['k_' + str(n_step+1) + '_' + str(j)] - big_M * model_max.variables['d' + str(j)],
                ub=0.0))
        model_max.add(
            Constraint(model_max.variables['k_' + str(n_step+1) + '_' + str(j)] - model_max.variables['v'],
                       ub=0.0))
        model_max.add(Constraint(model_max.variables['k_' + str(n_step+1) + '_' + str(j)] - (
                model_max.variables['v'] + big_M * model_max.variables['d' + str(j)]), lb=-1*big_M))

        model_max.add(Variable('m_' + str(n_step+1) + '_' + str(j), lb=0))
        model_max.add(
            Constraint(
                model_max.variables['m_' + str(n_step+1) + '_' + str(j)] - big_M * model_max.variables['d_a' + str(j)],
                ub=0.0))
        model_max.add(
            Constraint(model_max.variables['m_' + str(n_step+1) + '_' + str(j)] - model_max.variables['v'],
                       ub=0.0))
        model_max.add(Constraint(model_max.variables['m_' + str(n_step+1) + '_' + str(j)] - (
                model_max.variables['v'] + big_M * model_max.variables['d_a' + str(j)]), lb=-1*big_M))


    #todo update these
    "Constraints"
    c_mass_cons = Constraint(model_max.variables['E']  + model_max.variables['ES'] + model_max.variables['EP'], lb=1, ub=1)
    c_1 = Constraint(model_max.variables['z1']  - model_max.variables['ES'],  ub=0)  # k1f
    c_2 = Constraint(model_max.variables['z2']  - model_max.variables['EP'],  ub=0)  # k2f
    c_3 = Constraint(model_max.variables['z3']  - model_max.variables['E'],  ub=0)  # k1b
    c_4 = Constraint(model_max.variables['z4']  - model_max.variables['E'],  ub=0)  # k1B
    c_5 = Constraint(model_max.variables['z5']  - model_max.variables['ES'],  ub=0)  # k1f
    c_6 = Constraint(model_max.variables['z6']  - model_max.variables['EP'],  ub=0)  # k2f


    gamma_1= Variable('gamma_1', lb=gamma_ov, ub=1)
    model_max.add(gamma_1)

    #this is the multiplication variable
    a= Variable('a', lb=gamma_ov, ub=1)
    model_max.add(a)

    expression_gamma=0+model_max.variables['gamma_1'].lb
    multip_factor=(model_max.variables['gamma_1'].ub-model_max.variables['gamma_1'].lb)/N
    for i in range(binary_variables + 1):
        expression_gamma+=2**i*model_max.variables['d' + str(i)]*multip_factor


    expression_a=0+model_max.variables['a'].lb
    multip_factor_a=(model_max.variables['a'].ub-model_max.variables['a'].lb)/N
    for i in range(binary_variables + 1):
        expression_a+=2**i*model_max.variables['d_a' + str(i)]*multip_factor_a

    c_gamma=Constraint(gamma_1-expression_gamma,lb=0,ub=0)

    c_a=Constraint(a-expression_a,lb=0,ub=0)

    c_a_gamma=Constraint(a-gamma_1, ub=0)

    model_max.add([c_mass_cons,c_1,c_2,c_3,c_4,c_5,c_6,c_gamma,c_a,c_a_gamma])
    #write expressions in for loop

    expression={}
    expression_m={}
    for i in range(n_step):
        expression[i]=model_max.variables['z'+str(i+1)]*gamma_ov
        expression_m[i]=model_max.variables['z'+str(i+1)]*gamma_ov

        for j in range(binary_variables + 1):
            expression[i]+=multip_factor*2**j*(model_max.variables['k_'+str(i+1)+'_'+str(j)])
            expression_m[i]+=multip_factor_a*2**j*(model_max.variables['m_'+str(i+1)+'_'+str(j)])


    expression[n_step]=model_max.variables['v']*gamma_ov
    expression_m[n_step]=model_max.variables['v']*gamma_ov
    for j in range(binary_variables + 1):
        expression[n_step] += multip_factor * 2 ** j * (model_max.variables['k_' + str(n_step+1) + '_' + str(j)])
        expression_m[n_step] += multip_factor_a * 2 ** j * (model_max.variables['m_' + str(n_step+1) + '_' + str(j)])


    #todo to check
    c_rate_1=Constraint((model_max.variables['z1']-expression[0])-expression[n_step],lb=0, ub=0)
    c_rate_2=Constraint(expression_m[1]-model_max.variables['z2']*gamma_ov-model_max.variables['v']*gamma_ov,lb=0, ub=0)
    c_rate_3=Constraint(P*(expression[2]-expression_m[2])-expression_m[n_step],lb=0, ub=0)
    c_rate_4=Constraint(S*(model_max.variables['z4']-expression[3])-model_max.variables['v'],lb=0, ub=0)
    c_rate_5=Constraint((expression_m[4]-model_max.variables['z5']*gamma_ov)-expression_m[n_step],lb=0, ub=0)
    c_rate_6=Constraint((expression[5]-expression_m[5])-expression[n_step],lb=0, ub=0)


    model_max.add([c_rate_1,c_rate_2,c_rate_3,c_rate_4,c_rate_5,c_rate_6])
    obj = Objective(0, direction='max')
    model_max.objective = obj
    model_max.configuration.tolerances.feasibility = 1e-9
    model_max.configuration.tolerances.optimality = 1e-9
    model_max.configuration.tolerances.integrality = 1e-9

    status = model_max.optimize()
    print(status)

    "Post-Processing"
    dict_values = dict()
    for key in model_max.variables.keys():
        dict_values[key] = (model_max.variables[key].primal)
        #print(key, (model_max.variables[key].primal))

    # add also these values to have the scaled version
    dict_values['P'] = P
    dict_values['S'] = S

    # dict_values['gamma_1'] = gamma_1
    # dict_values['gamma_2'] = gamma_2

    dict_values['gamma_2'] = gamma_ov / a.primal
    dict_values['gamma_3'] = a.primal / gamma_1.primal
    dict_values['gamma_ov'] = gamma_ov

    dict_values['k1f'] = model_max.variables['z4'].primal /model_max.variables['E'].primal
    dict_values['k2f'] = model_max.variables['z5'].primal /model_max.variables['ES'].primal
    dict_values['k3f'] = model_max.variables['z6'].primal /model_max.variables['EP'].primal

    dict_values['k1b'] = model_max.variables['z1'].primal /model_max.variables['ES'].primal
    dict_values['k2b'] = model_max.variables['z2'].primal /model_max.variables['EP'].primal
    dict_values['k3b'] = model_max.variables['z3'].primal /model_max.variables['E'].primal


    #todo variability analysis
    if variability_analysis == True:
        df_var_analysis = pd.DataFrame()
        df_var_analysis= variability_analysis_linearized_feasibility(model_max)
        # print(df_var_analysis)
        #dict_values['unique'] = a
        # return (np.exp(model_max.variables.ln_v.primal)),(dict_values),df_var_analysis
    # for fitness you are directly returning the v value not log format

    return ((model_max.variables.v.primal)), (dict_values),(model_max), (df_var_analysis)

'for fixed v can also make slices or '
def linear_problem_for_3_step_linearized_feasibility(v_net_lb,v_net_ub,gamma_1, gamma_2, gamma_3, gamma_ov, q=2.0, S=1.0, variability_analysis=False):
    epsilon = -50
    big_M=1
    #k1b = Variable('k1b', lb=0, ub=1)
    #k2b = Variable('k2b', lb=0, ub=1)
    #k3b = Variable('k3b', lb=0, ub=1)

    #k1f = Variable('k1f', lb=0, ub=1)
    #k2f = Variable('k2f', lb=0, ub=1)
    #k3f = Variable('k3f', lb=0, ub=1)
    E = Variable('E', lb=0, ub=1)
    ES = Variable('ES', lb=0, ub=1)
    EP= Variable('EP', lb=0, ub=1)

    v = Variable('v', lb=v_net_lb, ub=v_net_ub)
    #v=v_net
    #Define binary variables for each multiplication

    z1 = Variable('z1', lb=0)# , ub=1)
    z2 = Variable('z2', lb=0)# , ub=1)
    z3 = Variable('z3', lb=0)# , ub=1)

    z4 = Variable('z4', lb=0) # , ub=1)
    z5 = Variable('z5', lb=0 )# , ub=1)
    z6 = Variable('z6', lb=0 )#, ub=1)

    P = S * q * gamma_ov

    #S = Variable('ln_S', lb=np.log(S), ub=np.log(S))
    #P = Variable('ln_P', lb=np.log(0.001), ub=np.log(10.0))

    "Constraints"
    # 1st order rate constants
    c1_max = Constraint(z1*(1 - gamma_1)/(gamma_1)-v , lb=0, ub=0) #z1
    c2_max = Constraint(z2*(1 - gamma_2)/(gamma_2)-v , lb=0, ub=0) #z2
    c3_max = Constraint(z3*P*(1 - gamma_3)/(gamma_3)-v , lb=0, ub=0) #z3

    c1f_max = Constraint(z4*S*(1 - gamma_1)-v , lb=0, ub=0) #z4
    c2f_max = Constraint(z5*(1 - gamma_2)-v , lb=0, ub=0) #z5
    c3f_max = Constraint(z6*(1 - gamma_3)-v , lb=0, ub=0) #z6

    #mass conservation for enzyme
    c_mass_cons = Constraint(E+ES+EP, lb=1, ub=1) #z6


    "binary constraints even no need for the remaining two"
    "similar to Petersen linearization but you dont have discrete variables so " \
    "no  need for additional 2 constraints"

    b1_1=Constraint(ES-z1,lb=0)
    #b1_2=Constraint(k1b-z1,lb=0)
    #b1_3=Constraint(z1-ES-k1b,lb=-1.0)

    b2_1 = Constraint(EP - z2, lb=0)
    #b2_2 = Constraint(k2b - z2, lb=0)
    #b2_3 = Constraint(z2 - EP - k2b, lb=-1.0)

    b3_1 = Constraint(E - z3, lb=0)
    #b3_2 = Constraint(k3b - z3, lb=0)
    #b3_3 = Constraint(z3 - E - k3b, lb=-1.0)

    b4_1 = Constraint(E - z4, lb=0)
    #b4_2 = Constraint(k1f - z4, lb=0)
    #b4_3 = Constraint(z4 - E - k1f, lb=-1.0)

    b5_1 = Constraint(ES - z5, lb=0)
    #b5_2 = Constraint(k2f - z5, lb=0)
    #b5_3 = Constraint(z5 - ES - k2f, lb=-1.0)

    b6_1 = Constraint(EP - z6, lb=0)
    #b6_2 = Constraint(k3f - z6, lb=0)
    #b6_3 = Constraint(z6 - EP - k3f, lb=-1.0)


    # P/AB constraint
    #cK_eq = Constraint(ln_P - ln_S - np.log(q * gamma_ov), lb=0, ub=0)

    # eq constant constraint
    #  cK_eq2 = Constraint(ln_k1f+ln_k3f+ln_k5f+ln_k6f-(ln_k1b+ln_k3b+ln_k5b+ln_k6b)-np.log(q) , lb=0, ub=0)
    #  cK_eq3 = Constraint(ln_k2f+ln_k4f+ln_k5f+ln_k6f-(ln_k2b+ln_k4b+ln_k5b+ln_k6b)-np.log(q) , lb=0, ub=0)

    obj = Objective(0, direction='max')

    model_max = Model(name='Max model')
    model_max.objective = obj
    model_max.configuration.tolerances.feasibility = 1e-9
   # model_max.add([c1_max, c2_max, c3_max, c1f_max, c2f_max, c3f_max, c_mass_cons, b1_1, b1_2, b1_3, b2_1, b2_2 ,b2_3, \
   #                b3_1, b3_2, b3_3, b4_1, b4_2, b4_3, b5_1, b5_2, b5_3,b6_1, b6_2, b6_3])  # ,cK_eq2, cK_eq3 ])


    model_max.add([c1_max, c2_max, c3_max, c1f_max, c2f_max, c3f_max, c_mass_cons, b1_1, b2_1, \
                   b3_1,b4_1, b5_1, b6_1])

    status = model_max.optimize()
    # print(status)

    ##until here too
    """
    Perform variability analysis as it is done in pyTFA it
    """

    # to check for A,B pairs if they are outside of the bounds I want to sample in just to assure but not useful
    if status == 'infeasible':
        print('NOT FEASIBLE')

    """
    store your k values 
    """
    # for 1 optimization problem variables are stored in a dict and returned along with the fitness value
    dict_values = dict()
    for key in model_max.variables.keys():
        dict_values[key] = (eval(key).primal)

    # add also these values to have the scaled version
    #dict_values['E'] = E
    #dict_values['ES'] = ES
    #dict_values['EP'] = EP
    dict_values['P']=P
    dict_values['S']=S

    dict_values['gamma_1'] = gamma_1
    dict_values['gamma_2'] = gamma_2
    dict_values['gamma_3'] = gamma_3
    dict_values['gamma_ov'] = gamma_ov



    dict_values['k1b']=z1.primal/ES.primal
    dict_values['k2b']=z2.primal/EP.primal
    dict_values['k3b']=z3.primal/E.primal

    dict_values['k1f']=z4.primal/E.primal
    dict_values['k2f']=z5.primal/ES.primal
    dict_values['k3f']=z6.primal/EP.primal
    dict_values['v']=v


#todo you can change this depending on which results you want to take
    if variability_analysis == True:
        df_var_analysis = pd.DataFrame()
        df_var_analysis= variability_analysis_linearized_feasibility(model_max)
        #print(df_var_analysis)
        #dict_values['unique'] = a

        # return (np.exp(model_max.variables.ln_v.primal)),(dict_values),df_var_analysis
    # for fitness you are directly returning the v value not log format
    return ((model_max.objective.value)), (dict_values), (model_max),(df_var_analysis)





'random order mechanism where alpha is a variable of the MILP'
def milp_problem_4step_biuni_random_split_ratio_feasibility_with_gammas(v_net_lb,v_net_ub,gamma_1,gamma_2, \
                                                                        gamma_3,gamma_4,gamma_5, gamma_6,gamma_ov, q=2.0, S=1.0, P=1.0,variability_analysis=False):
    "ideally this big M is e_TOT but here as it is scaled it is 1 E_tot=1"
    'irreversible number of steps'
    n_step=12
    model_max = Model(name='Max model')
    "E+S-->EP-->E+P"

    A=S
    B = P/(A*q*gamma_ov)

    epsilon=0
    model_max.add(Variable('E' ,lb=epsilon, ub=1))
    model_max.add(Variable('EA' ,lb=epsilon, ub=1))
    model_max.add(Variable('EAB' ,lb=epsilon, ub=1))
    model_max.add(Variable('EP' ,lb=epsilon, ub=1))
    #this enzyme state is for the random case
    model_max.add(Variable('EB' ,lb=epsilon, ub=1))


    # for z1 to z6
    for i in range(n_step):
        model_max.add(Variable('z' + str(i + 1), lb=0))


    v=Variable('v', lb=v_net_lb, ub=v_net_ub)
    model_max.add(v)
    'here you denote v_lower =v*(1-alpha)'
    'and v_upper=v*alpha'
    model_max.add(Variable('v_upper' ,lb=0))
    model_max.add(Variable('v_lower' ,lb=0))


    #todo update these
    "Constraints"
    c_mass_cons = Constraint(model_max.variables['E']  + model_max.variables['EA'] + model_max.variables['EB']+model_max.variables['EAB']+ model_max.variables['EP'], lb=1, ub=1)
    c_1 = Constraint(model_max.variables['z1']  - model_max.variables['EA'],  ub=0)  # k1b
    c_2 = Constraint(model_max.variables['z2']  - model_max.variables['EAB'],  ub=0)  # k2b
    c_3 = Constraint(model_max.variables['z3']  - model_max.variables['EP'],  ub=0)  # k3b
    c_4 = Constraint(model_max.variables['z4']  - model_max.variables['EAB'],  ub=0)  # k3f
    c_5 = Constraint(model_max.variables['z5']  - model_max.variables['EP'],  ub=0)  # k4f
    c_6 = Constraint(model_max.variables['z6']  - model_max.variables['E'],  ub=0)  # k1f
    c_7 = Constraint(model_max.variables['z7']  - model_max.variables['EA'],  ub=0)  # k2f
    c_8 = Constraint(model_max.variables['z8']  - model_max.variables['E'],  ub=0)  # k4b

    c_9 = Constraint(model_max.variables['z9']  - model_max.variables['E'],  ub=0)  # k5f

    c_10 = Constraint(model_max.variables['z10']  - model_max.variables['EB'],  ub=0)  # k5b
    c_11 = Constraint(model_max.variables['z11']  - model_max.variables['EB'],  ub=0)  # k6f
    c_12 = Constraint(model_max.variables['z12']  - model_max.variables['EAB'],  ub=0)  # k6b



    model_max.add([c_mass_cons,c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_10, c_11,c_12]) #,c_gamma3_a,c_gamma3_a_1,c_gamma3_a_2])
    #write expressions in for loop

    #todo to UPDATE
    c_rate_1=Constraint(model_max.variables['z1']*(1-gamma_1)/gamma_1-model_max.variables['v_upper'],lb=0, ub=0) #k1b
    c_rate_2=Constraint(model_max.variables['z2']*(1-gamma_2)/gamma_2-model_max.variables['v_upper'],lb=0, ub=0) #k2b
    c_rate_3=Constraint(model_max.variables['z3']*(1-gamma_3)/gamma_3-model_max.variables['v'],lb=0, ub=0) #k3b
    c_rate_4=Constraint(model_max.variables['z4']*(1-gamma_3)-model_max.variables['v'],lb=0, ub=0)  #k3f
    c_rate_5=Constraint(model_max.variables['z5']*(1-gamma_4)-model_max.variables['v'],lb=0, ub=0) #k4f

    c_rate_6=Constraint(A*(model_max.variables['z6']*(1-gamma_1))-model_max.variables['v_upper'],lb=0, ub=0) #k1f
    c_rate_7=Constraint(B*(model_max.variables['z7']*(1-gamma_2))-model_max.variables['v_upper'],lb=0, ub=0) #k2f

    c_rate_8=Constraint(P*(model_max.variables['z8']*(1-gamma_4)/gamma_4)-model_max.variables['v'],lb=0, ub=0) #k4b

    c_rate_9=Constraint(B*(model_max.variables['z9']*(1-gamma_5))-model_max.variables['v_lower'],lb=0, ub=0)  #k5f
    c_rate_10=Constraint((model_max.variables['z10']*(1-gamma_5)/gamma_5)-model_max.variables['v_lower'],lb=0, ub=0)#k5b
    c_rate_11=Constraint(A*(model_max.variables['z11']*(1-gamma_6))-model_max.variables['v_lower'],lb=0, ub=0) #k6f
    c_rate_12=Constraint(model_max.variables['z12']*(1-gamma_6)/gamma_6-model_max.variables['v_lower'],lb=0, ub=0) #k6b

    c_split_ratio = Constraint(model_max.variables['v'] - (model_max.variables['v_upper']+model_max.variables['v_lower']),lb=0, ub=0)

    model_max.add([c_rate_1,c_rate_2,c_rate_3,c_rate_4,c_rate_5,c_rate_6,c_rate_7,c_rate_8, c_rate_9, c_rate_10, c_rate_11, c_rate_12,c_split_ratio])


    obj = Objective(0, direction='max')
    model_max.objective = obj
    model_max.configuration.tolerances.feasibility = 1e-9
    status = model_max.optimize()
    print(status)

    "Post-Processing"
    dict_values = dict()
    for key in model_max.variables.keys():
        dict_values[key] = (model_max.variables[key].primal)
        #print(key, (model_max.variables[key].primal))

    # add also these values to have the scaled version
    dict_values['A'] = A
    dict_values['B'] = B
    dict_values['P'] = P
    #dict_values['alpha'] = alpha

    dict_values['gamma_1'] = gamma_1
    dict_values['gamma_2'] = gamma_2
    dict_values['gamma_3'] = gamma_3
    dict_values['gamma_4'] = gamma_4
    dict_values['gamma_5'] = gamma_5
    dict_values['gamma_6'] = gamma_6

    dict_values['gamma_ov'] = gamma_ov
    #
    # dict_values['k1f'] = model_max.variables['z6'].primal /model_max.variables['E'].primal
    # dict_values['k2f'] = model_max.variables['z7'].primal /model_max.variables['EA'].primal
    # dict_values['k3f'] = model_max.variables['z4'].primal /model_max.variables['EAB'].primal
    # dict_values['k4f'] = model_max.variables['z5'].primal /model_max.variables['EP'].primal
    # 'for random'
    # dict_values['k5f'] = model_max.variables['z9'].primal /model_max.variables['E'].primal
    # dict_values['k6f'] = model_max.variables['z11'].primal /model_max.variables['EB'].primal
    #
    # dict_values['k1b'] = model_max.variables['z1'].primal /model_max.variables['EA'].primal
    # dict_values['k2b'] = model_max.variables['z2'].primal /model_max.variables['EAB'].primal
    # dict_values['k3b'] = model_max.variables['z3'].primal /model_max.variables['EP'].primal
    # dict_values['k4b'] = model_max.variables['z8'].primal /model_max.variables['E'].primal
    # 'for random'
    # dict_values['k5b'] = model_max.variables['z10'].primal /model_max.variables['EB'].primal
    # dict_values['k6b'] = model_max.variables['z12'].primal /model_max.variables['EAB'].primal
    #
    # dict_values['alpha'] = model_max.variables['v_upper'].primal /model_max.variables['v'].primal

    #todo variability analysis
    if variability_analysis == True:
        df_var_analysis = pd.DataFrame()
        df_var_analysis= variability_analysis_linearized_feasibility(model_max)
        # print(df_var_analysis)
        #dict_values['unique'] = a
        # return (np.exp(model_max.variables.ln_v.primal)),(dict_values),df_var_analysis
    # for fitness you are directly returning the v value not log format

    return ((model_max.variables.v.primal)), (dict_values),(model_max), (df_var_analysis)



