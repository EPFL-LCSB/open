from optlang.cplex_interface import Model, Variable, Constraint, Objective
import numpy as np
from open.optim.variability_analysis import *
import math
import pandas as pd

def milp_problem_5step_bibi(gamma_ov, q=2.0, S1=1.0,S2=1.0,P1=1.0, variability_analysis=False):
    "ideally this big M is e_TOT but here as it is scaled it is 1 E_tot=1"
    big_M=1.0
    'irreversible number of steps'
    n_step=10
    model_max = Model(name='Max model')
    "E+S-->EP-->E+P"
    #N=math.ceil((1-gamma_ov)/0.005)
    #N=2000
    if (gamma_ov>=0.7 and gamma_ov<0.9):
         bin_resol=1e-3
    elif gamma_ov>=0.9:
        bin_resol=1e-4
    else:
        bin_resol=1e-3


    N=math.ceil((1-gamma_ov)/bin_resol)

    binary_variables=math.ceil(math.log2(N))
    print(binary_variables)

    #A=S
    P2=S1*S2*q*gamma_ov/P1
    epsilon=0

    #5 enzyme states
    model_max.add(Variable('E' ,lb=epsilon, ub=1))
    model_max.add(Variable('EA' ,lb=epsilon, ub=1))
    model_max.add(Variable('EAB' ,lb=epsilon, ub=1))
    model_max.add(Variable('EPQ' ,lb=epsilon, ub=1))
    model_max.add(Variable('EP' ,lb=epsilon, ub=1))


    "binary variables  ceil(log_2_N)"
    #although it is to be seen if you should add the same number of binary variables for different approximated variables
    #gamma1,gamma1-2,gamma3,gamma3-4
    for i in range(binary_variables + 1):
        'for gamma1'
        model_max.add(Variable('d' + str(i), type='binary'))
        'for a=gamma1*gamma2'
        model_max.add(Variable('d_a' + str(i), type='binary'))
        'for gamma3'
        model_max.add(Variable('d_3' + str(i), type='binary'))
        'for gamma3-4 b=gamma3*gamma4'
        model_max.add(Variable('d_b' + str(i), type='binary'))
        'for gamma3a multiplication'

    'you also need to add other constraints binary for this t based on d_a and d_b'

    for j in range(binary_variables + 1):
        for k in range(binary_variables + 1):
            model_max.add(Variable('t_' + str(j) + '_' +str(k), type='binary'))

            model_max.add(Constraint(
            model_max.variables['t_' + str(j) + '_' +str(k)] - model_max.variables['d_b' + str(j)] ,
            ub=0.0))
            model_max.add(Constraint(
            model_max.variables['t_' + str(j) + '_' +str(k)] - model_max.variables['d_a' + str(k)] ,
            ub=0.0))
            model_max.add(Constraint(
            model_max.variables['t_' + str(j) + '_' +str(k)] - model_max.variables['d_a' + str(k)] - model_max.variables['d_b' + str(j)] +1 ,
            lb=0.0))



    # for z1 to z6
    for i in range(n_step):
        model_max.add(Variable('z' + str(i + 1), lb=0))
        for j in range(binary_variables + 1):
            # todo you check if it can be 0 or not
            #Lets say k is fir gamma d
            model_max.add(Variable('k_' + str(i + 1) + '_' + str(j), lb=0))

            #lets m is for the multiplocation variable gamma1*gamma2 d_a
            model_max.add(Variable('m_' + str(i + 1) + '_' + str(j), lb=0))

            # k3 for gamma 3 d_3
            model_max.add(Variable('k3_' + str(i + 1) + '_' + str(j), lb=0))

            #lets p is for the multiplocation variable gamma3*gamma4 d_b
            model_max.add(Variable('p_' + str(i + 1) + '_' + str(j), lb=0))

            # for t binary multiplication variables and the continuous one z
            #model_max.add(Variable('mt_' + str(i + 1) + '_' + str(j), lb=0))

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


            #and three for multiplication variable
            model_max.add(Constraint(
                model_max.variables['k3_' + str(i + 1) + '_' + str(j)] - big_M * model_max.variables['d_3' + str(j)],
                ub=0.0))
            model_max.add(Constraint(
                model_max.variables['k3_' + str(i + 1) + '_' + str(j)] - model_max.variables['z' + str(i + 1)], ub=0.0))
            model_max.add(Constraint(model_max.variables['k3_' + str(i + 1) + '_' + str(j)] - (
                        model_max.variables['z' + str(i + 1)] + big_M * model_max.variables['d_3' + str(j)]), lb=-1*big_M))


            #and three for multiplication variable p and d_b
            model_max.add(Constraint(
                model_max.variables['p_' + str(i + 1) + '_' + str(j)] - big_M * model_max.variables['d_b' + str(j)],
                ub=0.0))
            model_max.add(Constraint(
                model_max.variables['p_' + str(i + 1) + '_' + str(j)] - model_max.variables['z' + str(i + 1)], ub=0.0))
            model_max.add(Constraint(model_max.variables['p_' + str(i + 1) + '_' + str(j)] - (
                        model_max.variables['z' + str(i + 1)] + big_M * model_max.variables['d_b' + str(j)]), lb=-1*big_M))




            #and three for multiplication variable
           # model_max.add(Constraint(
            #    model_max.variables['mt_' + str(i + 1) + '_' + str(j)] - big_M * model_max.variables['t' + str(j)],
             #   ub=0.0))
            #model_max.add(Constraint(
             #   model_max.variables['mt_' + str(i + 1) + '_' + str(j)] - model_max.variables['z' + str(i + 1)], ub=0.0))
            #model_max.add(Constraint(model_max.variables['mt_' + str(i + 1) + '_' + str(j)] - (
             #           model_max.variables['z' + str(i + 1)] + big_M * model_max.variables['t' + str(j)]), lb=-1*big_M))


    #z*t
    for i in range(n_step):
        for j in range(binary_variables + 1):
            for k in range(binary_variables + 1):
                model_max.add(Variable('mt_' + str(i + 1) + '_' + str(j)+  '_' + str(k), lb=0))

                model_max.add(Constraint(
                    model_max.variables['mt_' + str(i + 1) + '_' + str(j)+  '_' + str(k)] - big_M * model_max.variables['t_' + str(j) + '_' +str(k)],
                    ub=0.0))
                model_max.add(Constraint(
                    model_max.variables['mt_' + str(i + 1) + '_' + str(j)+  '_' + str(k)] - model_max.variables['z' + str(i + 1)],
                    ub=0.0))
                model_max.add(Constraint(model_max.variables['mt_' + str(i + 1) + '_' + str(j)+  '_' + str(k)] - (
                        model_max.variables['z' + str(i + 1)] + big_M * model_max.variables['t_' + str(j) + '_' +str(k)]),
                                         lb=-1 * big_M))

    v=Variable('v', lb=0)
    model_max.add(v)


    # for v related constraints
    for j in range(binary_variables + 1):
        "for k"
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
        "for m"
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

        "for k3"
        model_max.add(Variable('k3_' + str(n_step+1) + '_' + str(j), lb=0))
        model_max.add(
            Constraint(
                model_max.variables['k3_' + str(n_step+1) + '_' + str(j)] - big_M * model_max.variables['d_3' + str(j)],
                ub=0.0))
        model_max.add(
            Constraint(model_max.variables['k3_' + str(n_step+1) + '_' + str(j)] - model_max.variables['v'],
                       ub=0.0))
        model_max.add(Constraint(model_max.variables['k3_' + str(n_step+1) + '_' + str(j)] - (
                model_max.variables['v'] + big_M * model_max.variables['d_3' + str(j)]), lb=-1*big_M))


        "for p"
        model_max.add(Variable('p_' + str(n_step+1) + '_' + str(j), lb=0))
        model_max.add(
            Constraint(
                model_max.variables['p_' + str(n_step+1) + '_' + str(j)] - big_M * model_max.variables['d_b' + str(j)],
                ub=0.0))
        model_max.add(
            Constraint(model_max.variables['p_' + str(n_step+1) + '_' + str(j)] - model_max.variables['v'],
                       ub=0.0))
        model_max.add(Constraint(model_max.variables['p_' + str(n_step+1) + '_' + str(j)] - (
                model_max.variables['v'] + big_M * model_max.variables['d_b' + str(j)]), lb=-1*big_M))



    for j in range(binary_variables + 1):
        for k in range(binary_variables + 1):

            "for mt"
            model_max.add(Variable('mt_' + str(n_step+1) + '_' + str(j)+ '_' + str(k), lb=0))
            model_max.add(
            Constraint(
                model_max.variables['mt_' + str(n_step+1) + '_' + str(j)+ '_' + str(k)] - big_M * model_max.variables['t_' + str(j)+ '_' + str(k)],
                ub=0.0))
            model_max.add(
            Constraint(model_max.variables['mt_' + str(n_step+1) + '_' + str(j)+ '_' + str(k)] - model_max.variables['v'],
                       ub=0.0))
            model_max.add(Constraint(model_max.variables['mt_' + str(n_step+1) + '_' + str(j)+ '_' + str(k)] - (
                model_max.variables['v'] + big_M * model_max.variables['t_' + str(j)+ '_' + str(k)]), lb=-1*big_M))


    #todo update these
    "Constraints"
    c_mass_cons = Constraint(model_max.variables['E']  + model_max.variables['EA'] + model_max.variables['EAB']+ model_max.variables['EPQ']+ model_max.variables['EP'], lb=1, ub=1)
    c_1 = Constraint(model_max.variables['z1']  - model_max.variables['EA'],  ub=0)  # k1b
    c_2 = Constraint(model_max.variables['z2']  - model_max.variables['EAB'],  ub=0)  # k2b
    c_3 = Constraint(model_max.variables['z3']  - model_max.variables['EPQ'],  ub=0)  # k3b
    c_4 = Constraint(model_max.variables['z4']  - model_max.variables['EP'],  ub=0)  # k4b
    c_5 = Constraint(model_max.variables['z5']  - model_max.variables['E'],  ub=0)  # k5b

    c_6 = Constraint(model_max.variables['z6']  - model_max.variables['E'],  ub=0)  # k1f
    c_7 = Constraint(model_max.variables['z7']  - model_max.variables['EA'],  ub=0)  # k2f
    c_8 = Constraint(model_max.variables['z8']  - model_max.variables['EAB'],  ub=0)  # k3f
    c_9 = Constraint(model_max.variables['z9']  - model_max.variables['EPQ'],  ub=0)  # k3f
    c_10 = Constraint(model_max.variables['z10']  - model_max.variables['EP'],  ub=0)  # k3f



    gamma_1= Variable('gamma_1', lb=gamma_ov, ub=1)
    model_max.add(gamma_1)

    #this is the multiplication variable
    a= Variable('a', lb=gamma_ov, ub=1)
    model_max.add(a)

    gamma_3= Variable('gamma_3', lb=gamma_ov, ub=1)
    model_max.add(gamma_3)

    b= Variable('b', lb=gamma_ov, ub=1)
    model_max.add(b)

    b_a= Variable('b_a', lb=gamma_ov, ub=1)
    model_max.add(b_a)


    expression_gamma=0+model_max.variables['gamma_1'].lb
    multip_factor=(model_max.variables['gamma_1'].ub-model_max.variables['gamma_1'].lb)/N
    for i in range(binary_variables + 1):
        expression_gamma+=2**i*model_max.variables['d' + str(i)]*multip_factor


    expression_a=0+model_max.variables['a'].lb
    multip_factor_a=(model_max.variables['a'].ub-model_max.variables['a'].lb)/N
    for i in range(binary_variables + 1):
        expression_a+=2**i*model_max.variables['d_a' + str(i)]*multip_factor_a


    expression_gamma3_=0+model_max.variables['gamma_3'].lb
    multip_factor_gamma3=(model_max.variables['gamma_3'].ub-model_max.variables['gamma_3'].lb)/N
    for i in range(binary_variables + 1):
        expression_gamma3_+=2**i*model_max.variables['d_3' + str(i)]*multip_factor_gamma3

    #added for b
    expression_b = 0 + model_max.variables['b'].lb
    multip_factor_b = (model_max.variables['b'].ub - model_max.variables['b'].lb) / N
    for i in range(binary_variables + 1):
        expression_a += 2 ** i * model_max.variables['d_b' + str(i)] * multip_factor_b

    #convert this to expression a_b
    expression_b_a=0+model_max.variables['b_a'].lb**2
    for j in range(binary_variables + 1):
        expression_b_a += gamma_ov * multip_factor_a * 2 ** j * (model_max.variables['d_a' + str(j)]) + \
                           gamma_ov * multip_factor_b* 2 ** j * (model_max.variables['d_b' + str(j)])
        for k in range(binary_variables + 1):
            expression_b_a += multip_factor_b * multip_factor_a * 2 ** (k + j) * (model_max.variables['t_' + str(j) + '_' +str(k)])



    c_gamma=Constraint(gamma_1-expression_gamma,lb=0,ub=0)

    c_a=Constraint(a-expression_a,lb=0,ub=0)

    c_b=Constraint(b-expression_b,lb=0,ub=0)

    c_a_gamma=Constraint(a-gamma_1, ub=0)

    c_b_gamma=Constraint(b-gamma_3, ub=0)

    c_gamma3=Constraint(gamma_3-expression_gamma3_,lb=0,ub=0)

    #if this is required or not

    c_b_a=Constraint(b_a-expression_b_a,lb=0,ub=0)

    #c_gamma3_a_1=Constraint(gamma_3_a-gamma_3,ub=0)

    #c_gamma3_a_2=Constraint(gamma_3_a-a,ub=0)


    model_max.add([c_mass_cons,c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_10,c_gamma,c_a,c_b,c_a_gamma,c_b_gamma,c_gamma3, c_b_a ]) #,c_gamma3_a,c_gamma3_a_1,c_gamma3_a_2])
    #write expressions in for loop

    expression={}
    expression_m={}
    expression_p = {}
    expression_gamma3={}
    expression_t={}
    for i in range(n_step):
        expression[i]=model_max.variables['z'+str(i+1)]*gamma_ov
        expression_m[i]=model_max.variables['z'+str(i+1)]*gamma_ov
        expression_p[i]=model_max.variables['z'+str(i+1)]*gamma_ov
        expression_gamma3[i]=model_max.variables['z'+str(i+1)]*gamma_ov
        expression_t[i]=model_max.variables['z'+str(i+1)]*gamma_ov**2


        for j in range(binary_variables + 1):
            expression[i]+=multip_factor*2**j*(model_max.variables['k_'+str(i+1)+'_'+str(j)])
            expression_m[i]+=multip_factor_a*2**j*(model_max.variables['m_'+str(i+1)+'_'+str(j)])
            expression_gamma3[i]+=multip_factor_gamma3*2**j*(model_max.variables['k3_'+str(i+1)+'_'+str(j)])
            expression_p[i]+=multip_factor_b*2**j*(model_max.variables['p_'+str(i+1)+'_'+str(j)])

            expression_t[i]+=gamma_ov*multip_factor_a*2**j*(model_max.variables['m_'+str(i+1)+'_'+str(j)])+ \
                             gamma_ov * multip_factor_b* 2 ** j * ( model_max.variables['p_' + str(i + 1) + '_' + str(j)])
            for k in range(binary_variables + 1):
                expression_t[i] +=multip_factor_b*multip_factor_a*2**(k+j) * ( model_max.variables['mt_' + str(i + 1) + '_' + str(j)+ '_' + str(k)])




    expression[n_step]=model_max.variables['v']*gamma_ov
    expression_m[n_step]=model_max.variables['v']*gamma_ov
    expression_p[n_step]=model_max.variables['v']*gamma_ov

    expression_gamma3[n_step]=model_max.variables['v']*gamma_ov
    expression_t[n_step]=model_max.variables['v']*gamma_ov**2
    for j in range(binary_variables + 1):
        expression[n_step] += multip_factor * 2 ** j * (model_max.variables['k_' + str(n_step+1) + '_' + str(j)])
        expression_m[n_step] += multip_factor_a * 2 ** j * (model_max.variables['m_' + str(n_step+1) + '_' + str(j)])

        expression_gamma3[n_step] += multip_factor_gamma3 * 2 ** j * (model_max.variables['k3_' + str(n_step+1) + '_' + str(j)])
        expression_p[n_step] += multip_factor_b * 2 ** j * (model_max.variables['p_' + str(n_step+1) + '_' + str(j)])

        expression_t[n_step] += gamma_ov * multip_factor_a * 2 ** j * (
        model_max.variables['m_' + str(n_step + 1) + '_' + str(j)]) + \
                           gamma_ov * multip_factor_b * 2 ** j * (
                           model_max.variables['p_' + str(n_step + 1) + '_' + str(j)])
        for k in range(binary_variables + 1):
            expression_t[n_step] += multip_factor_b * multip_factor_a * 2 ** (k + j) * (
                #todo error was here you need to solve this make it automatic
            model_max.variables['mt_' + str(n_step + 1) + '_' + str(j)+ '_' + str(k)])




    #todo to UPDATE
    c_rate_1=Constraint((model_max.variables['z1']-expression[0])-expression[n_step],lb=0, ub=0) #k1b
    c_rate_2=Constraint(expression[1]-expression_m[1]-expression_m[n_step],lb=0, ub=0)    #k2b
    c_rate_3=Constraint(model_max.variables['z3']-expression_gamma3[2]-expression_gamma3[n_step],lb=0, ub=0) #k3b
    c_rate_4=Constraint((expression_gamma3[3]-expression_p[3])*P1-expression_p[n_step],lb=0, ub=0) #k4b
    c_rate_5=Constraint((expression_t[4]-model_max.variables['z5']*gamma_ov)*P2-model_max.variables['v']*gamma_ov,lb=0, ub=0) #k5b
    c_rate_6=Constraint(S1*(model_max.variables['z6']-expression[5])-model_max.variables['v'],lb=0, ub=0) #k1f
    c_rate_7=Constraint(S2*(expression[6]-expression_m[6])-expression[n_step],lb=0, ub=0) #k2f
    c_rate_8=Constraint((model_max.variables['z8']-expression_gamma3[7])-model_max.variables['v'],lb=0, ub=0) #k3f
    c_rate_9=Constraint((expression_gamma3[8]-expression_p[8])-expression_gamma3[n_step],lb=0, ub=0) #k4f
    c_rate_10=Constraint((expression_t[9]-model_max.variables['z10']*gamma_ov)-expression_t[n_step],lb=0, ub=0) #k4f

    model_max.add([c_rate_1,c_rate_2,c_rate_3,c_rate_4,c_rate_5,c_rate_6,c_rate_7,c_rate_8,c_rate_9,c_rate_10])
    obj = Objective(model_max.variables['v'], direction='max')
    model_max.objective = obj
    model_max.configuration.tolerances.feasibility = 1e-9
    model_max.configuration.tolerances.integrality = 1e-9
    model_max.configuration.tolerances.optimality = 1e-9


    status = model_max.optimize()
    print(status)

    "Post-Processing"
    dict_values = dict()
    for key in model_max.variables.keys():
        dict_values[key] = (model_max.variables[key].primal)
        #print(key, (model_max.variables[key].primal))

    # add also these values to have the scaled version
    dict_values['A'] = S1
    dict_values['B'] = S2
    dict_values['Q'] = P1
    dict_values['P'] = P2

    # dict_values['gamma_1'] = gamma_1
    # dict_values['gamma_2'] = gamma_2
    #dict_values['gamma_4'] = gamma_ov / (a.primal*gamma_3.primal)

    dict_values['gamma_5'] = gamma_ov / (a.primal*b.primal)
    dict_values['gamma_4'] = b.primal / gamma_3.primal
    dict_values['gamma_2'] = a.primal / gamma_1.primal
    dict_values['gamma_ov'] = gamma_ov

    dict_values['k1f'] = model_max.variables['z6'].primal /model_max.variables['E'].primal
    dict_values['k2f'] = model_max.variables['z7'].primal /model_max.variables['EA'].primal
    dict_values['k3f'] = model_max.variables['z8'].primal /model_max.variables['EAB'].primal
    dict_values['k4f'] = model_max.variables['z9'].primal /model_max.variables['EPQ'].primal
    dict_values['k5f'] = model_max.variables['z10'].primal /model_max.variables['EP'].primal



    dict_values['k1b'] = model_max.variables['z1'].primal /model_max.variables['EA'].primal
    dict_values['k2b'] = model_max.variables['z2'].primal /model_max.variables['EAB'].primal
    dict_values['k3b'] = model_max.variables['z3'].primal /model_max.variables['EPQ'].primal
    dict_values['k4b'] = model_max.variables['z4'].primal /model_max.variables['EP'].primal
    dict_values['k5b'] = model_max.variables['z5'].primal /model_max.variables['E'].primal



    #todo variability analysis
    if variability_analysis == True:
        df_var_analysis = pd.DataFrame()
        df_var_analysis= variability_analysis_linearized(model_max)
        # print(df_var_analysis)
        #dict_values['unique'] = a

        # return (np.exp(model_max.variables.ln_v.primal)),(dict_values),df_var_analysis
    # for fitness you are directly returning the v value not log format
    else:
        df_var_analysis=pd.DataFrame()

    return ((model_max.variables.v.primal)), (dict_values),(df_var_analysis)
