
from optlang.cplex_interface import Model, Variable, Constraint, Objective
import numpy as np
from open.optim.variability_analysis import *
import math
import pandas as pd

#TODO THIS WHOLE FUNCTION IS TO BE REWRITTEN IN A PROPER REPRODUCIBLE MANNER by 02/21!!!!!!!
'''random order mechanism where you fix splitting ratio alpha
alpha is a parameter here'''
def milp_problem_4step_biuni_random(gamma_ov, q=2.0, S=1.0, P=1.0, alpha=0.8,variability_analysis=False):
    "ideally this big M is e_TOT but here as it is scaled it is 1 E_tot=1"
    big_M=1.0
    'irreversible number of steps'
    n_step=12
    model_max = Model(name='Max model')
    "E+S-->EP-->E+P"
    #N=math.ceil((1-gamma_ov)/0.005)
    #N=2000

    N=math.ceil((1-gamma_ov)/0.09)

    binary_variables=math.ceil(math.log2(N))
    print(binary_variables)
    A=S
    B = P/(A*q*gamma_ov)
    epsilon=0.0
    model_max.add(Variable('E' ,lb=epsilon, ub=1))
    model_max.add(Variable('EA' ,lb=epsilon, ub=1))
    model_max.add(Variable('EAB' ,lb=epsilon, ub=1))
    model_max.add(Variable('EP' ,lb=epsilon, ub=1))
    #this enzyme state is for the random case
    model_max.add(Variable('EB' ,lb=epsilon, ub=1))


    "binary variables  ceil(log_2_N)"
    #although it is to be seen if you should add the same number of binary variables for different approximated variables
    for i in range(binary_variables + 1):
        'for gamma1'
        model_max.add(Variable('d' + str(i), type='binary'))
        'for a=gamma1*gamma2'
        model_max.add(Variable('d_a' + str(i), type='binary'))
        'for gamma3'
        model_max.add(Variable('d_2' + str(i), type='binary'))

        'for random ordered mechanism'
        'for gamma5'
        model_max.add(Variable('d_5' + str(i), type='binary'))

        'for gamma3a multiplication'
        'you also need to add other constraints for this t based on d_a and d_2'


    for j in range(binary_variables + 1):
        for k in range(binary_variables + 1):
            model_max.add(Variable('t_' + str(j) + '_' +str(k), type='binary'))

            model_max.add(Constraint(
            model_max.variables['t_' + str(j) + '_' +str(k)] - model_max.variables['d_2' + str(j)] ,
            ub=0.0))
            model_max.add(Constraint(
            model_max.variables['t_' + str(j) + '_' +str(k)] - model_max.variables['d_a' + str(k)] ,
            ub=0.0))
            model_max.add(Constraint(
            model_max.variables['t_' + str(j) + '_' +str(k)] - model_max.variables['d_a' + str(k)] - model_max.variables['d_2' + str(j)] +1 ,
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

            # k3 for gamma 3 d_2
            model_max.add(Variable('k3_' + str(i + 1) + '_' + str(j), lb=0))

            # k5 for gamma 5 d_5
            model_max.add(Variable('k5_' + str(i + 1) + '_' + str(j), lb=0))

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
                model_max.variables['k3_' + str(i + 1) + '_' + str(j)] - big_M * model_max.variables['d_2' + str(j)],
                ub=0.0))
            model_max.add(Constraint(
                model_max.variables['k3_' + str(i + 1) + '_' + str(j)] - model_max.variables['z' + str(i + 1)], ub=0.0))
            model_max.add(Constraint(model_max.variables['k3_' + str(i + 1) + '_' + str(j)] - (
                        model_max.variables['z' + str(i + 1)] + big_M * model_max.variables['d_2' + str(j)]), lb=-1*big_M))

            'for random gamma 5 constraints '
            model_max.add(Constraint(
                model_max.variables['k5_' + str(i + 1) + '_' + str(j)] - big_M * model_max.variables['d_5' + str(j)],
                ub=0.0))
            model_max.add(Constraint(
                model_max.variables['k5_' + str(i + 1) + '_' + str(j)] - model_max.variables['z' + str(i + 1)], ub=0.0))
            model_max.add(Constraint(model_max.variables['k5_' + str(i + 1) + '_' + str(j)] - (
                        model_max.variables['z' + str(i + 1)] + big_M * model_max.variables['d_5' + str(j)]), lb=-1*big_M))
            #and three for multiplication variable
           # model_max.add(Constraint(
            #    model_max.variables['mt_' + str(i + 1) + '_' + str(j)] - big_M * model_max.variables['t' + str(j)],
             #   ub=0.0))
            #model_max.add(Constraint(
             #   model_max.variables['mt_' + str(i + 1) + '_' + str(j)] - model_max.variables['z' + str(i + 1)], ub=0.0))
            #model_max.add(Constraint(model_max.variables['mt_' + str(i + 1) + '_' + str(j)] - (
             #           model_max.variables['z' + str(i + 1)] + big_M * model_max.variables['t' + str(j)]), lb=-1*big_M))



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
                model_max.variables['k3_' + str(n_step+1) + '_' + str(j)] - big_M * model_max.variables['d_2' + str(j)],
                ub=0.0))
        model_max.add(
            Constraint(model_max.variables['k3_' + str(n_step+1) + '_' + str(j)] - model_max.variables['v'],
                       ub=0.0))
        model_max.add(Constraint(model_max.variables['k3_' + str(n_step+1) + '_' + str(j)] - (
                model_max.variables['v'] + big_M * model_max.variables['d_2' + str(j)]), lb=-1*big_M))

        "for k5 random extra"
        model_max.add(Variable('k5_' + str(n_step+1) + '_' + str(j), lb=0))
        model_max.add(
            Constraint(
                model_max.variables['k5_' + str(n_step+1) + '_' + str(j)] - big_M * model_max.variables['d_5' + str(j)],
                ub=0.0))
        model_max.add(
            Constraint(model_max.variables['k5_' + str(n_step+1) + '_' + str(j)] - model_max.variables['v'],
                       ub=0.0))
        model_max.add(Constraint(model_max.variables['k5_' + str(n_step+1) + '_' + str(j)] - (
                model_max.variables['v'] + big_M * model_max.variables['d_5' + str(j)]), lb=-1*big_M))


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




    gamma_1= Variable('gamma_1', lb=gamma_ov, ub=1)
    model_max.add(gamma_1)

    #this is the multiplication variable
    a= Variable('a', lb=gamma_ov, ub=1)
    model_max.add(a)

    gamma_3= Variable('gamma_3', lb=gamma_ov, ub=1)
    model_max.add(gamma_3)

    gamma_5 = Variable('gamma_5', lb=gamma_ov, ub=1)
    model_max.add(gamma_5)

    # gamma_3_a= Variable('gamma_3_a', lb=gamma_ov, ub=1)
    # model_max.add(gamma_3_a)


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
        expression_gamma3_+=2**i*model_max.variables['d_2' + str(i)]*multip_factor_gamma3

    'for random gamma 5'
    expression_gamma5_=0+model_max.variables['gamma_5'].lb
    multip_factor_gamma5=(model_max.variables['gamma_5'].ub-model_max.variables['gamma_5'].lb)/N
    for i in range(binary_variables + 1):
        expression_gamma5_+=2**i*model_max.variables['d_5' + str(i)]*multip_factor_gamma5

    #
    # expression_gamma3_a=0+model_max.variables['gamma_3_a'].lb**2
    # for j in range(binary_variables + 1):
    #     expression_gamma3_a += gamma_ov * multip_factor_a * 2 ** j * (model_max.variables['d_a' + str(j)]) + \
    #                        gamma_ov * multip_factor_gamma3 * 2 ** j * (model_max.variables['d_2' + str(j)])
    #     for k in range(binary_variables + 1):
    #         expression_gamma3_a += multip_factor_gamma3 * multip_factor_a * 2 ** (k + j) * (model_max.variables['t_' + str(j) + '_' +str(k)])
    #


    c_gamma=Constraint(gamma_1-expression_gamma,lb=0,ub=0)

    c_a=Constraint(a-expression_a,lb=0,ub=0)

    c_gamma_5=Constraint(gamma_5-expression_gamma5_,lb=0,ub=0)


    c_a_gamma=Constraint(a-gamma_1, ub=0)

    'for random'
    c_a_gamma_5=Constraint(a-gamma_5, ub=0)


    c_gamma3=Constraint(gamma_3-expression_gamma3_,lb=0,ub=0)

    #if this is required or not

    #c_gamma3_a=Constraint(gamma_3_a-expression_gamma3_a,lb=0,ub=0)

    #c_gamma3_a_1=Constraint(gamma_3_a-gamma_3,ub=0)

    #c_gamma3_a_2=Constraint(gamma_3_a-a,ub=0)


    model_max.add([c_mass_cons,c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_10, c_11,c_12, c_gamma,c_a,c_a_gamma,c_a_gamma_5,c_gamma3,c_gamma_5]) #,c_gamma3_a,c_gamma3_a_1,c_gamma3_a_2])
    #write expressions in for loop

    expression={}
    expression_m={}
    expression_gamma3={}
    expression_t={}
    expression_gamma5={}

    for i in range(n_step):
        'this is the lower bound lower bound  + 2** terms'
        expression[i]=model_max.variables['z'+str(i+1)]*gamma_ov
        expression_m[i]=model_max.variables['z'+str(i+1)]*gamma_ov
        expression_gamma3[i]=model_max.variables['z'+str(i+1)]*gamma_ov

        'for random'
        expression_gamma5[i]=model_max.variables['z'+str(i+1)]*gamma_ov


        expression_t[i]=model_max.variables['z'+str(i+1)]*gamma_ov**2


        for j in range(binary_variables + 1):
            expression[i]+=multip_factor*2**j*(model_max.variables['k_'+str(i+1)+'_'+str(j)])
            expression_m[i]+=multip_factor_a*2**j*(model_max.variables['m_'+str(i+1)+'_'+str(j)])
            expression_gamma3[i]+=multip_factor_gamma3*2**j*(model_max.variables['k3_'+str(i+1)+'_'+str(j)])
            'for random'
            expression_gamma5[i]+=multip_factor_gamma5*2**j*(model_max.variables['k5_'+str(i+1)+'_'+str(j)])

            'this expression t is more complicated as it is also multiplication of 2 integers'
            expression_t[i]+=gamma_ov*multip_factor_a*2**j*(model_max.variables['m_'+str(i+1)+'_'+str(j)])+ \
                             gamma_ov * multip_factor_gamma3 * 2 ** j * ( model_max.variables['k3_' + str(i + 1) + '_' + str(j)])
            for k in range(binary_variables + 1):
                expression_t[i] +=multip_factor_gamma3*multip_factor_a*2**(k+j) * ( model_max.variables['mt_' + str(i + 1) + '_' + str(j)+ '_' + str(k)])




    expression[n_step]=model_max.variables['v']*gamma_ov
    expression_m[n_step]=model_max.variables['v']*gamma_ov
    expression_gamma3[n_step]=model_max.variables['v']*gamma_ov
    'for random'
    expression_gamma5[n_step]=model_max.variables['v']*gamma_ov

    expression_t[n_step]=model_max.variables['v']*gamma_ov**2
    for j in range(binary_variables + 1):
        expression[n_step] += multip_factor * 2 ** j * (model_max.variables['k_' + str(n_step+1) + '_' + str(j)])
        expression_m[n_step] += multip_factor_a * 2 ** j * (model_max.variables['m_' + str(n_step+1) + '_' + str(j)])
        expression_gamma3[n_step] += multip_factor_gamma3 * 2 ** j * (model_max.variables['k3_' + str(n_step+1) + '_' + str(j)])
        'for random'
        expression_gamma5[n_step] += multip_factor_gamma5 * 2 ** j * (model_max.variables['k5_' + str(n_step+1) + '_' + str(j)])

        expression_t[n_step] += gamma_ov * multip_factor_a * 2 ** j * (
        model_max.variables['m_' + str(n_step + 1) + '_' + str(j)]) + \
                           gamma_ov * multip_factor_gamma3 * 2 ** j * (
                           model_max.variables['k3_' + str(n_step + 1) + '_' + str(j)])
        for k in range(binary_variables + 1):
            expression_t[n_step] += multip_factor_gamma3 * multip_factor_a * 2 ** (k + j) * (
                #todo error was here you need to solve this make it automatic
            model_max.variables['mt_' + str(n_step + 1) + '_' + str(j)+ '_' + str(k)])

    #todo to UPDATE
    c_rate_1=Constraint((model_max.variables['z1']-expression[0])-expression[n_step]*alpha,lb=0, ub=0) #k1b
    c_rate_2=Constraint(expression[1]-expression_m[1]-expression_m[n_step]*alpha,lb=0, ub=0) #k2b
    c_rate_3=Constraint(model_max.variables['z3']-expression_gamma3[2]-expression_gamma3[n_step],lb=0, ub=0) #k3b
    c_rate_4=Constraint(model_max.variables['z4']-expression_gamma3[3]-model_max.variables['v'],lb=0, ub=0)  #k3f
    c_rate_5=Constraint(expression_t[4]-model_max.variables['z5']*gamma_ov-expression_t[n_step],lb=0, ub=0) #k4f
    c_rate_6=Constraint(A*(model_max.variables['z6']-expression[5])-model_max.variables['v']*alpha,lb=0, ub=0) #k1f
    c_rate_7=Constraint(B*(expression[6]-expression_m[6])-expression[n_step]*alpha,lb=0, ub=0) #k2f
    c_rate_8=Constraint(P*(expression_t[7]-model_max.variables['z8']*gamma_ov)-model_max.variables['v']*gamma_ov,lb=0, ub=0) #k4b

    c_rate_9=Constraint(B*(model_max.variables['z9']-expression_gamma5[8])-model_max.variables['v']*(1-alpha),lb=0, ub=0)  #k5f
    c_rate_10=Constraint((model_max.variables['z10']-expression_gamma5[9])-expression_gamma5[n_step]*(1-alpha),lb=0, ub=0)#k5b
    c_rate_11=Constraint(A*(expression_gamma5[10]-expression_m[10])-expression_gamma5[n_step]*(1-alpha),lb=0, ub=0) #k6f
    c_rate_12=Constraint(expression_gamma5[11]-expression_m[11]-expression_m[n_step]*(1-alpha),lb=0, ub=0) #k6b



    model_max.add([c_rate_1,c_rate_2,c_rate_3,c_rate_4,c_rate_5,c_rate_6,c_rate_7,c_rate_8, c_rate_9, c_rate_10, c_rate_11, c_rate_12])
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
    dict_values['A'] = A
    dict_values['B'] = B
    dict_values['P'] = P
    dict_values['alpha'] = alpha

    # dict_values['gamma_1'] = gamma_1
    # dict_values['gamma_2'] = gamma_2
    dict_values['gamma_4'] = gamma_ov / (a.primal*gamma_3.primal)
    dict_values['gamma_2'] = a.primal / gamma_1.primal
    dict_values['gamma_6'] = a.primal / gamma_5.primal

    dict_values['gamma_ov'] = gamma_ov

    dict_values['k1f'] = model_max.variables['z6'].primal /model_max.variables['E'].primal
    dict_values['k2f'] = model_max.variables['z7'].primal /model_max.variables['EA'].primal
    dict_values['k3f'] = model_max.variables['z4'].primal /model_max.variables['EAB'].primal
    dict_values['k4f'] = model_max.variables['z5'].primal /model_max.variables['EP'].primal
    'for random'
    dict_values['k5f'] = model_max.variables['z9'].primal /model_max.variables['E'].primal
    dict_values['k6f'] = model_max.variables['z11'].primal /model_max.variables['EB'].primal

    dict_values['k1b'] = model_max.variables['z1'].primal /model_max.variables['EA'].primal
    dict_values['k2b'] = model_max.variables['z2'].primal /model_max.variables['EAB'].primal
    dict_values['k3b'] = model_max.variables['z3'].primal /model_max.variables['EP'].primal
    dict_values['k4b'] = model_max.variables['z8'].primal /model_max.variables['E'].primal
    'for random'
    dict_values['k5b'] = model_max.variables['z10'].primal /model_max.variables['EB'].primal
    dict_values['k6b'] = model_max.variables['z12'].primal /model_max.variables['EAB'].primal



    if variability_analysis == True:
        df_var_analysis = pd.DataFrame()
        df_var_analysis= variability_analysis_general(model_max)
        # print(df_var_analysis)
        #dict_values['unique'] = a
    else:
        df_var_analysis= []

        # return (np.exp(model_max.variables.ln_v.primal)),(dict_values),df_var_analysis
    # for fitness you are directly returning the v value not log format

    return ((model_max.variables.v.primal)), (dict_values), (df_var_analysis)



'random order mechanism where alpha is a variable of the MILP'
def milp_problem_4step_biuni_random_split_ratio(gamma_ov, q=2.0, S=1.0, P=1.0,variability_analysis=False):
    "ideally this big M is e_TOT but here as it is scaled it is 1 E_tot=1"
    big_M=1.0
    'irreversible number of steps'
    n_step=12
    model_max = Model(name='Max model')
    "E+S-->EP-->E+P"
    #N=math.ceil((1-gamma_ov)/0.005)
    #N=2000
    A=S
    B = P/(A*q*gamma_ov)
    #todo check this
    # if (gamma_ov>1):
    #     gamma_ov=1/gamma_ov
    #     q=1/q
    #     B=gamma_ov*P*q/A
    'bin resolution of 0.001'
    if (gamma_ov>=0.1 and gamma_ov<0.9):
         bin_resol=0.001
    elif gamma_ov>=0.9:
        bin_resol=0.0001
    else:
        bin_resol=0.01
    N=math.ceil((1-gamma_ov)/bin_resol)
    #N=math.ceil((1-gamma_ov)/0.1)

    binary_variables=math.ceil(math.log2(N))
    print(binary_variables)
    # A=S
    # B = P/(A*q*gamma_ov)
    epsilon=1e-8
    model_max.add(Variable('E' ,lb=epsilon, ub=1))
    model_max.add(Variable('EA' ,lb=epsilon, ub=1))
    model_max.add(Variable('EAB' ,lb=epsilon, ub=1))
    model_max.add(Variable('EP' ,lb=epsilon, ub=1))
    #this enzyme state is for the random case
    model_max.add(Variable('EB' ,lb=epsilon, ub=1))



    "binary variables  ceil(log_2_N)"
    #although it is to be seen if you should add the same number of binary variables for different approximated variables
    for i in range(binary_variables + 1):
        'for gamma1'
        model_max.add(Variable('d' + str(i), type='binary'))
        'for a=gamma1*gamma2'
        model_max.add(Variable('d_a' + str(i), type='binary'))
        'for gamma3'
        model_max.add(Variable('d_2' + str(i), type='binary'))

        'for random ordered mechanism'
        'for gamma5'
        model_max.add(Variable('d_5' + str(i), type='binary'))

        'for gamma3a multiplication'
        'you also need to add other constraints for this t based on d_a and d_2'


    for j in range(binary_variables + 1):
        for k in range(binary_variables + 1):
            model_max.add(Variable('t_' + str(j) + '_' +str(k), type='binary'))

            model_max.add(Constraint(
            model_max.variables['t_' + str(j) + '_' +str(k)] - model_max.variables['d_2' + str(j)] ,
            ub=0.0))
            model_max.add(Constraint(
            model_max.variables['t_' + str(j) + '_' +str(k)] - model_max.variables['d_a' + str(k)] ,
            ub=0.0))
            model_max.add(Constraint(
            model_max.variables['t_' + str(j) + '_' +str(k)] - model_max.variables['d_a' + str(k)] - model_max.variables['d_2' + str(j)] +1 ,
            lb=0.0))



    # for z1 to z6
    for i in range(n_step):
        model_max.add(Variable('z' + str(i + 1), lb=epsilon))
        for j in range(binary_variables + 1):
            # todo you check if it can be 0 or not
            #Lets say k is fir gamma d
            model_max.add(Variable('k_' + str(i + 1) + '_' + str(j), lb=0))

            #lets m is for the multiplocation variable gamma1*gamma2 d_a
            model_max.add(Variable('m_' + str(i + 1) + '_' + str(j), lb=0))

            # k3 for gamma 3 d_2
            model_max.add(Variable('k3_' + str(i + 1) + '_' + str(j), lb=0))

            # k5 for gamma 5 d_5
            model_max.add(Variable('k5_' + str(i + 1) + '_' + str(j), lb=0))

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
                model_max.variables['k3_' + str(i + 1) + '_' + str(j)] - big_M * model_max.variables['d_2' + str(j)],
                ub=0.0))
            model_max.add(Constraint(
                model_max.variables['k3_' + str(i + 1) + '_' + str(j)] - model_max.variables['z' + str(i + 1)], ub=0.0))
            model_max.add(Constraint(model_max.variables['k3_' + str(i + 1) + '_' + str(j)] - (
                        model_max.variables['z' + str(i + 1)] + big_M * model_max.variables['d_2' + str(j)]), lb=-1*big_M))

            'for random gamma 5 constraints '
            model_max.add(Constraint(
                model_max.variables['k5_' + str(i + 1) + '_' + str(j)] - big_M * model_max.variables['d_5' + str(j)],
                ub=0.0))
            model_max.add(Constraint(
                model_max.variables['k5_' + str(i + 1) + '_' + str(j)] - model_max.variables['z' + str(i + 1)], ub=0.0))
            model_max.add(Constraint(model_max.variables['k5_' + str(i + 1) + '_' + str(j)] - (
                        model_max.variables['z' + str(i + 1)] + big_M * model_max.variables['d_5' + str(j)]), lb=-1*big_M))
            #and three for multiplication variable
           # model_max.add(Constraint(
            #    model_max.variables['mt_' + str(i + 1) + '_' + str(j)] - big_M * model_max.variables['t' + str(j)],
             #   ub=0.0))
            #model_max.add(Constraint(
             #   model_max.variables['mt_' + str(i + 1) + '_' + str(j)] - model_max.variables['z' + str(i + 1)], ub=0.0))
            #model_max.add(Constraint(model_max.variables['mt_' + str(i + 1) + '_' + str(j)] - (
             #           model_max.variables['z' + str(i + 1)] + big_M * model_max.variables['t' + str(j)]), lb=-1*big_M))



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
    'here you denote v_lower =v*(1-alpha)'
    'and v_upper=v*alpha'
    model_max.add(Variable('v_upper' ,lb=0))
    model_max.add(Variable('v_lower' ,lb=0))

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
                model_max.variables['v_upper'] + big_M * model_max.variables['d' + str(j)]), lb=-1*big_M))


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

        "for m_upper"
        model_max.add(Variable('m_upper' + str(n_step+1) + '_' + str(j), lb=0))
        model_max.add(
            Constraint(
                model_max.variables['m_upper' + str(n_step+1) + '_' + str(j)] - big_M * model_max.variables['d_a' + str(j)],
                ub=0.0))
        model_max.add(
            Constraint(model_max.variables['m_upper' + str(n_step+1) + '_' + str(j)] - model_max.variables['v_upper'],
                       ub=0.0))
        model_max.add(Constraint(model_max.variables['m_upper' + str(n_step+1) + '_' + str(j)] - (
                model_max.variables['v_upper'] + big_M * model_max.variables['d_a' + str(j)]), lb=-1*big_M))


        "for m_lower"
        model_max.add(Variable('m_lower' + str(n_step+1) + '_' + str(j), lb=0))
        model_max.add(
            Constraint(
                model_max.variables['m_lower' + str(n_step+1) + '_' + str(j)] - big_M * model_max.variables['d_a' + str(j)],
                ub=0.0))
        model_max.add(
            Constraint(model_max.variables['m_lower' + str(n_step+1) + '_' + str(j)] - model_max.variables['v_lower'],
                       ub=0.0))
        model_max.add(Constraint(model_max.variables['m_lower' + str(n_step+1) + '_' + str(j)] - (
                model_max.variables['v_lower'] + big_M * model_max.variables['d_a' + str(j)]), lb=-1*big_M))


        #todo for this you dont need but check
        "for k3"
        model_max.add(Variable('k3_' + str(n_step+1) + '_' + str(j), lb=0))
        model_max.add(
            Constraint(
                model_max.variables['k3_' + str(n_step+1) + '_' + str(j)] - big_M * model_max.variables['d_2' + str(j)],
                ub=0.0))
        model_max.add(
            Constraint(model_max.variables['k3_' + str(n_step+1) + '_' + str(j)] - model_max.variables['v'],
                       ub=0.0))
        model_max.add(Constraint(model_max.variables['k3_' + str(n_step+1) + '_' + str(j)] - (
                model_max.variables['v'] + big_M * model_max.variables['d_2' + str(j)]), lb=-1*big_M))

        "for k5 random extra"
        model_max.add(Variable('k5_' + str(n_step+1) + '_' + str(j), lb=0))
        model_max.add(
            Constraint(
                model_max.variables['k5_' + str(n_step+1) + '_' + str(j)] - big_M * model_max.variables['d_5' + str(j)],
                ub=0.0))
        model_max.add(
            Constraint(model_max.variables['k5_' + str(n_step+1) + '_' + str(j)] - model_max.variables['v'],
                       ub=0.0))
        model_max.add(Constraint(model_max.variables['k5_' + str(n_step+1) + '_' + str(j)] - (
                model_max.variables['v_lower'] + big_M * model_max.variables['d_5' + str(j)]), lb=-1*big_M))


    for j in range(binary_variables + 1):
        for k in range(binary_variables + 1):
            #todo this you shouldnt need but verify
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


    #indicate upper and lower bounds
    gamma_1= Variable('gamma_1', lb=gamma_ov, ub=1)
    model_max.add(gamma_1)

    #this is the multiplication variable
    a= Variable('a', lb=gamma_ov, ub=1)
    model_max.add(a)

    gamma_3= Variable('gamma_3', lb=gamma_ov, ub=1)
    model_max.add(gamma_3)

    gamma_5 = Variable('gamma_5', lb=gamma_ov, ub=1)
    model_max.add(gamma_5)

    # gamma_3_a= Variable('gamma_3_a', lb=gamma_ov, ub=1)
    # model_max.add(gamma_3_a)


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
        expression_gamma3_+=2**i*model_max.variables['d_2' + str(i)]*multip_factor_gamma3

    'for random gamma 5'
    expression_gamma5_=0+model_max.variables['gamma_5'].lb
    multip_factor_gamma5=(model_max.variables['gamma_5'].ub-model_max.variables['gamma_5'].lb)/N
    for i in range(binary_variables + 1):
        expression_gamma5_+=2**i*model_max.variables['d_5' + str(i)]*multip_factor_gamma5
    #
    #
    # expression_gamma3_a=0+model_max.variables['gamma_3_a'].lb**2
    # for j in range(binary_variables + 1):
    #     expression_gamma3_a += gamma_ov * multip_factor_a * 2 ** j * (model_max.variables['d_a' + str(j)]) + \
    #                        gamma_ov * multip_factor_gamma3 * 2 ** j * (model_max.variables['d_2' + str(j)])
    #     for k in range(binary_variables + 1):
    #         expression_gamma3_a += multip_factor_gamma3 * multip_factor_a * 2 ** (k + j) * (model_max.variables['t_' + str(j) + '_' +str(k)])



    c_gamma=Constraint(gamma_1-expression_gamma,lb=0,ub=0)

    c_a=Constraint(a-expression_a,lb=0,ub=0)

    'stating that gamma_1*gamma_2<gamma_1 & gamma_1*gamma_2<gamma_5'
    c_a_gamma=Constraint(a-gamma_1, ub=0)

    'for random'
    c_a_gamma_5=Constraint(a-gamma_5, ub=0)

    c_gamma_5=Constraint(gamma_5-expression_gamma5_,lb=0,ub=0)

    c_gamma3=Constraint(gamma_3-expression_gamma3_,lb=0,ub=0)

    #if this is required or not

    #c_gamma3_a=Constraint(gamma_3_a-expression_gamma3_a,lb=0,ub=0)

    #c_gamma3_a_1=Constraint(gamma_3_a-gamma_3,ub=0)

    #c_gamma3_a_2=Constraint(gamma_3_a-a,ub=0)


    model_max.add([c_mass_cons,c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_10, c_11,c_12, c_gamma,c_a,c_a_gamma,c_a_gamma_5,c_gamma3,c_gamma_5]) #,c_gamma3_a,c_gamma3_a_1,c_gamma3_a_2])
    #write expressions in for loop

    expression={}
    expression_m={}
    expression_gamma3={}
    expression_t={}
    expression_gamma5={}

    for i in range(n_step):
        'this is the lower bound lower bound  + 2** terms'
        expression[i]=model_max.variables['z'+str(i+1)]*gamma_ov
        expression_m[i]=model_max.variables['z'+str(i+1)]*gamma_ov
        expression_gamma3[i]=model_max.variables['z'+str(i+1)]*gamma_ov

        'for random'
        expression_gamma5[i]=model_max.variables['z'+str(i+1)]*gamma_ov


        expression_t[i]=model_max.variables['z'+str(i+1)]*gamma_ov**2


        for j in range(binary_variables + 1):
            expression[i]+=multip_factor*2**j*(model_max.variables['k_'+str(i+1)+'_'+str(j)])
            expression_m[i]+=multip_factor_a*2**j*(model_max.variables['m_'+str(i+1)+'_'+str(j)])
            expression_gamma3[i]+=multip_factor_gamma3*2**j*(model_max.variables['k3_'+str(i+1)+'_'+str(j)])
            'for random'
            expression_gamma5[i]+=multip_factor_gamma5*2**j*(model_max.variables['k5_'+str(i+1)+'_'+str(j)])

            'this expression t is more complicated as it is also multiplication of 2 integers'
            expression_t[i]+=gamma_ov*multip_factor_a*2**j*(model_max.variables['m_'+str(i+1)+'_'+str(j)])+ \
                             gamma_ov * multip_factor_gamma3 * 2 ** j * ( model_max.variables['k3_' + str(i + 1) + '_' + str(j)])
            for k in range(binary_variables + 1):
                expression_t[i] +=multip_factor_gamma3*multip_factor_a*2**(k+j) * ( model_max.variables['mt_' + str(i + 1) + '_' + str(j)+ '_' + str(k)])

    expression_m_lower = {}
    'for k'
    expression[n_step]=model_max.variables['v_upper']*gamma_ov
    expression_m[n_step]=model_max.variables['v_upper']*gamma_ov
    expression_m_lower[n_step]=model_max.variables['v_lower']*gamma_ov

    expression_gamma3[n_step]=model_max.variables['v']*gamma_ov
    'for random'
    expression_gamma5[n_step]=model_max.variables['v_lower']*gamma_ov

    expression_t[n_step]=model_max.variables['v']*gamma_ov**2
    for j in range(binary_variables + 1):
        expression[n_step] += multip_factor * 2 ** j * (model_max.variables['k_' + str(n_step+1) + '_' + str(j)])
        expression_m[n_step] += multip_factor_a * 2 ** j * (model_max.variables['m_upper' + str(n_step+1) + '_' + str(j)])

        expression_m_lower[n_step] += multip_factor_a * 2 ** j * (model_max.variables['m_lower' + str(n_step+1) + '_' + str(j)])

        expression_gamma3[n_step] += multip_factor_gamma3 * 2 ** j * (model_max.variables['k3_' + str(n_step+1) + '_' + str(j)])
        'for random'
        expression_gamma5[n_step] += multip_factor_gamma5 * 2 ** j * (model_max.variables['k5_' + str(n_step+1) + '_' + str(j)])

        expression_t[n_step] += gamma_ov * multip_factor_a * 2 ** j * (
        model_max.variables['m_' + str(n_step + 1) + '_' + str(j)]) + \
                           gamma_ov * multip_factor_gamma3 * 2 ** j * (
                           model_max.variables['k3_' + str(n_step + 1) + '_' + str(j)])
        for k in range(binary_variables + 1):
            expression_t[n_step] += multip_factor_gamma3 * multip_factor_a * 2 ** (k + j) * (
                #todo error was here you need to solve this make it automatic
            model_max.variables['mt_' + str(n_step + 1) + '_' + str(j)+ '_' + str(k)])

    #todo to UPDATE
    c_rate_1=Constraint((model_max.variables['z1']-expression[0])-expression[n_step],lb=0, ub=0) #k1b
    c_rate_2=Constraint(expression[1]-expression_m[1]-expression_m[n_step],lb=0, ub=0) #k2b
    c_rate_3=Constraint(model_max.variables['z3']-expression_gamma3[2]-expression_gamma3[n_step],lb=0, ub=0) #k3b
    c_rate_4=Constraint(model_max.variables['z4']-expression_gamma3[3]-model_max.variables['v'],lb=0, ub=0)  #k3f
    c_rate_5=Constraint(expression_t[4]-model_max.variables['z5']*gamma_ov-expression_t[n_step],lb=0, ub=0) #k4f
    c_rate_6=Constraint(A*(model_max.variables['z6']-expression[5])-model_max.variables['v_upper'],lb=0, ub=0) #k1f
    c_rate_7=Constraint(B*(expression[6]-expression_m[6])-expression[n_step],lb=0, ub=0) #k2f
    c_rate_8=Constraint(P*(expression_t[7]-model_max.variables['z8']*gamma_ov)-model_max.variables['v']*gamma_ov,lb=0, ub=0) #k4b

    c_rate_9=Constraint(B*(model_max.variables['z9']-expression_gamma5[8])-model_max.variables['v_lower'],lb=0, ub=0)  #k5f
    c_rate_10=Constraint((model_max.variables['z10']-expression_gamma5[9])-expression_gamma5[n_step],lb=0, ub=0)#k5b
    c_rate_11=Constraint(A*(expression_gamma5[10]-expression_m[10])-expression_gamma5[n_step],lb=0, ub=0) #k6f
    c_rate_12=Constraint(expression_gamma5[11]-expression_m[11]-expression_m_lower[n_step],lb=0, ub=0) #k6b

    c_split_ratio = Constraint(model_max.variables['v'] - (model_max.variables['v_upper']+model_max.variables['v_lower']),lb=0, ub=0)

    model_max.add([c_rate_1,c_rate_2,c_rate_3,c_rate_4,c_rate_5,c_rate_6,c_rate_7,c_rate_8, c_rate_9, c_rate_10, c_rate_11, c_rate_12,c_split_ratio])
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
    dict_values['A'] = A
    dict_values['B'] = B
    dict_values['P'] = P
    #dict_values['alpha'] = alpha

    # dict_values['gamma_1'] = gamma_1
    # dict_values['gamma_2'] = gamma_2
    dict_values['gamma_4'] = gamma_ov / (a.primal*gamma_3.primal)
    dict_values['gamma_2'] = a.primal / gamma_1.primal
    dict_values['gamma_6'] = a.primal / gamma_5.primal

    dict_values['gamma_ov'] = gamma_ov

    'do not forget to uncomment'
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

    dict_values['alpha'] = model_max.variables['v_upper'].primal /model_max.variables['v'].primal



    #todo variability analysis
    if variability_analysis == True:
        df_var_analysis = pd.DataFrame()
        df_var_analysis = variability_analysis_linearized(model_max)
        # print(df_var_analysis)
        #dict_values['unique'] = a

    else:
        df_var_analysis = []

        # return (np.exp(model_max.variables.ln_v.primal)),(dict_values),df_var_analysis
    # for fitness you are directly returning the v value not log format

    return ((model_max.variables.v.primal)), (dict_values),(df_var_analysis)




'FEASIBILITY PROBLEM'


'random order mechanism where alpha is a variable of the MILP'
def milp_problem_4step_biuni_random_split_ratio_feasibility(v_net_lb,v_net_ub,gamma_ov, q=2.0, S=1.0, P=1.0,variability_analysis=False):
    "ideally this big M is e_TOT but here as it is scaled it is 1 E_tot=1"
    big_M=1.0
    'irreversible number of steps'
    n_step=12
    model_max = Model(name='Max model')
    "E+S-->EP-->E+P"
    #N=math.ceil((1-gamma_ov)/0.005)
    #N=2000
    'bin resolution of 0.001'
    'bin resolution of 0.001'
    if (gamma_ov>=0.1 and gamma_ov<0.9):
         bin_resol=0.001
    elif gamma_ov>=0.9:
        bin_resol=0.0001
    else:
        bin_resol=0.01

    N=math.ceil((1-gamma_ov)/bin_resol)

    binary_variables=math.ceil(math.log2(N))
    print(binary_variables)
    A=S
    B = P/(A*q*gamma_ov)

    epsilon=0.0
    model_max.add(Variable('E' ,lb=epsilon, ub=1))
    model_max.add(Variable('EA' ,lb=epsilon, ub=1))
    model_max.add(Variable('EAB' ,lb=epsilon, ub=1))
    model_max.add(Variable('EP' ,lb=epsilon, ub=1))
    #this enzyme state is for the random case
    model_max.add(Variable('EB' ,lb=epsilon, ub=1))



    "binary variables  ceil(log_2_N)"
    #although it is to be seen if you should add the same number of binary variables for different approximated variables
    for i in range(binary_variables + 1):
        'for gamma1'
        model_max.add(Variable('d' + str(i), type='binary'))
        'for a=gamma1*gamma2'
        model_max.add(Variable('d_a' + str(i), type='binary'))
        'for gamma3'
        model_max.add(Variable('d_2' + str(i), type='binary'))

        'for random ordered mechanism'
        'for gamma5'
        model_max.add(Variable('d_5' + str(i), type='binary'))

        'for gamma3a multiplication'
        'you also need to add other constraints for this t based on d_a and d_2'


    for j in range(binary_variables + 1):
        for k in range(binary_variables + 1):
            model_max.add(Variable('t_' + str(j) + '_' +str(k), type='binary'))

            model_max.add(Constraint(
            model_max.variables['t_' + str(j) + '_' +str(k)] - model_max.variables['d_2' + str(j)] ,
            ub=0.0))
            model_max.add(Constraint(
            model_max.variables['t_' + str(j) + '_' +str(k)] - model_max.variables['d_a' + str(k)] ,
            ub=0.0))
            model_max.add(Constraint(
            model_max.variables['t_' + str(j) + '_' +str(k)] - model_max.variables['d_a' + str(k)] - model_max.variables['d_2' + str(j)] +1 ,
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

            # k3 for gamma 3 d_2
            model_max.add(Variable('k3_' + str(i + 1) + '_' + str(j), lb=0))

            # k5 for gamma 5 d_5
            model_max.add(Variable('k5_' + str(i + 1) + '_' + str(j), lb=0))

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
                model_max.variables['k3_' + str(i + 1) + '_' + str(j)] - big_M * model_max.variables['d_2' + str(j)],
                ub=0.0))
            model_max.add(Constraint(
                model_max.variables['k3_' + str(i + 1) + '_' + str(j)] - model_max.variables['z' + str(i + 1)], ub=0.0))
            model_max.add(Constraint(model_max.variables['k3_' + str(i + 1) + '_' + str(j)] - (
                        model_max.variables['z' + str(i + 1)] + big_M * model_max.variables['d_2' + str(j)]), lb=-1*big_M))

            'for random gamma 5 constraints '
            model_max.add(Constraint(
                model_max.variables['k5_' + str(i + 1) + '_' + str(j)] - big_M * model_max.variables['d_5' + str(j)],
                ub=0.0))
            model_max.add(Constraint(
                model_max.variables['k5_' + str(i + 1) + '_' + str(j)] - model_max.variables['z' + str(i + 1)], ub=0.0))
            model_max.add(Constraint(model_max.variables['k5_' + str(i + 1) + '_' + str(j)] - (
                        model_max.variables['z' + str(i + 1)] + big_M * model_max.variables['d_5' + str(j)]), lb=-1*big_M))
            #and three for multiplication variable
           # model_max.add(Constraint(
            #    model_max.variables['mt_' + str(i + 1) + '_' + str(j)] - big_M * model_max.variables['t' + str(j)],
             #   ub=0.0))
            #model_max.add(Constraint(
             #   model_max.variables['mt_' + str(i + 1) + '_' + str(j)] - model_max.variables['z' + str(i + 1)], ub=0.0))
            #model_max.add(Constraint(model_max.variables['mt_' + str(i + 1) + '_' + str(j)] - (
             #           model_max.variables['z' + str(i + 1)] + big_M * model_max.variables['t' + str(j)]), lb=-1*big_M))



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

    v=Variable('v', lb=v_net_lb,ub=v_net_ub)
    model_max.add(v)
    'here you denote v_lower =v*(1-alpha)'
    'and v_upper=v*alpha'
    model_max.add(Variable('v_upper' ,lb=0))
    model_max.add(Variable('v_lower' ,lb=0))

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
                model_max.variables['v_upper'] + big_M * model_max.variables['d' + str(j)]), lb=-1*big_M))


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

        "for m_upper"
        model_max.add(Variable('m_upper' + str(n_step+1) + '_' + str(j), lb=0))
        model_max.add(
            Constraint(
                model_max.variables['m_upper' + str(n_step+1) + '_' + str(j)] - big_M * model_max.variables['d_a' + str(j)],
                ub=0.0))
        model_max.add(
            Constraint(model_max.variables['m_upper' + str(n_step+1) + '_' + str(j)] - model_max.variables['v_upper'],
                       ub=0.0))
        model_max.add(Constraint(model_max.variables['m_upper' + str(n_step+1) + '_' + str(j)] - (
                model_max.variables['v_upper'] + big_M * model_max.variables['d_a' + str(j)]), lb=-1*big_M))


        "for m_lower"
        model_max.add(Variable('m_lower' + str(n_step+1) + '_' + str(j), lb=0))
        model_max.add(
            Constraint(
                model_max.variables['m_lower' + str(n_step+1) + '_' + str(j)] - big_M * model_max.variables['d_a' + str(j)],
                ub=0.0))
        model_max.add(
            Constraint(model_max.variables['m_lower' + str(n_step+1) + '_' + str(j)] - model_max.variables['v_lower'],
                       ub=0.0))
        model_max.add(Constraint(model_max.variables['m_lower' + str(n_step+1) + '_' + str(j)] - (
                model_max.variables['v_lower'] + big_M * model_max.variables['d_a' + str(j)]), lb=-1*big_M))


        #todo for this you dont need but check
        "for k3"
        model_max.add(Variable('k3_' + str(n_step+1) + '_' + str(j), lb=0))
        model_max.add(
            Constraint(
                model_max.variables['k3_' + str(n_step+1) + '_' + str(j)] - big_M * model_max.variables['d_2' + str(j)],
                ub=0.0))
        model_max.add(
            Constraint(model_max.variables['k3_' + str(n_step+1) + '_' + str(j)] - model_max.variables['v'],
                       ub=0.0))
        model_max.add(Constraint(model_max.variables['k3_' + str(n_step+1) + '_' + str(j)] - (
                model_max.variables['v'] + big_M * model_max.variables['d_2' + str(j)]), lb=-1*big_M))

        "for k5 random extra"
        model_max.add(Variable('k5_' + str(n_step+1) + '_' + str(j), lb=0))
        model_max.add(
            Constraint(
                model_max.variables['k5_' + str(n_step+1) + '_' + str(j)] - big_M * model_max.variables['d_5' + str(j)],
                ub=0.0))
        model_max.add(
            Constraint(model_max.variables['k5_' + str(n_step+1) + '_' + str(j)] - model_max.variables['v'],
                       ub=0.0))
        model_max.add(Constraint(model_max.variables['k5_' + str(n_step+1) + '_' + str(j)] - (
                model_max.variables['v_lower'] + big_M * model_max.variables['d_5' + str(j)]), lb=-1*big_M))


    for j in range(binary_variables + 1):
        for k in range(binary_variables + 1):
            #todo this you shouldnt need but verify
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


    #indicate upper and lower bounds
    gamma_1= Variable('gamma_1', lb=gamma_ov, ub=1)
    model_max.add(gamma_1)

    #this is the multiplication variable
    a= Variable('a', lb=gamma_ov, ub=1)
    model_max.add(a)

    gamma_3= Variable('gamma_3', lb=gamma_ov, ub=1)
    model_max.add(gamma_3)

    gamma_5 = Variable('gamma_5', lb=gamma_ov, ub=1)
    model_max.add(gamma_5)


    # gamma_3_a= Variable('gamma_3_a', lb=gamma_ov, ub=1)
    # model_max.add(gamma_3_a)


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
        expression_gamma3_+=2**i*model_max.variables['d_2' + str(i)]*multip_factor_gamma3

    'for random gamma 5'
    expression_gamma5_=0+model_max.variables['gamma_5'].lb
    multip_factor_gamma5=(model_max.variables['gamma_5'].ub-model_max.variables['gamma_5'].lb)/N
    for i in range(binary_variables + 1):
        expression_gamma5_+=2**i*model_max.variables['d_5' + str(i)]*multip_factor_gamma5


    # expression_gamma3_a=0+model_max.variables['gamma_3_a'].lb**2
    # for j in range(binary_variables + 1):
    #     expression_gamma3_a += gamma_ov * multip_factor_a * 2 ** j * (model_max.variables['d_a' + str(j)]) + \
    #                        gamma_ov * multip_factor_gamma3 * 2 ** j * (model_max.variables['d_2' + str(j)])
    #     for k in range(binary_variables + 1):
    #         expression_gamma3_a += multip_factor_gamma3 * multip_factor_a * 2 ** (k + j) * (model_max.variables['t_' + str(j) + '_' +str(k)])



    c_gamma=Constraint(gamma_1-expression_gamma,lb=0,ub=0)

    c_a=Constraint(a-expression_a,lb=0,ub=0)

    'stating that gamma_1*gamma_2<gamma_1 & gamma_1*gamma_2<gamma_5'
    c_a_gamma=Constraint(a-gamma_1, ub=0)

    'for random'
    c_a_gamma_5=Constraint(a-gamma_5, ub=0)


    c_gamma3=Constraint(gamma_3-expression_gamma3_,lb=0,ub=0)

    c_gamma_5=Constraint(gamma_5-expression_gamma5_,lb=0,ub=0)

    #if this is required or not

    #c_gamma3_a=Constraint(gamma_3_a-expression_gamma3_a,lb=0,ub=0)

    #c_gamma3_a_1=Constraint(gamma_3_a-gamma_3,ub=0)

    #c_gamma3_a_2=Constraint(gamma_3_a-a,ub=0)


    model_max.add([c_mass_cons,c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_10, c_11,c_12, c_gamma,c_a,c_a_gamma,c_a_gamma_5,c_gamma3,c_gamma_5]) #,c_gamma3_a,c_gamma3_a_1,c_gamma3_a_2])
    #write expressions in for loop

    expression={}
    expression_m={}
    expression_gamma3={}
    expression_t={}
    expression_gamma5={}

    for i in range(n_step):
        'this is the lower bound lower bound  + 2** terms'
        expression[i]=model_max.variables['z'+str(i+1)]*gamma_ov
        expression_m[i]=model_max.variables['z'+str(i+1)]*gamma_ov
        expression_gamma3[i]=model_max.variables['z'+str(i+1)]*gamma_ov

        'for random'
        expression_gamma5[i]=model_max.variables['z'+str(i+1)]*gamma_ov


        expression_t[i]=model_max.variables['z'+str(i+1)]*gamma_ov**2


        for j in range(binary_variables + 1):
            expression[i]+=multip_factor*2**j*(model_max.variables['k_'+str(i+1)+'_'+str(j)])
            expression_m[i]+=multip_factor_a*2**j*(model_max.variables['m_'+str(i+1)+'_'+str(j)])
            expression_gamma3[i]+=multip_factor_gamma3*2**j*(model_max.variables['k3_'+str(i+1)+'_'+str(j)])
            'for random'
            expression_gamma5[i]+=multip_factor_gamma5*2**j*(model_max.variables['k5_'+str(i+1)+'_'+str(j)])

            'this expression t is more complicated as it is also multiplication of 2 integers'
            expression_t[i]+=gamma_ov*multip_factor_a*2**j*(model_max.variables['m_'+str(i+1)+'_'+str(j)])+ \
                             gamma_ov * multip_factor_gamma3 * 2 ** j * ( model_max.variables['k3_' + str(i + 1) + '_' + str(j)])
            for k in range(binary_variables + 1):
                expression_t[i] +=multip_factor_gamma3*multip_factor_a*2**(k+j) * ( model_max.variables['mt_' + str(i + 1) + '_' + str(j)+ '_' + str(k)])

    expression_m_lower = {}
    'for k'
    expression[n_step]=model_max.variables['v_upper']*gamma_ov
    expression_m[n_step]=model_max.variables['v_upper']*gamma_ov
    expression_m_lower[n_step]=model_max.variables['v_lower']*gamma_ov

    expression_gamma3[n_step]=model_max.variables['v']*gamma_ov
    'for random'
    expression_gamma5[n_step]=model_max.variables['v_lower']*gamma_ov

    expression_t[n_step]=model_max.variables['v']*gamma_ov**2
    for j in range(binary_variables + 1):
        expression[n_step] += multip_factor * 2 ** j * (model_max.variables['k_' + str(n_step+1) + '_' + str(j)])
        expression_m[n_step] += multip_factor_a * 2 ** j * (model_max.variables['m_upper' + str(n_step+1) + '_' + str(j)])

        expression_m_lower[n_step] += multip_factor_a * 2 ** j * (model_max.variables['m_lower' + str(n_step+1) + '_' + str(j)])

        expression_gamma3[n_step] += multip_factor_gamma3 * 2 ** j * (model_max.variables['k3_' + str(n_step+1) + '_' + str(j)])
        'for random'
        expression_gamma5[n_step] += multip_factor_gamma5 * 2 ** j * (model_max.variables['k5_' + str(n_step+1) + '_' + str(j)])

        expression_t[n_step] += gamma_ov * multip_factor_a * 2 ** j * (
        model_max.variables['m_' + str(n_step + 1) + '_' + str(j)]) + \
                           gamma_ov * multip_factor_gamma3 * 2 ** j * (
                           model_max.variables['k3_' + str(n_step + 1) + '_' + str(j)])
        for k in range(binary_variables + 1):
            expression_t[n_step] += multip_factor_gamma3 * multip_factor_a * 2 ** (k + j) * (
                #todo error was here you need to solve this make it automatic
            model_max.variables['mt_' + str(n_step + 1) + '_' + str(j)+ '_' + str(k)])

    #todo to UPDATE
    c_rate_1=Constraint((model_max.variables['z1']-expression[0])-expression[n_step],lb=0, ub=0) #k1b
    c_rate_2=Constraint(expression[1]-expression_m[1]-expression_m[n_step],lb=0, ub=0) #k2b
    c_rate_3=Constraint(model_max.variables['z3']-expression_gamma3[2]-expression_gamma3[n_step],lb=0, ub=0) #k3b
    c_rate_4=Constraint(model_max.variables['z4']-expression_gamma3[3]-model_max.variables['v'],lb=0, ub=0)  #k3f
    c_rate_5=Constraint(expression_t[4]-model_max.variables['z5']*gamma_ov-expression_t[n_step],lb=0, ub=0) #k4f
    c_rate_6=Constraint(A*(model_max.variables['z6']-expression[5])-model_max.variables['v_upper'],lb=0, ub=0) #k1f
    c_rate_7=Constraint(B*(expression[6]-expression_m[6])-expression[n_step],lb=0, ub=0) #k2f
    c_rate_8=Constraint(P*(expression_t[7]-model_max.variables['z8']*gamma_ov)-model_max.variables['v']*gamma_ov,lb=0, ub=0) #k4b

    c_rate_9=Constraint(B*(model_max.variables['z9']-expression_gamma5[8])-model_max.variables['v_lower'],lb=0, ub=0)  #k5f
    c_rate_10=Constraint((model_max.variables['z10']-expression_gamma5[9])-expression_gamma5[n_step],lb=0, ub=0)#k5b
    c_rate_11=Constraint(A*(expression_gamma5[10]-expression_m[10])-expression_gamma5[n_step],lb=0, ub=0) #k6f
    c_rate_12=Constraint(expression_gamma5[11]-expression_m[11]-expression_m_lower[n_step],lb=0, ub=0) #k6b

    c_split_ratio = Constraint(model_max.variables['v'] - (model_max.variables['v_upper']+model_max.variables['v_lower']),lb=0, ub=0)

    model_max.add([c_rate_1,c_rate_2,c_rate_3,c_rate_4,c_rate_5,c_rate_6,c_rate_7,c_rate_8, c_rate_9, c_rate_10, c_rate_11, c_rate_12,c_split_ratio])


    obj = Objective(0, direction='max')
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
    dict_values['A'] = A
    dict_values['B'] = B
    dict_values['P'] = P
    #dict_values['alpha'] = alpha

    # dict_values['gamma_1'] = gamma_1
    # dict_values['gamma_2'] = gamma_2
    dict_values['gamma_4'] = gamma_ov / (a.primal*gamma_3.primal)
    dict_values['gamma_2'] = a.primal / gamma_1.primal
    dict_values['gamma_6'] = a.primal / gamma_5.primal

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

    dict_values['alpha'] = model_max.variables['v_upper'].primal /model_max.variables['v'].primal

    #todo variability analysis
    if variability_analysis == True:
        df_var_analysis = pd.DataFrame()
        df_var_analysis= variability_analysis_linearized_feasibility_biuni_random(model_max)
        # print(df_var_analysis)
        #dict_values['unique'] = a
        # return (np.exp(model_max.variables.ln_v.primal)),(dict_values),df_var_analysis
    # for fitness you are directly returning the v value not log format

    return ((model_max.variables.v.primal)), (dict_values),(model_max), (df_var_analysis)

