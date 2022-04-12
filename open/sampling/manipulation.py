from optlang.cplex_interface import Model, Variable, Constraint, Objective

'this is also from pytfa'
INTEGER_VARIABLE_TYPES = ('binary','integer')
def apply_generic_variability(tmodel,va, inplace = True):
    """
    Reactions a dealt with cobra, but the other variables added use pytfa's
    interface: the class GenericVariable. We use a different method to apply
    variability directly in the solver

    :param tmodel:
    :param va:
    :param inplace:
    :return:
    """
    if inplace:
        _tmodel = tmodel
    else:
        _tmodel = tmodel.copy()


    variables_list=tmodel.variables.keys()

    for varname in variables_list:
    #for varname in va.index:

        diff=abs((va.loc[varname].maximum - va.loc[varname].minimum))
        if (diff<=1e-8):
             the_min,the_max = va.loc[varname,['minimum','maximum']]
             _tmodel.variables[varname].lb = the_max
             _tmodel.variables[varname].ub = the_max
        else:
            the_min,the_max = va.loc[varname,['minimum','maximum']]
            _tmodel.variables[varname].lb = the_min
            _tmodel.variables[varname].ub = the_max

    return _tmodel

def apply_generic_variability_2(tmodel,va, inplace = True):
    """
    Reactions a dealt with cobra, but the other variables added use pytfa's
    interface: the class GenericVariable. We use a different method to apply
    variability directly in the solver

    :param tmodel:
    :param va:
    :param inplace:
    :return:
    """
    if inplace:
        _tmodel = tmodel
    else:
        _tmodel = tmodel.copy()


    variables_list=tmodel.variables.keys()

    for varname in variables_list:
        if varname=='v':
        #for varname in va.index:
            the_min,the_max = va.loc[varname,['minimum','maximum']]
            _tmodel.variables[varname].lb = the_min
            _tmodel.variables[varname].ub = the_max
        else:
            _tmodel.variables[varname].lb = 0
            _tmodel.variables[varname].ub = 1

    return _tmodel

#but me do i need this
def strip_from_integer_variables(tmodel):
    """
    Removes all integer and binary variables of a cobra_model, to make it sample-able
    :param tmodel:
    :return:
    """
    #changed .copy to clone
    #continuous_model = tmodel.copy()
    continuous_model = Model.clone(tmodel)
    #Model_optim(continuous_model)
    continuous_model.name = tmodel.name + ' - continuous'

    integer_variables = set()

    constraints_with_integer_variables = []

    # We go through all the constraint descriptors and check if at least one of
    # their variables is in the integer variable list
    for this_cons in continuous_model.constraints.values():
        has_integer_variable = False
        for this_var in this_cons.variables:
            if this_var.type in INTEGER_VARIABLE_TYPES:
                has_integer_variable += True
                this_var_descriptor = continuous_model.variables[this_var.name]
                integer_variables.add(this_var_descriptor)
        if has_integer_variable:
            constraints_with_integer_variables.append(this_cons)

    for this_cons in constraints_with_integer_variables:
        continuous_model.remove(this_cons)

    for this_var in integer_variables:
        continuous_model.remove(this_var)

    continuous_model.update()
    # This will update the values =
    print('Is the optimization model still integer ? {}' \
          .format(continuous_model.is_integer))

    return continuous_model
