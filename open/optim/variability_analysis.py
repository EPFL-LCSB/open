import re
from collections import namedtuple
from functools import partial
from types import ModuleType
from warnings import warn
import optlang
from optlang.cplex_interface import Model, Variable, Constraint, Objective

import pandas as pd
from optlang.interface import (
    FEASIBLE, INFEASIBLE, ITERATION_LIMIT, NUMERIC, OPTIMAL, SUBOPTIMAL,
    TIME_LIMIT)
from optlang.symbolics import Basic, Zero
from pandas import DataFrame
from numpy import zeros
import logging
import multiprocessing
from tqdm import tqdm


" Taken from pytfa "


def slim_optimize(self, error_value=float('nan'), message=None):
    """Optimize model without creating a solution object.

    Creating a full solution object implies fetching shadow prices and
    flux values for all reactions and metabolites from the solver
    object. This necessarily takes some time and in cases where only one
    or two values are of interest, it is recommended to instead use this
    function which does not create a solution object returning only the
    value of the objective. Note however that the `optimize()` function
    uses efficient means to fetch values so if you need fluxes/shadow
    prices for more than say 4 reactions/metabolites, then the total
    speed increase of `slim_optimize` versus `optimize` is  expected to
    be small or even negative depending on how you fetch the values
    after optimization.

    Parameters
    ----------
    error_value : float, None
       The value to return if optimization failed due to e.g.
       infeasibility. If None, raise `OptimizationError` if the
       optimization fails.
    message : string
       Error message to use if the model optimization did not succeed.

    Returns
    -------
    float
        The objective value.
    """
    self.optimize()
    if self.status == optlang.interface.OPTIMAL:
        return self.objective.value
    elif error_value is not None:
        return error_value
    else:
        assert_optimal(self, message)


def assert_optimal(model, message='optimization failed'):
    """Assert model solver status is optimal.

    Do nothing if model solver status is optimal, otherwise throw
    appropriate exception depending on the status.

    Parameters
    ----------
    model : cobra.Model
        The model to check the solver status for.
    message : str (optional)
        Message to for the exception if solver status was not optimal.
    """
    status = model.status
    if status != OPTIMAL:
        exception_cls = OPTLANG_TO_EXCEPTIONS_DICT.get(
            status, OptimizationError)
        raise exception_cls("{} ({})".format(message, status))

#epsilon=1e-20


def variability_analysis(_tmodel, kind='reactions'):

    tmodel= Model.clone(_tmodel)


    objective = _tmodel.objective

    tmodel.variables.v.lb=objective.value
    results = {'min':{}, 'max':{}}
    for sense in ['min','max']:
        for k,var in tqdm(tmodel.variables.items(), desc=sense+'imizing'):
            #tmodel.logger.debug(sense + '-' + k)
            results[sense][k] = _variability_analysis_element(tmodel,var,k,sense)
    #tmodel.objective = objective
    df = pd.DataFrame(results)
    df.rename(columns={'min': 'minimum', 'max': 'maximum'}, inplace=True)
    #
    if ((df.maximum-df.minimum)>1e-10).any():
        print('NOT UNIQUE')
        print(df)
          # raise ValueError('This is not a unique Solution')

    if ((df.maximum-df.minimum)<=1e-10).all():
         print('UNIQUE')

    return df


def _variability_analysis_element(tmodel, var,k ,sense):

    obj2 = Objective(var, direction=sense)
    tmodel.objective = obj2
    tmodel.objective.direction = sense
    sol = slim_optimize(tmodel)

    return sol



def variability_analysis_ln(_tmodel, kind='reactions'):

    tmodel= Model.clone(_tmodel)
    accuracy=tmodel.configuration.tolerances.feasibility

    objective = _tmodel.objective

    tmodel.variables.ln_v.lb=objective.value
    results = {'min':{}, 'max':{}}
    for sense in ['min','max']:
        for k,var in tqdm(tmodel.variables.items(), desc=sense+'imizing'):
            #tmodel.logger.debug(sense + '-' + k)
            results[sense][k] = _variability_analysis_element(tmodel,var,k,sense)
    #tmodel.objective = objective
    df = pd.DataFrame(results)
    df.rename(columns={'min': 'minimum', 'max': 'maximum'}, inplace=True)
    print(df)
    #
    if ((df.maximum-df.minimum)>accuracy).any():
        print('NOT UNIQUE PROBLEM!!!!')
        print(df)
          # raise ValueError('This is not a unique Solution')

    if ((df.maximum-df.minimum)<=accuracy).all():
         print('UNIQUE')

    return df

# todo check this

def variability_analysis_linearized(_tmodel, kind='reactions'):

    tmodel= Model.clone(_tmodel)
    accuracy=tmodel.configuration.tolerances.feasibility

    objective = _tmodel.objective
    print(objective.value)
    tmodel.variables.v.lb=objective.value
    #tmodel.variables.v.ub=objective.value

    results = {'min':{}, 'max':{}}
    for sense in ['min','max']:
        for k,var in tqdm(tmodel.variables.items(), desc=sense+'imizing'):
            'todo fix this this is only to run for EA and EB '
            # if (k=='v_upper')or (k=='v_lower'):
            if (k=='v_upper')or (k=='v_lower') :\
                    #or (k=='E') or (k=='EA') or (k=='EAB') or (k=='EP') or (k=='EPQ') or(k=='EB') or (k=='z1') \
                    #or (k=='z2') or (k=='z3') or (k=='z4') or (k=='z5') or (k=='z6') or (k=='z7') or (k=='z8') or (k=='z9') \
                    #or (k=='z10')or (k=='z11')or (k=='z12') or (k=='a')or (k=='b') or  (k=='b_a') or (k=='gamma_3')or (k=='gamma_1'):
            #tmodel.logger.debug(sense + '-' + k)
                results[sense][k] = _variability_analysis_element(tmodel,var,k,sense)
    #tmodel.objective = objective
    df = pd.DataFrame(results)
    df.rename(columns={'min': 'minimum', 'max': 'maximum'}, inplace=True)
    print(df)
    #
    if ((df.maximum-df.minimum)>accuracy).any():
        print('NOT UNIQUE PROBLEM!!!!')
        print(df)
          # raise ValueError('This is not a unique Solution')

    if ((df.maximum-df.minimum)<=accuracy).all():
         print('UNIQUE')

    return df

def variability_analysis_general(_tmodel, kind='reactions',specify=True):

    tmodel= Model.clone(_tmodel)
    accuracy=tmodel.configuration.tolerances.feasibility

    objective = _tmodel.objective
    print(objective.value)
    tmodel.variables.v.lb=objective.value
    #tmodel.variables.v.ub=objective.value

    results = {'min':{}, 'max':{}}
    for sense in ['min','max']:
        for k,var in tqdm(tmodel.variables.items(), desc=sense+'imizing'):
            if specify:
                # if (k=='v_upper')or (k=='v_lower'):
                if (k=='v_upper')or (k=='v_lower') \
                        or (k=='E') or (k=='EA') or (k=='EAB') or (k=='EP') or (k=='EPQ') or(k=='EB') or (k=='z1') \
                        or (k=='z2') or (k=='z3') or (k=='z4') or (k=='z5') or (k=='z6') or (k=='z7') or (k=='z8') or (k=='z9') \
                        or (k=='z10')or (k=='z11')or (k=='z12') or (k=='a')or (k=='b') or  (k=='b_a') or (k=='gamma_3')or (k=='gamma_1'):
                #tmodel.logger.debug(sense + '-' + k)
                    results[sense][k] = _variability_analysis_element(tmodel,var,k,sense)

            else:
                results[sense][k] = _variability_analysis_element(tmodel, var, k, sense)

    #tmodel.objective = objective
    df = pd.DataFrame(results)
    df.rename(columns={'min': 'minimum', 'max': 'maximum'}, inplace=True)
    print(df)
    #
    if ((df.maximum-df.minimum)>accuracy).any():
        print('NOT UNIQUE PROBLEM!!!!')
        print(df)
          # raise ValueError('This is not a unique Solution')

    if ((df.maximum-df.minimum)<=accuracy).all():
         print('UNIQUE')

    return df
def variability_analysis_linearized_feasibility(_tmodel, kind='reactions'):

    tmodel= Model.clone(_tmodel)
    accuracy=tmodel.configuration.tolerances.feasibility

    objective = _tmodel.objective

    #tmodel.variables.v.lb=objective.value
    #tmodel.variables.v.ub=objective.value

    results = {'min':{}, 'max':{}}
    for sense in ['min','max']:
        for k,var in tqdm(tmodel.variables.items(), desc=sense+'imizing'):
            #tmodel.logger.debug(sense + '-' + k)
            results[sense][k] = _variability_analysis_element(tmodel,var,k,sense)
    #tmodel.objective = objective
    df = pd.DataFrame(results)
    df.rename(columns={'min': 'minimum', 'max': 'maximum'}, inplace=True)
    print(df)
    #
    if ((df.maximum-df.minimum)>accuracy).any():
        print('NOT UNIQUE PROBLEM!!!!')
        print(df)
          # raise ValueError('This is not a unique Solution')

    if ((df.maximum-df.minimum)<=accuracy).all():
         print('UNIQUE')

    return df


def variability_analysis_linearized_feasibility_biuni_random(_tmodel, kind='reactions'):

    tmodel= Model.clone(_tmodel)
    accuracy=tmodel.configuration.tolerances.feasibility

    objective = _tmodel.objective

    #tmodel.variables.v.lb=objective.value
    #tmodel.variables.v.ub=objective.value

    results = {'min':{}, 'max':{}}
    for sense in ['min','max']:
        for k,var in tqdm(tmodel.variables.items(), desc=sense+'imizing'):
            if (k=='gamma_1')or (k=='gamma_3') or (k=='gamma_5') or (k=='a') or (k=='v_upper')  or (k=='v_lower'):
            #tmodel.logger.debug(sense + '-' + k)
                results[sense][k] = _variability_analysis_element(tmodel,var,k,sense)
    #tmodel.objective = objective
    df = pd.DataFrame(results)
    df.rename(columns={'min': 'minimum', 'max': 'maximum'}, inplace=True)
    print(df)
    #
    if ((df.maximum-df.minimum)>accuracy).any():
        print('NOT UNIQUE PROBLEM!!!!')
        print(df)
          # raise ValueError('This is not a unique Solution')

    if ((df.maximum-df.minimum)<=accuracy).all():
         print('UNIQUE')

    return df


def variability_analysis_dual(_tmodel, kind='reactions'):

    tmodel= Model.clone(_tmodel)
    accuracy=tmodel.configuration.tolerances.feasibility

    objective = _tmodel.objective

    tmodel.variables.P13.ub=objective.value
    results = {'min':{}, 'max':{}}
    for sense in ['min','max']:
        for k,var in tqdm(tmodel.variables.items(), desc=sense+'imizing'):
            #tmodel.logger.debug(sense + '-' + k)
            results[sense][k] = _variability_analysis_element(tmodel,var,k,sense)
    #tmodel.objective = objective
    df = pd.DataFrame(results)
    df.rename(columns={'min': 'minimum', 'max': 'maximum'}, inplace=True)
    print(df)
    #
    if ((df.maximum-df.minimum)>accuracy).any():
        print('NOT UNIQUE PROBLEM!!!!')
        print(df)
          # raise ValueError('This is not a unique Solution')

    if ((df.maximum-df.minimum)<=accuracy).all():
         print('UNIQUE')

    return df