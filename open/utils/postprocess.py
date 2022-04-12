import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#convert from dictionary to dataframe because you also store LP variables

"""Back calculate elementary rate constants and 
    macroscopic kinetic parameters"""

def calculate_rate_constants_uniuni(df):
    df['k1f']=df['z4']/df['E']
    df['k2f']=df['z5']/df['ES']
    df['k3f']=df['z6']/df['EP']


    df['k1b']=df['z1']/df['ES']
    df['k2b']=df['z2']/df['EP']
    df['k3b']=df['z3']/df['E']

    return df


def calculate_kinetic_params_uniuni(df):
    df['kcat_scaled'] = df.k2f * df.k3f / (
            df.k2f + df.k3f + df.k2b)


    df['KM_S'] = (df.k2f * df.k3f + df.k1b * df.k3f + df.k1b * df.k2b) / \
                    (df.k1f * (df.k2f + df.k3f + df.k2b))
    df['KM_P'] = (df.k2f * df.k3f + df.k1b * df.k3f + df.k1b * df.k2b) / \
                    (df.k3b * (df.k2f + df.k1b + df.k2b))
    df['kcat_kms'] = df.kcat_scaled / df.KM_S
    df['kcat_back'] = df.k1b * df.k2b / (df.k2f + df.k1b + df.k2b)

    return df



def calculate_rate_constants_random_biuni(df):
    df['k1f'] = df['z6'] /df['E']
    df['k2f'] = df['z7'] /df['EA']
    df['k3f'] = df['z4']/df['EAB']
    df['k4f'] = df['z5']/df['EP']
    'for random'
    df['k5f'] = df['z9']/df['E']
    df['k6f'] = df['z11']/df['EB']

    df['k1b'] = df['z1']/df['EA']
    df['k2b'] = df['z2']/df['EAB']
    df['k3b'] = df['z3']/df['EP']
    df['k4b'] = df['z8']/df['E']
    'for random'
    df['k5b'] = df['z10']/df['EB']
    df['k6b'] = df['z12']/df['EAB']


    return df