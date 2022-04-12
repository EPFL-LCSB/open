#!/bin/bash


#overall_gamma=np.array([0.02,0.05,0.08,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99])

#In order to run in parallel for isolines (thermodynamic displacements GAMMA) generate equally spaced P,S pairs on each isoline and run for them
#overall_gamma=np.array([0.02,0.021,0.023,0.025,0.027,0.03,0.035,0.04,0.045,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.99])

#alphas to be tried 0.8 0.5 0.2
alpha=0.5

#
python BiUni_random_parallel.py  0.02  $alpha
#python BiUni_random_parallel.py  0.021 $alpha
#python BiUni_random_parallel.py   0.023 $alpha
#python BiUni_random_parallel.py   0.025 $alpha
#python BiUni_random_parallel.py   0.027 $alpha
python BiUni_random_parallel.py   0.03  $alpha
#python BiUni_random_parallel.py   0.035 $alpha
python BiUni_random_parallel.py   0.04  $alpha
#python BiUni_random_parallel.py   0.045 $alpha
python BiUni_random_parallel.py   0.05  $alpha
python BiUni_random_parallel.py   0.06  $alpha
python BiUni_random_parallel.py   0.07  $alpha
python BiUni_random_parallel.py   0.08  $alpha
python BiUni_random_parallel.py   0.09  $alpha
python BiUni_random_parallel.py   0.1   $alpha
#python BiUni_random_parallel.py   0.12  $alpha
#python BiUni_random_parallel.py   0.15  $alpha
python BiUni_random_parallel.py   0.2   $alpha
python BiUni_random_parallel.py   0.3   $alpha
python BiUni_random_parallel.py   0.4   $alpha
python BiUni_random_parallel.py   0.5   $alpha
python BiUni_random_parallel.py   0.6   $alpha
python BiUni_random_parallel.py   0.7   $alpha
python BiUni_random_parallel.py   0.8   $alpha
python BiUni_random_parallel.py   0.9   $alpha
#python BiUni_random_parallel.py   0.99  $alpha &