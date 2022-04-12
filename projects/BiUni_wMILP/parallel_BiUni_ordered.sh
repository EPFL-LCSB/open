#!/bin/bash


#overall_gamma=np.array([0.02,0.05,0.08,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99])

#In order to run in parallel add &&
# for isolines (thermodynamic displacements GAMMA) generate equally spaced P,S pairs on each isoline and run for them


python BiUni_ordered.py  0.02
python BiUni_ordered.py  0.021
python BiUni_ordered.py  0.023
python BiUni_ordered.py  0.025
python BiUni_ordered.py  0.027
python BiUni_ordered.py  0.03
python BiUni_ordered.py  0.035
python BiUni_ordered.py  0.04
python BiUni_ordered.py  0.045
python BiUni_ordered.py  0.05
python BiUni_ordered.py  0.06
python BiUni_ordered.py  0.07
python BiUni_ordered.py  0.08
python BiUni_ordered.py  0.09
python BiUni_ordered.py  0.1
python BiUni_ordered.py  0.12
python BiUni_ordered.py  0.15
python BiUni_ordered.py  0.2
python BiUni_ordered.py  0.3
python BiUni_ordered.py  0.4
python BiUni_ordered.py  0.5
python BiUni_ordered.py  0.6
python BiUni_ordered.py  0.7
python BiUni_ordered.py  0.8
python BiUni_ordered.py  0.9
#python BiUni_ordered.py  0.99