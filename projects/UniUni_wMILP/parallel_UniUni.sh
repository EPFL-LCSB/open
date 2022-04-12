#!/bin/bash


#overall_gamma=np.array([0.02,0.05,0.08,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99])

#In order to run in parallel for isolines (thermodynamic displacements GAMMA) generate equally spaced P,S pairs on each isoline and run for them

#
#python UniUni2step.py 0.02 &
#python UniUni2step.py 0.05 &
#python UniUni2step.py 0.08 &
#python UniUni2step.py 0.1 &
#python UniUni2step.py 0.15 &
#python UniUni2step.py 0.2 &
#python UniUni2step.py 0.25 &
#python UniUni2step.py 0.3 &
#python UniUni2step.py 0.4 &
#python UniUni2step.py 0.5 &
#python UniUni2step.py 0.6 &
#python UniUni2step.py 0.7 &
#python UniUni2step.py 0.8 &
#python UniUni2step.py 0.9 &


#python UniUni2step.py 0.99 &

#
python UniUni3step.py 0.005
python UniUni3step.py 0.01
python UniUni3step.py 0.02
python UniUni3step.py 0.05
python UniUni3step.py 0.08
python UniUni3step.py 0.1
python UniUni3step.py 0.15
python UniUni3step.py 0.2
python UniUni3step.py 0.25
python UniUni3step.py 0.3
python UniUni3step.py 0.35
python UniUni3step.py 0.4
python UniUni3step.py 0.45
python UniUni3step.py 0.5
python UniUni3step.py 0.55
python UniUni3step.py 0.6
python UniUni3step.py 0.7
python UniUni3step.py 0.8
python UniUni3step.py 0.9
python UniUni3step.py 0.95
python UniUni3step.py 0.99
