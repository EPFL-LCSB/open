OpEn: Optimal Enzyme
==========================================
This repository contains the workflow to framework to explore the catalytically optimal modes of operations of enzymatic reactions.

More information can be found in the manuscript
`Optimal enzyme utilization suggests concentrations and thermodynamics favor condition-specific saturations and binding mechanisms <https://www.biorxiv.org/content/10.1101/2022.04.12.488028v1.full>`_
The package is developed using python 3.6 and run in Docker (20.10.6) containers.
Tested with solvers cplex and gurobi (default implemented on cplex_interface)

Recommended to be run in docker containers, with dedicated solver installed.
`Setting up the python API of cplex <https://www.ibm.com/docs/en/icos/12.8.0.0?topic=cplex-setting-up-python-api>`_  in case docker based installation is not used

Requirements
------------
You will need to have `Git-LFS <https://git-lfs.github.com/>`_ in order to properly download some binary files:

.. code:: bash

    git clone https://github.com/EPFL-LCSB/open.git /path/to/open
    cd /path/to/open
    git lfs install
    git lfs pull

Further the following pip-python packages are required (can be found in detail in requirements.txt

    - optlang
    - cobra==0.17.1
    - numpy<=1.17.4
    - pandas
    - matplotlib
    - tables
    - sklearn
    - ipython
    - jedi==0.17.2
    - tqdm
    - scipy



Container-based install
-----------------------

You might want to use this program inside of a container. The
|docker|_
subfolder has all the necessary information and source files to set it
up.

.. |docker| replace:: ``docker/``
.. _docker: https://github.com/EPFL-LCSB/open/tree/master/docker


.. code:: bash

    cd open/docker
    ./build.sh
    ./run.sh

Building the docker image takes approximately 5 mins.



Setup
=====
You can install this module from source using ``pip``:
*For Python 3, you might have to use* ``pip3`` *instead of* ``pip``

.. code:: bash

    git clone https://github.com/EPFL-LCSB/open.git /path/to/open
    pip3 install -e /path/to/open

The installation process should not exceed a minute if the requirements are installed. If they are not, it might take longer as the installer installs them first.


Quick start
===========
As a demo to generate optimal operating conditions for a 3-step MichaelisMenten mechanism please find below a simple example to get started:

.. code-block:: python


    from open.optim.LP_MILP_wpiecewise import *
    import time
    import pandas as pd

    t=time.time()
    S_concentration=3.0
    P_concentration=2.0
    q_equilibrium=2.0
    gamma_overall=P_concentration/S_concentration/q_equilibrium
    obj_primal_milp, variables_primal_milp, var_analysis = milp_problem_3step(gamma_overall, q=q_equilibrium,
                                                                          S=S_concentration, variability_analysis=False)
    df_milp_variables = pd.DataFrame(variables_primal_milp, index=[0])

    elapsed = time.time() - t
    print('time for optimization', elapsed)

Similarly for ordered multisubstrate mechanism A+B-->P

.. code-block:: python


    from open.optim.LP_MILP_wpiecewise import *
    import time
    import pandas as pd

    t=time.time()
    P_conc=1.0
    q_equilibrium=2.0
    A_concentration=3.0
    gamma_overall=0.5

    B_concentration = P_conc/(A_concentration*q_equilibrium*gamma_overall)

    obj_primal_milp_4step, variables_primal_milp_4step,var_analysis = milp_problem_4step_biuni(gamma_overall, q=q_equilibrium,
                                                                                      S=A_concentration, P=P_conc,
                                                                                      variability_analysis=False)
    df_milp_variables = pd.DataFrame(variables_primal_milp_4step,index=[0])

    elapsed = time.time() - t
    print('time for optim', elapsed)


And for random ordered multisubstrate mechanism A+B-->P

.. code-block:: python

    from open.optim.LP_MILP_random import *
    import time
    import pandas as pd
    P_concentration=1.0
    A_concentration=3.0
    B_concentration=3.0
    q_equilibrium=2.0
    df_st = pd.DataFrame(columns=['A', 'B', 'P', 'q', 'alpha_max', 'alpha_min', 'v_net', 'gamma_ov'])
    gamma_overall = P_concentration / A_concentration / q_equilibrium / B_concentration

    t = time.time()
    obj_primal_milp_4step_random_split, variables_primal_milp_4step_random_split, var_analysis_feasibility = milp_problem_4step_biuni_random_split_ratio(
        gamma_overall, q=q_equilibrium,
        S=A_concentration, P=P_concentration,
        variability_analysis=True)
    elapsed = time.time() - t
    print('time for optimization', elapsed)
    split_max = var_analysis_feasibility.loc['v_upper', 'maximum'] / obj_primal_milp_4step_random_split
    split_min = var_analysis_feasibility.loc['v_upper', 'minimum'] / obj_primal_milp_4step_random_split
    df_milp_variables_random_split = pd.DataFrame(variables_primal_milp_4step_random_split, index=[0])

    values_to_add = {'A': A_concentration, 'B': B_concentration, 'P': P_concentration, 'q': q_equilibrium,
                     'alpha_max': split_max, \
                     'alpha_min': split_min, 'v_net': obj_primal_milp_4step_random_split, 'gamma_ov': gamma_overall}
    row_to_add = pd.Series(values_to_add, name=str(0))
    df_st = df_st.append(row_to_add)


Generating optimal operating conditions for one data point should take around 2-10 seconds depending
on if variability analysis is performed or not


License
=======
The software in this repository is put under an APACHE licensing scheme - please see the `LICENSE <https://github.com/EPFL-LCSB/open/blob/master/LICENSE.txt>`_ file for more details.




