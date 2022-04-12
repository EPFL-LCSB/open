OpEn: Optimal Enzyme
==========================================
This repository contains the workflow to framework to explore the catalytically optimal modes of operations of enzymatic reactions.

More information can be found in the manuscript
"Optimal enzyme utilization suggests concentrations and thermodynamics favor condition-specific saturations and binding mechanisms"

The package is developed using python 3.6 and run in Docker (20.10.6) containers.
Tested with solvers cplex and gurobi (default implemented on cplex_interface)

Recommended to be run in docker containers, with dedicated solver installed

Requirements
------------
You will need to have Git-LFS in order to properly download some binary files:
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

Setup
=====
You can install this module from source using ``pip``:
*For Python 3, you might have to use* ``pip3`` *instead of* ``pip``

.. code:: bash

    git clone https://github.com/EPFL-LCSB/open.git /path/to/open
    pip3 install -e /path/to/open






