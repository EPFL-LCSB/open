OpEn: Optimal Enzyme
==========================================
This repository contains the workflow to framework to explore the catalytically optimal modes of operations of enzymatic reactions.

More information can be found in the manuscript
"Optimal enzyme utilization suggests concentrations and thermodynamics favor condition-specific saturations and binding mechanisms"

The package is developped using python 3.6 and run in Docker (20.10.6) containers.

Recommended to be run in docker containers

Requirements
------------
You will need to have Git-LFS in order to properly download some binary files:
.. code:: bash
    git clone https://github.com/EPFL-LCSB/open.git /path/to/open
    cd /path/to/open
    git lfs install
    git lfs pull

Further the following pip-python packages are required (can be found in detail in requirements.txt
    - sympy >= 1.1.
    - pytest
    - scipy
    - numpy
    - bokeh
    - pandas
    - matplotlib
    - optlang
    - cobra
    - tables

Container-based install
-----------------------

You might want to use this program inside of a container. The
|docker|_
subfolder has all the necessary information and source files to set it
up.




cd open/docker
and ru
./build.sh

./run.sh

