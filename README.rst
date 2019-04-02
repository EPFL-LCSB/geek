GEEK
=====

GEneralized Elementary Kinetics

Implements the GEEK analysis as described in : Daniel R. Weilandt and Vassily
Hatzimanikatis. "Particle-based simulation reveals macromolecular crowding effects on the Michaelis-Menten mechanism",
bioRxiv 429316; doi: `https://doi.org/10.1101/429316`_


Setup
=====

Container-based install
-----------------------

We recommend to use this package inside of the provided DOCKER container.
See |docker|_ for the setup and use of the container.

.. |docker| replace:: ``docker/``
.. _docker: https://github.com/EPFL-LCSB/geek/tree/master/docker



Source-based install
--------------------

*This step is not required if you're using the container, which bundles all this.*

If you wish not to use the container based install you can install the package from source

.. code:: bash

    git clone https://github.com/EPFL-LCSB/geek.git /path/to/geek
    pip3 install -e /path/to/geek


Requirements
------------

This module was developed in Python 3.5, and it is recommended to run Python 3.5.
The module also was tested in Python 3.6.

This module requires |OPENBREAD|_ to work properly.

.. |OPENBREAD| replace:: ``OPENBReaD``
.. _OPENBREAD: https://github.com/EPFL-LCSB/openbread/tree/master

Further the following pip-python packages are required:
    - sympy >= 1.1.
    - pytest
    - scipy
    - pandas
    - statsmodels


Examples and Validation
=======================

You find examples of the analysis performed in the paper:

::

    geek
    ├── data
    └── paper
        ├── geek_example_analysis_pgm.py
        ├── geek_qssa_example.py
        └── geek_validation_dissociation_association.py

   
License
========

The software in this repository is put under an APACHE-2.0 licensing scheme - please see the `LICENSE <https://github.com/EPFL-LCSB/geek/blob/master/LICENSE.txt>`_ file for more details.
