GEEK
=====
|PyPI| |Documentation Status| |Build Status| |Codecov| |Codacy branch grade| |Code climate| |license| 

GEneralized Elementary Kinetics

Implements the GEEK analysis as described in : Daniel R. Weilandt and Vassily
Hatzimanikatis. "Particle-based simulation reveals macromolecular crowding effects on the Michaelis-Menten mechanism",
bioRxiv 429316; doi: `https://doi.org/10.1101/429316`_

Requirements
------------

**This module was developed in Python 3.5, and it is recommended to run Python 3.5 


Setup
=====

*This step is not required if you're using the container, which bundles all this.*

You can install this module with ``pip``:

*For Python 3, you might have to use* ``pip3`` *instead of* ``pip``

.. code:: bash

    pip3 install geek

or from source

.. code:: bash

    git clone https://github.com/EPFL-LCSB/pytfa.git /path/to/geek
    pip3 install -e /path/to/geek

Quick start
===========

A simple example can be found at:

::

    geek
    ├── data
    └── examples
        └── geek_example_analysis_pgm.py

   
License
========

The software in this repository is put under an APACHE-2.0 licensing scheme - please see the `LICENSE <https://github.com/EPFL-LCSB/geek/blob/master/LICENSE.txt>`_ file for more details.
