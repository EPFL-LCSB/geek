# -*- coding: utf-8 -*-
"""
.. module:: geek
   :platform: Unix, Windows
   :synopsis: GEneralised Elementary Kinetics

.. moduleauthor:: geek team

[---------]

Copyright 2018 Laboratory of Computational Systems Biotechnology (LCSB),
Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""

from sympy.utilities.autowrap import ufuncify
import numpy as np


class OdeFun:
    def __init__(self, variables, parameters, expression_dict):
        """
        Class to create a ode function from symbolic expressions
        :param variables: list of symbols
        :param parameters: list of sybols
        :param expression_dict: dict of expressions for the rate of change of each variable
                                indexed by variable symbols {v1: p1*var1*var2, var2: ...}
        """

        self.variables = variables
        self.parameters = parameters

        self.input = variables + parameters
        self.expressions = [expression_dict[v] for v in variables]

        self.function = []
        for e in self.expressions:
            self.function.append(ufuncify(self.input, e, backend='Cython'))

    def __call__(self, t, x, parameter_dict):
        """
        Evaluates ODE expression create from expressions
        :param t: time scalar
        :param x: states numpy array
        :param parameter_dict: parameter dict indexed with symbols {sym_p1:value_p1, ....}
        :return: numpy array
        """

        parameters_values = [parameter_dict[str(p)] for p in self.parameters]
        x_list = list(x)
        input = x_list + parameters_values

        dxdt = np.array([f(*[np.array([i], dtype=np.float) for i in input])[0] for f in self.function]).T
        return dxdt


