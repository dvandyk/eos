# Copyright (c) 2018 Frederik Beaujean
# Copyright (c) 2017, 2018 Danny van Dyk
#
# This file is part of the EOS project. EOS is free software;
# you can redistribute it and/or modify it under the terms of the GNU General
# Public License version 2, as published by the Free Software Foundation.
#
# EOS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA  02111-1307  USA

from _eos import _Parameters

class Parameters(_Parameters):
    @staticmethod
    def FromWCxf(wc):

        # the following coefficients are treated as real-valued in EOS
        real_coeffs = [
            'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c8', 'c8\''
        ]

        parameters = _Parameters.Defaults()
        for name in wc.dict:
            qn     = eos.QualifiedName(name)
            prefix = qn.prefix_part()
            coeff  = qn.name_part()
            value  = wc.dict[name]
            #if (not value.imag == 0) and (coeff in real_coeffs):
            #    raise ValueError('WC {0} does not support non-zero imaginary part'.format(name))

            # Add values provided by WCxf to EOS central (SM) values
            if str(coeff) in real_coeffs:
                p = parameters[str(prefix) + '::' + str(coeff)]
                p.set(p.central() + value.real)
            else:
                pr = parameters[str(prefix) + '::Re{' + str(coeff) + '}']
                pi = parameters[str(prefix) + '::Im{' + str(coeff) + '}']
                pr.set(pr.central() + value.real)
                pi.set(pi.central() + value.imag)

        return parameters

from _eos import *

from .data import *
from .plot import *
from .analysis import Analysis, BestFitPoint
from .observable import Observables
from .constraint import Constraints

import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

__ipython__ = False
try:
    if __IPYTHON__:
        __ipython__ = True
        ip = get_ipython()
        html_formatter = ip.display_formatter.formatters['text/html']

        from .ipython import __format_Parameter, __format_KinematicVariable, __format_Observable, __format_ObservableEntry, __format_GoodnessOfFit
        html_formatter.for_type(Parameter, __format_Parameter)
        html_formatter.for_type(KinematicVariable, __format_KinematicVariable)
        html_formatter.for_type(Observable, __format_Observable)
        html_formatter.for_type(ObservableEntry, __format_ObservableEntry)
        html_formatter.for_type(GoodnessOfFit, __format_GoodnessOfFit)
except NameError as e:
    pass
