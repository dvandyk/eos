/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Danny van Dyk
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef EOS_GUARD_EOS_FORM_FACTORS_ANALYTIC_B_TO_GAMMA_HH
#define EOS_GUARD_EOS_FORM_FACTORS_ANALYTIC_B_TO_GAMMA_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>

namespace eos
{
    /*!
     * Implementation of the B -> gamma form factors within the framework of
     * QCD Factorization at the NLO level and with NLL accuracy.
     *
     * The formuals are taken from [BF2011].
     */
    class AnalyticFormFactorBToGammaQCDF :
        public FormFactors<PToGamma>,
        PrivateImplementationPattern<AnalyticFormFactorBToGammaQCDF>
    {
        public:
            AnalyticFormFactorBToGammaQCDF(const Parameters &, const Options &);

            ~AnalyticFormFactorBToGammaQCDF();

            static FormFactors<PToGamma> * make(const Parameters &, const Options & );

            /* Form factors */
            virtual double f_v(const double & Egamma) const;
            virtual double f_a(const double & Egamma) const;

            /* Diagnostics for unit tests */
            Diagnostics diagnostics() const;
    };
}

#endif
