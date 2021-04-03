/* vim: set sw=4 sts=4 et tw=140 foldmethod=syntax : */

/*
 * Copyright (c) 2020 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_FORM_FACTORS_PARAMETRIC_DGKVD2020_HH
#define EOS_GUARD_SRC_FORM_FACTORS_PARAMETRIC_DGKVD2020_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>
#include <eos/utils/qualified-name.hh>

namespace eos
{
    class DGKvD2020FormFactors :
        public FormFactors<VacuumToPGamma>
    {
        private:
            Parameter m_pi;
            Parameter f_pi;
            Parameter m_rho;
            Parameter f_rho;
            Parameter m_a1;
            Parameter f_a1;
            Parameter eF_pi_rho;
            Parameter eF_g;
            Parameter q_02;

        public:
            DGKvD2020FormFactors(const Parameters &, const Options &);

            static FormFactors<VacuumToPGamma> * make(const Parameters &, const Options &);

            virtual complex<double> a(const double & q2) const;
            virtual complex<double> v(const double & q2) const;
    };
}

#endif
