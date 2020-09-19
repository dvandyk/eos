/* vim: set sw=4 sts=4 et tw=120 foldmethod=syntax : */

/*
 * Copyright (c) 2020 Danny van Dyk
 * Copyright (c) 2020 Nico Gubernari
 * Copyright (c) 2020 Christoph Bobeth
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BGL1997_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BGL1997_HH 1

#include <eos/form-factors/form-factors-fwd.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>

namespace eos
{
    template <typename Process_> class BGL1997FormFactors;

    template <>
    class BGL1997FormFactors<BToDstar> :
        public FormFactors<PToV>
    {
        private:
            std::array<UsedParameter, 4> _a_g, _a_f, _a_F1, _a_F2;
            UsedParameter _t_0;

            const double _mB, _mB2, _mV, _mV2;
            const double _t_p, _t_m;

            static constexpr double chi_1m = 1.0;
            static constexpr double chi_1p = 1.0;

            double _z(const double & t, const double & t_0) const;

            double _phi(const double & s, const double & K,  const unsigned & a, const unsigned & b, const unsigned & c, const double & chi) const;

            static std::string _par_name(const std::string & ff_name);

        public:
            BGL1997FormFactors(const Parameters &, const Options &);

            ~BGL1997FormFactors();

            static FormFactors<PToV> * make(const Parameters & parameters, const Options & options);

            double g(const double & s) const;

            double f(const double & s) const;

            virtual double v(const double & s) const;

            virtual double a_0(const double & s) const;

            virtual double a_1(const double & s) const;

            virtual double a_2(const double & s) const;

            virtual double a_12(const double & s) const;

            virtual double t_1(const double & s) const;

            virtual double t_2(const double & s) const;

            virtual double t_3(const double & s) const;

            virtual double t_23(const double & s) const;
    };
}

#endif
