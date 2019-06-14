/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_B_DECAYS_B_TO_X_C_L_NU_HH
#define EOS_GUARD_EOS_B_DECAYS_B_TO_X_C_L_NU_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    /*
     * Decay: B -> X_c l nu
     */
    class BToXcLeptonNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<BToXcLeptonNeutrino>
    {
        public:
            BToXcLeptonNeutrino(const Parameters & parameters, const Options & options);
            ~BToXcLeptonNeutrino();

            // Integrated Moments
            // k_l_m: E_L^k q_0^l M_X^m
            double integrated_moment_0_0_0(const double & E_l_min, const double & E_l_max,const double & q_0_min, const double & q_0_max,const double & M_X_min, const double & M_X_max) const;

            double integrated_branching_ratio(const double & E_l_min, const double & E_l_max,const double & q_0_min, const double & q_0_max,const double & M_X_min, const double & M_X_max) const;
    };
}

#endif
