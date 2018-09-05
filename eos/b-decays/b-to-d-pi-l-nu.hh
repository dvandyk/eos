/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018 Danny van Dyk
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

#ifndef MASTER_GUARD_EOS_B_DECAYS_B_TO_D_PI_L_NU_HH
#define MASTER_GUARD_EOS_B_DECAYS_B_TO_D_PI_L_NU_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    class BToDPiLeptonNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<BToDPiLeptonNeutrino>
    {
        public:
            BToDPiLeptonNeutrino(const Parameters & parameters, const Options & options);
            ~BToDPiLeptonNeutrino();

            /*!
             * 1-dim. PDFs as functions of cos(theta_D), cost(theta_L), and phi
             */
            double differential_pdf_d(const double & c_d) const;
            double differential_pdf_l(const double & c_l) const;
            double differential_pdf_phi(const double & phi) const;

            /*!
             * Partially-integrated 1-dim. PDFs for cos(theta_D), cost(theta_L), and phi
             */
            double integrated_pdf_d(const double & c_d_min, const double & c_d_max) const;
            double integrated_pdf_l(const double & c_l_min, const double & c_l_max) const;
            double integrated_pdf_phi(const double & phi_min, const double & phi_max) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_c_d;
            static const std::string kinematics_description_c_l;
            static const std::string kinematics_description_phi;
    };
}

#endif
