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

#include <eos/form-factors/parametric-dgkvd2020.hh>

namespace eos
{
    DGKvD2020FormFactors::DGKvD2020FormFactors(const Parameters & p, const Options &) :
        m_pi(p["mass::pi^+"]),
        f_pi(p["decay-constant::pi"]),
        m_rho(p["mass::rho^+"]),
        f_rho(p["decay-constant::rho"]),
        m_a1(p["mass::a_1(1260)^+"]),
        f_a1(p["decay-constant::a_1"]),
        eF_pi_rho(p["rho->pi::eV(0)"]),
        eF_g(p["a_1->pi::eF_g(0)"]),
        q_02(p["0->pigamma::q_0^2"])

    {

    }

    FormFactors<VacuumToPGamma> *
    DGKvD2020FormFactors::make(const Parameters & p, const Options & o)
    {
        return new DGKvD2020FormFactors(p,o);
    }

    complex<double>
    DGKvD2020FormFactors::a(const double & q2) const
    {
        const double F_pi=1.0; // electromagnetic pion form factor is normalized to 1.0 for an on-shell photon    
        const double a_q0=1.0; // axial vector form factor at subtraction point  

        return a_q0 + (q2-q_02)*f_a1/(m_pi*(m_a1*m_a1-q2)*(m_a1*m_a1-q_02))*
                      (eF_g/(1.602177e-19)*m_a1*m_a1);
    }

    complex<double>
    DGKvD2020FormFactors::v(const double & q2) const
    {
        const double v_q0=1.0; // vector form factor at subtraction point
     
        return v_q0 + (q2-q_02)*f_rho*m_pi*eF_pi_rho/((m_rho*m_rho-q2)*(m_rho*m_rho-q_02)*(1.602177e-19));
    }
}
