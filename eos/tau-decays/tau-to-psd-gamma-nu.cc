/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020 Danny van Dyk
 * Copyright (c) 2020 Katarina Dugic
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

#include <eos/form-factors/mesonic.hh>
#include <eos/tau-decays/tau-to-psd-gamma-nu.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    using std::norm;

    /*
     * Decay: tau -> P gamma nu, cf. [DGvD:2020A]
     */
    template <>
    struct Implementation<TauToPseudoscalarGammaNeutrino>
    {
        std::shared_ptr<Model> model;

        UsedParameter hbar;

        UsedParameter g_fermi;

        UsedParameter m_tau;

        UsedParameter m_pi;

        UsedParameter m_pi0;

        UsedParameter f_pi;

        UsedParameter tau_tau;

        SwitchOption opt_P;

        /*UsedParameter f_P;

        UsedParameter m_P;*/

        std::shared_ptr<FormFactors<VacuumToPGamma>> ff;


        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            hbar(p["hbar"], u),
            g_fermi(p["G_Fermi"], u),
            m_tau(p["mass::tau"], u),
            tau_tau(p["life_time::tau"], u),
            opt_P(o, "P", {"pi", "K"}, "pi"),
            /*f_P(p["decay-constant::" + opt_P.value() + "^-"], u),
            m_P(p["mass::" + opt_P.value() + "^-"], u),*/
            f_pi(p["decay-constant::pi"], u),
            m_pi0(p["mass::pi^0"], u),
            m_pi(p["mass::pi^+"], u),
            ff(FormFactorFactory<VacuumToPGamma>::create("0->" + opt_P.value() + "gamma::DGKvD2020", p, o))
        {
            u.uses(*model);
            u.uses(*ff);
        }

/*TODO -use V_us if P = K
       -replace delta
       -define C,D
       -apply form factor correction
       -alpha_e at tau scale!
       -check jacobi determinant xy -> E_gamma E_pi
*/
        double differential_ratio(const double & E_gamma, const double & E_pi) const
        {
            const double q2 = 2.0 * E_gamma * m_tau + 2.0 * E_pi * m_tau - m_tau * m_tau;
            const double delta = pow( m_pi / m_tau, 2);
            const double x = 2.0 * E_gamma / m_tau;
            const double y = 2.0 * E_pi / m_tau;
            const double F_BB = -2.0 * delta * (2.0 - x - y) / (pow(x + y - 1.0 - delta, 2)) - 2.0 * (1.0 - x - delta) / (pow(x, 2)) + 2.0 * y * (2.0 - x -y) / (x * (x + y - 1.0 - delta)) 
                                - (1.0 - delta - x) / x - y * (1.0 + delta -y) / (x * (x + y - 1.0 - delta)) + (1.0 - delta - x) / (x + y - 1.0 - delta);
            const double F_VV = (x + y - 1.0 - delta) * (x * (1.0 - delta - x) + y * (1.0 + delta - y)) - 2.0 * delta * x * (1.0 + delta - y);
            const double F_AA = (x + y - 1.0 - delta) * (x * (1.0 - delta - x) + y * (1.0 + delta - y)) - 2.0 * delta * x * (1.0 + delta - y);
            const double F_VB = (x + y - 1.0 - delta) * (1.0 + delta - y) / x;
            const double F_AB = y * (1.0 + delta - y) / x - (x + y - 1.0 - delta) * (1.0 - delta - x) / x + (1.0 - delta - x) - 2.0 * delta * (1.0 + delta - y) / (x + y - 1.0 - delta);
            const double F_AV = (x + y - 1.0 - delta) * (x * (1.0 - delta - x) - y * (1.0 + delta - y));
            const complex<double> v_q2 = ff->v(q2) / m_pi();
            const complex<double> a_q2 = ff->a(q2) * 2.0 * m_pi() / (q2 - (m_pi * m_pi));
            const double alpha_e = 1.0 / 137.0; 
            /*const double decay_width_pi0 = 7.73e-9;*/
            const double C = pow(m_tau, 4) * alpha_e / (8.0 * M_PI * pow(f_pi * (1.0 - delta), 2));
            const double D = pow(m_tau, 2) * alpha_e / (4.0 * M_PI * f_pi * pow(1.0 - delta, 2));
            
            return alpha_e / ( 2.0 * M_PI * (1.0-delta)*(1.0-delta)) * F_BB 
                   + C * std::norm(v_q2) * F_VV 
                   + C * std::norm(a_q2) * F_AA 
                   - D * 2.0 * std::real(v_q2) * F_VB 
                   + D * 2.0 * std::real(a_q2) * F_AB 
                   + C * 2.0 * std::real(a_q2 * v_q2) * F_AV;
        }

  /*      double branching_ratio(const double & q2) const
        {
            return decay_width(q2) * tau_tau / hbar;
        }*/
    };

    TauToPseudoscalarGammaNeutrino::TauToPseudoscalarGammaNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<TauToPseudoscalarGammaNeutrino>(new Implementation<TauToPseudoscalarGammaNeutrino>(parameters, options, *this))
    {
    }

    TauToPseudoscalarGammaNeutrino::~TauToPseudoscalarGammaNeutrino()
    {
    }

   /* double
    TauToPseudoscalarGammaNeutrino::branching_ratio(const double & q2) const
    {
        return _imp->branching_ratio(q2);
    }
*/
    double
    TauToPseudoscalarGammaNeutrino::differential_ratio(const double & E_gamma, const double & E_pi) const
    {
        // return _imp->differential_ratio(E_gamma, E_pi);
        return std::exp(abs(E_gamma-E_pi));
    }

    double
    TauToPseudoscalarGammaNeutrino::dummy(const double & ) const
    {
        return 1.0;
    }

    const std::set<ReferenceName>
    TauToPseudoscalarGammaNeutrino::references
    {
        "DGvD:2020A"_rn
    };

    const std::string
    TauToPseudoscalarGammaNeutrino::description = "Decay tau -> P gamma neutrino with P = {pi, K}";

    const std::string
    TauToPseudoscalarGammaNeutrino::kinematics_description_E_gamma = "Energy of the photon in the tau rest frame";

    const std::string
    TauToPseudoscalarGammaNeutrino::kinematics_description_E_pi = "Energy of the pion in the tau rest frame";
    
}
