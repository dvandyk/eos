/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016 Danny van Dyk
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

#include <eos/k-decays/k-to-l-nu.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/polylog.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    using std::norm;

    /*
     * Decay: K^- -> l^- nubar, cf. [CR:2007A]
     */
    template <>
    struct Implementation<KToLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        UsedParameter hbar;

        UsedParameter g_fermi;

        UsedParameter alpha_e_0;

        UsedParameter mu;

        UsedParameter m_Z;

        UsedParameter m_K;

        UsedParameter f_K;

        UsedParameter tau_K;

        SwitchOption opt_l;

        UsedParameter m_l;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            hbar(p["hbar"], u),
            g_fermi(p["G_Fermi"], u),
            alpha_e_0(p["QED::alpha_e(0)"], u),
            mu(p["K->lnu::mu@ChPT"], u),
            m_Z(p["mass::Z"], u),
            m_K(p["mass::K_u"], u),
            f_K(p["decay-constant::K_u"], u),
            tau_K(p["life_time::K_u"], u),
            opt_l(o, "l", { "e", "mu" }),
            m_l(p["mass::" + opt_l.value()], u)
        {
            u.uses(*model);
        }

        /*
         * Decay width to lowest order in ChPT, cf. [CR:2007A], eq. (5), p. 4.
         */
        inline double decay_width_zero() const
        {
            return power_of<2>(g_fermi * std::abs(model->ckm_us()) * f_K * m_l * (1.0 - power_of<2>(m_l() / m_K())))
                * m_K / (4.0 * M_PI);
        }

        /*
         * Decay width for the decay K -> l nu (gamma), including corrections to order e^2 p^4
         * as shown in [CR:2007A], eq. (114), p. 26.
         */
        double decay_width(const double & omega) const
        {
            using std::log;

            // cf. [CF:2007A], eq. (114), p. 26
            const double gamma_0   = decay_width_zero();
            const double alpha     = alpha_e_0;
            const double sirlin    = 1.0 + alpha * 2.0 / M_PI * log(m_Z() / mu());
            const double z_l       = m_l * m_l / m_K * m_K;
            const double log_z_l   = log(z_l);
            const double li2_1mz_l = std::real(dilog(complex<double>(1.0 - z_l, 0.0)));
            const double rho2      = pow(mu / m_K, 2); // [CR:2007A] uses mu = m_rho; we leave them uncorrelated
            const double log_rho2  = log(rho2);

            // cf. [CR:2007A], eq. (118), p. 27
            const double F_soft    = - 3.0 / 4.0 + 3.0 / 4.0 * log_z_l - 2.0 * z_l / (1.0 - z_l) * log_z_l
                                   - (1.0 + z_l) / (1.0 - z_l) * li2_1mz_l
                                   - (2.0 + (1.0 + z_l) / (1.0 - z_l) * log_z_l) * log(2.0 * omega / m_K);

            // cf. [CR:2007A], eq. (116), p. 27
            const double c_1       = -1.98;
            // cf. [CR:2007A], table 1, p. 23
            const double ctilde_2  = +0.0784;
            const double c_2       = +4.3;
            const double c_3       = -4.73;
            const double c_4       = +0.22;

            // cf. [CF:2007A], eq. (114), p. 26
            const double chi_pt    = 3.0 / 4.0 * log_rho2 + c_1
                                   + z_l / rho2 * (c_2 * (log_rho2 - log_z_l) + c_3 + c_4)
                                   - 1.0 / rho2 * ctilde_2 * (log_rho2 - log_z_l);

            return gamma_0 * sirlin * (1 + alpha / M_PI * F_soft) * (1 - alpha / M_PI * chi_pt);
        }

        double branching_ratio(const double & omega) const
        {
            return decay_width(omega) * tau_K / hbar;
        }
    };

    KToLeptonNeutrino::KToLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<KToLeptonNeutrino>(new Implementation<KToLeptonNeutrino>(parameters, options, *this))
    {
    }

    KToLeptonNeutrino::~KToLeptonNeutrino()
    {
    }

    double
    KToLeptonNeutrino::branching_ratio(const double & omega) const
    {
        return _imp->branching_ratio(omega);
    }

    double
    KToLeptonNeutrino::decay_width(const double & omega) const
    {
        return _imp->decay_width(omega);
    }
}
