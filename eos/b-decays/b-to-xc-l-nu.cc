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

#include <eos/b-decays/b-to-xc-l-nu.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    template <>
    struct Implementation<BToXcLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter g_fermi;

        UsedParameter hbar;

        //UsedParameter mu_b_kin;
        //UsedParameter mu_c_MSbar;

        UsedParameter _mu2pi;
        UsedParameter _mu2G;

        SwitchOption opt_l;

        SwitchOption opt_q;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            m_B(p["mass::B_" + o.get("q", "d")], u),
            tau_B(p["life_time::B_" + o.get("q", "d")], u),
            g_fermi(p["G_Fermi"], u),
            hbar(p["hbar"], u),
            //mu_b_kin(p["B->X_clnu::mu_b@Kinetic"], u),
            //mu_c_MSbar(p["B->X_clnu::mu_c@MSbar"], u),
            _mu2pi(p["B->B::mu_pi^2@1GeV"], u),
            _mu2G(p["B->B::mu_G^2@1GeV"], u),
            opt_l(o, "l", {"e", "mu"}, "mu"),
            opt_q(o, "q", {"u", "d"}, "d")
        {
            u.uses(*model);
        }

        inline double m_b() const { return 4.2; } // TODO: kinetic scheme: model->m_b_kinetic(...);
        inline double m_c() const { return 1.5; } // TODO: what scheme?

        inline double mu2pi() const
        {
            return _mu2pi(); // TODO: RGE running to scale mu!
        }

        inline double mu2G() const
        {
            return _mu2G(); // TODO: RGE running to scale mu!
        }

        // {{{ integrated moments

        double integrated_moment_0_0_0(const double & E_l_min, const double & E_l_max,const double & q_0_min, const double & q_0_max,const double & M_X_min, const double & M_X_max) const
        {
            const double a_s    = model->alpha_s(4.2) / M_PI;
            const double m_b    = this->m_b();
            const double m_c    = this->m_c();
            const double r      = m_b * m_b / (m_c * m_c);
            const double mu2pi  = this->mu2pi();
            const double mu2G   = this->mu2G();

            // {{{ i = V_L, j = V_L
            const double z0_VL_VL        = 1.0 - 8.0 * r * 8.0 * pow(r, 3) - pow(r, 4) - 12.0 * r * r * log(r);
            const double A3_VL_VL        = -0.94 * a_s;
            const double A5_VL_VL        = 0.0; // ???
            const double res_VL_VL_0     = z0_VL_VL * (1.0 + A3_VL_VL);
            const double res_VL_VL_mu2pi = -0.5 * z0_VL_VL * (1.0 + A3_VL_VL);
            const double res_VL_VL_mu2G  = +0.5 * z0_VL_VL * (1.0 + A3_VL_VL)
                                           -2.0 * pow(1.0 - r, 4) * (1.0 + A5_VL_VL);

            const double res_VL_VL =
                      res_VL_VL_0
                    + res_VL_VL_mu2pi * mu2pi
                    + res_VL_VL_mu2G * mu2G;
            // }}}

            const auto wc = model->wilson_coefficients_b_to_c(opt_l.value(), false);
            const double result =
                      norm(wc.cvl()) * res_VL_VL
                    + 0.0;

            return result;
        }

        // }}}

        inline double normalization()
        {
            const double result = power_of<2>(g_fermi() * abs(model->ckm_cb()))
                    / (192.0 * power_of<3>(M_PI)) * pow(m_b(), 5) * tau_B() / hbar;

            return result;
        }
    };

    BToXcLeptonNeutrino::BToXcLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToXcLeptonNeutrino>(new Implementation<BToXcLeptonNeutrino>(parameters, options, *this))
    {
    }

    BToXcLeptonNeutrino::~BToXcLeptonNeutrino() = default;

    double
    BToXcLeptonNeutrino::integrated_branching_ratio(const double & E_l_min, const double & E_l_max,const double & q_0_min, const double & q_0_max, const double & M_X_min, const double & M_X_max) const
    {
        return _imp->integrated_moment_0_0_0(E_l_min, E_l_max, q_0_min, q_0_max, M_X_min, M_X_max) * _imp->normalization();
    }
}
