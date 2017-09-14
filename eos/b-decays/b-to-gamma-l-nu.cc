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

#include <eos/form-factors/mesonic.hh>
#include <eos/b-decays/b-to-gamma-l-nu.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    template <>
    struct Implementation<BToGammaLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        std::shared_ptr<FormFactors<PToGamma>> form_factors;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter m_l;

        UsedParameter alpha_e;

        UsedParameter g_fermi;

        UsedParameter hbar;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            m_B(p["mass::B_u"], u),
            tau_B(p["life_time::B_u"], u),
            m_l(p["mass::" + o.get("l", "mu")], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            g_fermi(p["G_Fermi"], u),
        hbar(p["hbar"], u)
    {
            form_factors = FormFactorFactory<PToGamma>::create("B->gamma@" + o.get("form-factors", "QCDF"), p, o);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            u.uses(*form_factors);
            u.uses(*model);
        }

        double normalized_differential_decay_width(const double & Egamma) const
        {
            const double xgamma = 2.0 * Egamma / m_B;
            const double f_a = form_factors->f_a(Egamma);
            const double f_v = form_factors->f_v(Egamma);

            return (1.0 - xgamma) * pow(xgamma, 3) * (f_a * f_a + f_v * f_v);
        }

        double normalized_forward_backward_asymmetry(const double & Egamma) const
        {
            const double xgamma = 2.0 * Egamma / m_B;
            const double f_a = form_factors->f_a(Egamma);
            const double f_v = form_factors->f_v(Egamma);

            return 0.5 * (1.0 - xgamma) * pow(xgamma, 3) * f_a * f_v;
        }

        double integrand_moments(const double Egamma, const unsigned idx) const
        {
            const double xgamma = 2.0 * Egamma / m_B;

            return pow(xgamma, idx) * normalized_differential_decay_width(Egamma);
        }

        double differential_branching_ratio(const double & s) const
        {
            const double norm = alpha_e * pow(g_fermi() * abs(model->ckm_ub()), 2)
                              * pow(m_B, 4) / (48.0 * pow(M_PI, 2)) * tau_B / hbar;

            return normalized_differential_decay_width(s) * norm;
        }
    };

    BToGammaLeptonNeutrino::BToGammaLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToGammaLeptonNeutrino>(new Implementation<BToGammaLeptonNeutrino>(parameters, options, *this))
    {
    }

    BToGammaLeptonNeutrino::~BToGammaLeptonNeutrino()
    {
    }

    double
    BToGammaLeptonNeutrino::differential_branching_ratio(const double & Egamma) const
    {
        return _imp->differential_branching_ratio(Egamma);
    }

    double
    BToGammaLeptonNeutrino::integrated_branching_ratio(const double & Egamma_min, const double & Egamma_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToGammaLeptonNeutrino>::differential_branching_ratio,
                _imp.get(), std::placeholders::_1);

        return integrate(f, 128, Egamma_min, Egamma_max);
    }

    double
    BToGammaLeptonNeutrino::integrated_forward_backward_asymmetry(const double & Egamma_min, const double & Egamma_max) const
    {
        std::function<double (const double &)> num   = std::bind(&Implementation<BToGammaLeptonNeutrino>::normalized_forward_backward_asymmetry,
                _imp.get(), std::placeholders::_1);

        std::function<double (const double &)> denom = std::bind(&Implementation<BToGammaLeptonNeutrino>::normalized_differential_decay_width,
                _imp.get(), std::placeholders::_1);

        return integrate(num, 128, Egamma_min, Egamma_max) / integrate(denom, 128, Egamma_min, Egamma_max);
    }

    double
    BToGammaLeptonNeutrino::integrated_photon_energy_moment_1(const double & Egamma_min, const double & Egamma_max) const
    {
        std::function<double (const double &)> num   = std::bind(&Implementation<BToGammaLeptonNeutrino>::integrand_moments,
                _imp.get(), std::placeholders::_1, 1u);

        std::function<double (const double &)> denom = std::bind(&Implementation<BToGammaLeptonNeutrino>::normalized_differential_decay_width,
                _imp.get(), std::placeholders::_1);

        return integrate(num, 128, Egamma_min, Egamma_max) / integrate(denom, 128, Egamma_min, Egamma_max);
    }

    double
    BToGammaLeptonNeutrino::integrated_photon_energy_moment_2(const double & Egamma_min, const double & Egamma_max) const
    {
        std::function<double (const double &)> num   = std::bind(&Implementation<BToGammaLeptonNeutrino>::integrand_moments,
                _imp.get(), std::placeholders::_1, 2u);

        std::function<double (const double &)> denom = std::bind(&Implementation<BToGammaLeptonNeutrino>::normalized_differential_decay_width,
                _imp.get(), std::placeholders::_1);

        return integrate(num, 128, Egamma_min, Egamma_max) / integrate(denom, 128, Egamma_min, Egamma_max);
    }

    const std::string
    BToGammaLeptonNeutrino::description = "\
The decay B->gamma l nu, where l=e,mu is a lepton.";

    const std::string
    BToGammaLeptonNeutrino::kinematics_description_Egamma = "\
The energy of the photon in the B-meson rest frame.";

}
