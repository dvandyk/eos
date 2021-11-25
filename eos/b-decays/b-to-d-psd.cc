/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Danny van Dyk
 * Copyright (c) 2021 Stefan Meiser
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

#include <eos/b-decays/b-to-d-psd.hh>
#include <eos/form-factors/k-lcdas.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/polylog.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/options-impl.hh>

#include <string>

namespace eos
{
    using std::norm;
    using std::real;
    using std::imag;

    /*
     * Decay: B_q -> D_q P, cf. [BBNS:2000A]
     */
    template <>
    struct Implementation<BToDPseudoscalar>
    {
        std::shared_ptr<Model> model;

        SwitchOption opt_q;

        UsedParameter hbar;

        UsedParameter g_fermi;

        UsedParameter m_b;

        UsedParameter m_c;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter m_D;

        UsedParameter m_P;

        UsedParameter f_P;

        GSL::QAGS::Config int_config;

        bool cp_conjugate;

        UsedParameter mu;

        KaonLCDAs lcdas;

        std::shared_ptr<FormFactors<PToP>> form_factors;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            opt_q(o, "q", { "d", "s" }),
            hbar(p["QM::hbar"], u),
            g_fermi(p["WET::G_Fermi"], u),
            m_b(p["mass::b(MSbar)"], u),
            m_c(p["mass::c"], u),
            m_B(p["mass::B_" + opt_q.value()], u),
            tau_B(p["life_time::B_" + opt_q.value()], u),
            m_D(p["mass::D_" + opt_q.value()], u),
            m_P(p["mass::K_u"], u), // TODO: adapt to the pi case
            f_P(p["decay-constant::K_u"], u), // TODO: adapt to the pi case
            cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false"))),
            mu(p["cbsu::mu"], u),
            lcdas(p, o)
        {
            form_factors = FormFactorFactory<PToP>::create("B->D::" + o.get("form-factors", "BSZ2015"), p, o);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");
        }

        // building block of the one-loop functions
        inline static complex<double> fV(const double & u, const double & z)
        {
            // Auxiliary quantities
            double ubar = 1.0 - u;
            double A    = u    * (1.0 - z * z);
            double Abar = ubar * (1.0 - z * z);

            return -A * (3.0 * (1.0 - A) + z) / power_of<2>(1.0 - A) * std::log(A) - z / (1.0 - A)
                    + 2.0 * (
                        +std::log(A)    / (1.0 - A)    - power_of<2>(A)    - dilog(1.0 - A)
                        -std::log(Abar) / (1.0 - Abar) + power_of<2>(Abar) + dilog(1.0 - Abar)
                    );
        }

        inline static complex<double> FVLL(const double & u, const double & z)
        {
            return (3.0 + 2.0 * std::log(u / (1.0 - u))) * std::log(z * z) - 7.0 + fV(u, z) + fV(1 - u, 1.0 / z);
        }

        complex<double> a1() const
        {
            // alpha_s at low scale
            const double alpha_s = model->alpha_s(mu);

            // Auxiliary quantities
            const auto wc = model->wet_cbsu(cp_conjugate);
            const double Nc = 3.0;
            const double B = 11.0;   //NDR
            const double B_p = +(Nc - 1.0) / (2.0 * Nc) * B;
            const double B_m = -(Nc + 1.0) / (2.0 * Nc) * B;
            const double Cf = (Nc * Nc - 1.0) / (2.0 * Nc);
            const double z = m_c / m_b;

            // Wilson coefficients
            std::array<complex<double>, 20> wc_MvD = wc.convertToMvD();         //Convert from Bern Basis to our basis
            const complex<double> C_p = wc_MvD[4] + wc_MvD[14] / 3.0;                    //Convert to Wilson coefficients as defined in [BBNS:2000A] eqs. (45) - (47)
            const complex<double> C_m = wc_MvD[4] - 2.0 * wc_MvD[14] / 3.0;
            const complex<double> barC_p = (1.0 - alpha_s * B_p / (4.0 * M_PI)) * C_p;
            const complex<double> barC_m = (1.0 - alpha_s * B_m / (4.0 * M_PI)) * C_m;

            // Twist 2, LO contribution
            complex<double> tw2_LO = (Nc + 1.0) / (2.0 * Nc) * barC_p + (Nc - 1.0) / (2.0  * Nc) * barC_m;

            // Twist 2, NLO contribution: requires convolution of LCDA and loop function
            const auto & tw2_NLO_integrand_re = [&] (const double & u) -> double {
                return real(wc_MvD[14] * FVLL(u, z)) * lcdas.phi(u, mu);
            };
            const auto & tw2_NLO_integrand_im = [&] (const double & u) -> double {
                return imag(wc_MvD[14] * FVLL(u, z)) * lcdas.phi(u, mu);
            };
            complex<double> tw2_NLO_integral = complex<double>{
                integrate<GSL::QAGS>(tw2_NLO_integrand_re, 0.0, 1.0, this->int_config),
                integrate<GSL::QAGS>(tw2_NLO_integrand_im, 0.0, 1.0, this->int_config)
            };
            complex<double> tw2_NLO = alpha_s * Cf  / (8.0 * M_PI * Nc)
                    * (
                        -1.0 * std::log(mu * mu / (m_b * m_b))
                        + tw2_NLO_integral
                    );

            return tw2_LO + tw2_NLO;
        }

        double decay_width() const
        {
            // masses
            const double m_B = this->m_B(), m_B2 = m_B * m_B;
            const double m_D = this->m_D(), m_D2 = m_D * m_D;
            const double m_P = this->m_P(), m_P2 = m_P * m_P;

            // auxiliary quantities
            const double F0 = form_factors->f_0(m_P2);
            const double lambda = eos::lambda(m_B2, m_D2, m_P2);
            const double qvec = std::sqrt(lambda) / (2.0 * m_B);
            const double norm = power_of<2>(g_fermi * std::abs(model->ckm_cb() * std::conj(model->ckm_us())) * (m_B2 - m_D2))
                    * qvec / (16.0 * M_PI * m_B2);

            return norm * std::norm(this->a1()) * power_of<2>(f_P * F0);
        }

        double branching_ratio() const
        {
            return decay_width() * tau_B / hbar;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BToDPseudoscalar>::options
    {
	    { "q", { "d", "s" }, ""}
    };

    BToDPseudoscalar::BToDPseudoscalar(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToDPseudoscalar>(new Implementation<BToDPseudoscalar>(parameters, options, *this))
    {
    }

    BToDPseudoscalar::~BToDPseudoscalar()
    {
    }

    double
    BToDPseudoscalar::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    const std::set<ReferenceName>
    BToDPseudoscalar::references
    {
        "BBNS:2000A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    BToDPseudoscalar::begin_options()
    {
        return Implementation<BToDPseudoscalar>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToDPseudoscalar::end_options()
    {
        return Implementation<BToDPseudoscalar>::options.cend();
    }
}
