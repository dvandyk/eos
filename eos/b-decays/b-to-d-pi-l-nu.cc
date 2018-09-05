/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018, 2019 Danny van Dyk
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

#include <eos/b-decays/b-to-d-pi-l-nu.hh>
#include <eos/b-decays/b-to-v-l-nu.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/utils/integrate-impl.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    using namespace eos::btovlnu;
    using std::conj;
    using std::norm;
    using std::real;
    using std::sqrt;

    template <>
    struct Implementation<BToDPiLeptonNeutrino>
    {
        // model
        SwitchOption opt_model;

        // meson masses
        UsedParameter m_B, m_Dstar;

        // lepton mass
        SwitchOption opt_l;
        UsedParameter m_l;

        // form factors
        std::shared_ptr<FormFactors<PToV>> ff;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_model(o, "model", { "SM", "CKMScan" }, "SM"),
            m_B(p["mass::B_d"], u),
            m_Dstar(p["mass::D^*_d"], u),
            opt_l(o, "l", { "e", "mu", "tau" }, "mu"),
            m_l(p["mass::" + opt_l.value()], u),
            ff(FormFactorFactory<PToV>::create("B->D^*::" + o.get("form-factors", "HQET"), p, o))
        {
            if (! ff.get())
                throw InternalError("Form factors not found!");

            u.uses(*ff);
        }

        inline double pdf_normalization(const double & q2) const
        {
            // we only require q2-dependent terms of the normalization, since
            // the constants parts drop out in the PDF.

            const double m_B     = this->m_B(),     m_B2 = m_B * m_B;
            const double m_Dstar = this->m_Dstar(), m_Dstar2 = m_Dstar * m_Dstar;

            const double pDstar = sqrt(eos::lambda(m_B2, m_Dstar2, q2)) / (2.0 * m_B);
            const double beta   = 1.0 - m_l() * m_l() / q2;

            //const double N = 1.0 / (192.0 * pow(M_PI, 3) * m_B2);

            return pDstar * q2 * beta * beta;// * N;
        }

        inline double a_long(const double & q2) const
        {
            const double m_B     = this->m_B(),     m_B2 = m_B * m_B;
            const double m_Dstar = this->m_Dstar(), m_Dstar2 = m_Dstar * m_Dstar;
            const double lambda  = eos::lambda(m_B2, m_Dstar2, q2);

            double result = (m_B + m_Dstar) * (m_B2 - m_Dstar2 - q2) * ff->a_1(q2);
            result -= lambda / (m_B + m_Dstar) * ff->a_2(q2);
            result /= (2.0 * m_Dstar * sqrt(q2));

            // cf. [DDS2014], eq. (22), p. 17
            return result;
        }

        inline double a_para(const double & q2) const
        {
            const double m_B     = this->m_B();
            const double m_Dstar = this->m_Dstar();

            // cf. [DDS2014], eq. (22), p. 17
            return sqrt(2.0) * (m_B + m_Dstar) * ff->a_1(q2);
        }

        inline double a_perp(const double & q2) const
        {
            const double m_B     = this->m_B(),     m_B2 = m_B * m_B;
            const double m_Dstar = this->m_Dstar(), m_Dstar2 = m_Dstar * m_Dstar;
            const double lambda  = eos::lambda(m_B2, m_Dstar2, q2);

            // cf. [DDS2014], eq. (22), p. 17
            return -sqrt(2.0) * sqrt(lambda) / (m_B + m_Dstar) * ff->v(q2);
        }

        inline double a_time(const double & q2) const
        {
            const double m_B     = this->m_B(),     m_B2 = m_B * m_B;
            const double m_Dstar = this->m_Dstar(), m_Dstar2 = m_Dstar * m_Dstar;
            const double lambda  = eos::lambda(m_B2, m_Dstar2, q2);

            // cf. [DDS2014], eq. (22), p. 17
            return sqrt(lambda / q2) * ff->a_0(q2);
        }

        double lepton_polarization_numerator(const double & q2) const
        {
            const double nf = pdf_normalization(q2);

            const double m_l2     = m_l() * m_l();
            const double m_B      = this->m_B(),     m_B2     = m_B     * m_B;
            const double m_Dstar  = this->m_Dstar(), m_Dstar2 = m_Dstar * m_Dstar;
            const double p_Dstar  = sqrt(eos::lambda(m_B2, m_Dstar2, q2)) / (2.0 * m_B),
                         p_Dstar2 = p_Dstar * p_Dstar;
            const double sqrt_q2  = sqrt(q2);

            const double a_0 = ff->a_0(q2);
            const double a_1 = ff->a_1(q2);
            const double a_2 = ff->a_2(q2);
            const double v   = ff->v(q2);

            // cf. [CJLP2012]
            const double H_pp = (m_B + m_Dstar) * a_1 - 2.0 * m_B / (m_B + m_Dstar) * p_Dstar * v;
            const double H_mm = (m_B + m_Dstar) * a_1 + 2.0 * m_B / (m_B + m_Dstar) * p_Dstar * v;
            const double H_00 = ((m_B2 - m_Dstar2 - q2) * (m_B + m_Dstar) * a_1 - 4.0 * m_B2 * p_Dstar2 * a_2 / (m_B + m_Dstar))
                              / (2.0 * m_Dstar * sqrt_q2);
            const double H_0t = 2.0 * m_B * p_Dstar / sqrt_q2 * a_0;

            const double H_pp2 = H_pp * H_pp;
            const double H_mm2 = H_mm * H_mm;
            const double H_002 = H_00 * H_00;
            const double H_0t2 = H_0t * H_0t;

            // cf. [CJLP2012]], eq. (22), p. 17
            const double num   = (H_pp2 + H_mm2 + H_002) * (1.0 - m_l2 / (2.0 * q2)) - 3.0 * m_l2 / (2.0 * q2) * H_0t2;

            return nf * num;
        }

        double lepton_polarization_denominator(const double & q2) const
        {
            const double nf = pdf_normalization(q2);

            const double m_l2     = m_l() * m_l();
            const double m_B      = this->m_B(),     m_B2     = m_B     * m_B;
            const double m_Dstar  = this->m_Dstar(), m_Dstar2 = m_Dstar * m_Dstar;
            const double p_Dstar  = sqrt(eos::lambda(m_B2, m_Dstar2, q2)) / (2.0 * m_B),
                         p_Dstar2 = p_Dstar * p_Dstar;
            const double sqrt_q2  = sqrt(q2);

            const double a_0 = ff->a_0(q2);
            const double a_1 = ff->a_1(q2);
            const double a_2 = ff->a_2(q2);
            const double v   = ff->v(q2);

            // cf. [CJLP2012]
            const double H_pp = (m_B + m_Dstar) * a_1 - 2.0 * m_B / (m_B + m_Dstar) * p_Dstar * v;
            const double H_mm = (m_B + m_Dstar) * a_1 + 2.0 * m_B / (m_B + m_Dstar) * p_Dstar * v;
            const double H_00 = ((m_B2 - m_Dstar2 - q2) * (m_B + m_Dstar) * a_1 - 4.0 * m_B2 * p_Dstar2 * a_2 / (m_B + m_Dstar))
                              / (2.0 * m_Dstar * sqrt_q2);
            const double H_0t = 2.0 * m_B * p_Dstar / sqrt_q2 * a_0;

            const double H_pp2 = H_pp * H_pp;
            const double H_mm2 = H_mm * H_mm;
            const double H_002 = H_00 * H_00;
            const double H_0t2 = H_0t * H_0t;

            // cf. [CJLP2012]], eq. (22), p. 17
            const double denom = (H_pp2 + H_mm2 + H_002) * (1.0 + m_l2 / (2.0 * q2)) + 3.0 * m_l2 / (2.0 * q2) * H_0t2;

            return nf * denom;
        }

        double lepton_polarization() const
        {
            const double q2_min = power_of<2>(m_l());
            const double q2_max = 10.68;

            std::function<double (const double &)> integrand_num   = std::bind(&Implementation<BToDPiLeptonNeutrino>::lepton_polarization_numerator,   this, std::placeholders::_1);
            std::function<double (const double &)> integrand_denom = std::bind(&Implementation<BToDPiLeptonNeutrino>::lepton_polarization_denominator, this, std::placeholders::_1);
            const auto num   = integrate<GSL::QAGS>(integrand_num,   q2_min, q2_max);
            const auto denom = integrate<GSL::QAGS>(integrand_denom, q2_min, q2_max);

            return num / denom;
        }

        double longitudinal_polarization_L(const double & q2) const
        {
            const double nf = pdf_normalization(q2);

            const double m_l2    = m_l() * m_l();
            const double a_long2 = norm(a_long(q2));
            const double a_time2 = norm(a_time(q2));

            // cf. [DDS2014], eq. (11), p. p
            const double result = 2.0       * ((1.0 + m_l2 / q2) * a_long2 + 2.0 * m_l2 / q2 * a_time2)
                                + 2.0 / 3.0 * ((1.0 - m_l2 / q2) * a_long2);

            return nf * result;
        }

        double longitudinal_polarization_T(const double & q2) const
        {
            const double nf = pdf_normalization(q2);

            const double m_l2    = m_l() * m_l();
            const double a_para2 = norm(a_para(q2));
            const double a_perp2 = norm(a_perp(q2));

            // cf. [DDS2014], eq. (11), p. p
            const double result =             (3.0 + m_l2 / q2) * (a_para2 + a_perp2)
                                - 1.0 / 3.0 * (1.0 - m_l2 / q2) * (a_para2 + a_perp2);

            return nf * result;
        }

        double longitudinal_polarization() const
        {
            const double q2_min = power_of<2>(m_l());
            const double q2_max = 10.68;

            std::function<double (const double &)> integrand_L = std::bind(&Implementation<BToDPiLeptonNeutrino>::longitudinal_polarization_L, this, std::placeholders::_1);
            std::function<double (const double &)> integrand_T = std::bind(&Implementation<BToDPiLeptonNeutrino>::longitudinal_polarization_T, this, std::placeholders::_1);
            const auto L = integrate<GSL::QAGS>(integrand_L, q2_min, q2_max);
            const auto T = integrate<GSL::QAGS>(integrand_T, q2_min, q2_max);

            return L / (L + T);
        }

        double dist_q2(const double & q2) const
        {
            const double nf = pdf_normalization(q2);

            const double m_l2    = m_l() * m_l();
            const double a_long2 = norm(a_long(q2));
            const double a_para2 = norm(a_para(q2));
            const double a_perp2 = norm(a_perp(q2));
            const double a_time2 = norm(a_time(q2));

            const double a = 2.0 * (a_long2 + a_para2 + a_perp2) * (1.0 + m_l2 / (2.0 * q2))
                           + 3.0 * a_time2 * m_l2 / q2;

            return nf * a;
        }

        double pdf_q2(const double & q2) const
        {
            const double q2_min = power_of<2>(m_l());
            const double q2_max = 10.68;

            std::function<double (const double &)> f = std::bind(&Implementation<BToDPiLeptonNeutrino>::dist_q2, this, std::placeholders::_1);
            const double num   = dist_q2(q2);
            const double denom = integrate1D(f, 32, q2_min, q2_max);

            return num / denom;
        }

        double pdf_w(const double & w) const
        {
            const double m_B     = this->m_B(),     m_B2 = m_B * m_B;
            const double m_Dstar = this->m_Dstar(), m_Dstar2 = m_Dstar * m_Dstar;
            const double q2      = m_B2 + m_Dstar2 - 2.0 * m_B * m_Dstar * w;

            return 2.0 * m_B * m_Dstar * pdf_q2(q2);
        }

        double integrated_pdf_w(const double & w_min, const double & w_max) const
        {
            std::function<double (const double &)> f = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_w, this, std::placeholders::_1);

            return integrate<GSL::QAGS>(f, w_min, w_max) / (w_max - w_min);
        }

        double pdf_coefficients_q2d_a(const double & q2) const
        {
            const double nf = pdf_normalization(q2);

            const double m_l2    = m_l() * m_l();

            const double a_para2 = norm(a_para(q2));
            const double a_perp2 = norm(a_perp(q2));

            const double a = (a_para2 + a_perp2) * (1.0 + m_l2 / (2.0 * q2));

            return nf * a;
        }

        double pdf_coefficients_q2d_b(const double & q2) const
        {
            const double nf = pdf_normalization(q2);

            const double m_l2    = m_l() * m_l();

            const double a_long2 = norm(a_long(q2));
            const double a_para2 = norm(a_para(q2));
            const double a_perp2 = norm(a_perp(q2));
            const double a_time2 = norm(a_time(q2));

            const double b = (2.0 * a_long2 - a_para2 - a_perp2) * (1.0 + m_l2 / (2.0 * q2))
                           + 3.0 * m_l2 / q2 * a_time2;

            return nf * b;
        }

        double pdf_q2d(const double & q2, const double & c_d) const
        {
            auto a = pdf_coefficients_q2d_a(q2);
            auto b = pdf_coefficients_q2d_b(q2);

            return 3.0 / 2.0 * (a + b * c_d * c_d);
        }

        double pdf_d(const double & c_d) const
        {
            const double q2_min = power_of<2>(m_l());
            const double q2_max = 10.68;

            std::function<double (const double &)> integrand_a = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_coefficients_q2d_a, this, std::placeholders::_1);
            std::function<double (const double &)> integrand_b = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_coefficients_q2d_b, this, std::placeholders::_1);
            const auto a = integrate<GSL::QAGS>(integrand_a, q2_min, q2_max);
            const auto b = integrate<GSL::QAGS>(integrand_b, q2_min, q2_max);

            const double num   = 3.0 / 2.0 * (a + b * c_d * c_d);
            const double denom = 3.0 * a + b;

            return num / denom;
        }

        double pdf_d(const double & c_d_min, const double & c_d_max) const
        {
            const double q2_min = power_of<2>(m_l());
            const double q2_max = 10.68;

            const double c_d_max3 = pow(c_d_max, 3);
            const double c_d_min3 = pow(c_d_min, 3);

            std::function<double (const double &)> integrand_a = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_coefficients_q2d_a, this, std::placeholders::_1);
            std::function<double (const double &)> integrand_b = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_coefficients_q2d_b, this, std::placeholders::_1);
            const auto a = integrate<GSL::QAGS>(integrand_a, q2_min, q2_max);
            const auto b = integrate<GSL::QAGS>(integrand_b, q2_min, q2_max);

            const double num   = 3.0 / 2.0 * (a * (c_d_max - c_d_min) + b * (c_d_max3 - c_d_min3) / 3.0);
            const double denom = 3.0 * a + b;

            return num / denom;
        }

        double pdf_coefficients_q2l_a(const double & q2) const
        {
            const double nf = pdf_normalization(q2);

            const double m_l2         = m_l() * m_l();

            const double a_long       = this->a_long(q2);
            const double a_para       = this->a_para(q2);
            const double a_perp       = this->a_perp(q2);
            const double a_time       = this->a_time(q2);

            const double a_long2      = norm(a_long);
            const double a_para2      = norm(a_para);
            const double a_perp2      = norm(a_perp);
            const double a_time2      = norm(a_time);

            const double a = 2.0 * a_long2 + (a_para2 + a_perp2) * (1.0 + m_l2 / q2)
                           + 2.0 * a_time2 * m_l2 / q2;

            return nf * a;
        }

        double pdf_coefficients_q2l_b(const double & q2) const
        {
            const double nf = pdf_normalization(q2);

            const double m_l2         = m_l() * m_l();

            const double a_long       = this->a_long(q2);
            const double a_para       = this->a_para(q2);
            const double a_perp       = this->a_perp(q2);
            const double a_time       = this->a_time(q2);

            const double re_para_perp = real(a_para * conj(a_perp));
            const double re_time_long = real(a_time * conj(a_long));

            const double b = -4.0 * (re_para_perp + re_time_long * m_l2 / q2);

            return nf * b;
        }

        double pdf_coefficients_q2l_c(const double & q2) const
        {
            const double nf = pdf_normalization(q2);

            const double m_l2         = m_l() * m_l();

            const double a_long       = this->a_long(q2);
            const double a_para       = this->a_para(q2);
            const double a_perp       = this->a_perp(q2);

            const double a_long2      = norm(a_long);
            const double a_para2      = norm(a_para);
            const double a_perp2      = norm(a_perp);

            const double c = -(2.0 * a_long2 - a_para2 - a_perp2) * (1.0 - m_l2 / q2);

            return nf * c;
        }

        double pdf_l(const double & c_l) const
        {
            const double q2_min = power_of<2>(m_l());
            const double q2_max = 10.68;

            std::function<double (const double &)> integrand_a = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_coefficients_q2l_a, this, std::placeholders::_1);
            std::function<double (const double &)> integrand_b = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_coefficients_q2l_b, this, std::placeholders::_1);
            std::function<double (const double &)> integrand_c = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_coefficients_q2l_c, this, std::placeholders::_1);
            const double a = integrate<GSL::QAGS>(integrand_a, q2_min, q2_max);
            const double b = integrate<GSL::QAGS>(integrand_b, q2_min, q2_max);
            const double c = integrate<GSL::QAGS>(integrand_c, q2_min, q2_max);

            const double num   = 3.0 / 4.0 * (a + b * c_l + c * c_l * c_l);
            const double denom = (3.0 * a + c) / 2.0;

            return num / denom;
        }

        double pdf_l(const double & c_l_min, const double & c_l_max) const
        {
            const double q2_min = power_of<2>(m_l());
            const double q2_max = 10.68;

            std::function<double (const double &)> integrand_a = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_coefficients_q2l_a, this, std::placeholders::_1);
            std::function<double (const double &)> integrand_b = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_coefficients_q2l_b, this, std::placeholders::_1);
            std::function<double (const double &)> integrand_c = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_coefficients_q2l_c, this, std::placeholders::_1);
            const double a = integrate<GSL::QAGS>(integrand_a, q2_min, q2_max);
            const double b = integrate<GSL::QAGS>(integrand_b, q2_min, q2_max);
            const double c = integrate<GSL::QAGS>(integrand_c, q2_min, q2_max);

            double num         = a * (c_l_max - c_l_min);
            num               += b * (c_l_max * c_l_max - c_l_min * c_l_min) / 2.0;
            num               += c * (c_l_max * c_l_max * c_l_max - c_l_min * c_l_min * c_l_min) / 3.0;
            num               *= 3.0 / 4.0;
            const double denom = (3.0 * a + c) / 2.0;

            return num / denom;
        }

        double pdf_coefficients_q2chi_a(const double & q2) const
        {
            const double nf = pdf_normalization(q2);

            const double m_l2         = m_l() * m_l();

            const double a_long       = this->a_long(q2);
            const double a_para       = this->a_para(q2);
            const double a_perp       = this->a_perp(q2);
            const double a_time       = this->a_time(q2);

            const double a_long2      = norm(a_long);
            const double a_para2      = norm(a_para);
            const double a_perp2      = norm(a_perp);
            const double a_time2      = norm(a_time);

            const double a = (2.0 +       m_l2 / q2) * a_long2
                           +  3.0                    * a_para2
                           + (1.0 + 2.0 * m_l2 / q2) * a_perp2
                           +        3.0 * m_l2 / q2  * a_time2;

            return nf * a;
        }

        double pdf_coefficients_q2chi_b(const double & q2) const
        {
            const double nf = pdf_normalization(q2);

            const double m_l2         = m_l() * m_l();

            const double a_para       = this->a_para(q2);
            const double a_perp       = this->a_perp(q2);

            const double a_para2      = norm(a_para);
            const double a_perp2      = norm(a_perp);

            const double b = -2.0 * (a_para2 - a_perp2) * (1.0 - m_l2 / q2);

            return nf * b;
        }

        double pdf_chi(const double & chi) const
        {
            const double q2_min = power_of<2>(m_l());
            const double q2_max = 10.68;

            std::function<double (const double &)> integrand_a = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_coefficients_q2chi_a, this, std::placeholders::_1);
            std::function<double (const double &)> integrand_b = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_coefficients_q2chi_b, this, std::placeholders::_1);

            const double a = integrate<GSL::QAGS>(integrand_a, q2_min, q2_max, GSL::QAGS::Config().epsrel(1e-5));
            const double b = integrate<GSL::QAGS>(integrand_b, q2_min, q2_max, GSL::QAGS::Config().epsrel(1e-5));

            const double c_chi = cos(chi);

            const double num   = (a + b * c_chi * c_chi) / (2.0 * M_PI);
            const double denom = a + b / 2.0;

            return num / denom;
        }

        double pdf_chi(const double & chi_min, const double & chi_max) const
        {
            const double q2_min = power_of<2>(m_l());
            const double q2_max = 10.68;

            std::function<double (const double &)> integrand_a = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_coefficients_q2chi_a, this, std::placeholders::_1);
            std::function<double (const double &)> integrand_b = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_coefficients_q2chi_b, this, std::placeholders::_1);
            const double a = integrate<GSL::QAGS>(integrand_a, q2_min, q2_max);
            const double b = integrate<GSL::QAGS>(integrand_b, q2_min, q2_max);

            const double c_chi_min = cos(chi_min), c_chi_max = cos(chi_max);
            const double s_chi_min = sin(chi_min), s_chi_max = sin(chi_max);

            double num         = a * (chi_max - chi_min);
            num               += b * (chi_max - chi_min + s_chi_max * c_chi_max - s_chi_min * c_chi_min) / 2.0;
            num               /= 2.0 * M_PI;
            const double denom = a + b / 2.0;

            return num / denom;
        }

        Diagnostics diagnostics() const
        {
            Diagnostics results;

            const double q20 = opt_l.value() == "tau" ? 4.0 : 1.0;

            // reference point
            {
                results.add({ q20, "q2^0" });
            }

            // amplitudes
            {
                // q2 =  q2^0
                results.add({ a_time( q20), "A_time(q2 = q2^0)" });
                results.add({ a_perp( q20), "A_perp(q2 = q2^0)" });
                results.add({ a_para( q20), "A_para(q2 = q2^0)" });
                results.add({ a_long( q20), "A_long(q2 = q2^0)" });

                // q2 = 10.0
                results.add({ a_time(10.0), "A_time(q2 = 10.0)" });
                results.add({ a_perp(10.0), "A_perp(q2 = 10.0)" });
                results.add({ a_para(10.0), "A_para(q2 = 10.0)" });
                results.add({ a_long(10.0), "A_long(q2 = 10.0)" });
            }

            // normalization
            {
                // q2 =  1.0
                results.add({ pdf_normalization( q20), "N(q2 = q2^0)" });

                // q2 = 10.0
                results.add({ pdf_normalization(10.0), "N(q2 = 10.0)" });
            }

            // (unintegrated!) coefficients of the PDFs
            {
                results.add({ pdf_coefficients_q2d_a( q20),   "P(c_d) term ~          1 @ q2 = q2^0"});
                results.add({ pdf_coefficients_q2d_b( q20),   "P(c_d) term ~        c_d @ q2 = q2^0"});

                results.add({ pdf_coefficients_q2d_a(10.0),   "P(c_d) term ~          1 @ q2 = 10.0"});
                results.add({ pdf_coefficients_q2d_b(10.0),   "P(c_d) term ~        c_d @ q2 = 10.0"});

                results.add({ pdf_coefficients_q2l_a( q20),   "P(c_l) term ~          1 @ q2 = q2^0"});
                results.add({ pdf_coefficients_q2l_b( q20),   "P(c_l) term ~        c_l @ q2 = q2^0"});
                results.add({ pdf_coefficients_q2l_c( q20),   "P(c_l) term ~      c_l^2 @ q2 = q2^0"});

                results.add({ pdf_coefficients_q2l_a(10.0),   "P(c_l) term ~          1 @ q2 = 10.0"});
                results.add({ pdf_coefficients_q2l_b(10.0),   "P(c_l) term ~        c_l @ q2 = 10.0"});
                results.add({ pdf_coefficients_q2l_c(10.0),   "P(c_l) term ~      c_l^2 @ q2 = 10.0"});

                results.add({ pdf_coefficients_q2chi_a( q20), "P(chi) term ~          1 @ q2 = q2^0"});
                results.add({ pdf_coefficients_q2chi_b( q20), "P(chi) term ~ cos(chi)^2 @ q2 = q2^0"});

                results.add({ pdf_coefficients_q2chi_a(10.0), "P(chi) term ~          1 @ q2 = 10.0"});
                results.add({ pdf_coefficients_q2chi_b(10.0), "P(chi) term ~ cos(chi)^2 @ q2 = 10.0"});
            }

            return results;
        }
    };

    BToDPiLeptonNeutrino::BToDPiLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToDPiLeptonNeutrino>(new Implementation<BToDPiLeptonNeutrino>(parameters, options, *this))
    {
    }

    BToDPiLeptonNeutrino::~BToDPiLeptonNeutrino()
    {
    }

    double
    BToDPiLeptonNeutrino::lepton_polarization() const
    {
        return _imp->lepton_polarization();
    }

    double
    BToDPiLeptonNeutrino::longitudinal_polarization() const
    {
        return _imp->longitudinal_polarization();
    }

    double
    BToDPiLeptonNeutrino::differential_pdf_d(const double & c_d) const
    {
        return _imp->pdf_d(c_d);
    }

    double
    BToDPiLeptonNeutrino::differential_pdf_l(const double & c_l) const
    {
        return _imp->pdf_l(c_l);
    }

    double
    BToDPiLeptonNeutrino::differential_pdf_chi(const double & chi) const
    {
        return _imp->pdf_chi(chi);
    }

    double
    BToDPiLeptonNeutrino::differential_pdf_q2(const double & q2) const
    {
        return _imp->pdf_q2(q2);
    }

    double
    BToDPiLeptonNeutrino::differential_pdf_w(const double & w) const
    {
        return _imp->pdf_w(w);
    }

    double
    BToDPiLeptonNeutrino::integrated_pdf_d(const double & c_d_min, const double & c_d_max) const
    {
        return _imp->pdf_d(c_d_min, c_d_max);
    }

    double
    BToDPiLeptonNeutrino::integrated_pdf_l(const double & c_l_min, const double & c_l_max) const
    {
        return _imp->pdf_l(c_l_min, c_l_max);
    }

    double
    BToDPiLeptonNeutrino::integrated_pdf_chi(const double & chi_min, const double & chi_max) const
    {
        return _imp->pdf_chi(chi_min, chi_max);
    }

    double
    BToDPiLeptonNeutrino::integrated_pdf_w(const double & w_min, const double & w_max) const
    {
        return _imp->integrated_pdf_w(w_min, w_max);
    }

    Diagnostics
    BToDPiLeptonNeutrino::diagnostics() const
    {
        return _imp->diagnostics();
    }

    const std::string
    BToDPiLeptonNeutrino::description = "\
The decay B->D pi l nu, where l is a massless lepton.";

    const std::string
    BToDPiLeptonNeutrino::kinematics_description_c_d = "\
The cosine of the helicity angle theta_D in the D-pi rest frame.";

    const std::string
    BToDPiLeptonNeutrino::kinematics_description_c_l = "\
The cosine of the helicity angle theta_L in the l-nu rest frame.";

    const std::string
    BToDPiLeptonNeutrino::kinematics_description_chi = "\
The azimuthal angle between the decay planes.";

    const std::string
    BToDPiLeptonNeutrino::kinematics_description_q2 = "\
The squared mass of the dilepton pair.";
}
