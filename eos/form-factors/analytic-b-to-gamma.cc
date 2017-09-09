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

#include <eos/form-factors/analytic-b-to-gamma.hh>
#include <eos/utils/model.hh>
#include <eos/utils/polylog.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>

#include <iostream>

namespace eos
{
    template <>
    struct Implementation<AnalyticFormFactorBToGammaQCDF>
    {
        std::shared_ptr<Model> model;

        // hadronic parameters
        UsedParameter f_B;
        UsedParameter m_B;

        // renormalization scale
        UsedParameter mu;

        // B-meson LCDA moments
        UsedParameter lambda_B_p_mu0;
        UsedParameter sigma_1_B_p_mu0;
        UsedParameter sigma_2_B_p_mu0;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            f_B(p["decay-constant::B_u"], u),
            m_B(p["mass::B_u"], u),
            mu(p["B->gamma::mu@QCDF"], u),
            lambda_B_p_mu0(p["lambda_B_p"], u),
            sigma_1_B_p_mu0(p["sigma_1_B_p"], u),
            sigma_2_B_p_mu0(p["sigma_2_B_p"], u)
        {
            //std::cout << "[lambda_B(mu_0) / lambda_B(mu) x R(Egamma = 2.0, mu)]          = " << lambda_B_p(1.0) / lambda_B_p(mu()) * R(2.0) << std::endl;
            //std::cout << "mu                           = " << mu() << std::endl;
            //std::cout << "lambda_B_p(1.0)              = " << lambda_B_p(1.0) << std::endl;
            //std::cout << "lambda_B_p(mu)               = " << lambda_B_p(mu()) << std::endl;
            //std::cout << "lambda_B_p ratio             = " << lambda_B_p(1.0) / lambda_B_p(mu()) << std::endl;
            //std::cout << "alpha_s(mu = 1.0)            = " << model->alpha_s(1.0) << std::endl;
            //std::cout << "sigma_1_B_p(1.5)             = " << sigma_1_B_p(1.5) << std::endl;
            //std::cout << "sigma_2_B_p(1.5)             = " << sigma_2_B_p(1.5) << std::endl;
            //std::cout << "alpha_s(mu = 4.2)            = " << model->alpha_s(4.2) << std::endl;
            //std::cout << "C(Egamma = 2.0, mu_h1 = 4.2) = " << C(2.0, 4.2) << std::endl;
            //std::cout << "K(mu_h2 = 4.2)               = " << K(4.2) << std::endl;
            //std::cout << "K(mu_h2 = 1.5)               = " << K(1.5) << std::endl;
            //std::cout << "J(Egamma = 2.0)              = " << J(2.0) << std::endl;
        }

        // cf. [BR2013], eq. (3.1), p. 6
        double lambda_B_p(const double & mu) const
        {
            static const double C_F = 4.0 / 3.0;
            static const double mu_0 = 1.0;

            const double alpha_s = model->alpha_s(mu_0);
            const double L       = std::log(mu / mu_0);

            double result = 1.0;
            result += alpha_s * C_F / (4.0 * M_PI) * L * (2.0 - 2.0 * L - 4.0 * sigma_1_B_p_mu0);
            result = lambda_B_p_mu0 / result;

            return result;
            // DvD 2017-09-13: checked
        }

        // cf. [BFWY2013], eq. (2.28), p. 9 in the limit g -> 0
        double sigma_1_B_p(const double & mu) const
        {
            static const double mu_0 = 1.0;

            const double L = std::log(mu_0 / mu);

            return sigma_1_B_p_mu0 - L;
        }

        // cf. [BFWY2013], eq. (2.28), p. 9 in the limit g -> 0
        double sigma_2_B_p(const double & mu) const
        {
            static const double mu_0 = 1.0;

            const double L = std::log(mu_0 / mu);

            return sigma_2_B_p_mu0 - 2.0 * L * sigma_1_B_p_mu0 + L * L;
        }

        double m_b_msbar(const double & mu_) const
        {
            return 4.2;
            return model->m_b_msbar(mu_);
        }

        // matching between HQET and QCD definitions of the B-meson decay constant,
        // cf. [BR2011], eq. (2.15), p. 4
        double K(const double & mu_h2) const
        {
            static const double C_F = 4.0 / 3.0;

            const double alpha_s = model->alpha_s(mu_h2);
            const double m_b     = m_b_msbar(mu_h2);
            const double L       = std::log(m_b / mu_h2);

            return 1.0 + alpha_s * C_F / (4.0 * M_PI) * (3.0 * L - 2.0);
            // DvD 2017-09-13: checked
        }

        // hard SCET/QCD matching coefficient, cf. [BR2011], eq. (2.12), p. 4
        double C(const double & Egamma, const double & mu_h1) const
        {
            static const double C_F = 4.0 / 3.0;

            const double alpha_s  = model->alpha_s(mu_h1);
            const double m_b      = m_b_msbar(mu_h1);
            const double x        = 2.0 * Egamma / m_b;
            const double L        = std::log(2.0 * Egamma / mu_h1);
            const double logx     = std::log(x);
            const double dilog1mx = real(dilog(1.0 - x));

            double result = (-2.0 * L + 5.0) * L - (3.0 - 2.0 * x) / (1.0 - x) * logx - 2.0 * dilog1mx - 6.0 - M_PI * M_PI / 12.0;
            result *= alpha_s * C_F / (4.0 * M_PI);
            result += 1.0;

            return result;
            // DvD 2017-09-13: checked
        }

        // hard-collinear (Jet-) function, cf. [BR2011], eq. (2.14), p. 4
        double J(const double & Egamma) const
        {
            static const double C_F = 4.0 / 3.0;
            static const double mu_0 = 1.0;

            const double alpha_s = model->alpha_s(mu());
            const double L       = std::log((2.0 * Egamma * mu_0) / (mu * mu));

            double result = L * L - 2.0 * sigma_1_B_p(mu()) * L - 1.0 - M_PI * M_PI / 6.0 + sigma_2_B_p(mu());
            result *= alpha_s * C_F / (4.0 * M_PI);
            result += 1.0;

            return result;
            // DvD 2017-09-13: checked
        }

        // factor resumming logarithms of ratio from hard (SCET) to hard-collinear scales to NLL accuracy
        double U_1(const double & Egamma, const double & mu_h1) const
        {
            static const double pi2 = M_PI * M_PI;
            static const double C_F = 4.0 / 3.0;
            static const double n_l = 4.0;
            static const double zeta_3 = 1.20205690315959428540;
            // coefficents of the cusp anomalous dimension, cf. [BR2011], eq. (A.4), p. 12
            static const double Gamma_0 = C_F * 4.0;
            static const double Gamma_1 = C_F * (268.0 / 3.0 - 4.0 * pi2 - 40.0 / 9.0 * n_l);
            static const double Gamma_2 = C_F * (1470.0 - 536.0 / 3.0 * pi2 + 44.0 / 5.0 * pi2 * pi2
                                        + 264.0 * zeta_3 - 16.0 / 27.0 * n_l * n_l
                                        + n_l * (-1276.0 / 9.0 + 80.0 / 9.0 * pi2 - 208.0 / 3.0 * zeta_3));
            // coefficents of the SCET heavy-light current's anomalous dimension,
            // cf. [BR2011], eq. (A.5), p. 13
            static const double gamma_0 = C_F * (-5.0);
            static const double gamma_1 = C_F * (-1585.0 / 18.0 - 5.0 * pi2 / 6.0 + 34.0 * zeta_3
                                        + (125.0 / 27.0 + pi2 / 3.0) * n_l);
            static const auto   beta    = QCD::beta_function_nf_4;

            const double a_s_h = model->alpha_s(mu_h1) / (4.0 * M_PI);
            const double a_s   = model->alpha_s(mu())  / (4.0 * M_PI);
            const double r     = a_s / a_s_h;
            const double lnr   = std::log(r);

            // exponential prefactors of eq. (A.3)
            const double fact_cusp1 = std::exp(-Gamma_0 / (4.0 * beta[0] * beta[0]) * (
                    (lnr - 1.0 + 1.0 / r) / a_s_h - beta[1] / (2.0 * beta[0]) * lnr * lnr
                    + (Gamma_1 / Gamma_0 - beta[1] / beta[0]) * (r - 1.0 - lnr)));
            const double fact_cusp2 = std::pow(2.0 * Egamma / mu_h1, -Gamma_0 * lnr / (2.0 * beta[0]));
            const double fact_ll    = std::pow(r, -gamma_0 / (2.0 * beta[0]));

            const double nlo_nll   = (gamma_1 / (2.0 * beta[0]) - gamma_0 * beta[1] / (2.0 * beta[0] * beta[0])) * (1 - r);
            const double nlo_cusp1 = std::log(2.0 * Egamma / mu_h1) * (Gamma_1 / (2.0 * beta[0]) - Gamma_0 * beta[1] / (2.0 * beta[0] * beta[0])) * (1.0 - r);
            const double nlo_cusp2 = Gamma_0 / (4.0 * beta[0] * beta[0]) * (Gamma_2 / (2.0 * Gamma_0) * pow(1.0 - r, 2)
                                   + beta[2] / (2.0 * beta[0]) * (1 - r * r + 2.0 * lnr)
                                   - Gamma_1 * beta[1] / (2.0 * Gamma_0 * beta[0]) * (3.0 - 4.0 * r + r * r + 2.0 * r * lnr)
                                   + beta[1] * beta[1] / (2.0 * beta[0] * beta[0]) * (1.0 - r) * (1.0 - r - 2.0 * lnr)
                                    );

            return fact_cusp1 * fact_cusp2 * fact_ll * (1.0 + a_s_h * (nlo_nll + nlo_cusp1 - nlo_cusp2));
        }

        // factor resumming logarithms of ratio from hard (HQET) to hard-collinear scales to NLL accuracy
        double U_2(const double & mu_h2) const
        {
            static const double pi2 = M_PI * M_PI;
            static const double C_F = 4.0 / 3.0;
            static const double n_l = 4.0;
            // coefficents of the HQET heavy-light current's anomalous dimension,
            // cf. [BR2011], eq. (A.6), p. 13
            static const double gamma_0 = C_F * (-3.0);
            static const double gamma_1 = C_F * (-127.0 / 6.0 - 14.0 * pi2 / 9.0 + 5.0 / 3.0 * n_l);
            static const auto   beta    = QCD::beta_function_nf_4;

            const double a_s_h = model->alpha_s(mu_h2) / (4.0 * M_PI);
            const double a_s   = model->alpha_s(mu())  / (4.0 * M_PI);
            const double r     = a_s / a_s_h;

            const double fact_ll = std::pow(r, -gamma_0 / (2.0 * beta[0]));
            const double nlo_nll = (gamma_1 / (2.0 * beta[0]) - gamma_0 * beta[1] / (2.0 * beta[0] * beta[0])) * (1.0 - r);

            return fact_ll * (1.0 + a_s_h * nlo_nll);
        }
#if 1
        double R(const double & Egamma) const
        {
            const double mu_h1 = 4.2;
            const double mu_h2 = 4.2;

            return C(Egamma, mu_h1) / K(mu_h2) * U_1(Egamma, mu_h1) / U_2(mu_h2) * J(Egamma);
        }
#endif

#if 0
        double R(const double Egamma) const
        {
            const double mu_h1 = 4.2;
            const double mu_h2 = 4.2;

            double deltaC = 0.0;
            {
                static const double C_F = 4.0 / 3.0;

                const double alpha_s  = model->alpha_s(mu_h1);
                const double m_b      = m_b_msbar(mu_h1);
                const double x        = 2.0 * Egamma / m_b;
                const double L        = std::log(2.0 * Egamma / mu_h1);
                const double logx     = std::log(x);
                const double dilog1mx = real(dilog(1.0 - x));

                deltaC = (-2.0 * L + 5.0) * L - (3.0 - 2.0 * x) / (1.0 - x) * logx - 2.0 * dilog1mx - 6.0 - M_PI * M_PI / 12.0;
                deltaC *= alpha_s * C_F / (4.0 * M_PI);
            }

            double deltaK = 0.0;
            {
                static const double C_F = 4.0 / 3.0;

                const double alpha_s = model->alpha_s(mu_h2);
                const double m_b     = m_b_msbar(mu_h2);
                const double L       = std::log(m_b / mu_h2);

                deltaK = alpha_s * C_F / (4.0 * M_PI) * (3.0 * L - 2.0);
            }

            double deltaJ = 0.0;
            {
                static const double C_F = 4.0 / 3.0;
                static const double mu_0 = 1.0;

                const double alpha_s = model->alpha_s(mu());
                const double L       = std::log((2.0 * Egamma * mu_0) / (mu * mu));

                deltaJ  = L * L - 2.0 * sigma_1_B_p(mu()) * L - 1.0 - M_PI * M_PI / 6.0 + sigma_2_B_p(mu());
                deltaJ *= alpha_s * C_F / (4.0 * M_PI);
            }


            double rU1 = 1.0;
            double deltaU1 = 0.0;
            {
                static const double pi2 = M_PI * M_PI;
                static const double C_F = 4.0 / 3.0;
                static const double n_l = 4.0;
                static const double zeta_3 = 1.20205690315959428539;
                // coefficents of the cusp anomalous dimension, cf. [BR2011], eq. (A.4), p. 12
                static const double Gamma_0 = C_F * 4.0;
                static const double Gamma_1 = C_F * (268.0 / 3.0 - 4.0 * pi2 - 40.0 / 9.0 * n_l);
                static const double Gamma_2 = C_F * (1470.0 - 536.0 / 3.0 * pi2 + 44.0 / 5.0 * pi2 * pi2
                                            + 264.0 * zeta_3 - 16.0 / 27.0 * n_l * n_l
                                            + n_l * (-1276.0 / 9.0 + 80.0 / 9.0 * pi2 - 208.0 / 3.0 * zeta_3));
                // coefficents of the SCET heavy-light current's anomalous dimension,
                // cf. [BR2011], eq. (A.5), p. 13
                static const double gamma_0 = C_F * (-5.0);
                static const double gamma_1 = C_F * (-1585.0 / 18.0 - 5.0 * pi2 / 6.0 + 34.0 * zeta_3
                                            + (125.0 / 27.0 + pi2 / 3.0) * n_l);
                static const auto   beta    = QCD::beta_function_nf_4;

                const double a_s_h = model->alpha_s(mu_h1) / (4.0 * M_PI);
                const double a_s   = model->alpha_s(mu())  / (4.0 * M_PI);
                const double r     = a_s / a_s_h;
                const double lnr   = std::log(r);

                // exponential prefactors of eq. (A.3)
                const double fact_cusp1 = std::exp(-Gamma_0 / (4.0 * beta[0] * beta[0]) * (
                        (lnr - 1.0 + 1.0 / r) / a_s_h - beta[1] / (2.0 * beta[0]) * lnr * lnr
                        + (Gamma_1 / Gamma_0 - beta[1] / beta[0]) * (r - 1.0 - lnr)));
                const double fact_cusp2 = std::pow(2.0 * Egamma / mu_h1, -Gamma_0 * lnr / (2.0 * beta[0]));
                const double fact_ll    = std::pow(r, -gamma_0 / (2.0 * beta[0]));

                const double nlo_nll   = (gamma_1 / (2.0 * beta[0]) - gamma_0 * beta[1] / (2.0 * beta[0] * beta[0])) * (1 - r);
                const double nlo_cusp1 = std::log(2.0 * Egamma / mu_h1) * (Gamma_1 / (2.0 * beta[0]) - Gamma_0 * beta[1] / (2.0 * beta[0] * beta[0])) * (1.0 - r);
                const double nlo_cusp2 = Gamma_0 / (4.0 * beta[0] * beta[0]) * (Gamma_2 / (2.0 * Gamma_0) * pow(1.0 - r, 2)
                                       + beta[2] / (2.0 * beta[0]) * (1 - r * r + 2.0 * lnr)
                                       - Gamma_1 * beta[1] / (2.0 * Gamma_0 * beta[0]) * (3.0 - 4.0 * r + r * r + 2.0 * r * lnr)
                                       + beta[1] * beta[1] / (2.0 * beta[0] * beta[0]) * (1.0 - r) * (1.0 - r - 2.0 * lnr)
                                        );

                rU1 = fact_cusp1 * fact_cusp2 * fact_ll;
                deltaU1 = a_s_h * (nlo_nll + nlo_cusp1 - nlo_cusp2);
            }

            double rU2 = 1.0;
            double deltaU2 = 0.0;
            {
                static const double pi2 = M_PI * M_PI;
                static const double C_F = 4.0 / 3.0;
                static const double n_l = 4.0;
                // coefficents of the HQET heavy-light current's anomalous dimension,
                // cf. [BR2011], eq. (A.6), p. 13
                static const double gamma_0 = C_F * (-3.0);
                static const double gamma_1 = C_F * (-127.0 / 6.0 - 14.0 * pi2 / 9.0 + 5.0 / 3.0 * n_l);
                static const auto   beta    = QCD::beta_function_nf_4;

                const double a_s_h = model->alpha_s(mu_h2) / (4.0 * M_PI);
                const double a_s   = model->alpha_s(mu())  / (4.0 * M_PI);
                const double r     = a_s / a_s_h;

                const double fact_ll = std::pow(r, -gamma_0 / (2.0 * beta[0]));
                const double nlo_nll = (gamma_1 / (2.0 * beta[0]) - gamma_0 * beta[1] / (2.0 * beta[0] * beta[0])) * (1.0 - r);

                rU2 = fact_ll;
                deltaU2 = a_s_h * nlo_nll;
            }

            return (1.0 + deltaC - deltaK + deltaJ + deltaU1 - deltaU2) * rU1 / rU2;
        }
#endif

        // cf. [BR2011], eq. (2.9), p. 5
        double f_v(const double & Egamma) const
        {
            static const double Q_b = -1.0 / 3.0;
            static const double Q_u = +2.0 / 3.0;
            const double R   = this->R(Egamma);
            const double m_b = model->m_b_pole();

            const double f_symm =  Q_u * m_B * f_B / (2.0 * lambda_B_p(mu()) * Egamma) * R;
            const double f_asym = +Q_b * m_B * f_B / (2.0 * m_b              * Egamma)
                                +  Q_u * m_B * f_B / pow(2.0 * Egamma, 2);

            return f_symm + f_asym;
        }

        // cf. [BR2011], eq. (2.9), p. 5
        double f_a(const double & Egamma) const
        {
            static const double Q_b = -1.0 / 3.0;
            static const double Q_l = -1.0;
            static const double Q_u = +2.0 / 3.0;
            const double R   = this->R(Egamma);
            const double m_b = model->m_b_pole();

            const double f_symm =  Q_u * m_B * f_B / (2.0 * lambda_B_p(mu()) * Egamma) * R;
            const double f_asym = -Q_b * m_B * f_B / (2.0 * m_b              * Egamma)
                                -  Q_u * m_B * f_B / pow(2.0 * Egamma, 2)
                                +  Q_l *       f_B / Egamma;

            return f_symm + f_asym;
        }

        Diagnostics diagnostics() const
        {
            Diagnostics results;

            return results;
        }
    };

    AnalyticFormFactorBToGammaQCDF::AnalyticFormFactorBToGammaQCDF(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<AnalyticFormFactorBToGammaQCDF>(new Implementation<AnalyticFormFactorBToGammaQCDF>(p, o, *this))
    {
    }

    AnalyticFormFactorBToGammaQCDF::~AnalyticFormFactorBToGammaQCDF()
    {
    }

    FormFactors<PToGamma> *
    AnalyticFormFactorBToGammaQCDF::make(const Parameters & p, const Options & o)
    {
        return new AnalyticFormFactorBToGammaQCDF(p, o);
    }

    double
    AnalyticFormFactorBToGammaQCDF::f_v(const double & Egamma) const
    {
        return _imp->f_v(Egamma);
    }

    double
    AnalyticFormFactorBToGammaQCDF::f_a(const double & Egamma) const
    {
        return _imp->f_a(Egamma);
    }

    Diagnostics
    AnalyticFormFactorBToGammaQCDF::diagnostics() const
    {
        return _imp->diagnostics();
    }
}

