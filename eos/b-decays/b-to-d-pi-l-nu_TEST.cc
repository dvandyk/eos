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

#include <test/test.hh>
#include <eos/observable.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/b-decays/b-to-d-pi-l-nu.hh>

#include <eos/utils/integrate.hh>

using namespace test;
using namespace eos;

class BToDPiLeptonNeutrinoTest :
    public TestCase
{
    public:
        BToDPiLeptonNeutrinoTest() :
            TestCase("b_to_d_pi_l_nu_test")
        {
        }

        virtual void run() const
        {
            Parameters p = Parameters::Defaults();

            // reference values from Martin Jung as for the HQET form factor tests, for an l=mu lepton
            {
                p["mass::B_d"]                   =  5.279;
                p["mass::D^*_d"]                 =  2.0103;
                p["mass::mu"]                    =  0.10566;
                p["B(*)->D(*)::xi'(1)@HQET"]     = -1.5;
                p["B(*)->D(*)::xi''(1)@HQET"]    = +3.0;
                p["B(*)->D(*)::xi'''(1)@HQET"]   = +6.0;
                p["B(*)->D(*)::xi''''(1)@HQET"]  = -9.0;
                p["B(*)->D(*)::chi_2(1)@HQET"]   = +0.5;
                p["B(*)->D(*)::chi_2'(1)@HQET"]  = -1.0;
                p["B(*)->D(*)::chi_2''(1)@HQET"] = +2.0;
                p["B(*)->D(*)::chi_3'(1)@HQET"]  = -1.5;
                p["B(*)->D(*)::chi_3''(1)@HQET"] = +2.5;
                p["B(*)->D(*)::eta(1)@HQET"]     = +0.25;
                p["B(*)->D(*)::eta'(1)@HQET"]    = -1.25;
                p["B(*)->D(*)::eta''(1)@HQET"]   = +1.75;
                p["B(*)->D(*)::l_1(1)@HQET"]     = +0.5;
                p["B(*)->D(*)::l_2(1)@HQET"]     = -2.0;
                p["B(*)->D(*)::l_3(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_4(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_5(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_6(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::a@HQET"]          =  1.0000e+00;

                Options o{ { "z-order-lp", "3" }, { "z-order-slp", "1" }, { "l", "mu" } };

                const double eps = 1e-5;

                // Ensure that the form factors are constant across q^2
                auto ff = FormFactorFactory<PToV>::create("B->D^*::HQET", p, o);
                TEST_CHECK(ff != nullptr);

                BToDPiLeptonNeutrino d(p, o);
                auto diag = d.diagnostics();

                for (auto & d : diag)
                {
                    std::cerr << d.description << ": " << d.value << std::endl;
                }
                std::vector<std::pair<double, double>> ref
                {
                    // reference point
                    std::make_pair(  1.0,     1e-5), // q2^0 = 1.0

                    // amplitudes
                    std::make_pair( 17.35381, 1e-5), // A_time(q2 =  1.0)
                    std::make_pair(- 4.44119, 1e-5), // A_perp(q2 =  1.0)
                    std::make_pair(  7.91009, 1e-5), // A_para(q2 =  1.0)
                    std::make_pair( 18.00975, 1e-5), // A_long(q2 =  1.0)

                    std::make_pair(  1.90900, 1e-5), // A_time(q2 = 10.0)
                    std::make_pair(- 1.37388, 1e-5), // A_perp(q2 = 10.0)
                    std::make_pair(  8.09716, 1e-5), // A_para(q2 = 10.0)
                    std::make_pair(  6.01952, 1e-5), // A_long(q2 = 10.0)

                    // normalization
                    std::make_pair(  2.08096, 1e-5), // N(q2 =  1.0)
                    std::make_pair(  5.13467, 1e-5), // N(q2 = 10.0)

                    // unintegrated coefficients P(c_d)
                    std::make_pair(  172.205, 1e-3), // a(q2 =  1.0)
                    std::make_pair( 1206.240, 1e-3), // b(q2 =  1.0)
                    std::make_pair(  346.535, 1e-3), // a(q2 = 10.0)
                    std::make_pair(   25.842, 1e-3), // b(q2 = 10.0)

                    // unintegrated coefficients P(c_l)
                    std::make_pair( 1537.076, 1e-3), // a(q2 =  1.0)
                    std::make_pair(  263.374, 1e-3), // b(q2 =  1.0)
                    std::make_pair(-1165.513, 1e-3), // c(q2 =  1.0)
                    std::make_pair(  718.876, 1e-3), // a(q2 = 10.0)
                    std::make_pair(  228.220, 1e-3), // b(q2 = 10.0)
                    std::make_pair(-  25.736, 1e-3), // c(q2 = 10.0)

                    // unintegrated coefficients P(chi)
                    std::make_pair( 1811.022, 1e-3), // a(q2 =  1.0)
                    std::make_pair(- 176.327, 1e-3), // b(q2 =  1.0)
                    std::make_pair( 1392.038, 1e-3), // a(q2 = 10.0)
                    std::make_pair(- 653.185, 1e-3), // b(q2 = 10.0)
                };
                TEST_CHECK_DIAGNOSTICS(diag, ref);

                // distribution in cos(theta_D)
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-1.00), 0.804958, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.80), 0.640281, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.60), 0.512198, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.40), 0.420711, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.20), 0.365819, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d( 0.00), 0.347521, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+0.20), 0.365819, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+0.40), 0.420711, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+0.60), 0.512198, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+0.80), 0.640281, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+1.00), 0.804958, eps);

                // distribution in cos(theta_l)
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-1.00), 0.127192, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.80), 0.2517  , eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.60), 0.358581, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.40), 0.447833, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.20), 0.519456, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l( 0.00), 0.573451, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+0.20), 0.609818, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+0.40), 0.628557, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+0.60), 0.629667, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+0.80), 0.613149, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+1.00), 0.579003, eps);

                // distribution in chi
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(0.00000), 0.132799, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(0.62832), 0.151011, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(1.25664), 0.180477, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(1.88496), 0.180477, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(2.51327), 0.151011, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(3.14159), 0.132799, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(3.76991), 0.151011, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(4.39823), 0.180477, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(5.02655), 0.180477, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(5.65487), 0.15101 , eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(6.28319), 0.132799, eps);

                // normalization of the integrated distributions
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_d(-1.0,  0.0),         0.50000, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_d(-1.0, +1.0),         1.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_l(-1.0,  0.0),         0.38705, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_l(-1.0, +1.0),         1.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_chi( 0.0,       M_PI), 0.50000, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_chi( 0.0, 2.0 * M_PI), 1.0,     eps);

                // polarizations
                TEST_CHECK_NEARLY_EQUAL(d.lepton_polarization(),       0.98299, eps);
                TEST_CHECK_NEARLY_EQUAL(d.longitudinal_polarization(), 0.53664, eps);
            }

            // reference values from Martin Jung as for the HQET form factor tests, for an l=tau lepton
            {
                p["mass::B_d"]                   =  5.279;
                p["mass::D^*_d"]                 =  2.0103;
                p["mass::tau"]                   =  1.7768;
                p["B(*)->D(*)::xi'(1)@HQET"]     = -1.5;
                p["B(*)->D(*)::xi''(1)@HQET"]    = +3.0;
                p["B(*)->D(*)::xi'''(1)@HQET"]   = +6.0;
                p["B(*)->D(*)::xi''''(1)@HQET"]  = -9.0;
                p["B(*)->D(*)::chi_2(1)@HQET"]   = +0.5;
                p["B(*)->D(*)::chi_2'(1)@HQET"]  = -1.0;
                p["B(*)->D(*)::chi_2''(1)@HQET"] = +2.0;
                p["B(*)->D(*)::chi_3'(1)@HQET"]  = -1.5;
                p["B(*)->D(*)::chi_3''(1)@HQET"] = +2.5;
                p["B(*)->D(*)::eta(1)@HQET"]     = +0.25;
                p["B(*)->D(*)::eta'(1)@HQET"]    = -1.25;
                p["B(*)->D(*)::eta''(1)@HQET"]   = +1.75;
                p["B(*)->D(*)::l_1(1)@HQET"]     = +0.5;
                p["B(*)->D(*)::l_2(1)@HQET"]     = -2.0;
                p["B(*)->D(*)::l_3(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_4(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_5(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_6(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::a@HQET"]          =  1.0000e+00;

                Options o{ { "z-order-lp", "3" }, { "z-order-slp", "1" }, { "l", "tau" } };

                const double eps = 1e-5;

                // Ensure that the form factors are constant across q^2
                auto ff = FormFactorFactory<PToV>::create("B->D^*::HQET", p, o);
                TEST_CHECK(ff != nullptr);

                BToDPiLeptonNeutrino d(p, o);
                auto diag = d.diagnostics();

                for (auto & d : diag)
                {
                    std::cerr << d.description << ": " << d.value << std::endl;
                }
                std::vector<std::pair<double, double>> ref
                {
                    // reference point
                    std::make_pair(  4.0,     1e-5), // q2^0 = 1.0

                    // amplitudes
                    std::make_pair(  7.30668, 1e-5), // A_time(q2 =  4.0)
                    std::make_pair(- 3.52987, 1e-5), // A_perp(q2 =  4.0)
                    std::make_pair(  7.25697, 1e-5), // A_para(q2 =  4.0)
                    std::make_pair(  8.70693, 1e-5), // A_long(q2 =  4.0)

                    std::make_pair(  1.90900, 1e-5), // A_time(q2 = 10.0)
                    std::make_pair(- 1.37388, 1e-5), // A_perp(q2 = 10.0)
                    std::make_pair(  8.09716, 1e-5), // A_para(q2 = 10.0)
                    std::make_pair(  6.01952, 1e-5), // A_long(q2 = 10.0)

                    // normalization
                    std::make_pair(  0.30494, 1e-5), // N(q2 =  4.0)
                    std::make_pair(  2.40976, 1e-5), // N(q2 = 10.0)

                    // unintegrated coefficients P(c_d)
                    std::make_pair(   27.696, 1e-3), // a(q2 =  4.0)
                    std::make_pair(   75.333, 1e-3), // b(q2 =  4.0)
                    std::make_pair(  188.199, 1e-3), // a(q2 = 10.0)
                    std::make_pair(   22.318, 1e-3), // b(q2 = 10.0)

                    // unintegrated coefficients P(c_l)
                    std::make_pair(  107.467, 1e-3), // a(q2 =  4.0)
                    std::make_pair(-  30.001, 1e-3), // b(q2 =  4.0)
                    std::make_pair(-   5.559, 1e-3), // c(q2 =  4.0)
                    std::make_pair(  394.035, 1e-3), // a(q2 = 10.0)
                    std::make_pair(   72.261, 1e-3), // b(q2 = 10.0)
                    std::make_pair(-   8.274, 1e-3), // c(q2 = 10.0)

                    // unintegrated coefficients P(chi)
                    std::make_pair(  161.004, 1e-3), // a(q2 =  4.0)
                    std::make_pair(-   5.167, 1e-3), // b(q2 =  4.0)
                    std::make_pair(  691.918, 1e-3), // a(q2 = 10.0)
                    std::make_pair(- 210.004, 1e-3), // b(q2 = 10.0)
                };
                TEST_CHECK_DIAGNOSTICS(diag, ref);

                // distribution in cos(theta_D)
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-1.00), 0.693699, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.80), 0.589101, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.60), 0.507748, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.40), 0.449638, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.20), 0.414773, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d( 0.00), 0.403151, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+0.20), 0.414773, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+0.40), 0.449638, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+0.60), 0.507748, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+0.80), 0.589101, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+1.00), 0.693699, eps);

                // distribution in cos(theta_l)
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-1.00), 0.420344, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.80), 0.443628, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.60), 0.464316, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.40), 0.48241 , eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.20), 0.497909, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l( 0.00), 0.510812, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+0.20), 0.521121, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+0.40), 0.528835, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+0.60), 0.533954, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+0.80), 0.536477, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+1.00), 0.536406, eps);

                // distribution in chi
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(0.00000), 0.142349, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(0.62832), 0.153962, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(1.25664), 0.172751, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(1.88496), 0.172751, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(2.51327), 0.153962, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(3.14159), 0.142349, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(3.76991), 0.153962, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(4.39823), 0.172751, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(5.02655), 0.172751, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(5.65487), 0.153962, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi(6.28319), 0.142349, eps);

                // normalization of the integrated distributions
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_d(-1.0,  0.0),         0.50000, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_d(-1.0, +1.0),         1.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_l(-1.0,  0.0),         0.47098, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_l(-1.0, +1.0),         1.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_chi( 0.0,       M_PI), 0.50000, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_chi( 0.0, 2.0 * M_PI), 1.0,     eps);

                // polarizations
                TEST_CHECK_NEARLY_EQUAL(d.lepton_polarization(),       0.50376, eps);
                TEST_CHECK_NEARLY_EQUAL(d.longitudinal_polarization(), 0.46247, eps);
            }
        }
} b_to_d_pi_l_nu_test;
