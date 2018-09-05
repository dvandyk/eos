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

#include <test/test.hh>
#include <eos/observable.hh>
#include <eos/b-decays/b-to-d-pi-l-nu.hh>

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

            // compared to Danny's Mathematica notebook
            {
                p["B->Dpilnu::A_FB" ] = 0.120;
                p["B->Dpilnu::A_T^2"] =-0.230;
                p["B->Dpilnu::A_im" ] = 0.340;
                p["B->Dpilnu::A_S"  ] =-0.045;
                p["B->Dpilnu::F_S"  ] = 0.130;
                p["B->Dpilnu::F_L"  ] = 0.570;

                Options o;

                BToDPiLeptonNeutrino d(p, o);

                const double eps = 1e-5;

                // distribution in cos(theta_D)
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.70), 0.60408, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.50), 0.48389, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d( 0.00), 0.34558, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d( 0.30), 0.37377, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d( 0.90), 0.68033, eps);

                // distribution in cos(theta_L)
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.80), 0.31555, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.40), 0.51529, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l( 0.00), 0.60971, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l( 0.20), 0.61743, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l( 0.70), 0.52152, eps);

                // distribution in phi
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_phi( 1.5708),  0.16600, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_phi( 2.35619), 0.11208, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_phi( 3.14159), 0.15231, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_phi( 3.45575), 0.18129, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_phi( 5.49779), 0.11208, eps);

                // normalization of the integrated distributions
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_d(-1.0,  0.0),     0.52250, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_d(-1.0, +1.0),     1.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_l(-1.0,  0.0),     0.44780, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_l(-1.0, +1.0),     1.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_phi( 0.0,  +M_PI), 0.50000, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_phi(-M_PI, +M_PI), 1.0,     eps);
            }
        }
} b_to_d_pi_l_nu_test;
