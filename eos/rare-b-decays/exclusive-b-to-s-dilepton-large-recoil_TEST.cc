/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013, 2014 Danny van Dyk
 * Copyright (c) 2014 Frederik Beaujean
 * Copyright (c) 2014 Christoph Bobeth
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
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton-large-recoil.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using namespace test;
using namespace eos;

class BToKDileptonLargeRecoilBobethCompatibilityTest :
    public TestCase
{
    public:
        BToKDileptonLargeRecoilBobethCompatibilityTest() :
            TestCase("b_to_k_dilepton_large_recoil_bobeth_compatibility_test")
        {
        }

        virtual void run() const
        {
            // Christoph uses \Delta C instead of C for C9, C10
            // important to agree on alpha_s, can change values by 1%
            Parameters p = Parameters::Defaults();
            p["c1"] = -0.3231323312;
            p["c2"] = 1.009301831;
            p["c3"] = -0.005233499106;
            p["c4"] = -0.08829686414;
            p["c5"] = 0.0003601965805;
            p["c6"] = 0.001020749573;
            p["Re{c7}"] = -0.3370422989 + 0.1;
            p["Im{c7}"] = 0.2;
            p["Re{c7'}"] = 0.3;
            p["Im{c7'}"] = 0.4;
            p["c8"] = -0.1827530948;
            p["Re{c9}"] = 4.294489364 + 1;
            p["Im{c9}"] = 0.5;
            p["Re{c9'}"] = 2;
            p["Im{c9'}"] = 1.5;
            p["Re{c10}"] = -4.196294696 + 3;
            p["Im{c10}"] = 2.5;
            p["Re{c10'}"] = 4;
            p["Im{c10'}"] = 3.5;
            p["Re{cS}"] = 0.5;
            p["Im{cS}"] = 1;
            p["Re{cS'}"] = 0.6;
            p["Im{cS'}"] = 1.1;
            p["Re{cP}"] = 0.7;
            p["Im{cP}"] = 1.2;
            p["Re{cP'}"] = 0.8;
            p["Im{cP'}"] = 1.3;
            p["Re{cT}"] = 0.9;
            p["Im{cT}"] = 1.4;
            p["Re{cT5}"] = 1.0;
            p["Im{cT5}"] = 1.5;

            Options oo;
            oo.set("model", "WilsonScan");
            oo.set("scan-mode", "cartesian");
            oo.set("form-factors", "KMPW2010");
            oo.set("l", "mu");
            oo.set("q", "u");

            BToKDilepton<LargeRecoil> d(p, oo);
            double eps = 1e-3;
            static const double s = 6.0;

            TEST_CHECK_RELATIVE_ERROR_C(d.F_A(s), complex<double>(2.803705304, 6), 1e-14);
            TEST_CHECK_RELATIVE_ERROR_C(d.F_S(s), complex<double>(3.277235546, 6.256540588), eps);
            TEST_CHECK_RELATIVE_ERROR_C(d.F_T(s), complex<double>(7.695315895, 11.97049139), eps);
            TEST_CHECK_RELATIVE_ERROR_C(d.F_T5(s), complex<double>(8.550350995, 12.82552649), eps);
            TEST_CHECK_RELATIVE_ERROR_C(d.F_P(s), complex<double>(4.010492477, 6.467135768), eps);

            /* difference comes from cal_T, F_V affects everything below */
            TEST_CHECK_RELATIVE_ERROR(std::real(d.F_V(s)), 7.756362368, eps);
            TEST_CHECK_RELATIVE_ERROR(std::imag(d.F_V(s)), 3.191642172, 6 * eps);

            eps *= 2.5;
            TEST_CHECK_RELATIVE_ERROR(d.a_l(s),  3.92053702e-20, eps);
            TEST_CHECK_RELATIVE_ERROR(d.b_l(s),  9.694697008e-21, eps);
            TEST_CHECK_RELATIVE_ERROR(d.c_l(s), -2.756810607e-20, eps);

            const double tau_over_hbar = p["life_time::B_u"] / p["hbar"];
            TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(1, 6),
                                      2.898727023e-19 * tau_over_hbar, eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio_cp_averaged(1, 6),
                                      2.8855929e-19 * tau_over_hbar, eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry(1, 6), 0.1097985735, eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_flat_term(1, 6), 0.2788261376, eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_ratio_muons_electrons(1, 6), 1.073039657, eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_cp_asymmetry(1, 6), 0.00455162022, 5 * eps);
        }
} b_to_k_dilepton_large_recoil_bobeth_compatibility_test;
