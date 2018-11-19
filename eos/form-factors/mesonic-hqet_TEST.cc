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
#include <eos/form-factors/form-factors.hh>
#include <eos/form-factors/mesonic-hqet.hh>
#include <eos/form-factors/mesonic-impl.hh>

#include <vector>

using namespace test;
using namespace eos;

class BToDHQETFormFactorsTest :
    public TestCase
{
    public:
        BToDHQETFormFactorsTest() :
            TestCase("b_to_d_hqet_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.0e-6;

            Parameters p = Parameters::Defaults();
            HQETFormFactors<BToD, PToP> ff(p, Options{ });

            Diagnostics diag = ff.diagnostics();
            for (auto d : diag)
            {
                std::cout << d.description << ": " << d.value << std::endl;
            }
            static const std::vector<std::pair<double, double>> ref
            {
                /* Inputs */
                std::make_pair(+0.217926, eps),
                std::make_pair(+2.403316, eps),

                /* xi(w) */
                std::make_pair(+0.881073, eps), // w = 1.10
                std::make_pair(+0.937683, eps), // w = 1.05
                std::make_pair(+1.000000, eps), // w = 1.00

                /* chi2(w) */
                std::make_pair(-0.0581,   eps), // w = 1.10
                std::make_pair(-0.05805,  eps), // w = 1.05
                std::make_pair(-0.058,    eps), // w = 1.00

                /* chi3(w) */
                std::make_pair(+0.0035,   eps), // w = 1.10
                std::make_pair(+0.00175,  eps), // w = 1.05
                std::make_pair(+0.0,      eps), // w = 1.00

                /* eta(w) */
                std::make_pair(+0.3624,   eps), // w = 1.10
                std::make_pair(+0.3602,   eps), // w = 1.05
                std::make_pair(+0.358,    eps), // w = 1.00

                /* r(w) */
                std::make_pair(+0.967945, eps), // w = 1.1
                std::make_pair(+0.999767, eps), // w = 1.0007
                std::make_pair(+0.999967, eps), // w = 1.0001
                std::make_pair(+0.999983, eps), // w = 1.00005
                std::make_pair(+1.0,      eps), // w = 1.0

                /* Omega(w, z = 0.25) */
                std::make_pair(+1.294026, eps), // w = 1.1
                std::make_pair(+1.310389, eps), // w = 1.0007
                std::make_pair(+1.310476, eps), // w = 1.0001
                std::make_pair(+1.310483, eps), // w = 1.00005
                std::make_pair(+1.310491, eps), // w = 1.0

                /* Omega(w, z = 0.20) */
                std::make_pair(+1.403808, eps), // w = 1.1
                std::make_pair(+1.414099, eps), // w = 1.0007
                std::make_pair(+1.414149, eps), // w = 1.0001
                std::make_pair(+1.414153, eps), // w = 1.00005
                std::make_pair(+1.414157, eps), // w = 1.0

                /* WCs at (w = 1.2, z = 0.20) */
                std::make_pair(+1.123905, eps), // C_{V_1}
                std::make_pair(-0.454499, eps), // C_{V_2}
                std::make_pair(-0.162046, eps), // C_{V_3}
                std::make_pair(-0.127091, eps), // C_{A_1}
                std::make_pair(-1.247185, eps), // C_{A_2}
                std::make_pair( 0.316106, eps), // C_{A_3}

                /* WCs at (w = 1.0, z = 0.25) */
                std::make_pair(+0.977157, eps), // C_{V_1}
                std::make_pair(-0.478135, eps), // C_{V_2}
                std::make_pair(-0.188532, eps), // C_{V_3}
                std::make_pair(-0.356176, eps), // C_{A_1}
                std::make_pair(-1.250411, eps), // C_{A_2}
                std::make_pair( 0.381601, eps), // C_{A_3}

                /* HQET form factors at w = 1.4 */
                std::make_pair(+0.706675, eps), // h_{p}
                std::make_pair(-0.033688, eps), // h_{m}

                /* HQET form factors at w = 1.2 */
                std::make_pair(+0.837003, eps), // h_{p}
                std::make_pair(-0.043229, eps), // h_{m}

                /* HQET form factors at w = 1.0 */
                std::make_pair(+1.036225, eps), // h_{p}
                std::make_pair(-0.057905, eps), // h_{m}
            };

            TEST_CHECK_DIAGNOSTICS(diag, ref);
        }
} b_to_d_hqet_form_factors_test;

class BToDstarHQETFormFactorsTest :
    public TestCase
{
    public:
        BToDstarHQETFormFactorsTest() :
            TestCase("b_to_dstar_hqet_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.0e-6;

            Parameters p = Parameters::Defaults();
            HQETFormFactors<BToDstar, PToV> ff(p, Options{ });

            Diagnostics diag = ff.diagnostics();
            for (auto d : diag)
            {
                std::cout << d.description << ": " << d.value << std::endl;
            }
            static const std::vector<std::pair<double, double>> ref
            {
                /* Inputs */
                std::make_pair(+0.217926, eps),
                std::make_pair(+2.403316, eps),

                /* xi(w) */
                std::make_pair(+0.642745, eps), // w = 1.40
                std::make_pair(+0.783551, eps), // w = 1.20
                std::make_pair(+0.881073, eps), // w = 1.10
                std::make_pair(+0.937683, eps), // w = 1.05
                std::make_pair(+1.000000, eps), // w = 1.00

                /* chi2(w) */
                std::make_pair(-0.0581,   eps), // w = 1.10
                std::make_pair(-0.05805,  eps), // w = 1.05
                std::make_pair(-0.058,    eps), // w = 1.00

                /* chi3(w) */
                std::make_pair(+0.0035,   eps), // w = 1.10
                std::make_pair(+0.00175,  eps), // w = 1.05
                std::make_pair(+0.0,      eps), // w = 1.00

                /* eta(w) */
                std::make_pair(+0.3624,   eps), // w = 1.10
                std::make_pair(+0.3602,   eps), // w = 1.05
                std::make_pair(+0.358,    eps), // w = 1.00

                /* r(w) */
                std::make_pair(+0.967945, eps), // w = 1.1
                std::make_pair(+0.999767, eps), // w = 1.0007
                std::make_pair(+0.999967, eps), // w = 1.0001
                std::make_pair(+0.999983, eps), // w = 1.00005
                std::make_pair(+1.0,      eps), // w = 1.0

                /* Omega(w, z = 0.25) */
                std::make_pair(+1.294026, eps), // w = 1.1
                std::make_pair(+1.310389, eps), // w = 1.0007
                std::make_pair(+1.310476, eps), // w = 1.0001
                std::make_pair(+1.310483, eps), // w = 1.00005
                std::make_pair(+1.310491, eps), // w = 1.0

                /* Omega(w, z = 0.20) */
                std::make_pair(+1.403808, eps), // w = 1.1
                std::make_pair(+1.414099, eps), // w = 1.0007
                std::make_pair(+1.414149, eps), // w = 1.0001
                std::make_pair(+1.414153, eps), // w = 1.00005
                std::make_pair(+1.414157, eps), // w = 1.0

                /* WCs at (w = 1.2, z = 0.20) */
                std::make_pair(+1.123905, eps), // C_{V_1}
                std::make_pair(-0.454499, eps), // C_{V_2}
                std::make_pair(-0.162046, eps), // C_{V_3}
                std::make_pair(-0.127091, eps), // C_{A_1}
                std::make_pair(-1.247185, eps), // C_{A_2}
                std::make_pair( 0.316106, eps), // C_{A_3}

                /* WCs at (w = 1.0, z = 0.25) */
                std::make_pair(+0.977157, eps), // C_{V_1}
                std::make_pair(-0.478135, eps), // C_{V_2}
                std::make_pair(-0.188532, eps), // C_{V_3}
                std::make_pair(-0.356176, eps), // C_{A_1}
                std::make_pair(-1.250411, eps), // C_{A_2}
                std::make_pair( 0.381601, eps), // C_{A_3}

                /* HQET form factors at w = 1.4 */
                std::make_pair(+0.605958, eps), // h_{A_1}
                std::make_pair(-0.217704, eps), // h_{A_2}
                std::make_pair(+0.672385, eps), // h_{A_3}
                std::make_pair(+0.821481, eps), // h_{V}

                /* HQET form factors at w = 1.2 */
                std::make_pair(+0.719796, eps), // h_{A_1}
                std::make_pair(-0.282823, eps), // h_{A_2}
                std::make_pair(+0.800482, eps), // h_{A_3}
                std::make_pair(+0.999360, eps), // h_{V}

                /* HQET form factors at w = 1.0 */
                std::make_pair(+0.889539, eps), // h_{A_1}
                std::make_pair(-0.387433, eps), // h_{A_2}
                std::make_pair(+0.992008, eps), // h_{A_3}
                std::make_pair(+1.271847, eps), // h_{V}
            };

            TEST_CHECK_DIAGNOSTICS(diag, ref);
        }
} b_to_dstar_hqet_form_factors_test;

class BstarToDHQETFormFactorsTest :
    public TestCase
{
    public:
        BstarToDHQETFormFactorsTest() :
            TestCase("bstar_to_d_hqet_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.0e-6;

            Parameters p = Parameters::Defaults();
            HQETFormFactors<BstarToD, VToP> ff(p, Options{ });

            Diagnostics diag = ff.diagnostics();
            for (auto d : diag)
            {
                std::cout << d.description << ": " << d.value << std::endl;
            }
            static const std::vector<std::pair<double, double>> ref
            {
                /* Inputs */
                std::make_pair(+0.217926, eps),
                std::make_pair(+2.403316, eps),

                /* xi(w) */
                std::make_pair(+0.642745, eps), // w = 1.40
                std::make_pair(+0.783551, eps), // w = 1.20
                std::make_pair(+0.881073, eps), // w = 1.10
                std::make_pair(+0.937683, eps), // w = 1.05
                std::make_pair(+1.000000, eps), // w = 1.00

                /* chi2(w) */
                std::make_pair(-0.0581,   eps), // w = 1.10
                std::make_pair(-0.05805,  eps), // w = 1.05
                std::make_pair(-0.058,    eps), // w = 1.00

                /* chi3(w) */
                std::make_pair(+0.0035,   eps), // w = 1.10
                std::make_pair(+0.00175,  eps), // w = 1.05
                std::make_pair(+0.0,      eps), // w = 1.00

                /* eta(w) */
                std::make_pair(+0.3624,   eps), // w = 1.10
                std::make_pair(+0.3602,   eps), // w = 1.05
                std::make_pair(+0.358,    eps), // w = 1.00

                /* r(w) */
                std::make_pair(+0.967945, eps), // w = 1.1
                std::make_pair(+0.999767, eps), // w = 1.0007
                std::make_pair(+0.999967, eps), // w = 1.0001
                std::make_pair(+0.999983, eps), // w = 1.00005
                std::make_pair(+1.0,      eps), // w = 1.0

                /* Omega(w, z = 0.25) */
                std::make_pair(+1.294026, eps), // w = 1.1
                std::make_pair(+1.310389, eps), // w = 1.0007
                std::make_pair(+1.310476, eps), // w = 1.0001
                std::make_pair(+1.310483, eps), // w = 1.00005
                std::make_pair(+1.310491, eps), // w = 1.0

                /* Omega(w, z = 0.20) */
                std::make_pair(+1.403808, eps), // w = 1.1
                std::make_pair(+1.414099, eps), // w = 1.0007
                std::make_pair(+1.414149, eps), // w = 1.0001
                std::make_pair(+1.414153, eps), // w = 1.00005
                std::make_pair(+1.414157, eps), // w = 1.0

                /* WCs at (w = 1.2, z = 0.20) */
                std::make_pair(+1.123905, eps), // C_{V_1}
                std::make_pair(-0.454499, eps), // C_{V_2}
                std::make_pair(-0.162046, eps), // C_{V_3}
                std::make_pair(-0.127091, eps), // C_{A_1}
                std::make_pair(-1.247185, eps), // C_{A_2}
                std::make_pair( 0.316106, eps), // C_{A_3}

                /* WCs at (w = 1.0, z = 0.25) */
                std::make_pair(+0.977157, eps), // C_{V_1}
                std::make_pair(-0.478135, eps), // C_{V_2}
                std::make_pair(-0.188532, eps), // C_{V_3}
                std::make_pair(-0.356176, eps), // C_{A_1}
                std::make_pair(-1.250411, eps), // C_{A_2}
                std::make_pair( 0.381601, eps), // C_{A_3}

                /* HQET form factors at w = 1.4 */
                std::make_pair(+0.670369, eps), // h_{Abar_1}
                std::make_pair(-0.067144, eps), // h_{Abar_2}
                std::make_pair(+0.748717, eps), // h_{Abar_3}
                std::make_pair(+0.785744, eps), // h_{Vbar}

                /* HQET form factors at w = 1.2 */
                std::make_pair(+0.788237, eps), // h_{Abar_1}
                std::make_pair(-0.087100, eps), // h_{Abar_2}
                std::make_pair(+0.891686, eps), // h_{Abar_3}
                std::make_pair(+0.941604, eps), // h_{Vbar}

                /* HQET form factors at w = 1.0 */
                std::make_pair(+0.965109, eps), // h_{Abar_1}
                std::make_pair(-0.119119, eps), // h_{Abar_2}
                std::make_pair(+1.109067, eps), // h_{Abar_3}
                std::make_pair(+1.180068, eps), // h_{Vbar}
            };

            TEST_CHECK_DIAGNOSTICS(diag, ref);
        }
} bstar_to_d_hqet_form_factors_test;
