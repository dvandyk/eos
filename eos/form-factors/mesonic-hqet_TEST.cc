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

            // using z_* with a = 1.0
            {
                Parameters p = Parameters::Defaults();
                p["B(*)->D(*)::xi'(1)@HQET"]    = -1.5;
                p["B(*)->D(*)::xi''(1)@HQET"]   = +3.0;
                p["B(*)->D(*)::xi'''(1)@HQET"]  = +6.0;
                p["B(*)->D(*)::xi''''(1)@HQET"] = +0.0;
                p["B(*)->D(*)::a@HQET"]         =  1.0;

                HQETFormFactors<BToD, PToP> ff(p, Options{ });

                Diagnostics diag = ff.diagnostics();
                //for (auto d : diag)
                //{
                //    std::cout << d.description << ": " << d.value << std::endl;
                //}
                static const std::vector<std::pair<double, double>> ref
                {
                    /* Inputs */
                    std::make_pair(+0.217926, eps),
                    std::make_pair(+2.403316, eps),

                    /* z(w) */
                    std::make_pair(0.01219690, eps), // w = 1.10
                    std::make_pair(0.00617307, eps), // w = 1.05
                    std::make_pair(0.0,        eps), // w = 1.00

                    /* xi(w) */
                    std::make_pair(+1.706810, eps), // w = 2.10
                    std::make_pair(+0.764544, eps), // w = 1.60
                    std::make_pair(+0.865908, eps), // w = 1.10
                    std::make_pair(+0.928869, eps), // w = 1.05
                    std::make_pair(+1.000000, eps), // w = 1.00
                };

                // TEST_CHECK_DIAGNOSTICS(diag, ref);
            }

            {
                Parameters p = Parameters::Defaults();
                p["B(*)->D(*)::xi'(1)@HQET"]    = -1.5;
                p["B(*)->D(*)::xi''(1)@HQET"]   = +3.0;
                p["B(*)->D(*)::xi'''(1)@HQET"]  = +6.0;
                p["B(*)->D(*)::xi''''(1)@HQET"] = +0.0;
                p["B(*)->D(*)::a@HQET"]         =  1.25;

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

                    /* z(w) */
                    std::make_pair(-0.0354427, eps), // w = 1.10
                    std::make_pair(-0.0414582, eps), // w = 1.05
                    std::make_pair(-0.0476190, eps), // w = 1.00

                    /* xi(w) */
                    std::make_pair(+1.706810, eps), // w = 2.10
                    std::make_pair(+0.769847, eps), // w = 1.60
                    std::make_pair(+0.865911, eps), // w = 1.10
                    std::make_pair(+0.928869, eps), // w = 1.05
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
                    std::make_pair(+0.694295, eps), // C_{T_1}
                    std::make_pair(-0.931381, eps), // C_{T_2}
                    std::make_pair( 0.319615, eps), // C_{T_3}

                    /* WCs at (w = 1.0, z = 0.25) */
                    std::make_pair(+0.977157, eps), // C_{V_1}
                    std::make_pair(-0.478135, eps), // C_{V_2}
                    std::make_pair(-0.188532, eps), // C_{V_3}
                    std::make_pair(-0.356176, eps), // C_{A_1}
                    std::make_pair(-1.250411, eps), // C_{A_2}
                    std::make_pair( 0.381601, eps), // C_{A_3}
                    std::make_pair(+0.413987, eps), // C_{T_1}
                    std::make_pair(-0.956270, eps), // C_{T_2}
                    std::make_pair( 0.377063, eps), // C_{T_3}

                    /* HQET form factors at w = 1.4 */
                    std::make_pair(+0.706675, eps), // h_{p}
                    std::make_pair(-0.033688, eps), // h_{m}
                    std::make_pair(+0.822889, eps), // h_{T}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.837003, eps), // h_{p}
                    std::make_pair(-0.043229, eps), // h_{m}
                    std::make_pair(+0.979887, eps), // h_{T}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+1.036225, eps), // h_{p}
                    std::make_pair(-0.057905, eps), // h_{m}
                    std::make_pair(+1.219910, eps), // h_{T}
                };

                TEST_CHECK(false);
                //TEST_CHECK_DIAGNOSTICS(diag, ref);
            }
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
                std::make_pair(+0.623812, eps), // h_{A_3}
                std::make_pair(+0.772908, eps), // h_{V}
                std::make_pair(+0.614276, eps), // h_{T_1}
                std::make_pair(-0.136754, eps), // h_{T_2}
                std::make_pair(-0.149709, eps), // h_{T_3}

                /* HQET form factors at w = 1.2 */
                std::make_pair(+0.719796, eps), // h_{A_1}
                std::make_pair(-0.282823, eps), // h_{A_2}
                std::make_pair(+0.741269, eps), // h_{A_3}
                std::make_pair(+0.940146, eps), // h_{V}
                std::make_pair(+0.747488, eps), // h_{T_1}
                std::make_pair(-0.165357, eps), // h_{T_2}
                std::make_pair(-0.199636, eps), // h_{T_3}

                /* HQET form factors at w = 1.0 */
                std::make_pair(+0.889539, eps), // h_{A_1}
                std::make_pair(-0.387433, eps), // h_{A_2}
                std::make_pair(+0.916438, eps), // h_{A_3}
                std::make_pair(+1.196276, eps), // h_{V}
                std::make_pair(+0.951473, eps), // h_{T_1}
                std::make_pair(-0.209220, eps), // h_{T_2}
                std::make_pair(-0.280823, eps), // h_{T_3}
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
                std::make_pair(+0.680616, eps), // h_{Abar1}
                std::make_pair(-0.067144, eps), // h_{Abar2}
                std::make_pair(+0.758964, eps), // h_{Abar3}
                std::make_pair(+0.795991, eps), // h_{Vbar}

                /* HQET form factors at w = 1.2 */
                std::make_pair(+0.800729, eps), // h_{1}
                std::make_pair(-0.087100, eps), // h_{2}
                std::make_pair(+0.904178, eps), // h_{3}
                std::make_pair(+0.954095, eps), // h_{4}

                /* HQET form factors at w = 1.0 */
                std::make_pair(+0.981052, eps), // h_{Abar1}
                std::make_pair(-0.119119, eps), // h_{Abar2}
                std::make_pair(+1.125010, eps), // h_{Abar3}
                std::make_pair(+1.196010, eps), // h_{Vbar}
            };

            TEST_CHECK_DIAGNOSTICS(diag, ref);
        }
} bstar_to_d_hqet_form_factors_test;

class BstarToDstarHQETFormFactorsTest :
    public TestCase
{
    public:
        BstarToDstarHQETFormFactorsTest() :
            TestCase("bstar_to_dstar_hqet_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.0e-6;

            Parameters p = Parameters::Defaults();
            HQETFormFactors<BstarToDstar, VToV> ff(p, Options{ });

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
                std::make_pair(+0.600025, eps), // h_{1}
                std::make_pair(-0.095579, eps), // h_{2}
                std::make_pair(+0.864368, eps), // h_{3}
                std::make_pair(+0.810176, eps), // h_{4}
                std::make_pair(+0.127283, eps), // h_{5}
                std::make_pair(+0.041092, eps), // h_{6}
                std::make_pair(+0.558896, eps), // h_{7}
                std::make_pair(-0.126632, eps), // h_{8}
                std::make_pair(+0.163261, eps), // h_{9}
                std::make_pair(+0.051192, eps), // h_{10}

                /* HQET form factors at w = 1.2 */
                std::make_pair(+0.736173, eps), // h_{1}
                std::make_pair(-0.116910, eps), // h_{2}
                std::make_pair(+1.061682, eps), // h_{3}
                std::make_pair(+0.992678, eps), // h_{4}
                std::make_pair(+0.171034, eps), // h_{5}
                std::make_pair(+0.054643, eps), // h_{6}
                std::make_pair(+0.689381, eps), // h_{7}
                std::make_pair(-0.152735, eps), // h_{8}
                std::make_pair(+0.216679, eps), // h_{9}
                std::make_pair(+0.067720, eps), // h_{10}

                /* HQET form factors at w = 1.0 */
                std::make_pair(+0.944712, eps), // h_{1}
                std::make_pair(-0.149684, eps), // h_{2}
                std::make_pair(+1.364124, eps), // h_{3}
                std::make_pair(+1.272345, eps), // h_{4}
                std::make_pair(+0.242485, eps), // h_{5}
                std::make_pair(+0.076650, eps), // h_{6}
                std::make_pair(+0.889539, eps), // h_{7}
                std::make_pair(-0.192686, eps), // h_{8}
                std::make_pair(+0.303308, eps), // h_{9}
                std::make_pair(+0.094471, eps), // h_{10}
            };

            TEST_CHECK_DIAGNOSTICS(diag, ref);
        }
} bstar_to_dstar_hqet_form_factors_test;
