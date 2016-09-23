/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2015, 2016 Danny van Dyk
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
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton-low-recoil.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

// Uncomment the following #define to generate new test data for the Bobeth compatibility tests
//#define EOS_GENERATE_TEST_DATA

using namespace test;
using namespace eos;

class BToKDileptonLowRecoilTest :
    public TestCase
{
    public:
        BToKDileptonLowRecoilTest() :
            TestCase("b_to_k_dilepton_low_recoil_test")
        {
        }

        virtual void run() const
        {
            /* Low Recoil (SM) */
            {
                Parameters p = Parameters::Defaults();
                p["life_time::B_d"] = 1.530e-12;
                p["b->s::c1"] = -0.32300000;
                p["b->s::c2"] = +1.00931000;
                p["b->s::c3"] = -0.00522869;
                p["b->s::c4"] = -0.08794730;
                p["b->s::c5"] = +0.00037476;
                p["b->s::c6"] = +0.00105859;
                p["b->s::Re{c7}"] = -0.331;
                p["b->s::c8"] = -0.181;
                p["b->smumu::Re{c9}"] = 4.27;
                p["b->smumu::Re{c10}"] = -4.17;
                // PDG 2008 CKM parameters
                p["CKM::A"] = 0.814;
                p["CKM::lambda"] = 0.2257;
                p["CKM::rhobar"] = 0.135;
                p["CKM::etabar"] = 0.349;
                // Kaon mass
                p["mass::K_d"] = 0.49761;
                // B mass
                p["mass::B_d"] = 5.27953;
                p["mass::b(MSbar)"] = 4.2;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BZ2004v2");

                BToKDilepton<LowRecoil> d(p, oo);

                /* q^2 = [14.18, 22.8] */
                {
                    const double eps = 1e-5;
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(14.18  ), 2.498607492e-08, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(15.2575), 2.454634936e-08, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(16.335 ), 2.374070832e-08, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(17.4125), 2.244476570e-08, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(18.49  ), 2.046333084e-08, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(19.5675), 1.749827352e-08, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(20.645 ), 1.312607210e-08, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(21.7225), 6.929281810e-09, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(22.8   ), 1.579971652e-10, eps);

                    TEST_CHECK_RELATIVE_ERROR(d.differential_flat_term(15.0), 0.005562348378, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_flat_term(22.0), 0.008213626852, eps);

                    TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(14.18, 22.8),       1.5267386e-07, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_flat_term(14.18, 22.8),             5.4236817e-03, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_ratio_muons_electrons(14.18, 22.8), 1.0015589,     eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_cp_asymmetry(14.18, 22.8),          2.2706273e-05, eps);
                }
            }

            /* Benchmark Point */
            {
                Parameters p = Parameters::Defaults();
                p["life_time::B_d"] = 1.530e-12;
                // PDG 2008 CKM parameters
                p["CKM::A"] = 0.814;
                p["CKM::lambda"] = 0.2257;
                p["CKM::rhobar"] = 0.135;
                p["CKM::etabar"] = 0.349;
                // B mass
                p["mass::B_d"] = 5.27953;
                // Kaon mass
                p["mass::K_d"] = 0.49761;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;
                p["b->s::c1"] = -0.32300000;
                p["b->s::c2"] = +1.00931000;
                p["b->s::c3"] = -0.00522869;
                p["b->s::c4"] = -0.08794730;
                p["b->s::c5"] = +0.00037476;
                p["b->s::c6"] = +0.00105859;
                p["b->s::Re{c7}"] = 0.0;
                p["b->s::Im{c7}"] = -0.331;
                p["b->s::c8"] = -0.181;
                p["b->smumu::Re{c9}"] = 0.0;
                p["b->smumu::Im{c9}"] = +4.27;
                p["b->smumu::Re{c10}"] = 0.0;
                p["b->smumu::Im{c10}"] = -4.17;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("l", "mu");
                oo.set("form-factors", "BZ2004v2");

                BToKDilepton<LowRecoil> d(p, oo);

                /* q^2 = [14.18, 22.8] */
                {
                    const double eps = 1e-5;

                    TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(14.18, 22.8),              1.5520940e-07, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio_cp_averaged(14.18, 22.8),  1.4629637e-07, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_flat_term(14.18, 22.8),                    5.3935506e-03, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_ratio_muons_electrons(14.18, 22.8),        1.0015315,     eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_cp_asymmetry(14.18, 22.8),                 0.0609245,     eps);
                }
            }
        }
} b_to_k_dilepton_low_recoil_test;

class BToKDileptonLowRecoilBobethCompatibilityTest :
    public TestCase
{
    public:
        BToKDileptonLowRecoilBobethCompatibilityTest() :
            TestCase("b_to_k_dilepton_low_recoil_bobeth_compatibility_test")
        {
        }

        virtual void run() const
        {
            static const std::vector<std::string> variation_names
            {
                "b->s::Re{c7}",      "b->s::Im{c7}",      "b->s::Re{c7'}",      "b->s::Im{c7'}",
                "b->smumu::Re{c9}",  "b->smumu::Im{c9}",  "b->smumu::Re{c9'}",  "b->smumu::Im{c9'}",
                "b->smumu::Re{c10}", "b->smumu::Im{c10}", "b->smumu::Re{c10'}", "b->smumu::Im{c10'}",
            };

            Parameters p = Parameters::Defaults();
            // old test data generated for K^+ mass set to K0 mass
            p["mass::K_u"] = 0.497614;
            Options o;
            o.set("model", "WilsonScan");
            o.set("form-factors", "KMPW2010");

            std::vector<Parameter> variations;

            for (auto n = variation_names.cbegin(), n_end = variation_names.cend() ; n != n_end ; ++n)
            {
                variations.push_back(p[*n]);
            }

            Kinematics k;
            k.declare("s_min"); k.set("s_min", 14.18);
            k.declare("s_max"); k.set("s_max", 22.86);

            std::vector<ObservablePtr> observables;
            std::vector<std::string> observable_names = {
                    "B->Kll::BR@LowRecoil,q=u,l=mu",
                    "B->Kll::F_H@LowRecoil,q=u,l=mu",
            };
            for (auto & s : observable_names)
            {
                observables.push_back(Observable::make(s, p, k, o));
                TEST_CHECK_MSG(observables.back(), "Could not create '" + s + "'");
            }

            std::string filename(EOS_SRCDIR "/eos/rare-b-decays/exclusive-b-to-s-dilepton-low-recoil_TEST-btokll.data");
#ifdef EOS_GENERATE_TEST_DATA
            {
                std::cout << "-- GENERATING test case data for B->Kll@LowRecoil --" << std::endl;
                RandomNumberGenerator rng;
                std::fstream file(filename.c_str(), std::fstream::out);
                file.precision(17);

                for (int i = 0 ; i < 1000 ; ++i)
                {
                    for (auto v = variations.begin(), v_end = variations.end() ; v != v_end ; ++v)
                    {
                        *v = v->min() + (v->max() - v->min()) * rng();

                        file << *v << '\t';
                    }

                    for (auto o = observables.cbegin(), o_end = observables.cend() ; o != o_end ; ++o)
                    {
                        file << (*o)->evaluate() << '\t';
                    }
                    file << std::endl;
                }
            }
#else
            // Verify the test case data
            {
                std::cout << "-- Verifying test case data for B->Kll@LowRecoil --" << std::endl;
                std::fstream file(filename.c_str(), std::fstream::in);
                TEST_CHECK_MSG(file, "'" + filename + "' does not exist");
                std::string line;
                while (file)
                {
                    std::getline(file, line);
                    if (line.empty())
                        break;

                    std::stringstream ss(line);

                    for (auto v = variations.begin(), v_end = variations.end() ; v != v_end ; ++v)
                    {
                        double value;
                        ss >> value;
                        *v = value;
                    }

                    for (auto o = observables.cbegin(), o_end = observables.cend() ; o != o_end ; ++o)
                    {
                        double reference;
                        ss >> reference;

                        TEST_CHECK_RELATIVE_ERROR(reference, (*o)->evaluate(), 1e-3);
                    }
                }
            }
#endif
        }
} b_to_k_dilepton_low_recoil_bobeth_compatibility_test;
