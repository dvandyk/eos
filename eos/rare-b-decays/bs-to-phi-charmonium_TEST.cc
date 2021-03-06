/*
 * Copyright (c) 2021 Méril Reboud
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
#include <eos/rare-b-decays/bs-to-phi-charmonium.hh>
#include <eos/rare-b-decays/nonlocal-formfactors.hh>

using namespace test;
using namespace eos;

class BsToPhiCharmoniumGvDV2020Test :
    public TestCase
{
    public:
    BsToPhiCharmoniumGvDV2020Test() :
            TestCase("bs_to_phi_charmonium_GvDV2020_test")
        {
        }

        virtual void run() const
        {

            Parameters p = Parameters::Defaults();
            p["mass::B_s"]                             = 5.36688;
            p["mass::phi"]                             = 1.019461;
            p["mass::J/psi"]                           = 3.0969;
            p["mass::psi(2S)"]                         = 3.6860;
            p["mass::D^0"]                             = 1.86723;
            p["b->sccbar::t_0"]                        = 9.0;
            p["b->sccbar::t_s"]                        = -17.4724;
            p["b->sccbar::chiOPE@GvDV2020"]            = 1.81e-4;

            p["B_s->phiccbar::Re{alpha_0^perp}@GvDV2020"]  = 2.0;
            p["B_s->phiccbar::Im{alpha_0^perp}@GvDV2020"]  = 3.0;
            p["B_s->phiccbar::Re{alpha_1^perp}@GvDV2020"]  = 4.0;
            p["B_s->phiccbar::Im{alpha_1^perp}@GvDV2020"]  = 5.0;
            p["B_s->phiccbar::Re{alpha_2^perp}@GvDV2020"]  = 6.0;
            p["B_s->phiccbar::Im{alpha_2^perp}@GvDV2020"]  = 7.0;
            p["B_s->phiccbar::Re{alpha_0^para}@GvDV2020"]  = 8.0;
            p["B_s->phiccbar::Im{alpha_0^para}@GvDV2020"]  = 9.0;
            p["B_s->phiccbar::Re{alpha_1^para}@GvDV2020"]  = 10.0;
            p["B_s->phiccbar::Im{alpha_1^para}@GvDV2020"]  = 11.0;
            p["B_s->phiccbar::Re{alpha_2^para}@GvDV2020"]  = 12.0;
            p["B_s->phiccbar::Im{alpha_2^para}@GvDV2020"]  = 13.0;
            p["B_s->phiccbar::Re{alpha_0^long}@GvDV2020"]  = 14.0;
            p["B_s->phiccbar::Im{alpha_0^long}@GvDV2020"]  = 15.0;
            p["B_s->phiccbar::Re{alpha_1^long}@GvDV2020"]  = 16.0;
            p["B_s->phiccbar::Im{alpha_1^long}@GvDV2020"]  = 17.0;
            p["B_s->phiccbar::Re{alpha_2^long}@GvDV2020"]  = 18.0;
            p["B_s->phiccbar::Im{alpha_2^long}@GvDV2020"]  = 19.0;


            Options oo;
            oo.set("model",               "WilsonScan");
            oo.set("nonlocal-formfactor", "GvDV2020");
            oo.set("psi",                 "J/psi");

            BsToPhiCharmonium c(p, oo);

            TEST_CHECK_NEARLY_EQUAL(c.branching_ratio(),  13598663., 1.);

        }
} bs_to_phi_charmonium_GvDV2020_test;


class BsToPhiCharmoniumGRvDV2021Test :
    public TestCase
{
    public:
    BsToPhiCharmoniumGRvDV2021Test() :
            TestCase("bs_to_phi_charmonium_GRvDV2021_test")
        {
        }

        virtual void run() const
        {

            Parameters p = Parameters::Defaults();
            p["mass::B_s"]                              = 5.36688;
            p["mass::phi"]                              = 1.019461;
            p["mass::J/psi"]                            = 3.0969;
            p["mass::psi(2S)"]                          = 3.6860;
            p["mass::D^0"]                              = 1.86723;
            p["b->sccbar::t_0"]                         = 9.0;
            p["b->sccbar::t_s"]                         = -17.4724;
            p["b->sccbar::chiOPE@GRvDV2021"]            = 1.81e-4;

            p["B_s->phiccbar::Re{alpha_0^perp}@GRvDV2021"]  = 2.0;
            p["B_s->phiccbar::Im{alpha_0^perp}@GRvDV2021"]  = 3.0;
            p["B_s->phiccbar::Re{alpha_1^perp}@GRvDV2021"]  = 4.0;
            p["B_s->phiccbar::Im{alpha_1^perp}@GRvDV2021"]  = 5.0;
            p["B_s->phiccbar::Re{alpha_2^perp}@GRvDV2021"]  = 6.0;
            p["B_s->phiccbar::Im{alpha_2^perp}@GRvDV2021"]  = 7.0;
            p["B_s->phiccbar::Re{alpha_0^para}@GRvDV2021"]  = 8.0;
            p["B_s->phiccbar::Im{alpha_0^para}@GRvDV2021"]  = 9.0;
            p["B_s->phiccbar::Re{alpha_1^para}@GRvDV2021"]  = 10.0;
            p["B_s->phiccbar::Im{alpha_1^para}@GRvDV2021"]  = 11.0;
            p["B_s->phiccbar::Re{alpha_2^para}@GRvDV2021"]  = 12.0;
            p["B_s->phiccbar::Im{alpha_2^para}@GRvDV2021"]  = 13.0;
            p["B_s->phiccbar::Re{alpha_0^long}@GRvDV2021"]  = 14.0;
            p["B_s->phiccbar::Im{alpha_0^long}@GRvDV2021"]  = 15.0;
            p["B_s->phiccbar::Re{alpha_1^long}@GRvDV2021"]  = 16.0;
            p["B_s->phiccbar::Im{alpha_1^long}@GRvDV2021"]  = 17.0;
            p["B_s->phiccbar::Re{alpha_2^long}@GRvDV2021"]  = 18.0;
            p["B_s->phiccbar::Im{alpha_2^long}@GRvDV2021"]  = 19.0;

            Options oo;
            oo.set("model",               "WilsonScan");
            oo.set("nonlocal-formfactor", "GRvDV2021");
            oo.set("psi",                 "J/psi");

            BsToPhiCharmonium c(p, oo);

            TEST_CHECK_NEARLY_EQUAL(c.branching_ratio(),  973155.8, .1);

        }
} bs_to_phi_charmonium_GRvDV2021_test;