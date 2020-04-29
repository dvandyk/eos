/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020 Danny van Dyk
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
#include <eos/utils/observable-arithmetic.hh>

using namespace test;
using namespace eos;

class WeightedObservableTest :
    public TestCase
{
    public:
        WeightedObservableTest() :
            TestCase("weighted_observable_test")
        {
        }

        virtual void run() const
        {
            {
                Parameters p = Parameters::Defaults();
                p.declare("test::param_1", 17.0);

                auto o = Observable::make("test::param_1", p, Kinematics{ }, Options{ });
                TEST_CHECK_EQUAL(+17.0, o->evaluate());

                ObservablePtr wo{ new WeightedObservable("test::param_1[*-1.0]", -1.0, o) };
                TEST_CHECK_EQUAL(-17.0, wo->evaluate());

                p["test::param_1"] = -9.0;
                TEST_CHECK_EQUAL(+9.0, wo->evaluate());
            }
        }
} weighted_observable_test;
