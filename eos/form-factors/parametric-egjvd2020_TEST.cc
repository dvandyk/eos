/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Danny van Dyk
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
#include <eos/form-factors/parametric-egjvd2020.hh>
#include <eos/form-factors/parametric-egjvd2020-impl.hh> 

#include <eos/utils/model.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class ParametricEGJvD2020Test :
    public TestCase
{
    public:
        ParametricEGJvD2020Test() :
            TestCase("parametric_egjvd2020_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            p["0->pipi::a_+^0@EGJvD2020"] = 0.5;
            p["0->pipi::a_+^1@EGJvD2020"] = 0.25;
            p["0->pipi::a_+^2@EGJvD2020"] = 0.125;
            p["0->pipi::a_+^3@EGJvD2020"] = 0.0625;
            p["0->pipi::a_+^4@EGJvD2020"] = 0.03125;
            p["0->pipi::a_+^5@EGJvD2020"] = 0.015625;
            p["0->pipi::a_+^6@EGJvD2020"] = 7.8125e-3;
            p["0->pipi::a_+^7@EGJvD2020"] = 3.90625e-3;
            p["0->pipi::a_+^8@EGJvD2020"] = 1.953125e-3;
            p["0->pipi::a_+^9@EGJvD2020"] = 9,765625e-4;

            /* f_+ */
            {
                std::shared_ptr<FormFactors<VacuumToPP>> ff = FormFactorFactory<VacuumToPP>::create("0->pipi::EGJvD2020", p, Options{ });

                // Formfactor (@ q2 = -1)
                TEST_CHECK_NEARLY_EQUAL( ,     real(ff->f_p(-1.0)),   eps); 
                TEST_CHECK_NEARLY_EQUAL( 0.0,                   imag(ff->f_p(-1.0)),   eps); 
                // Formfactor (@ q2 = -0.5)
                TEST_CHECK_NEARLY_EQUAL( ,      real(ff->f_p(-0.5)),   eps); 
                TEST_CHECK_NEARLY_EQUAL( 0.0,                   imag(ff->f_p(-0.8)),   eps); 
                // Formfactor (@ q2 =  0.0)
                TEST_CHECK_NEARLY_EQUAL( ,    real(ff->f_p(0.0)),   eps); 
                TEST_CHECK_NEARLY_EQUAL( 0.0,                   imag(ff->f_p(-0.6)),   eps); 
                // Formfactor (@ q2 =  0.5)
                TEST_CHECK_NEARLY_EQUAL( ,     real(ff->f_p(0.5)),   eps); 
                TEST_CHECK_NEARLY_EQUAL( 0.0,                   imag(ff->f_p(-0.4)),   eps); 
                // Formfactor (@ q2 =  1.0)
                TEST_CHECK_NEARLY_EQUAL( ,    real(ff->f_p(1.0)),   eps); 
                TEST_CHECK_NEARLY_EQUAL( 0.0,                   imag(ff->f_p(-0.2)),   eps); 
            }
        }
} parametric_egjvd2020_test;
