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
            p["0->pipi::a_+^0@EGJvD2020"] = 0.5; // enter four agreed upon parameters of a test formfactor
            p["0->pipi::a_+^1@EGJvD2020"] = 0.25;
            p["0->pipi::a_+^2@EGJvD2020"] = 0.125;
            p["0->pipi::a_+^3@EGJvD2020"] = 0.0625;

            /* f_+ */
            {
                std::shared_ptr<FormFactors<VacuumToPP>> ff = FormFactorFactory<VacuumToPP>::create("0->pipi::EGJvD2020", p, Options{ });

                
                //TEST_CHECK_NEARLY_EQUAL( 0.17,      real(ff->f_p(0.2)),   eps); //TODO: enter proper values
                //TEST_CHECK_NEARLY_EQUAL( 0.17,      imag(ff->f_p(0.2)),   eps); //TODO: enter proper values
                
                TEST_CHECK_NEARLY_EQUAL( -0.13104714470397733,      real(ff->f_p(0.5)),   eps); //TODO: enter proper values
                TEST_CHECK_NEARLY_EQUAL( +0.7154340963332151,      imag(ff->f_p(0.5)),   eps); //TODO: enter proper values
                
                //TEST_CHECK_NEARLY_EQUAL( 0.17,      real(ff->f_p(0.7)),   eps); //TODO: enter proper values
                //TEST_CHECK_NEARLY_EQUAL( 0.17,      imag(ff->f_p(0.7)),   eps); //TODO: enter proper values
                
            }
        }
        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            p["0->pipi::a_+^0@EGJvD2020"] = 1.0;
            p["0->pipi::a_+^1@EGJvD2020"] = 0.0;
            p["0->pipi::a_+^2@EGJvD2020"] = 0.0;
            p["0->pipi::a_+^3@EGJvD2020"] = 0.0;

            /* f_+ */
            {
                std::shared_ptr<FormFactors<VacuumToPP>> ff = FormFactorFactory<VacuumToPP>::create("0->pipi::EGJvD2020", p, Options{ });

                // Formfactor (@ q2 = -1): f(q2) = (8.176250218394324+0j)
                TEST_CHECK_NEARLY_EQUAL( 8.176250218394324,     real(ff->f_p(-1.0)),   eps); 
                TEST_CHECK_NEARLY_EQUAL( 0.0,                   imag(ff->f_p(-1.0)),   eps); 
                // Formfactor (@ q2 = -0.8): f(q2) = (9.18334561548279+0j)
                TEST_CHECK_NEARLY_EQUAL( 9.18334561548279,      real(ff->f_p(-0.8)),   eps); 
                TEST_CHECK_NEARLY_EQUAL( 0.0,                   imag(ff->f_p(-0.8)),   eps); 
                // Formfactor (@ q2 = -0.6): f(q2) = (10.661459069642154+0j)
                TEST_CHECK_NEARLY_EQUAL( 10.661459069642154,    real(ff->f_p(-0.6)),   eps); 
                TEST_CHECK_NEARLY_EQUAL( 0.0,                   imag(ff->f_p(-0.6)),   eps); 
                // Formfactor (@ q2 = -0.4): f(q2) = (13.09240490482987+0j)
                TEST_CHECK_NEARLY_EQUAL( 13.09240490482987,     real(ff->f_p(-0.4)),   eps); 
                TEST_CHECK_NEARLY_EQUAL( 0.0,                   imag(ff->f_p(-0.4)),   eps); 
                // Formfactor (@ q2 = -0.2): f(q2) = (18.052699030274198+0j)
                TEST_CHECK_NEARLY_EQUAL( 18.052699030274198,    real(ff->f_p(-0.2)),   eps); 
                TEST_CHECK_NEARLY_EQUAL( 0.0,                   imag(ff->f_p(-0.2)),   eps); 
                // Formfactor (@ q2 = 0): f(q2) = (37.09615019059266+0j)
                TEST_CHECK_NEARLY_EQUAL( 37.09615019059266,     real(ff->f_p(0.0)),   eps); 
                TEST_CHECK_NEARLY_EQUAL( 0.0,                   imag(ff->f_p(0.0)),   eps); 
                // Formfactor (@ q2 = 0.2): f(q2) = (-18.879501240331354-39.955515781364156j)
                TEST_CHECK_NEARLY_EQUAL( -18.879501240331354,   real(ff->f_p(0.2)),   eps); 
                TEST_CHECK_NEARLY_EQUAL( -39.955515781364156,   imag(ff->f_p(0.2)),   eps);
                // Formfactor (@ q2 = 0.4): f(q2) = (-11.01256109617313-13.987296979144421j)
                TEST_CHECK_NEARLY_EQUAL( -11.01256109617313,    real(ff->f_p(0.4)),   eps); 
                TEST_CHECK_NEARLY_EQUAL( -13.987296979144421,   imag(ff->f_p(0.4)),   eps);
                // Formfactor (@ q2 = 0.6): f(q2) = (-6.689600624311777-8.639330287012191j)
                TEST_CHECK_NEARLY_EQUAL( -6.689600624311777,    real(ff->f_p(0.6)),   eps); 
                TEST_CHECK_NEARLY_EQUAL( -8.639330287012191,    imag(ff->f_p(0.6)),   eps);
                // Formfactor (@ q2 = 0.8): f(q2) = (-4.5182490211790265-6.518648082523539j)
                TEST_CHECK_NEARLY_EQUAL( -4.5182490211790265,   real(ff->f_p(0.8)),   eps); 
                TEST_CHECK_NEARLY_EQUAL( -6.518648082523539,    imag(ff->f_p(0.8)),   eps);
                // Formfactor (@ q2 = 1): f(q2) = (-3.2518312169642507-5.388943873463791j)
                TEST_CHECK_NEARLY_EQUAL( -3.2518312169642507,   real(ff->f_p(1.0)),   eps); 
                TEST_CHECK_NEARLY_EQUAL( -5.388943873463791,    imag(ff->f_p(1.0)),   eps);
            }
        }
} parametric_egjvd2020_test;
