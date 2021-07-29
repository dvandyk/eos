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

#include <eos/utils/expression.hh>
#include <eos/utils/expression-evaluator.cc>
#include <eos/utils/expression-fwd.hh>
#include <eos/utils/expression-maker.cc>
#include <eos/utils/expression-parser-impl.hh>
#include <eos/utils/expression-printer.cc>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>

#include <iostream>

using namespace test;
using namespace eos::exp;
using namespace eos;

class ExpressionTest
{
    private:
        const std::string _input;

    public:
        Expression e;
        bool completed;

        using It = std::string::const_iterator;
        ExpressionParser<It> p;

        ExpressionTest(const std::string input) :
            _input(input)
        {
            It f(_input.begin()), l(_input.end());
            completed = qi::phrase_parse(f, l, p, ascii::space, e);
        }
};


class ExpressionParserTest :
    public TestCase
{
public:
    ExpressionParserTest() :
    TestCase("expression_parser_test")
    {
    }

    virtual void run() const
    {
        // testing basic parsing of constants and binary expressions
        {
            ExpressionTest test("1+2*3");

            std::stringstream out;
            ExpressionPrinter printer(out);
            test.e.accept(printer);

            ExpressionEvaluator evaluator;

            TEST_CHECK(test.completed);
            TEST_CHECK_EQUAL(test.e.accept_returning<double>(evaluator), 7.0);
            TEST_CHECK_EQUAL("BinaryExpression(ConstantExpression(1) + BinaryExpression(ConstantExpression(2) * ConstantExpression(3)))", out.str());
        }

        // testing parsing and evaluation of observables
        {
            ExpressionTest test2("{B->Dlnu::BR;l=tau}[q2_min=>q2_min_tau] / {B->Dlnu::BR;l=mu}[q2_min=0.0]");

            std::stringstream out;
            ExpressionPrinter printer(out);
            test2.e.accept(printer);

            ExpressionEvaluator evaluator;

            TEST_CHECK(test2.completed);
            TEST_CHECK_EQUAL_STR(
                "BinaryExpression(ObservableNameExpression(B->Dlnu::BR;l=tau, aliases=[q2_min=>q2_min_tau]) / ObservableNameExpression(B->Dlnu::BR;l=mu, values=[q2_min=0]))",
                out.str()
            );

            // Cannot evaluate expression with ObservableNameExpression objects
            TEST_CHECK_THROWS(InternalError, test2.e.accept_returning<double>(evaluator));

            Parameters p = Parameters::Defaults();
            Kinematics k = Kinematics({{"q2_min", 0.011164}, {"q2_max", 11.62}});
            Options o;
            ExpressionMaker maker(p, k, o);
            Expression evaluable_test2 = test2.e.accept_returning<Expression>(maker);

            // This test is to be improved with correct kinematics and default parameters
            TEST_CHECK_RELATIVE_ERROR(evaluable_test2.accept_returning<double>(evaluator), 0.295627626, 1e-5);
        }

    }
} expression_parser_test;
