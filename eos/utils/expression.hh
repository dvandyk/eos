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

#ifndef EXPRESSION_HH
#define EXPRESSION_HH 1

#include <eos/observable-fwd.hh>
#include <eos/utils/expression-fwd.hh>

#include <cassert>
#include <memory>
#include <iostream>
#include <map>

namespace eos::exp
{
    class BinaryExpression
    {
        public:
            char op;
            Expression lhs, rhs;

            using func = double(*)(const double &, const double &);

            static double sum(const double &, const double &);
            static double difference(const double &, const double &);
            static double product(const double &, const double &);
            static double ratio(const double &, const double &);

            static BinaryExpression::func Method(char op);

            BinaryExpression() {}

            BinaryExpression(char op, const Expression & l, const Expression & r) :
                op(op), lhs(l), rhs(r)
            {
            }
    };

    class ConstantExpression
    {
        public:
            double value;

            ConstantExpression(const double& v = 0) :
                value(v)
            {
            }
    };

    class ObservableNameExpression
    {
        public:
            std::string observable_name;

            ObservableNameExpression(const std::string & input = "") :
                 observable_name(input)
            {
            }

    };

    class ObservableExpression
    {
        public:
            ObservablePtr observable;

            ObservableExpression(ObservablePtr observable) :
                 observable(observable)
            {
            }
    };
}

#endif
