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

#include <eos/observable.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/expression.hh>
#include <eos/utils/expression-fwd.hh>

namespace eos::exp
{
    class ExpressionEvaluator
    {
        public:
            double visit(BinaryExpression & e)
            {
                BinaryExpression::func f = BinaryExpression::Method(e.op);

                return f(e.lhs.accept_returning<double>(*this), e.rhs.accept_returning<double>(*this));
            }

            double visit(ConstantExpression & e)
            {
                return e.value;
            }

            double visit(ObservableExpression & e)
            {
                return e.observable->evaluate();
            }

            double visit(ObservableNameExpression &)
            {
                throw InternalError("Encountered ObserableNameExpression in ExpressionEvaluator::visit");

                return 0.0;
            }
    };
}