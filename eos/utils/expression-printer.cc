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

#include <eos/utils/exception.hh>
#include <eos/utils/expression.hh>
#include <eos/utils/expression-fwd.hh>

#include <iostream>

namespace eos::exp
{
    class ExpressionPrinter
    {
        private:
            std::ostream & _os;

        public:
            ExpressionPrinter(std::ostream& ostream) : _os(ostream)
            {
            }

            void visit(BinaryExpression & e)
            {
                _os << "BinaryExpression(";
                e.lhs.accept(*this);
                _os << " " << e.op << " ";
                e.rhs.accept(*this);
                _os << ")";
            }

            void visit(ConstantExpression & e)
            {
                _os << "ConstantExpression(" << e.value << ")";
            }

            void visit(ObservableNameExpression & e)
            {
                _os << "ObservableNameExpression(" << e.observable_name << ")";
            }

            void visit(ObservableExpression &)
            {
                throw InternalError("Encountered ObserableExpression in ExpressionPrinter::visit");
            }
    };
}