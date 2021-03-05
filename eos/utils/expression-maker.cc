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
#include <eos/utils/expression.hh>
#include <eos/utils/expression-fwd.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/qualified-name.hh>

#include <iostream>

namespace eos::exp
{

    // copies one expression tree to a new tree, replacing ObservableNameExpression objects with ObservableExpression objects
    class ExpressionMaker
    {
        private:
            Parameters _parameters;
            Kinematics _kinematics;
            Options    _options;

        public:
            ExpressionMaker(const Parameters & parameters, const Kinematics & kinematics, const Options & options) :
                _parameters(parameters),
                _kinematics(kinematics),
                _options(options)
            {
            }

            Expression visit(BinaryExpression & e)
            {
                return BinaryExpression(e.op, e.lhs.accept_returning<Expression>(*this), e.rhs.accept_returning<Expression>(*this));
            }

            Expression visit(ConstantExpression & e)
            {
                return e;
            }

            Expression visit(ObservableNameExpression & e)
            {
                return ObservableExpression(Observable::make(
                    QualifiedName(e.observable_name),
                    this->_parameters,
                    this->_kinematics,
                    this->_options
                    )
                );
            }

            Expression visit(ObservableExpression &)
            {
                throw InternalError("Encountered ObserableExpression in ExpressionMaker::visit");
            }
    };
}
