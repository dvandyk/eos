/* vim: set sw=4 sts=4 et foldmethod=marker : */

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

#include <eos/observable.hh>
#include <eos/utils/observable-arithmetic.hh>

namespace eos
{
    // {{{
    WeightedObservable::WeightedObservable(const QualifiedName & name, const double & weight, const ObservablePtr & observable) :
        _name(name),
        _weight(weight),
        _observable(observable)
    {
    }

    WeightedObservable::~WeightedObservable() = default;

    const QualifiedName &
    WeightedObservable::name() const
    {
        return _name;
    }

    double
    WeightedObservable::evaluate() const
    {
        return _weight * _observable->evaluate();
    }

    Kinematics
    WeightedObservable::kinematics()
    {
        return _observable->kinematics();
    }

    Parameters
    WeightedObservable::parameters()
    {
        return _observable->parameters();
    }

    Options
    WeightedObservable::options()
    {
        return _observable->options();
    }

    ObservablePtr
    WeightedObservable::clone() const
    {
        return ObservablePtr{ new WeightedObservable(_name, _weight, _observable->clone()) };
    }

    ObservablePtr
    WeightedObservable::clone(const Parameters & parameters) const
    {
        return ObservablePtr{ new WeightedObservable(_name, _weight, _observable->clone(parameters)) };
    }
    // }}}
}
