/* vim: set sw=4 sts=4 et tw=150 foldmethod=syntax : */

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

#ifndef EOS_GUARD_EOS_UTILS_OBSERVABLE_ARITHMETIC_HH
#define EOS_GUARD_EOS_UTILS_OBSERVABLE_ARITHMETIC_HH 1

#include <eos/observable.hh>

namespace eos
{
    /* Observable weighted with a constant number */
    class WeightedObservable :
        public Observable
    {
        private:
            const QualifiedName _name;
            const double _weight;
            const ObservablePtr _observable;

        public:
            /// Constructor.
            WeightedObservable(const QualifiedName & name, const double & weight, const ObservablePtr & observable);

            /// Destructor.
            virtual ~WeightedObservable() final;

            virtual const QualifiedName & name() const;

            virtual double evaluate() const;

            virtual Kinematics kinematics();

            virtual Parameters parameters();

            virtual Options options();

            virtual ObservablePtr clone() const;

            virtual ObservablePtr clone(const Parameters & parameters) const;
    };
}

#endif
