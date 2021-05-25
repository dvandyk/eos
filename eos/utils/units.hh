/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_UNITS_HH
#define EOS_GUARD_EOS_UTILS_UNITS_HH 1

#include <eos/utils/instantiation_policy.hh>

#include <string>

namespace eos
{
    class Unit :
        public InstantiationPolicy<NonCopyable>
    {
        private:
            const std::string _latex;

        protected:
            Unit(const std::string & latex);
            virtual ~Unit();

        public:
            const std::string & latex() const { return _latex; }
    };

    namespace units
    {
        class Undefined :
            public Unit
        {
            public:
                Undefined();
                ~Undefined() = default;
        }

        class None :
            public Unit
        {
            public:
                None();
                ~None() = default;
        }

        class GeV :
            public Unit
        {
            public:
                GeV();
                ~GeV() = default;
        }

        class InversePicoSecond :
            public Unit
        {
            public:
                InversePicoSecond();
                ~InversePicoSecond() = default;
        }
    }
}

#endif