
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

#include <eos/utils/instantiation_policy-impl.hh>

#include <string>

namespace eos
{
    Unit::Unit(const std::string & latex) :
        _latex(latex)
    {
    }

    Unit::~Unit() = default;

    namespace units
    {
        Undefined::Undefined() :
            Unit(R"(\textrm{undefined})")
        {
        }

        None::None() :
            Unit("1")
        {
        }

        GeV::Gev() :
            Unit(R"($\textrm{GeV}$)")
        {
        }

        InversePicoSecond::InversePicoSecond() :
            Unit(R"(\textrm{ps}^{-1})")
        {
        }
    }
}