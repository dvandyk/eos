/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2022 Danny van Dyk
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

#include <eos/observable-impl.hh>
#include <eos/tau-decays/tau-to-3-leptons.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    // Leptonic tau decays
    // {{{
    ObservableGroup
    make_tau_to_3_l_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $tau^-\to \ell_1^-\ell_2^+\ell_3^-$ decays)",
            R"()",
            {
                make_observable("tau->mumumu::BR", R"(\mathcal{B}(\tau^- \to \mu^-\mu^+\mu^-))",
                        Unit::None(),
                        &TauToThreeLeptons::branching_ratio,
                        std::make_tuple(),
                        Options{ { "l1", "mu" }, { "l2", "mu" }, { "l3", "mu" } })
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    ObservableSection
    make_tau_decays_section()
    {
        auto imp = new Implementation<ObservableSection>(
            "Observables in $\tau$ decays",
            "",
            {
                // tau^- -> l1^- l2^+ l3^-
                make_tau_to_3_l_group()
            }
        );

        return ObservableSection(imp);
    }
}
