/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#include <eos/tau-decays/tau-to-3-leptons.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <string>

namespace eos
{
    using std::norm;

    /*
     * Decay: tau -> l1 l2 l3bar, cf. [FMvD2015]
     */
    template <>
    struct Implementation<TauToThreeLeptons>
    {
        std::shared_ptr<Model> model;

        UsedParameter hbar;

        UsedParameter g_fermi;

        UsedParameter m_tau;

        UsedParameter tau_tau;

        LeptonFlavorOption opt_l1;
        LeptonFlavorOption opt_l2;
        LeptonFlavorOption opt_l3;

        UsedParameter m_l1;
        UsedParameter m_l2;
        UsedParameter m_l3;

        UsedParameter mu;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            hbar(p["QM::hbar"], u),
            g_fermi(p["WET::G_Fermi"], u),
            m_tau(p["mass::tau"], u),
            tau_tau(p["life_time::tau"], u),
            opt_l1(o, options, "l1"),
            opt_l2(o, options, "l2"),
            opt_l3(o, options, "l3"),
            m_l1(p["mass::" + opt_l1.str()], u),
            m_l2(p["mass::" + opt_l2.str()], u),
            m_l3(p["mass::" + opt_l3.str()], u),
            mu(p[opt_l1.str() + "tau" + opt_l2.str() + opt_l3.str() + "::mu"], u)
        {
            u.uses(*model);
        }

        inline auto wc() const
        {
            return model->wet_ltaull(std::array<LeptonFlavor, 3>{{ opt_l1.value(), opt_l2.value(), opt_l3.value() }});
        }

        double decay_width() const
        {
            const WilsonCoefficients<wc::LLLL> wc = this->wc();

            // cf. [CNR:2013A], eq. (???)
            return 0.0;
        }

        double branching_ratio() const
        {
            return decay_width() * tau_tau / hbar;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<TauToThreeLeptons>::options
    {
        { "l1", { "e", "mu" }, "mu" },
        { "l2", { "e", "mu" }, "mu" },
        { "l3", { "e", "mu" }, "mu" }
    };

    TauToThreeLeptons::TauToThreeLeptons(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<TauToThreeLeptons>(new Implementation<TauToThreeLeptons>(parameters, options, *this))
    {
    }

    TauToThreeLeptons::~TauToThreeLeptons()
    {
    }

    double
    TauToThreeLeptons::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    const std::set<ReferenceName>
    TauToThreeLeptons::references
    {
        "CNR:2013A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    TauToThreeLeptons::begin_options()
    {
        return Implementation<TauToThreeLeptons>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    TauToThreeLeptons::end_options()
    {
        return Implementation<TauToThreeLeptons>::options.cend();
    }
}
