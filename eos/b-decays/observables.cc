/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2019 Danny van Dyk
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
#include <eos/b-decays/b-to-l-nu.hh>
#include <eos/b-decays/b-to-pi-l-nu.hh>
#include <eos/b-decays/b-to-pi-pi-l-nu.hh>
#include <eos/b-decays/b-to-d-l-nu.hh>
#include <eos/b-decays/b-to-dstar-l-nu.hh>
#include <eos/b-decays/bs-to-kstar-l-nu.hh>
#include <eos/b-decays/lambdab-to-lambdac-l-nu.hh>
#include <eos/b-decays/lambdab-to-lambdac2595-l-nu.hh>
#include <eos/b-decays/lambdab-to-lambdac2625-l-nu.hh>
#include <eos/b-decays/inclusive-b-to-u.hh>
#include <eos/b-decays/properties.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    // Leptonic B decays
    // {{{
    ObservableGroup
    make_b_to_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B^-\to \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavour.)",
            {
                make_observable("B_u->lnu::BR", R"(\mathcal{B}(B^- \to \ell^-\bar\nu))",
                        &BToLeptonNeutrino::branching_ratio),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // Semileptonic B -> P(seudoscalar) decays
    // {{{

    // B -> pi l nu
    // {{{
    ObservableGroup
    make_b_to_pi_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B\to \pi \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavour. The option "q" selects the spectator quark flavour. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("B->pilnu::dBR/dq2", R"(d\mathcal{B}(B\to\pi\ell^-\bar\nu)/dq^2)",
                        &BToPiLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("B->pilnu::BR", R"(\mathcal{B}(B\to\pi\ell^-\bar\nu))",
                        &BToPiLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->pilnu::zeta",
                        &BToPiLeptonNeutrino::integrated_zeta,
                        std::make_tuple("q2_min", "q2_max")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> D l nu
    // {{{
    ObservableGroup
    make_b_to_d_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B\to \bar{D} \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavour. The option "q" selects the spectator quark flavour. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("B->Dlnu::dBR/dq2", R"(d\mathcal{B}(B\to \bar{D}\ell^-\bar\nu)/dq^2)",
                        &BToDLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("B->Dlnu::BR", R"(\mathcal{B}(B\to \bar{D}\ell^-\bar\nu))",
                        &BToDLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Dlnu::normdBR/ds",
                        &BToDLeptonNeutrino::normalized_differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("B->Dlnu::normBR",
                        &BToDLeptonNeutrino::normalized_integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Dlnu::R_D(q2)",
                        &BToDLeptonNeutrino::differential_r_d,
                        std::make_tuple("q2")),

                make_observable("B->Dlnu::R_D",
                        &BToDLeptonNeutrino::integrated_r_d,
                        std::make_tuple("q2_mu_min", "q2_tau_min", "q2_mu_max", "q2_tau_max")),

                make_observable("B->Dlnu::A_FB(q2)", R"(A_{\text{FB}}(B\to \bar{D}\ell^-\bar\nu)(q^2))",
                        &BToDLeptonNeutrino::differential_a_fb_leptonic,
                        std::make_tuple("q2")),

                make_observable("B->Dlnu::A_FB", R"(A_{\text{FB}}(B\to \bar{D}\ell^-\bar\nu))",
                        &BToDLeptonNeutrino::integrated_a_fb_leptonic,
                        std::make_tuple("q2_min", "q2_max")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // Semileptonic B -> V(seudoscalar) decays
    // {{{

    // B -> D^* l nu
    // {{{
    ObservableGroup
    make_b_to_dstar_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B\to \bar{D}^* \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavour. The option "q" selects the spectator quark flavour. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                // B -> D^* l nu
                make_observable("B->D^*lnu::dBR/dq2", R"(d\mathcal{B}(B\to \bar{D}^*\ell^-\bar\nu)/dq^2)",
                                &BToDstarLeptonNeutrino::differential_branching_ratio,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::normdBR/dq2",
                                &BToDstarLeptonNeutrino::normalized_differential_branching_ratio,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::BR", R"(\mathcal{B}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToDstarLeptonNeutrino::integrated_branching_ratio,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::normBR",
                                &BToDstarLeptonNeutrino::normalized_integrated_branching_ratio,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::R_D^*(q2)",
                                &BToDstarLeptonNeutrino::differential_ratio_tau_mu,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::R_D^*",
                                &BToDstarLeptonNeutrino::integrated_ratio_tau_mu,
                                std::make_tuple("q2_mu_min", "q2_tau_min", "q2_mu_max", "q2_tau_max")),

                make_observable("B->D^*lnu::A_FB(q2)", R"(A_{\text{FB}}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                &BToDstarLeptonNeutrino::differential_a_fb_leptonic,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::A_FB", R"(A_{\text{FB}}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToDstarLeptonNeutrino::integrated_a_fb_leptonic,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::A_L",
                                &BToDstarLeptonNeutrino::integrated_amplitude_polarization_L,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::A_T",
                                &BToDstarLeptonNeutrino::integrated_amplitude_polarization_T,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::F_L", R"(F_{\text{L}}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToDstarLeptonNeutrino::integrated_f_L,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::A_C^1", R"(A_{\text{C}}^1(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToDstarLeptonNeutrino::integrated_a_c_1,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::A_C^2", R"(A_{\text{C}}^2(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToDstarLeptonNeutrino::integrated_a_c_2,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::A_C^3", R"(A_{\text{C}}^3(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToDstarLeptonNeutrino::integrated_a_c_3,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::A_T^1", R"(A_{\text{T}}^1(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToDstarLeptonNeutrino::integrated_a_t_1,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::A_T^2", R"(A_{\text{T}}^2(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToDstarLeptonNeutrino::integrated_a_t_2,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::A_T^3", R"(A_{\text{T}}^3(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToDstarLeptonNeutrino::integrated_a_t_3,
                                std::make_tuple("q2_min", "q2_max")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B_s -> K^* l nu
    // {{{
    ObservableGroup
    make_bs_to_kstar_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_s\to \bar{K}^* \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavour. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                // B_s -> K^* l nubar
                make_observable("B_s->K^*lnu::F_perp(q2)",
                        &BsToKstarLeptonNeutrino::Fperp,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::F_para(q2)",
                        &BsToKstarLeptonNeutrino::Fpara,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::F_long(q2)",
                        &BsToKstarLeptonNeutrino::Flong,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::d^4Gamma",
                        &BsToKstarLeptonNeutrino::four_differential_decay_width,
                        std::make_tuple("q2", "cos(theta_l)", "cos(theta_k)", "phi")),

                make_observable("B_s->K^*lnu::dBR/ds", R"(d\mathcal{B}(B_s\to \bar{K}^*\ell^-\bar\nu)/dq^2)",
                        &BsToKstarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::A_FB(q2)",
                        &BsToKstarLeptonNeutrino::differential_forward_backward_asymmetry,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::BR",
                        &BsToKstarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::A_FB",
                        &BsToKstarLeptonNeutrino::integrated_forward_backward_asymmetry,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::Shat_1s",
                        &BsToKstarLeptonNeutrino::integrated_s_1s,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::Shat_1c",
                        &BsToKstarLeptonNeutrino::integrated_s_1c,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::Shat_2s",
                        &BsToKstarLeptonNeutrino::integrated_s_2s,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::Shat_2c",
                        &BsToKstarLeptonNeutrino::integrated_s_2c,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::Shat_3",
                        &BsToKstarLeptonNeutrino::integrated_s_3,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::Shat_4",
                        &BsToKstarLeptonNeutrino::integrated_s_4,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::Shat_5",
                        &BsToKstarLeptonNeutrino::integrated_s_5,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::Shat_6s",
                        &BsToKstarLeptonNeutrino::integrated_s_6s,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::A_FB(q2)",
                        &BsToKstarLeptonNeutrino::differential_forward_backward_asymmetry,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::A_T^2(q2)",
                        &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_2,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::A_T^3(q2)",
                        &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_3,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::A_T^4(q2)",
                        &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_4,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::A_T^5(q2)",
                        &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_5,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::A_T^re(q2)",
                        &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_re,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::A_T^im(q2)",
                        &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_im,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::F_L(q2)",
                        &BsToKstarLeptonNeutrino::differential_longitudinal_polarisation,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::F_T(q2)",
                        &BsToKstarLeptonNeutrino::differential_transversal_polarisation,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::H_T^1(q2)",
                        &BsToKstarLeptonNeutrino::differential_h_1,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::H_T^2(q2)",
                        &BsToKstarLeptonNeutrino::differential_h_2,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::H_T^3(q2)",
                        &BsToKstarLeptonNeutrino::differential_h_3,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::H_T^4(q2)",
                        &BsToKstarLeptonNeutrino::differential_h_4,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::H_T^5(q2)",
                        &BsToKstarLeptonNeutrino::differential_h_5,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::A_FB",
                        &BsToKstarLeptonNeutrino::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::BR",
                        &BsToKstarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::F_L",
                        &BsToKstarLeptonNeutrino::integrated_longitudinal_polarisation,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::F_T",
                        &BsToKstarLeptonNeutrino::integrated_transversal_polarisation,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::A_T^2",
                        &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_2,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::A_T^3",
                        &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_3,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::A_T^4",
                        &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_4,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::A_T^5",
                        &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_5,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::A_T^re",
                        &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_re,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::A_T^im",
                        &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_im,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::H_T^1",
                        &BsToKstarLeptonNeutrino::integrated_h_1,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::H_T^2",
                        &BsToKstarLeptonNeutrino::integrated_h_2,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::H_T^3",
                        &BsToKstarLeptonNeutrino::integrated_h_3,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::H_T^4",
                        &BsToKstarLeptonNeutrino::integrated_h_4,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::H_T^5",
                        &BsToKstarLeptonNeutrino::integrated_h_5,
                        std::make_tuple("q2_min", "q2_max")),

                // B_s -> K^* l nubar Ratios
                make_observable("B_s->K^*lnu::R_long",
                        &BsToKstarLeptonNeutrinoRatios::ratio_long),

                make_observable("B_s->K^*lnu::R_para",
                        &BsToKstarLeptonNeutrinoRatios::ratio_para),

                make_observable("B_s->K^*lnu::R_perp",
                        &BsToKstarLeptonNeutrinoRatios::ratio_perp),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // Semileptonic B -> P(seudoscalar) P(seudoscalar) decays
    // {{{

    // B -> pi pi l nu
    // {{{
    ObservableGroup
    make_b_to_pi_pi_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B\to \pi\pi \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavour. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("B->pipilnu::BR(q2,k2)", R"(d^2\mathcal{B}(B\to \pi\pi \ell^-\bar\nu)/(dq^2\,dk^2))",
                        &BToPiPiLeptonNeutrino::double_differential_branching_ratio,
                        std::make_tuple("q2", "k2")),

                make_observable("B->pipilnu::BR(q2,k2,cos(theta_pi))",
                        &BToPiPiLeptonNeutrino::triple_differential_branching_ratio,
                        std::make_tuple("q2", "k2", "cos(theta_pi)")),

                make_observable("B->pipilnu::A_FB(q2,k2)", R"(A_{\text{FB}}(B\to \pi\pi \ell^-\bar\nu)(q^2,k^2))",
                        &BToPiPiLeptonNeutrino::double_differential_forward_backward_asymmetry,
                        std::make_tuple("q2", "k2")),

                make_observable("B->pipilnu::P(cos(theta_pi))",
                        &BToPiPiLeptonNeutrino::partial_waves,
                        std::make_tuple("q2", "k2", "cos(theta_pi)")),

                make_observable("B->pipilnu::BR", R"(\mathcal{B}(B\to \pi\pi \ell^-\bar\nu))",
                        &BToPiPiLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max", "k2_min", "k2_max", "z_min", "z_max")),

                make_observable("B->pipilnu::A_FB", R"(A_{\text{FB}}(B\to \pi\pi \ell^-\bar\nu))",
                        &BToPiPiLeptonNeutrino::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max", "k2_min", "k2_max")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // Semileptonic Lambda_b decays
    // {{{

    // Lambda_b -> Lambda_c l nu
    // {{{
    ObservableGroup
    make_lambdab_to_lambdac_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $\Lambda_b\to \Lambda_c \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavour. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                // Lambda_b -> Lambda_c l nu
                make_observable("Lambda_b->Lambda_clnu::dBR/dq2", R"(d\mathcal{B}(\Lambda_b\to\Lambda_c \ell^-\bar\nu)/dq^2)",
                        &LambdaBToLambdaCLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::A_FB^l(q2)", R"(A_{\text{FB}}^\ell(\Lambda_b\to\Lambda_c \ell^-\bar\nu)(q^2))",
                        &LambdaBToLambdaCLeptonNeutrino::differential_a_fb_leptonic,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::A_FB^h(q2)", R"(A_{\text{FB}}^h(\Lambda_b\to\Lambda_c \ell^-\bar\nu)(q^2))",
                        &LambdaBToLambdaCLeptonNeutrino::differential_a_fb_hadronic,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::A_FB^c(q2)", R"(A_{\text{FB}}^{h\ell}(\Lambda_b\to\Lambda_c \ell^-\bar\nu)(q^2))",
                        &LambdaBToLambdaCLeptonNeutrino::differential_a_fb_combined,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::F_0(q2)", R"(F_0(\Lambda_b\to\Lambda_c \ell^-\bar\nu)(q^2))",
                        &LambdaBToLambdaCLeptonNeutrino::differential_fzero,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::BR", R"(\mathcal{B}(\Lambda_b\to\Lambda_c \ell^-\bar\nu))",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::R_Lambda_c(q2)",
                        &LambdaBToLambdaCLeptonNeutrino::differential_ratio_tau_mu,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::R_Lambda_c",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_ratio_tau_mu,
                        std::make_tuple("q2_min_mu", "q2_min_tau", "q2_max_mu", "q2_max_tau")),

                make_observable("Lambda_b->Lambda_clnu::A_FB^l",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_a_fb_leptonic,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::A_FB^h",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_a_fb_hadronic,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::R_Lambda_c_A_FB^h(q2)",
                        &LambdaBToLambdaCLeptonNeutrino::differential_ratio_a_fb_hadronic_tau_mu,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::R_Lambda_c_A_FB^h",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_ratio_a_fb_hadronic_tau_mu,
                        std::make_tuple("q2_min_mu", "q2_min_tau", "q2_max_mu", "q2_max_tau")),

                make_observable("Lambda_b->Lambda_clnu::A_FB^c",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_a_fb_combined,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::F_0",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_fzero,
                        std::make_tuple("q2_min", "q2_max")),

                // Lambda_b -> Lambda_c(2595) l nubar
                make_observable("Lambda_b->Lambda_c(2595)lnu::dBR/ds",
                        &LambdaBToLambdaC2595LeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_c(2595)lnu::dBR/dsdtheta_l",
                        &LambdaBToLambdaC2595LeptonNeutrino::double_differential_branching_ratio,
                        std::make_tuple("q2", "theta_l")),

                make_observable("Lambda_b->Lambda_c(2595)lnu::BR", R"(\mathcal{B}(\Lambda_b\to\Lambda_c(2595) \ell^-\bar\nu))",
                        &LambdaBToLambdaC2595LeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_c(2595)lnu::A_FB",
                        &LambdaBToLambdaC2595LeptonNeutrino::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_c(2595)lnu::Gamma_normalized(q2_min,q2_max)",
                        &LambdaBToLambdaC2595LeptonNeutrino::normalized_integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_c(2595)lnu::R_Lambda_c(2595)(q2)",
                        &LambdaBToLambdaC2595LeptonNeutrino::differential_r_lambdac2595,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_c(2595)lnu::R_Lambda_c(2595)",
                        &LambdaBToLambdaC2595LeptonNeutrino::integrated_r_lambdac2595),

                // Lambda_b -> Lambda_c(2625) l nubar
                make_observable("Lambda_b->Lambda_c(2625)lnu::dBR/ds",
                        &LambdaBToLambdaC2625LeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::A_FB(q2)",
                        &LambdaBToLambdaC2625LeptonNeutrino::differential_forward_backward_asymmetry,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::dBR/dsdtheta_l",
                        &LambdaBToLambdaC2625LeptonNeutrino::double_differential_branching_ratio,
                        std::make_tuple("q2", "theta_l")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::BR", R"(\mathcal{B}(\Lambda_b\to\Lambda_c(2625) \ell^-\bar\nu))",
                        &LambdaBToLambdaC2625LeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::A_FB",
                        &LambdaBToLambdaC2625LeptonNeutrino::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::Gamma_normalized(q2_min,q2_max)",
                        &LambdaBToLambdaC2625LeptonNeutrino::normalized_integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::R_Lambda_c(2625)(q2)",
                        &LambdaBToLambdaC2625LeptonNeutrino::differential_r_lambdac2625,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::R_Lambda_c(2625)",
                        &LambdaBToLambdaC2625LeptonNeutrino::integrated_r_lambdac2625),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // Misc.
    // {{{
    ObservableGroup
    make_b_to_xu_semileptonic_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Miscellaneous matrix elements)",
            R"()",
            {
                /* B Meson Properties */
                make_observable("B::M_B^*-M_B",
                        &BMesonProperties::mass_splitting_j1_j0),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    ObservableSection
    make_b_decays_section()
    {
        auto imp = new Implementation<ObservableSection>(
            "Observables in (semi)leptonic $b$-hadron decays",
            "",
            {
                // B^- -> l^- nubar
                make_b_to_l_nu_group(),

                // B_{u,d} -> P l^- nubar
                make_b_to_pi_l_nu_group(),
                make_b_to_d_l_nu_group(),

                // B_{u,d} -> V l^- nubar
                make_b_to_dstar_l_nu_group(),

                // B_s -> V l^- nubar
                make_bs_to_kstar_l_nu_group(),

                // Lambda_b
                make_lambdab_to_lambdac_l_nu_group(),

                // B -> X_u l^- nubar
                make_b_to_xu_semileptonic_group(),
            }
        );

        return ObservableSection(imp);
    }
}
