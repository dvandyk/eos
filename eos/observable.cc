/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Danny van Dyk
 * Copyright (c) 2011 Christian Wacker
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
#include <eos/form-factors/form-factor-adapter.hh>
#include <eos/form-factors/analytic-b-to-pi-pi.hh>
#include <eos/form-factors/baryonic-impl.hh>
#include <eos/form-factors/mesonic-impl.hh>
#include <eos/form-factors/zero-recoil-sum-rule.hh>
#include <eos/b-decays/b-to-l-nu.hh>
#include <eos/b-decays/b-to-pi-l-nu.hh>
#include <eos/b-decays/b-to-pi-pi-l-nu.hh>
#include <eos/b-decays/b-to-d-l-nu.hh>
#include <eos/b-decays/bs-to-kstar-l-nu.hh>
#include <eos/b-decays/inclusive-b-to-u.hh>
#include <eos/b-decays/properties.hh>
#include <eos/rare-b-decays/exclusive-b-to-dilepton.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton-large-recoil.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton-low-recoil.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-gamma.hh>
#include <eos/rare-b-decays/inclusive-b-to-s-dilepton.hh>
#include <eos/rare-b-decays/inclusive-b-to-s-gamma.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda-dilepton.hh>
#include <eos/utils/concrete_observable.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/observable_stub.hh>

#include <algorithm>
#include <map>

namespace eos
{
    ObservableNameError::ObservableNameError(const std::string & name) :
        Exception("Observable name '" + name + "' is malformed")
    {
    }

    ObservableFactory::ObservableFactory()
    {
    }

    ObservableFactory::~ObservableFactory()
    {
    }

    /* regular observables */
    template <typename Decay_, typename ... Args_>
    std::pair<std::string, ObservableFactory *> make_observable(const char * name,
            double (Decay_::* function)(const Args_ & ...) const)
    {
        std::string sname(name);

        return std::make_pair(sname, make_concrete_observable_factory(sname, function, std::make_tuple()));
    }

    template <typename Decay_, typename Tuple_, typename ... Args_>
    std::pair<std::string, ObservableFactory *> make_observable(const char * name,
            double (Decay_::* function)(const Args_ & ...) const,
            const Tuple_ & kinematics_names)
    {
        std::string sname(name);

        return std::make_pair(sname, make_concrete_observable_factory(sname, function, kinematics_names));
    }

    /* form factors as observables */
    template <typename Transition_>
    std::pair<std::string, ObservableFactory *> make_observable(const char * name,
            const char * process,
            double (FormFactors<Transition_>::* function)(const double &) const)
    {
        std::string sname(name), sprocess(process);

        return std::make_pair(sname, new FormFactorAdapterFactory<Transition_>(sname, sprocess, function));
    }

    template <typename Transition_>
    std::pair<std::string, ObservableFactory *> make_observable(const char * name,
            const char * process,
            double (FormFactors<Transition_>::* numerator)(const double &) const,
            double (FormFactors<Transition_>::* denominator)(const double &) const)
    {
        std::string sname(name), sprocess(process);

        return std::make_pair(sname, new FormFactorRatioAdapterFactory<Transition_>(sname, sprocess, numerator, denominator));
    }

    ObservablePtr
    Observable::make(const std::string & _name, const Parameters & parameters, const Kinematics & kinematics, const Options & _options)
    {
        static const std::map<std::string, ObservableFactory *> simple_observables
        {
            /* B Meson Properties */
            make_observable("B::M_B^*-M_B",
                    &BMesonProperties::mass_splitting_j1_j0),

            /* Form Factor for the Exclusive Decays */

            // B -> pi Form Factors
            make_observable("B->pi::f_+(s)", "B->pi",
                    &FormFactors<PToP>::f_p),

            make_observable("B->pi::f_+'(s)", "B->pi",
                    &FormFactors<PToP>::f_p_d1),

            make_observable("B->pi::f_+''(s)", "B->pi",
                    &FormFactors<PToP>::f_p_d2),

            make_observable("B->pi::f_T(s)", "B->pi",
                    &FormFactors<PToP>::f_t),

            make_observable("B->pi::f_0(s)", "B->pi",
                    &FormFactors<PToP>::f_0),

            // B -> pi Form Factors (auxiliary variables, e.g. for determining the
            // LCSR/SVZ threshold parameters)
            make_observable("B->pi::M_B(LCSR)@DKMMO2008",
                    &AnalyticFormFactorBToPiDKMMO2008::MB_lcsr,
                    std::make_tuple("s")),

            make_observable("B->pi::M_B(SVZ)@DKMMO2008",
                    &AnalyticFormFactorBToPiDKMMO2008::MB_svz),

            make_observable("B->pi::f_B@DKMMO2008",
                    &AnalyticFormFactorBToPiDKMMO2008::decay_constant),

            // B -> K Form Factors
            make_observable("B->K::f_+(s)", "B->K",
                    &FormFactors<PToP>::f_p),

            make_observable("B->K::f_T(s)", "B->K",
                    &FormFactors<PToP>::f_t),

            make_observable("B->K::f_0(s)", "B->K",
                    &FormFactors<PToP>::f_0),

            // B -> K^* Form Factors
            make_observable("B->K^*::V(s)", "B->K^*",
                    &FormFactors<PToV>::v),

            make_observable("B->K^*::A_0(s)", "B->K^*",
                    &FormFactors<PToV>::a_0),

            make_observable("B->K^*::A_1(s)", "B->K^*",
                    &FormFactors<PToV>::a_1),

            make_observable("B->K^*::A_2(s)", "B->K^*",
                    &FormFactors<PToV>::a_2),

            make_observable("B->K^*::A_12(s)", "B->K^*",
                    &FormFactors<PToV>::a_12),

            make_observable("B->K^*::T_1(s)", "B->K^*",
                    &FormFactors<PToV>::t_1),

            make_observable("B->K^*::T_2(s)", "B->K^*",
                    &FormFactors<PToV>::t_2),

            make_observable("B->K^*::T_3(s)", "B->K^*",
                    &FormFactors<PToV>::t_3),

            make_observable("B->K^*::T_23(s)", "B->K^*",
                    &FormFactors<PToV>::t_23),

            make_observable("B->K^*::V(s)/A_1(s)", "B->K^*",
                    &FormFactors<PToV>::v, &FormFactors<PToV>::a_1),

            make_observable("B->K^*::A_2(s)/A_1(s)", "B->K^*",
                    &FormFactors<PToV>::a_2, &FormFactors<PToV>::a_1),

            make_observable("B->K^*::A_12(s)/A_1(s)", "B->K^*",
                    &FormFactors<PToV>::a_12, &FormFactors<PToV>::a_1),

            make_observable("B->K^*::T_23(s)/T_2(s)", "B->K^*",
                    &FormFactors<PToV>::t_23, &FormFactors<PToV>::t_2),

            // B -> D Form Factors
            make_observable("B->D::f_+(s)", "B->D",
                    &FormFactors<PToP>::f_p),

            make_observable("B->D::f_0(s)", "B->D",
                    &FormFactors<PToP>::f_0),

            // B_s -> K^* Form Factors
            make_observable("B_s->K^*::V(s)", "B_s->K^*",
                    &FormFactors<PToV>::v),

            make_observable("B_s->K^*::A_0(s)", "B_s->K^*",
                    &FormFactors<PToV>::a_0),

            make_observable("B_s->K^*::A_1(s)", "B_s->K^*",
                    &FormFactors<PToV>::a_1),

            make_observable("B_s->K^*::A_2(s)", "B_s->K^*",
                    &FormFactors<PToV>::a_2),

            make_observable("B_s->K^*::A_12(s)", "B_s->K^*",
                    &FormFactors<PToV>::a_12),

            make_observable("B->pipi::Im{F_perp}",
                    &AnalyticFormFactorBToPiPiBFvD2016::im_f_perp,
                    std::make_tuple("q2", "k2", "z")),

            make_observable("B->pipi::Im{F_para}",
                    &AnalyticFormFactorBToPiPiBFvD2016::im_f_para,
                    std::make_tuple("q2", "k2", "z")),

            make_observable("B->pipi::Im{F_long}",
                    &AnalyticFormFactorBToPiPiBFvD2016::im_f_long,
                    std::make_tuple("q2", "k2", "z")),

            make_observable("B->pipi::Im{F_time}",
                    &AnalyticFormFactorBToPiPiBFvD2016::im_f_time,
                    std::make_tuple("q2", "k2", "z")),

            // Lambda_b -> Lambda Form Factors
            make_observable("Lambda_b->Lambda::f_time^V(s)", "Lambda_b->Lambda",
                    &FormFactors<OneHalfPlusToOneHalfPlus>::f_time_v),

            make_observable("Lambda_b->Lambda::f_long^V(s)", "Lambda_b->Lambda",
                    &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_v),

            make_observable("Lambda_b->Lambda::f_perp^V(s)", "Lambda_b->Lambda",
                    &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_v),

            make_observable("Lambda_b->Lambda::f_time^A(s)", "Lambda_b->Lambda",
                    &FormFactors<OneHalfPlusToOneHalfPlus>::f_time_a),

            make_observable("Lambda_b->Lambda::f_long^A(s)", "Lambda_b->Lambda",
                    &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_a),

            make_observable("Lambda_b->Lambda::f_perp^A(s)", "Lambda_b->Lambda",
                    &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_a),

            make_observable("Lambda_b->Lambda::f_long^T(s)", "Lambda_b->Lambda",
                    &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_t),

            make_observable("Lambda_b->Lambda::f_perp^T(s)", "Lambda_b->Lambda",
                    &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_t),

            make_observable("Lambda_b->Lambda::f_long^T5(s)", "Lambda_b->Lambda",
                    &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_t5),

            make_observable("Lambda_b->Lambda::f_perp^T5(s)", "Lambda_b->Lambda",
                    &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_t5),

            // Zero-Recoil Sum Rule for the Lambda_b -> Lambda_c Form Factors
            make_observable("Lambda_b->Lambda_c::F(1)",
                    &ZeroRecoilSumRule<LambdaBToC>::vector_current),

            make_observable("Lambda_b->Lambda_c::G(1)",
                    &ZeroRecoilSumRule<LambdaBToC>::axialvector_current),

            /* Exclusive Decays */

            /* Exclusive B Decays */

            // B_q -> l nubar
            make_observable("B_u->lnu::BR",
                    &BToLeptonNeutrino::branching_ratio),

            // B -> pi l nu
            make_observable("B->pilnu::dBR/ds",
                    &BToPiLeptonNeutrino::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B->pilnu::BR",
                    &BToPiLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->pilnu::zeta",
                    &BToPiLeptonNeutrino::integrated_zeta,
                    std::make_tuple("s_min", "s_max")),

            // B -> pi pi l nu
            make_observable("B->pipilnu::BR(2)",
                    &BToPiPiLeptonNeutrino::double_differential_branching_ratio,
                    std::make_tuple("q2", "k2")),

            make_observable("B->pipilnu::BR(3)",
                    &BToPiPiLeptonNeutrino::triple_differential_branching_ratio,
                    std::make_tuple("q2", "k2", "cos(theta_pi)")),

            make_observable("B->pipilnu::A_FB(2)",
                    &BToPiPiLeptonNeutrino::double_differential_forward_backward_asymmetry,
                    std::make_tuple("q2", "k2")),

            make_observable("B->pipilnu::P(z)",
                    &BToPiPiLeptonNeutrino::partial_waves,
                    std::make_tuple("q2", "k2", "cos(theta_pi)")),

            make_observable("B->pipilnu::BR",
                    &BToPiPiLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple("q2_min", "q2_max", "k2_min", "k2_max", "z_min", "z_max")),

            make_observable("B->pipilnu::A_FB",
                    &BToPiPiLeptonNeutrino::integrated_forward_backward_asymmetry,
                    std::make_tuple("q2_min", "q2_max", "k2_min", "k2_max")),

            // B -> D l nu
            make_observable("B->Dlnu::dBR/ds",
                            &BToDLeptonNeutrino::differential_branching_ratio,
                            std::make_tuple("s")),

            make_observable("B->Dlnu::BR",
                            &BToDLeptonNeutrino::integrated_branching_ratio,
                            std::make_tuple("s_min", "s_max")),

            make_observable("B->Dlnu::R_D(s)",
                            &BToDLeptonNeutrino::differential_r_d,
                            std::make_tuple("s")),

            make_observable("B->Dlnu::R_D",
                            &BToDLeptonNeutrino::integrated_r_d),

            // B_s -> K^* l nubar
            make_observable("B_s->K^*lnu::F_perp(s)",
                    &BsToKstarLeptonNeutrino::Fperp,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::F_para(s)",
                    &BsToKstarLeptonNeutrino::Fpara,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::F_long(s)",
                    &BsToKstarLeptonNeutrino::Flong,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::d^4Gamma",
                    &BsToKstarLeptonNeutrino::four_differential_decay_width,
                    std::make_tuple("s", "cos(theta_l)", "cos(theta_k)", "phi")),

            make_observable("B_s->K^*lnu::dBR/ds",
                    &BsToKstarLeptonNeutrino::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::A_FB(s)",
                    &BsToKstarLeptonNeutrino::differential_forward_backward_asymmetry,
                    std::make_tuple("s")),

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

            make_observable("B_s->K^*lnu::A_FB(s)",
                    &BsToKstarLeptonNeutrino::differential_forward_backward_asymmetry,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::A_T^2(s)",
                    &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_2,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::A_T^3(s)",
                    &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_3,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::A_T^4(s)",
                    &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_4,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::A_T^5(s)",
                    &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_5,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::A_T^re(s)",
                    &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_re,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::A_T^im(s)",
                    &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_im,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::F_L(s)",
                    &BsToKstarLeptonNeutrino::differential_longitudinal_polarisation,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::F_T(s)",
                    &BsToKstarLeptonNeutrino::differential_transversal_polarisation,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::H_T^1(s)",
                    &BsToKstarLeptonNeutrino::differential_h_1,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::H_T^2(s)",
                    &BsToKstarLeptonNeutrino::differential_h_2,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::H_T^3(s)",
                    &BsToKstarLeptonNeutrino::differential_h_3,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::H_T^4(s)",
                    &BsToKstarLeptonNeutrino::differential_h_4,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::H_T^5(s)",
                    &BsToKstarLeptonNeutrino::differential_h_5,
                    std::make_tuple("s")),

            make_observable("B_s->K^*lnu::A_FB",
                    &BsToKstarLeptonNeutrino::integrated_forward_backward_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B_s->K^*lnu::BR",
                    &BsToKstarLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B_s->K^*lnu::F_L",
                    &BsToKstarLeptonNeutrino::integrated_longitudinal_polarisation,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B_s->K^*lnu::F_T",
                    &BsToKstarLeptonNeutrino::integrated_transversal_polarisation,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B_s->K^*lnu::A_T^2",
                    &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_2,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B_s->K^*lnu::A_T^3",
                    &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B_s->K^*lnu::A_T^4",
                    &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B_s->K^*lnu::A_T^5",
                    &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_5,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B_s->K^*lnu::A_T^re",
                    &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_re,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B_s->K^*lnu::A_T^im",
                    &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_im,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B_s->K^*lnu::H_T^1",
                    &BsToKstarLeptonNeutrino::integrated_h_1,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B_s->K^*lnu::H_T^2",
                    &BsToKstarLeptonNeutrino::integrated_h_2,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B_s->K^*lnu::H_T^3",
                    &BsToKstarLeptonNeutrino::integrated_h_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B_s->K^*lnu::H_T^4",
                    &BsToKstarLeptonNeutrino::integrated_h_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B_s->K^*lnu::H_T^5",
                    &BsToKstarLeptonNeutrino::integrated_h_5,
                    std::make_tuple("s_min", "s_max")),

            // B_s -> K^* l nubar Ratios
            make_observable("B_s->K^*lnu::R_long",
                    &BsToKstarLeptonNeutrinoRatios::ratio_long),

            make_observable("B_s->K^*lnu::R_para",
                    &BsToKstarLeptonNeutrinoRatios::ratio_para),

            make_observable("B_s->K^*lnu::R_perp",
                    &BsToKstarLeptonNeutrinoRatios::ratio_perp),

            /* Exclusive Rare B Decays */

            // B_q -> ll
            make_observable("B_q->ll::BR",
                    &BToDilepton::branching_ratio_time_zero),

            make_observable("B_q->ll::BR@Untagged",
                    &BToDilepton::branching_ratio_untagged_integrated),

            make_observable("B_q->ll::A_DeltaGamma",
                    &BToDilepton::cp_asymmetry_del_gamma),

            make_observable("B_q->ll::S",
                    &BToDilepton::cp_asymmetry_mixing_S),

            make_observable("B_q->ll::eff_lifetime",
                    &BToDilepton::effective_lifetime),

            // B -> K^* gamma
            make_observable("B->K^*gamma::BR",
                    &BToKstarGamma::branching_ratio),

            make_observable("B->K^*gamma::BRavg",
                    &BToKstarGamma::branching_ratio_cp_averaged),

            make_observable("B->K^*gamma::A_CP",
                    &BToKstarGamma::cp_asymmetry),

            make_observable("B->K^*gamma::S_K^*gamma",
                    &BToKstarGamma::s_kstar_gamma),

            make_observable("B->K^*gamma::C_K^*gamma",
                    &BToKstarGamma::c_kstar_gamma),

            make_observable("B->K^*gamma::A_I",
                    &BToKstarGamma::isospin_asymmetry),

            // B -> K ll, Large Recoil
            make_observable("B->Kll::d^2Gamma@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::two_differential_decay_width,
                    std::make_tuple("s", "cos(theta_l)")),

            make_observable("B->Kll::dBR/ds@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B->Kll::F_H(s)@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::differential_flat_term,
                    std::make_tuple("s")),

            make_observable("B->Kll::A_FB(s)@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::differential_forward_backward_asymmetry,
                    std::make_tuple("s")),

            make_observable("B->Kll::R_K(s)@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::differential_ratio_muons_electrons,
                    std::make_tuple("s")),

            make_observable("B->Kll::BR@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::BRavg@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_branching_ratio_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::A_CP@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_cp_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::Gamma@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_decay_width,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::F_H@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_flat_term,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::F_Havg@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_flat_term_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::A_FB@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_forward_backward_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::A_FBavg@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_forward_backward_asymmetry_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::R_K@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_ratio_muons_electrons,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::a_l@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::a_l,
                    std::make_tuple("s")),

            make_observable("B->Kll::b_l@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::b_l,
                    std::make_tuple("s")),

            make_observable("B->Kll::c_l@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::c_l,
                    std::make_tuple("s")),

            // B -> K ll, Low Recoil
            make_observable("B->Kll::d^2Gamma@LowRecoil",
                    &BToKDilepton<LowRecoil>::two_differential_decay_width,
                    std::make_tuple("s", "cos(theta_l)")),

            make_observable("B->Kll::dBR/ds@LowRecoil",
                    &BToKDilepton<LowRecoil>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B->Kll::F_H(s)@LowRecoil",
                    &BToKDilepton<LowRecoil>::differential_flat_term,
                    std::make_tuple("s")),

            make_observable("B->Kll::A_FB(s)@LowRecoil",
                    &BToKDilepton<LowRecoil>::differential_forward_backward_asymmetry,
                    std::make_tuple("s")),

            make_observable("B->Kll::R_K(s)@LowRecoil",
                    &BToKDilepton<LowRecoil>::differential_ratio_muons_electrons,
                    std::make_tuple("s")),

            make_observable("B->Kll::BR@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::BRavg@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_branching_ratio_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::A_CP@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_cp_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::Gamma@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_decay_width,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::F_H@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_flat_term,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::F_Havg@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_flat_term_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::A_FB@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_forward_backward_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::A_FBavg@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_forward_backward_asymmetry_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::R_K@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_ratio_muons_electrons,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::a_l@LowRecoil",
                    &BToKDilepton<LowRecoil>::a_l,
                    std::make_tuple("s")),

            make_observable("B->Kll::b_l@LowRecoil",
                    &BToKDilepton<LowRecoil>::b_l,
                    std::make_tuple("s")),

            make_observable("B->Kll::c_l@LowRecoil",
                    &BToKDilepton<LowRecoil>::c_l,
                    std::make_tuple("s")),

            make_observable("B->Kll::Re{c9eff}@LowRecoil",
                    &BToKDilepton<LowRecoil>::real_c9eff,
                    std::make_tuple("s")),

            make_observable("B->Kll::Im{c9eff}@LowRecoil",
                    &BToKDilepton<LowRecoil>::imag_c9eff,
                    std::make_tuple("s")),

            make_observable("B->Kll::Re{c7eff}@LowRecoil",
                    &BToKDilepton<LowRecoil>::real_c7eff,
                    std::make_tuple("s")),

            make_observable("B->Kll::Im{c7eff}@LowRecoil",
                    &BToKDilepton<LowRecoil>::imag_c7eff,
                    std::make_tuple("s")),

            // B -> K^* ll
            make_observable("B->K^*ll::s_0(A_FB)",
                    &BToKstarDilepton::a_fb_zero_crossing),

            make_observable("B->K^*ll::d^4Gamma",
                    &BToKstarDilepton::decay_width,
                    std::make_tuple("s", "cos(theta_l)", "cos(theta_k)", "phi")),

            make_observable("B->K^*ll::dBR/ds",
                    &BToKstarDilepton::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_FB(s)",
                    &BToKstarDilepton::differential_forward_backward_asymmetry,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::F_L(s)",
                    &BToKstarDilepton::differential_longitudinal_polarisation,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::F_T(s)",
                    &BToKstarDilepton::differential_transversal_polarisation,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^2(s)",
                    &BToKstarDilepton::differential_transverse_asymmetry_2,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^3(s)",
                    &BToKstarDilepton::differential_transverse_asymmetry_3,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^4(s)",
                    &BToKstarDilepton::differential_transverse_asymmetry_4,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^5(s)",
                    &BToKstarDilepton::differential_transverse_asymmetry_5,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^re(s)",
                    &BToKstarDilepton::differential_transverse_asymmetry_re,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^im(s)",
                    &BToKstarDilepton::differential_transverse_asymmetry_im,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::P'_4(s)",
                    &BToKstarDilepton::differential_p_prime_4,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::P'_5(s)",
                    &BToKstarDilepton::differential_p_prime_5,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::P'_6(s)",
                    &BToKstarDilepton::differential_p_prime_6,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_1s(s)",
                    &BToKstarDilepton::differential_j_1s,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_1c(s)",
                    &BToKstarDilepton::differential_j_1c,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_2s(s)",
                    &BToKstarDilepton::differential_j_2s,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_2c(s)",
                    &BToKstarDilepton::differential_j_2c,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_3(s)",
                    &BToKstarDilepton::differential_j_3,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_4(s)",
                    &BToKstarDilepton::differential_j_4,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_5(s)",
                    &BToKstarDilepton::differential_j_5,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_6s(s)",
                    &BToKstarDilepton::differential_j_6s,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_6c(s)",
                    &BToKstarDilepton::differential_j_6c,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_7(s)",
                    &BToKstarDilepton::differential_j_7,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_8(s)",
                    &BToKstarDilepton::differential_j_8,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_9(s)",
                    &BToKstarDilepton::differential_j_9,
                    std::make_tuple("s")),

#if 0
            make_observable("B->K^*ll::D_4(s)",
                    &BToKstarDilepton::differential_d_4,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::D_5(s)",
                    &BToKstarDilepton::differential_d_5,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::D_6s(s)",
                    &BToKstarDilepton::differential_d_6s,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::R_K^*(s)",
                    &BToKstarDilepton::differential_ratio_muons_electrons,
                    std::make_tuple("s")),
#endif
            make_observable("B->K^*ll::Gamma",
                    &BToKstarDilepton::integrated_decay_width,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::BR",
                    &BToKstarDilepton::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_FB",
                    &BToKstarDilepton::integrated_forward_backward_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::F_L",
                    &BToKstarDilepton::integrated_longitudinal_polarisation,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::F_T",
                    &BToKstarDilepton::integrated_transversal_polarisation,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^2",
                    &BToKstarDilepton::integrated_transverse_asymmetry_2,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^3",
                    &BToKstarDilepton::integrated_transverse_asymmetry_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^4",
                    &BToKstarDilepton::integrated_transverse_asymmetry_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^5",
                    &BToKstarDilepton::integrated_transverse_asymmetry_5,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^re",
                    &BToKstarDilepton::integrated_transverse_asymmetry_re,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^im",
                    &BToKstarDilepton::integrated_transverse_asymmetry_im,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::P'_4",
                    &BToKstarDilepton::integrated_p_prime_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::P'_5",
                    &BToKstarDilepton::integrated_p_prime_5,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::P'_6",
                    &BToKstarDilepton::integrated_p_prime_6,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^1(s)",
                    &BToKstarDilepton::differential_h_1,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::H_T^2(s)",
                    &BToKstarDilepton::differential_h_2,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::H_T^3(s)",
                    &BToKstarDilepton::differential_h_3,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::H_T^4(s)",
                    &BToKstarDilepton::differential_h_4,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::H_T^5(s)",
                    &BToKstarDilepton::differential_h_5,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::H_T^1",
                    &BToKstarDilepton::integrated_h_1,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^2",
                    &BToKstarDilepton::integrated_h_2,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^3",
                    &BToKstarDilepton::integrated_h_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^4",
                    &BToKstarDilepton::integrated_h_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^5",
                    &BToKstarDilepton::integrated_h_5,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_1s",
                    &BToKstarDilepton::integrated_j_1s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_1c",
                    &BToKstarDilepton::integrated_j_1c,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_2s",
                    &BToKstarDilepton::integrated_j_2s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_2c",
                    &BToKstarDilepton::integrated_j_2c,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_3",
                    &BToKstarDilepton::integrated_j_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_4",
                    &BToKstarDilepton::integrated_j_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_5",
                    &BToKstarDilepton::integrated_j_5,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_6s",
                    &BToKstarDilepton::integrated_j_6s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_6c",
                    &BToKstarDilepton::integrated_j_6c,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_7",
                    &BToKstarDilepton::integrated_j_7,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_8",
                    &BToKstarDilepton::integrated_j_8,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_9",
                    &BToKstarDilepton::integrated_j_9,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_CP",
                    &BToKstarDilepton::integrated_cp_asymmetry,
                    std::make_tuple("s_min", "s_max")),
#if 0
            make_observable("B->K^*ll::A_I(s)",
                    &BToKstarDilepton::differential_isospin_asymmetry,
                    std::make_tuple("s")),
#endif
            make_observable("B->K^*ll::S_1s",
                    &BToKstarDilepton::integrated_s_1s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::S_1c",
                    &BToKstarDilepton::integrated_s_1c,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::S_2s",
                    &BToKstarDilepton::integrated_s_2s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::S_2c",
                    &BToKstarDilepton::integrated_s_2c,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::S_3",
                    &BToKstarDilepton::integrated_s_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::S_4",
                    &BToKstarDilepton::integrated_s_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::S_5",
                    &BToKstarDilepton::integrated_s_5,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::S_6s",
                    &BToKstarDilepton::integrated_s_6s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::S_7",
                    &BToKstarDilepton::integrated_s_7,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::S_8",
                    &BToKstarDilepton::integrated_s_8,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::S_9",
                    &BToKstarDilepton::integrated_s_9,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_1s",
                    &BToKstarDilepton::integrated_a_1s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_1c",
                    &BToKstarDilepton::integrated_a_1c,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_2s",
                    &BToKstarDilepton::integrated_a_2s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_2c",
                    &BToKstarDilepton::integrated_a_2c,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_3",
                    &BToKstarDilepton::integrated_a_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_4",
                    &BToKstarDilepton::integrated_a_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_5",
                    &BToKstarDilepton::integrated_a_5,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_6s",
                    &BToKstarDilepton::integrated_a_6s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_7",
                    &BToKstarDilepton::integrated_a_7,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_8",
                    &BToKstarDilepton::integrated_a_8,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_9",
                    &BToKstarDilepton::integrated_a_9,
                    std::make_tuple("s_min", "s_max")),

#if 0
            make_observable("B->K^*ll::D_4",
                    &BToKstarDilepton::integrated_d_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::D_5",
                    &BToKstarDilepton::integrated_d_5,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::D_6s",
                    &BToKstarDilepton::integrated_d_6s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::R_K^*",
                    &BToKstarDilepton::integrated_ratio_muons_electrons,
                    std::make_tuple("s_min", "s_max")),
#endif
            // Lambda_b -> Lambda l^+ l^-, Large Recoil
            make_observable("Lambda_b->Lambdall::dBR/ds@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::A_FB^l(s)@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_leptonic,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::A_FB^h(s)@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_hadronic,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::A_FB^c(s)@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_combined,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::F_0(s)@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::differential_fzero,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::BR@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::A_FB^l@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_leptonic,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::A_FB^h@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_hadronic,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::A_FB^c@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_combined,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::F_0@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::integrated_fzero,
                    std::make_tuple("s_min", "s_max")),

            // Lambda_b -> Lambda l^+ l^-, Low Recoil
            make_observable("Lambda_b->Lambdall::dBR/ds@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::A_FB^l(s)@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_leptonic,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::A_FB^h(s)@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_hadronic,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::A_FB^c(s)@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_combined,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::F_0(s)@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::differential_fzero,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::BR@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::A_FB^l@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_leptonic,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::A_FB^h@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_hadronic,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::A_FB^c@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_combined,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::F_0@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_fzero,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::K_1ss@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_k1ss,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::K_1cc@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_k1cc,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::K_1c@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_k1c,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::K_2ss@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_k2ss,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::K_2cc@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_k2cc,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::K_2c@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_k2c,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::K_3sc@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_k3sc,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::K_3s@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_k3s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::K_4sc@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_k4sc,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::K_4s@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_k4s,
                    std::make_tuple("s_min", "s_max")),

            /* Inclusive Decays */

            // B->X_u l nu (naive)
            make_observable("B->X_ulnu::|V_ub|@Naive",
                    &BToXuLeptonNeutrino<Naive>::v_ub),

            // B->X_s ll, HLMW2005
            make_observable("B->X_sll::dBR/ds@HLMW2005",
                    &BToXsDilepton<HLMW2005>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B->X_sll::BR@HLMW2005",
                    &BToXsDilepton<HLMW2005>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            // B->X_s gamma
            make_observable("B->X_sgamma::BR@Minimal",
                    &BToXsGamma<Minimal>::integrated_branching_ratio),

            // B->X_s gamma, NLO implementation
            make_observable("B->X_sgamma::BR(E_min)@NLO",
                    &BToXsGamma<NLO>::integrated_branching_ratio,
                    std::make_tuple("E_min")),

            make_observable("B->X_sgamma::E_1(E_min)@NLO",
                    &BToXsGamma<NLO>::photon_energy_moment_1,
                    std::make_tuple("E_min")),

            make_observable("B->X_sgamma::E_2(E_min)@NLO",
                    &BToXsGamma<NLO>::photon_energy_moment_2,
                    std::make_tuple("E_min")),
        };

        Options options;
        std::string name(_name);

        std::string::size_type pos;
        while (std::string::npos != (pos = name.rfind(',')))
        {
            std::string::size_type sep(name.find('=', pos + 1));
            if (std::string::npos == sep)
                throw ObservableNameError(_name);

            std::string key(name.substr(pos + 1, sep - pos - 1));
            std::string value(name.substr(sep + 1));

            options.set(key, value);
            name.erase(pos);
        }

        // check if 'name' matches a simple observable
        {
            auto i = simple_observables.find(name);
            if (simple_observables.end() != i)
                return i->second->make(parameters, kinematics, options + _options);
        }

        // check if 'name' matches a parameter
        {
            auto i = std::find_if(parameters.begin(), parameters.end(), [&] (const Parameter & p) { return p.name() == name; });
            if (options.empty() && (parameters.end() != i))
                return ObservablePtr(new ObservableStub(parameters, name));
        }

        return ObservablePtr();
    }
}
