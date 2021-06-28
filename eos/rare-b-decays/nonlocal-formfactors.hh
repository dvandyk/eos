/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017-2020 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_RARE_B_DECAYS_NONLOCAL_FORMFACTORS_HH
#define EOS_GUARD_EOS_RARE_B_DECAYS_NONLOCAL_FORMFACTORS_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/baryonic.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/rare-b-decays/nonlocal-formfactors-fwd.hh>

#include <memory>
#include <string>

namespace eos
{
    /**
     * Provides the hadronic matrix element of the non-local operator T{ cbar gamma^mu c(x), C_1 O_1 + C_2 O_2 }.
     * We decompose this matrix element as in [BCvDV:2017A], eq. (4).
     */
    template <typename Transition_>
    class NonlocalFormFactor;

    template <typename Transition_>
    using NonlocalFormFactorPtr = std::shared_ptr<NonlocalFormFactor<Transition_>>;

    /**
     * Pseudoobservable in order to expose the nonlocal formfactor, see NonlocalFormFactor.
     */
    template <typename Process_, typename Transition_>
    class NonlocalFormFactorObservable;

    // P -> P

    template <>
    class NonlocalFormFactor<nff::PToP> :
        public ParameterUser
    {
        protected:
            ///@name Stubs to throw InternalError whenever an implementation without residues is called.
            ///@{

            complex<double> jpsi_residues_not_implemented() const;
            complex<double> psi2s_residues_not_implemented() const;
            complex<double> moments_not_implemented() const;

            ///@}

        public:
            ///@name Basic operations
            ///@{

            virtual ~NonlocalFormFactor();

            ///@}

            ///@name Evaluate the formfactor at arbitrary q2 values.
            ///@{

            virtual complex<double> H_plus(const double & q2) const = 0;
            virtual complex<double> Hhat_plus(const double & q2) const = 0;

            ///@}

            ///@name Evaluate the first normalized moment of the formfactor.
            ///@{

            virtual complex<double> normalized_moment_A(const double &) const { return moments_not_implemented(); };

            ///@}

            ///@name Ratio between non local and local formfactors
            ///@{

            virtual complex<double> ratio_plus(const double & q2) const = 0;

            ///@}

            ///@name Evaluate the residue of the formfactor on the J/psi pole.
            ///@{

            virtual complex<double> H_plus_residue_jpsi() const { return jpsi_residues_not_implemented(); };

            ///@}

            ///@name Evaluate the residue of the formfactor on the psi(2S) pole.
            ///@{

            virtual complex<double> H_plus_residue_psi2s() const { return psi2s_residues_not_implemented(); };

            ///@}

            /// Factory method.
            static NonlocalFormFactorPtr<nff::PToP> make(const QualifiedName & name, const Parameters & p, const Options & o);

            ///@name Internal diagnostics for unit tests
            ///@{
            virtual Diagnostics diagnostics() const = 0;
            ///@}
    };

    template <typename Process_>
    class NonlocalFormFactorObservable<Process_, nff::PToP> :
        public ParameterUser,
        public PrivateImplementationPattern<NonlocalFormFactorObservable<Process_, nff::PToP>>
    {
        public:
            ///@name Basic operations
            ///@{

            NonlocalFormFactorObservable(const Parameters &, const Options &);
            ~NonlocalFormFactorObservable();

            ///@}

            ///@name FormFactor as observable
            ///@{

            double re_H_plus(const double & q2) const;
            double im_H_plus(const double & q2) const;
            double abs_H_plus(const double & q2) const;
            double arg_H_plus(const double & q2) const;

            double re_Hhat_plus(const double & q2) const;
            double im_Hhat_plus(const double & q2) const;
            double abs_Hhat_plus(const double & q2) const;

            ///@}

            ///@name First moment of the formfactor as observable
            ///@{

            double re_normalized_moment_A(const double & q2) const;
            double im_normalized_moment_A(const double & q2) const;

            ///@}

            ///@name Ratio between non local and local formfactors
            ///@{

            double re_ratio_plus(const double & q2) const;
            double im_ratio_plus(const double & q2) const;
            double abs_ratio_plus(const double & q2) const;

            ///@}

    };
    extern template class NonlocalFormFactorObservable<nff::BToK, nff::PToP>;

    // P -> V

    template <>
    class NonlocalFormFactor<nff::PToV> :
        public ParameterUser
    {
        protected:
            ///@name Stubs to throw InternalError whenever an implementation without residues is called.
            ///@{

            complex<double> jpsi_residues_not_implemented() const;
            complex<double> psi2s_residues_not_implemented() const;
            complex<double> moments_not_implemented() const;

            ///@}

        public:
            ///@name Basic operations
            ///@{

            virtual ~NonlocalFormFactor();

            ///@}

            ///@name Evaluate the formfactor at arbitrary q2 values.
            ///@{

            virtual complex<double> H_perp(const double & q2) const = 0;
            virtual complex<double> Hhat_perp(const double & q2) const = 0;
            virtual complex<double> H_para(const double & q2) const = 0;
            virtual complex<double> Hhat_para(const double & q2) const = 0;
            virtual complex<double> H_long(const double & q2) const = 0;
            virtual complex<double> Hhat_long(const double & q2) const = 0;

            ///@}

            ///@name Evaluate the first normalized moment of the formfactor.
            ///@{

            virtual complex<double> normalized_moment_V1(const double & /*q2*/) const { return moments_not_implemented(); };
            virtual complex<double> normalized_moment_V2(const double & /*q2*/) const { return moments_not_implemented(); };
            virtual complex<double> normalized_moment_V23(const double & /*q2*/) const { return moments_not_implemented(); };

            ///@}


            ///@name Evaluate the ratio between non-local and local formfactors.
            ///@{

            virtual complex<double> ratio_perp(const double & q2) const = 0;
            virtual complex<double> ratio_para(const double & q2) const = 0;
            virtual complex<double> ratio_long(const double & q2) const = 0;

            ///@}


            ///@name Evaluate the residue of the formfactor on the J/psi pole.
            ///@{

            virtual complex<double> H_perp_residue_jpsi() const { return jpsi_residues_not_implemented(); };
            virtual complex<double> H_para_residue_jpsi() const { return jpsi_residues_not_implemented(); };
            virtual complex<double> H_long_residue_jpsi() const { return jpsi_residues_not_implemented(); };

            ///@}

            ///@name Evaluate the residue of the formfactor on the psi(2S) pole.
            ///@{

            virtual complex<double> H_perp_residue_psi2s() const { return psi2s_residues_not_implemented(); };
            virtual complex<double> H_para_residue_psi2s() const { return psi2s_residues_not_implemented(); };
            virtual complex<double> H_long_residue_psi2s() const { return psi2s_residues_not_implemented(); };

            ///@}

            /// Factory method.
            static NonlocalFormFactorPtr<nff::PToV> make(const QualifiedName & name, const Parameters & p, const Options & o);

            ///@name Internal diagnostics for unit tests
            ///@{
            virtual Diagnostics diagnostics() const = 0;
            ///@}
    };

    template <typename Process_>
    class NonlocalFormFactorObservable<Process_, nff::PToV> :
        public ParameterUser,
        public PrivateImplementationPattern<NonlocalFormFactorObservable<Process_, nff::PToV>>
    {
        public:
            ///@name Basic operations
            ///@{

            NonlocalFormFactorObservable(const Parameters &, const Options &);
            ~NonlocalFormFactorObservable();

            ///@}

            ///@name FormFactor as observable
            ///@{

            double re_H_perp(const double & q2) const;
            double im_H_perp(const double & q2) const;
            double abs_H_perp(const double & q2) const;
            double arg_H_perp(const double & q2) const;
            double re_H_para(const double & q2) const;
            double im_H_para(const double & q2) const;
            double abs_H_para(const double & q2) const;
            double arg_H_para(const double & q2) const;
            double re_H_long(const double & q2) const;
            double im_H_long(const double & q2) const;
            double abs_H_long(const double & q2) const;
            double arg_H_long(const double & q2) const;

            double re_Hhat_perp(const double & q2) const;
            double im_Hhat_perp(const double & q2) const;
            double abs_Hhat_perp(const double & q2) const;
            double re_Hhat_para(const double & q2) const;
            double im_Hhat_para(const double & q2) const;
            double abs_Hhat_para(const double & q2) const;
            double re_Hhat_long(const double & q2) const;
            double im_Hhat_long(const double & q2) const;
            double abs_Hhat_long(const double & q2) const;

            ///@}

            ///@name First moment of the formfactor as observable
            ///@{

            double re_normalized_moment_V1(const double & q2) const;
            double im_normalized_moment_V1(const double & q2) const;
            double re_normalized_moment_V2(const double & q2) const;
            double im_normalized_moment_V2(const double & q2) const;
            double re_normalized_moment_V23(const double & q2) const;
            double im_normalized_moment_V23(const double & q2) const;

            ///@}

            ///@name Ratio between non local and local formfactors
            ///@{

            double re_ratio_perp(const double & q2) const;
            double im_ratio_perp(const double & q2) const;
            double abs_ratio_perp(const double & q2) const;
            double re_ratio_para(const double & q2) const;
            double im_ratio_para(const double & q2) const;
            double abs_ratio_para(const double & q2) const;
            double re_ratio_long(const double & q2) const;
            double im_ratio_long(const double & q2) const;
            double abs_ratio_long(const double & q2) const;
            ///@}
    };
    extern template class NonlocalFormFactorObservable<nff::BToKstar, nff::PToV>;
    extern template class NonlocalFormFactorObservable<nff::BsToPhi, nff::PToV>;

    //  OneHalfPlus -> OneHalfPlus

    template <>
    class NonlocalFormFactor<nff::OneHalfPlusToOneHalfPlus> :
        public ParameterUser
    {
        public:

            virtual ~NonlocalFormFactor();

            virtual complex<double> H_V_perp(const double & q2) const = 0;
            virtual complex<double> H_V_long(const double & q2) const = 0;

            virtual complex<double> H_V_perp_residue_jpsi() const = 0;
            virtual complex<double> H_V_long_residue_jpsi() const = 0;

            virtual complex<double> H_V_perp_residue_psi2s() const = 0;
            virtual complex<double> H_V_long_residue_psi2s() const = 0;

            virtual complex<double> H_A_perp(const double & q2) const = 0;
            virtual complex<double> H_A_long(const double & q2) const = 0;

            virtual complex<double> H_A_perp_residue_jpsi() const = 0;
            virtual complex<double> H_A_long_residue_jpsi() const = 0;

            virtual complex<double> H_A_perp_residue_psi2s() const = 0;
            virtual complex<double> H_A_long_residue_psi2s() const = 0;

            virtual complex<double> ratio_H_V_long(const double & q2) const = 0;
            virtual complex<double> ratio_H_V_perp(const double & q2) const = 0;
            virtual complex<double> ratio_H_A_long(const double & q2) const = 0;
            virtual complex<double> ratio_H_A_perp(const double & q2) const = 0;


            /// Factory method.
            static NonlocalFormFactorPtr<nff::OneHalfPlusToOneHalfPlus> make(const QualifiedName & name, const Parameters & p, const Options & o);

            ///Internal diagnostics for unit tests
            virtual Diagnostics diagnostics() const = 0;

    };

    template <typename Process_>
    class NonlocalFormFactorObservable<Process_, nff::OneHalfPlusToOneHalfPlus> :
        public ParameterUser,
        public PrivateImplementationPattern<NonlocalFormFactorObservable<Process_, nff::OneHalfPlusToOneHalfPlus>>
    {
        public:
            NonlocalFormFactorObservable(const Parameters &, const Options &);
            ~NonlocalFormFactorObservable();

            double re_H_V_perp(const double & q2) const;
            double im_H_V_perp(const double & q2) const;

            double re_H_V_long(const double & q2) const;
            double im_H_V_long(const double & q2) const;

            double re_H_A_perp(const double & q2) const;
            double im_H_A_perp(const double & q2) const;

            double re_H_A_long(const double & q2) const;
            double im_H_A_long(const double & q2) const;


            double re_ratio_H_V_perp(const double & q2) const;
            double im_ratio_H_V_perp(const double & q2) const;
            double re_ratio_H_V_long(const double & q2) const;
            double im_ratio_H_V_long(const double & q2) const;

            double re_ratio_H_A_perp(const double & q2) const;
            double im_ratio_H_A_perp(const double & q2) const;
            double re_ratio_H_A_long(const double & q2) const;
            double im_ratio_H_A_long(const double & q2) const;
    };
    extern template class NonlocalFormFactorObservable<nff::LambdabToLambda, nff::OneHalfPlusToOneHalfPlus>;


    namespace nff_utils
    {

        complex<double> z(const double & q2, complex<double> s_plus, complex<double> s_0);
        complex<double> blaschke_cc(const complex<double> z, const complex<double> z_Jpsi, const complex<double> z_psi2S);
        complex<double> P(complex<double> z,
            const complex<double> & alpha_0, const complex<double> & alpha_1, const complex<double> & alpha_2);
        complex<double> PGvDV2020(complex<double> z, const complex<double> zXY,
            const complex<double> & alpha_0, const complex<double> & alpha_1, const complex<double> & alpha_2);

    }

}
#endif
