/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013, 2014, 2015, 2016 Danny van Dyk
 * Copyright (c) 2015 Christoph Bobeth
 * Copyright (c) 2010 Christian Wacker
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

#ifndef EOS_GUARD_SRC_FORM_FACTORS_MESONIC_HH
#define EOS_GUARD_SRC_FORM_FACTORS_MESONIC_HH 1

#include <eos/form-factors/form-factors-fwd.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>
#include <eos/utils/qualified-name.hh>

#include <memory>
#include <string>

namespace eos
{
    /* Mesonic Tags */

    /*
     * P -> V transitions
     *
     * P: (Heavy) pseudoscalar meson.
     * V: Light vector meson.
     */
    struct PToV { };

    struct PToP { };

    struct PToPP { };

    struct VToP { };

    struct VToV { };

    struct VacuumToPGamma { };

    template <>
    class FormFactors<PToV> :
        public virtual ParameterUser
    {
        public:
            virtual ~FormFactors();

            virtual double v(const double & s) const = 0;

            virtual double a_0(const double & s) const = 0;
            virtual double a_1(const double & s) const = 0;
            virtual double a_2(const double & s) const = 0;
            virtual double a_12(const double & s) const = 0;

            virtual double t_1(const double & s) const = 0;
            virtual double t_2(const double & s) const = 0;
            virtual double t_3(const double & s) const = 0;
            virtual double t_23(const double & s) const = 0;
    };

    template <>
    class FormFactorFactory<PToV>
    {
        public:
            static std::shared_ptr<FormFactors<PToV>> create(const QualifiedName & name, const Parameters & parameters, const Options & options = Options{ });
    };

    template <>
    class FormFactors<PToP> :
        public virtual ParameterUser
    {
        public:
            virtual ~FormFactors();

            virtual double f_p(const double & s) const = 0;
            virtual double f_0(const double & s) const = 0;
            virtual double f_t(const double & s) const = 0;
            virtual double f_m(const double & s) const;

            virtual double f_p_d1(const double & s) const;
            virtual double f_p_d2(const double & s) const;
    };

    template <>
    class FormFactorFactory<PToP>
    {
        public:
            static std::shared_ptr<FormFactors<PToP>> create(const QualifiedName & label, const Parameters & parameters, const Options & options = Options{ });
    };

    template <>
    class FormFactors<PToPP> :
        public virtual ParameterUser
    {
        public:
            virtual ~FormFactors();

            // form factors
            virtual complex<double> f_perp(const double & q2, const double & k2, const double & z) const = 0;
            virtual complex<double> f_para(const double & q2, const double & k2, const double & z) const = 0;
            virtual complex<double> f_long(const double & q2, const double & k2, const double & z) const = 0;
            virtual complex<double> f_time(const double & q2, const double & k2, const double & z) const = 0;

            double im_f_perp(const double & q2, const double & k2, const double & z) const { return std::imag(f_perp(q2, k2, z)); }
            double im_f_para(const double & q2, const double & k2, const double & z) const { return std::imag(f_para(q2, k2, z)); }
            double im_f_long(const double & q2, const double & k2, const double & z) const { return std::imag(f_long(q2, k2, z)); }
            double im_f_time(const double & q2, const double & k2, const double & z) const { return std::imag(f_time(q2, k2, z)); }

            // residues
            virtual double f_perp_im_res_qhat2(const double & q2, const double & k2) const = 0;
            virtual double f_para_im_res_qhat2(const double & q2, const double & k2) const = 0;
            virtual double f_long_im_res_qhat2(const double & q2, const double & k2) const = 0;
            virtual double f_time_im_res_qhat2(const double & q2, const double & k2) const = 0;
    };

    template <>
    class FormFactorFactory<PToPP>
    {
        public:
            static std::shared_ptr<FormFactors<PToPP>> create(const QualifiedName & name, const Parameters & parameters, const Options & options = Options{ });
    };

    template <>
    class FormFactors<VToP> :
        public virtual ParameterUser
    {
        public:
            virtual ~FormFactors();

            virtual double h_vbar(const double & s) const = 0;

            virtual double h_abar_1(const double & s) const = 0;
            virtual double h_abar_2(const double & s) const = 0;
            virtual double h_abar_3(const double & s) const = 0;
    };

    template <>
    class FormFactorFactory<VToP>
    {
        public:
            static std::shared_ptr<FormFactors<VToP>> create(const QualifiedName & label, const Parameters & parameters, const Options & options = Options{ });
    };

    template <>
    class FormFactors<VToV> :
        public virtual ParameterUser
    {
        public:
            virtual ~FormFactors();

            // vector current
            virtual double h_1(const double & s) const = 0;
            virtual double h_2(const double & s) const = 0;
            virtual double h_3(const double & s) const = 0;
            virtual double h_4(const double & s) const = 0;
            virtual double h_5(const double & s) const = 0;
            virtual double h_6(const double & s) const = 0;

            // axial current
            virtual double h_7(const double & s) const = 0;
            virtual double h_8(const double & s) const = 0;
            virtual double h_9(const double & s) const = 0;
            virtual double h_10(const double & s) const = 0;
    };

    template <>
    class FormFactorFactory<VToV>
    {
        public:
            static std::shared_ptr<FormFactors<VToV>> create(const QualifiedName & label, const Parameters & parameters, const Options & options = Options{ });
    };

    /* Vacuum -> P Gamma */

    template <>
    class FormFactors<VacuumToPGamma> :
        public virtual ParameterUser
    {
        public:
            virtual ~FormFactors();

            virtual complex<double> a(const double & s) const = 0;
            virtual complex<double> v(const double & s) const = 0;
    };

    template <>
    class FormFactorFactory<VacuumToPGamma>
    {
        public:
            static std::shared_ptr<FormFactors<VacuumToPGamma>> create(const QualifiedName & name, const Parameters & parameters, const Options & options = Options{ });
    };
}


#endif
