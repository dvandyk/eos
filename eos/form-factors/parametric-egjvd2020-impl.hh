/* vim: set sw=4 sts=4 et tw=120 foldmethod=syntax : */

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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_EGJVD2020_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_EGJVD2020_IMPL_HH 1

#include <eos/form-factors/parametric-egjvd2020.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/stringify.hh>

namespace eos
{
    struct VacuumToPiPi {
        typedef VacuumToPP Transition;
        static constexpr const char * label = "0->pipi";
        static constexpr const double m_1 = 0.137; //GeV
        static constexpr const double m_2 = 0.137; //GeV
        static constexpr const double t_p = (m_1 + m_2) * (m_1 + m_2);
        static constexpr const double t_m = (m_1 - m_2) * (m_1 - m_2);
        static constexpr const double Q2  = 2.0;
        static constexpr const bool has_scalar_form_factor = false;
    };

    struct PiToPi {
        typedef PToP Transition;
        static constexpr const char * label = "0->pipi"; // TODO
        static constexpr const double m_1 = 0.137; //GeV
        static constexpr const double m_2 = 0.137; //GeV
        static constexpr const double t_p = (m_1 + m_2) * (m_1 + m_2); // = 1.87e-8
        static constexpr const double t_m = (m_1 - m_2) * (m_1 - m_2);
        static constexpr const double Q2  = 2.0;
        static constexpr const bool has_scalar_form_factor = false;
    };

    template <typename Process_, typename Transition_, bool has_scalar_form_factor_> class EGJvD2020FormFactorBase;

    template <typename Process_> class EGJvD2020FormFactorBase<Process_, VacuumToPP, false> :
        public FormFactors<VacuumToPP>
    {
        private:
            // parameters for form factors f_+ and f_T
            std::array<UsedParameter, 4u> _a_fp, _a_ft;

            // parameter for zero point of z
            UsedParameter _t_0;

            std::string _par_name(const std::string & ff, const std::string & index) const
            {
                return stringify(Process_::label) + "::" + "a_" + ff + "^" + index + "@EGJvD2020";
            }

            complex<double> _z(const double & q2) const
            {   
                const double t_p = Process_::t_p;
                const double t_0 = this->_t_0;
                if (q2 > t_p) {
                    // assumes that Re(q2) > t_p and Im(q2) < 0 such that Im(z) > 0.
                    const double re = (q2 + t_0 - 2.0 * t_p) / (q2 - t_0);
                    const double im = 2.0 * sqrt((q2 - t_p) * (t_p - t_0)) / (q2 - t_0);
                    return complex<double>{ re, im };
                } else {
                    const double a = sqrt(t_p - t_0);
                    const double re = (sqrt(t_p - q2) - a) / (sqrt(t_p - q2) + a);
                    const double im = 0.0;
                    return complex<double>{ re, im };
                }
            }

        public:
            EGJvD2020FormFactorBase(const Parameters & p, const Options &) :
                _a_fp{{
                    UsedParameter(p[_par_name("+", "0")], *this),
                    UsedParameter(p[_par_name("+", "1")], *this),
                    UsedParameter(p[_par_name("+", "2")], *this),
                    UsedParameter(p[_par_name("+", "3")], *this),
                    UsedParameter(p[_par_name("+", "4")], *this),

                    UsedParameter(p[_par_name("+", "5")], *this),
                    UsedParameter(p[_par_name("+", "6")], *this),
                    UsedParameter(p[_par_name("+", "7")], *this),
                    UsedParameter(p[_par_name("+", "8")], *this),
                    UsedParameter(p[_par_name("+", "9")], *this)
                }},
                _a_ft{{
                    UsedParameter(p[_par_name("T", "0")], *this),
                    UsedParameter(p[_par_name("T", "1")], *this),
                    UsedParameter(p[_par_name("T", "2")], *this),
                    UsedParameter(p[_par_name("T", "3")], *this),
                    UsedParameter(p[_par_name("T", "4")], *this),

                    UsedParameter(p[_par_name("T", "5")], *this),
                    UsedParameter(p[_par_name("T", "6")], *this),
                    UsedParameter(p[_par_name("T", "7")], *this),
                    UsedParameter(p[_par_name("T", "8")], *this),
                    UsedParameter(p[_par_name("T", "9")], *this)
                }},
                _t_0(p[stringify(Process_::label) + "::t_0@EGJvD2020"], *this)
            {
            }

            /* f_+ */

            complex<double> phi_p(const complex<double> & z) const
            {
                // TODO (->EE): implement phi_+
                // TODO implement proper pi and chi(Q2)
                // TODO Is q2 implemented properly here? apperantly not, but how to fix this
                const double t_p = Process_::t_p;           // = 1.87e-8
                const double t_0 = this->_t_0;              // = -2.0
                const double tfactor = 1.0 - t_0 / t_p;     // = 1.0 for t_0 = 0
                const double chi = 0.00405;                 //GeV^-2 as a makeshift value
                const double Q2 = Process_::Q2;             // = 2.0

                const double part0 = 1.0 / sqrt(12.0 * M_PI * t_p * chi);                                               // = 18714.84031
                const complex<double> part1 = (1.0 + z) * (1.0 + z) * sqrt(1.0 - z) * pow(tfactor, 1.25);               // -> 0 (z->1)
                const complex<double> part2 = pow(sqrt(tfactor) * (1.0 + z) + (1.0 - z), -0.5);
                const complex<double> part3 = pow(sqrt(1.0 + Q2 / t_p) * (1.0 - z) + sqrt(tfactor) * (1.0 + z), -3.0);

                return part0 * part1 * part2 * part3;
            }

            complex<double> blaschke_p(const complex<double> & /*z*/) const
            {
                return 1.0;
            }

            complex<double> series_p(const complex<double> & z) const
            {
                complex<double> result = 0.0;

                // TODO(->EE): implement truncated series. How are the coefficients supposed to be fitted.
                // static const std::array<double, 4u> norm
                // {{
                //     1.0,
                //     std::sqrt((1.0 - (3.0 / 7.0) * (3.0 / 7.0))),
                //     std::sqrt((1.0 - (3.0 / 7.0) * (3.0 / 7.0)) * (1.0 - (5.0 / 9.0) * (5.0 / 9.0))),
                //     std::sqrt((1.0 - (3.0 / 7.0) * (3.0 / 7.0)) * (1.0 - (5.0 / 9.0) * (5.0 / 9.0)) * (1.0 - (3.0 / 11.0) * (3.0 / 11.0)))
                // }}; 

                // result += _a_fp[0]() * (1.0) / norm[0];
                // result += _a_fp[1]() * (- 3.0 / 7.0 + 1.0 * z) / norm[1];
                // result += _a_fp[2]() * (5.0 / 9.0 - 2.0 / 3.0 * z + 1.0 * z * z) / norm[2];
                // result += _a_fp[3]() * (- 3.0 / 11.0 + 73.0 / 99.0 * z - 9.0 / 11.0 * z * z + 1.0 * z * z * z) / norm[3];

                result += _a_fp[0]() * (1.0);
                result += _a_fp[1]() * (-0.4743416490252569 + 1.1067971810589328 * z);
                result += _a_fp[2]() * (0.7395099728874526  + -0.8874119674649428 * z + 1.3311179511974138 * z  * z);
                result += _a_fp[3]() * (-0.3773364712030903 + 1.0202060147342804 * z + -1.1320094136092695 * z  * z + 1.383567061077995 * z  * z  * z);
                result += _a_fp[4]() * (0.5764862754491666  + -0.8804517661405453 * z + 1.5303090221014224 * z  * z + -1.3835670610779975 * z  * z  * z + 1.4988643161678288 * z  * z  * z  * z);
                result += _a_fp[5]() * (-0.30595439735640384+ 0.8707932847836088 * z + -1.2109803419840846 * z  * z + 1.7415865695672144 * z  * z  * z + -1.5297719867820112 * z  * z  * z  * z + 1.5297719867820083 * z  * z  * z  * z  * z);
                result += _a_fp[6]() * (0.4707547867785388  + -0.7908680417879451 * z + 1.4470277907438744 * z  * z + -1.6396751342563207 * z  * z  * z + 2.0901512532967006 * z  * z  * z  * z + -1.6947172324027273 * z  * z  * z  * z  * z + 1.6005662750470162 * z  * z  * z  * z  * z  * z);
                result += _a_fp[7]() * (-0.25593140731099673+ 0.7477211703791818 * z + -1.1351310653675895 * z  * z + 1.7275949023825377 * z  * z  * z + -1.8918851089459727 * z  * z  * z  * z + 2.2431635111375248 * z  * z  * z  * z  * z + -1.7915198511769457 * z  * z  * z  * z  * z  * z + 1.6208989129696116 * z  * z  * z  * z  * z  * z  * z);
                result += _a_fp[8]() * (0.39735553782661237 + -0.7026918984723276 * z + 1.319761674861889 * z  * z + -1.6325284610446684 * z  * z  * z + 2.2022600297080053 * z  * z  * z  * z + -2.2261751741518077 * z  * z  * z  * z  * z + 2.4928831636279853 * z  * z  * z  * z  * z  * z + -1.907306581567715 * z  * z  * z  * z  * z  * z  * z + 1.6688932588717427 * z  * z  * z  * z  * z  * z  * z  * z);
                result += _a_fp[9]() * (-0.2195574324091433 + 0.6517022200080864 * z + -1.0367072748341433 * z  * z + 1.624006412919271 * z  * z  * z + -1.9363219782303405 * z  * z  * z  * z + 2.436009619378891 * z  * z  * z  * z  * z + -2.418983641279637 * z  * z  * z  * z  * z  * z + 2.606808880032302 * z  * z  * z  * z  * z  * z  * z + -1.9760168916822334 * z  * z  * z  * z  * z  * z  * z  * z + 1.6832736484700404 * z  * z  * z  * z  * z  * z  * z  * z  * z);

                return result;
            }

            virtual complex<double> f_p(const double & q2) const
            {
                const auto z = _z(q2);

                const auto phi      = this->phi_p(z);
                const auto blaschke = this->blaschke_p(z);
                const auto series   = this->series_p(z);
                const auto asymptotics = (1.0 + z) * (1.0 + z) * sqrt(1.0 - z);

                return 1.0 / (phi * blaschke) * series * asymptotics;
            }

            /* f_T */
            virtual complex<double> f_t(const double & q2) const
            {
                return 1.0;
            }

            /* f_0 */
            virtual complex<double> f_0(const double & /*q2*/) const
            {
                return 0.0; // vanishes exactly
            }
    };
    //TODO: Check this
    template <typename Process_> class EGJvD2020FormFactorBase<Process_, PToP, false> :
        public FormFactors<PToP>
    {
        private:
            // parameters for form factors f_+ and f_T
            std::array<UsedParameter, 4u> _a_fp, _a_ft;

            // parameter for zero point of z
            UsedParameter _t_0;

            std::string _par_name(const std::string & ff, const std::string & index) const
            {
                return stringify(Process_::label) + "::" + "a_" + ff + "^" + index + "@EGJvD2020";
            }

            double _z(const double & q2) const
            {
                const double t_p = Process_::t_p;
                const double t_0 = this->_t_0;
                const double a = sqrt(t_p - t_0);
                const double z = (sqrt(t_p - q2) - a) / (sqrt(t_p - q2) + a);

                return z;
            }

        public:
            EGJvD2020FormFactorBase(const Parameters & p, const Options &) :
                _a_fp{{
                    UsedParameter(p[_par_name("+", "0")], *this),
                    UsedParameter(p[_par_name("+", "1")], *this),
                    UsedParameter(p[_par_name("+", "2")], *this),
                    UsedParameter(p[_par_name("+", "3")], *this),
                    UsedParameter(p[_par_name("+", "4")], *this),

                    UsedParameter(p[_par_name("+", "5")], *this),
                    UsedParameter(p[_par_name("+", "6")], *this),
                    UsedParameter(p[_par_name("+", "7")], *this),
                    UsedParameter(p[_par_name("+", "8")], *this),
                    UsedParameter(p[_par_name("+", "9")], *this)
                }},
                _a_ft{{
                    UsedParameter(p[_par_name("T", "0")], *this),
                    UsedParameter(p[_par_name("T", "1")], *this),
                    UsedParameter(p[_par_name("T", "2")], *this),
                    UsedParameter(p[_par_name("T", "3")], *this),
                    UsedParameter(p[_par_name("T", "4")], *this),

                    UsedParameter(p[_par_name("T", "5")], *this),
                    UsedParameter(p[_par_name("T", "6")], *this),
                    UsedParameter(p[_par_name("T", "7")], *this),
                    UsedParameter(p[_par_name("T", "8")], *this),
                    UsedParameter(p[_par_name("T", "9")], *this)
                }},
                _t_0(p[stringify(Process_::label) + "::t_0@EGJvD2020"], *this)
            {
            }

            /* f_+ */

            double phi_p(const double & z) const
            {
                const double t_p = Process_::t_p;
                const double t_0 = this->_t_0;
                const double tfactor = 1.0 - t_0 / t_p;
                const double chi = 0.00405; //GeV^-2 as a makeshift value
                const double Q2 = Process_::Q2;

                const double part0 = 1.0 / sqrt(12.0 * M_PI * t_p * chi); // = 
                const double part1 = (1.0 + z) * (1.0 + z) * sqrt(1.0 - z) * pow(tfactor, 1.25);
                const double part2 = pow(sqrt(tfactor) * (1.0 + z) + (1.0 - z), -0.5);
                const double part3 = pow(sqrt(1.0 + Q2 / t_p) * (1.0 - z) + sqrt(tfactor) * (1.0 + z), -3.0);

                return part0 * part1 * part2 * part3;
            }

            double blaschke_p(const double & /*z*/) const
            {
                return 1.0;
            }

            double series_p(const double & z) const
            {
                double result = 0.0;

                // static const std::array<double, 4u> norm
                // {{
                //     1.0,
                //     std::sqrt((1.0 - (3.0 / 7.0) * (3.0 / 7.0))),
                //     std::sqrt((1.0 - (3.0 / 7.0) * (3.0 / 7.0)) * (1.0 - (5.0 / 9.0) * (5.0 / 9.0))),
                //     std::sqrt((1.0 - (3.0 / 7.0) * (3.0 / 7.0)) * (1.0 - (5.0 / 9.0) * (5.0 / 9.0)) * (1.0 - (3.0 / 11.0) * (3.0 / 11.0)))
                // }}; 

                // result += _a_fp[0]() * (1.0) / norm[0];
                // result += _a_fp[1]() * (- 3.0 / 7.0 + 1.0 * z) / norm[1];
                // result += _a_fp[2]() * (5.0 / 9.0 - 2.0 / 3.0 * z + 1.0 * z * z) / norm[2];
                // result += _a_fp[3]() * (- 3.0 / 11.0 + 73.0 / 99.0 * z - 9.0 / 11.0 * z * z + 1.0 * z * z * z) / norm[3];

                result += _a_fp[0]() * (1.0);
                result += _a_fp[1]() * (-0.4743416490252569 + 1.1067971810589328 * z);
                result += _a_fp[2]() * (0.7395099728874526  + -0.8874119674649428 * z + 1.3311179511974138 * z  * z);
                result += _a_fp[3]() * (-0.3773364712030903 + 1.0202060147342804 * z + -1.1320094136092695 * z  * z + 1.383567061077995 * z  * z  * z);
                result += _a_fp[4]() * (0.5764862754491666  + -0.8804517661405453 * z + 1.5303090221014224 * z  * z + -1.3835670610779975 * z  * z  * z + 1.4988643161678288 * z  * z  * z  * z);
                result += _a_fp[5]() * (-0.30595439735640384+ 0.8707932847836088 * z + -1.2109803419840846 * z  * z + 1.7415865695672144 * z  * z  * z + -1.5297719867820112 * z  * z  * z  * z + 1.5297719867820083 * z  * z  * z  * z  * z);
                result += _a_fp[6]() * (0.4707547867785388  + -0.7908680417879451 * z + 1.4470277907438744 * z  * z + -1.6396751342563207 * z  * z  * z + 2.0901512532967006 * z  * z  * z  * z + -1.6947172324027273 * z  * z  * z  * z  * z + 1.6005662750470162 * z  * z  * z  * z  * z  * z);
                result += _a_fp[7]() * (-0.25593140731099673+ 0.7477211703791818 * z + -1.1351310653675895 * z  * z + 1.7275949023825377 * z  * z  * z + -1.8918851089459727 * z  * z  * z  * z + 2.2431635111375248 * z  * z  * z  * z  * z + -1.7915198511769457 * z  * z  * z  * z  * z  * z + 1.6208989129696116 * z  * z  * z  * z  * z  * z  * z);
                result += _a_fp[8]() * (0.39735553782661237 + -0.7026918984723276 * z + 1.319761674861889 * z  * z + -1.6325284610446684 * z  * z  * z + 2.2022600297080053 * z  * z  * z  * z + -2.2261751741518077 * z  * z  * z  * z  * z + 2.4928831636279853 * z  * z  * z  * z  * z  * z + -1.907306581567715 * z  * z  * z  * z  * z  * z  * z + 1.6688932588717427 * z  * z  * z  * z  * z  * z  * z  * z);
                result += _a_fp[9]() * (-0.2195574324091433 + 0.6517022200080864 * z + -1.0367072748341433 * z  * z + 1.624006412919271 * z  * z  * z + -1.9363219782303405 * z  * z  * z  * z + 2.436009619378891 * z  * z  * z  * z  * z + -2.418983641279637 * z  * z  * z  * z  * z  * z + 2.606808880032302 * z  * z  * z  * z  * z  * z  * z + -1.9760168916822334 * z  * z  * z  * z  * z  * z  * z  * z + 1.6832736484700404 * z  * z  * z  * z  * z  * z  * z  * z  * z);


                return result;
            }

            virtual double f_p(const double & q2) const
            {
                const auto z = _z(q2);

                const auto phi      = this->phi_p(z);
                const auto blaschke = this->blaschke_p(z);
                const auto series   = this->series_p(z);
                const auto asymptotics = (1.0 + z) * (1.0 + z) * sqrt(1.0 - z);

                return 1.0 / (phi * blaschke) * series * asymptotics;
            }

            /* f_T */
            virtual double f_t(const double & q2) const
            {
                return 1.0;
            }

            /* f_0 */
            virtual double f_0(const double & /*q2*/) const
            {
                return 0.0; // vanishes exactly
            }
    };

    template <typename Process_> class EGJvD2020FormFactors :
        public EGJvD2020FormFactorBase<Process_, typename Process_::Transition, Process_::has_scalar_form_factor>
    {
        public:
            EGJvD2020FormFactors(const Parameters & p, const Options & o) :
                EGJvD2020FormFactorBase<Process_, typename Process_::Transition, Process_::has_scalar_form_factor>(p, o)
            {
            }

            ~EGJvD2020FormFactors() = default;

            static FormFactors<typename Process_::Transition> * make(const Parameters & parameters, const Options & options)
            {
                return new EGJvD2020FormFactors(parameters, options);
            }
    };
}

#endif
