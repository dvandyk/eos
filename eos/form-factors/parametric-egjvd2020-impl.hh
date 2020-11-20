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
        static constexpr const double m_1 = 0.000137;
        static constexpr const double m_2 = 0.000137;
        static constexpr const double t_p = (m_1 + m_2) * (m_1 + m_2);
        static constexpr const double t_m = (m_1 - m_2) * (m_1 - m_2);
        static constexpr const double Q2  = 2.0;
        static constexpr const bool has_scalar_form_factor = false;
    };

    struct PiToPi {
        typedef PToP Transition;
        static constexpr const char * label = "pi->pi"; //Is this label a good choice?
        static constexpr const double m_1 = 0.000137;
        static constexpr const double m_2 = 0.000137;
        static constexpr const double t_p = (m_1 + m_2) * (m_1 + m_2);
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

                // assumes that Re(q2) > t_p and Im(q2) < 0 such that Im(z) > 0.
                const double re = (q2 + t_0 - 2.0 * t_p) / (q2 - t_0);
                const double im = 2.0 * sqrt((q2 - t_p) * (t_p - t_0)) / (q2 - t_0);

                return complex<double>{ re, im };
            }

        public:
            EGJvD2020FormFactorBase(const Parameters & p, const Options &) :
                _a_fp{{
                    UsedParameter(p[_par_name("+", "0")], *this),
                    UsedParameter(p[_par_name("+", "1")], *this),
                    UsedParameter(p[_par_name("+", "2")], *this),
                    UsedParameter(p[_par_name("+", "3")], *this)
                }},
                _a_ft{{
                    UsedParameter(p[_par_name("T", "0")], *this),
                    UsedParameter(p[_par_name("T", "1")], *this),
                    UsedParameter(p[_par_name("T", "2")], *this),
                    UsedParameter(p[_par_name("T", "3")], *this)
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
                const double t_p = Process_::t_p;
                const double t_0 = this->_t_0;
                const double tfactor = 1.0 - t_0 / t_p;
                const double chi = 0.00405; //GeV^-2 as a makeshift value
                const double Q2 = Process_::Q2;

                const double part0 = 1.0 / sqrt(12.0 * M_PI * t_p * chi);
                const complex<double> part1 = (1.0 + z) * (1.0 + z) * sqrt(1.0 - z) * pow(tfactor, 1.25);
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
                static const std::array<double, 4u> norm
                {{
                    1.0,
                    std::sqrt((1.0 - (3.0 / 7.0) * (3.0 / 7.0))),
                    std::sqrt((1.0 - (3.0 / 7.0) * (3.0 / 7.0)) * (1.0 - (5.0 / 9.0) * (5.0 / 9.0))),
                    std::sqrt((1.0 - (3.0 / 7.0) * (3.0 / 7.0)) * (1.0 - (5.0 / 9.0) * (5.0 / 9.0)) * (1.0 - (3.0 / 11.0) * (3.0 / 11.0)))
                }}; 

                result += _a_fp[0]() * (1.0) / norm[0];
                result += _a_fp[1]() * (- 3.0 / 7.0 + 1.0 * z) / norm[1];
                result += _a_fp[2]() * (5.0 / 9.0 - 2.0 / 3.0 * z + 1.0 * z * z) / norm[2];
                result += _a_fp[3]() * (- 3.0 / 11.0 + 73.0 / 99.0 * z - 9.0 / 11.0 * z * z + 1.0 * z * z * z) / norm[3];

                return result;
            }

            virtual complex<double> f_p(const double & q2) const
            {
                const auto z = _z(q2);

                const auto phi      = this->phi_p(z);
                const auto blaschke = this->blaschke_p(z);
                const auto series   = this->series_p(z);
                const auto asymptotics = (1.0 + z) * (1.0 + z) * sqrt(1.0 - z);

                return phi * blaschke * series * asymptotics;
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

            complex<double> _z(const double & q2) const
            {   
                const double t_p = Process_::t_p;
                const double t_0 = this->_t_0;
                const double a = sqrt(t_p - t_0);
                const double re = (sqrt(t_p - q2) - a) / (sqrt(t_p - q2) + a);
                const double im = 0.0;

                return complex<double>{ re, im };
            }

        public:
            EGJvD2020FormFactorBase(const Parameters & p, const Options &) :
                _a_fp{{
                    UsedParameter(p[_par_name("+", "0")], *this),
                    UsedParameter(p[_par_name("+", "1")], *this),
                    UsedParameter(p[_par_name("+", "2")], *this),
                    UsedParameter(p[_par_name("+", "3")], *this)
                }},
                _a_ft{{
                    UsedParameter(p[_par_name("T", "0")], *this),
                    UsedParameter(p[_par_name("T", "1")], *this),
                    UsedParameter(p[_par_name("T", "2")], *this),
                    UsedParameter(p[_par_name("T", "3")], *this)
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
                const double t_p = Process_::t_p;
                const double t_0 = this->_t_0;
                const double tfactor = 1.0 - t_0 / t_p;
                const double chi = 0.00405; //GeV^-2 as a makeshift value
                const double Q2 = Process_::Q2;

                const double part0 = 1.0 / sqrt(12.0 * M_PI * t_p * chi);
                const complex<double> part1 = (1.0 + z) * (1.0 + z) * sqrt(1.0 - z) * pow(tfactor, 1.25);
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
                static const std::array<double, 4u> norm
                {{
                    1.0,
                    std::sqrt((1.0 - (3.0 / 7.0) * (3.0 / 7.0))),
                    std::sqrt((1.0 - (3.0 / 7.0) * (3.0 / 7.0)) * (1.0 - (5.0 / 9.0) * (5.0 / 9.0))),
                    std::sqrt((1.0 - (3.0 / 7.0) * (3.0 / 7.0)) * (1.0 - (5.0 / 9.0) * (5.0 / 9.0)) * (1.0 - (3.0 / 11.0) * (3.0 / 11.0)))
                }}; 

                result += _a_fp[0]() * (1.0) / norm[0];
                result += _a_fp[1]() * (- 3.0 / 7.0 + 1.0 * z) / norm[1];
                result += _a_fp[2]() * (5.0 / 9.0 - 2.0 / 3.0 * z + 1.0 * z * z) / norm[2];
                result += _a_fp[3]() * (- 3.0 / 11.0 + 73.0 / 99.0 * z - 9.0 / 11.0 * z * z + 1.0 * z * z * z) / norm[3];

                return result;
            }

            virtual complex<double> f_p(const double & q2) const
            {
                const auto z = _z(q2);

                const auto phi      = this->phi_p(z);
                const auto blaschke = this->blaschke_p(z);
                const auto series   = this->series_p(z);
                const auto asymptotics = (1.0 + z) * (1.0 + z) * sqrt(1.0 - z);

                return phi * blaschke * series * asymptotics;
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

    template <typename Process_> class EGJvD2020FormFactors :
        public EGJvD2020FormFactorBase<Process_, typename Process_::Transition, Process_::has_scalar_form_factor>
    {
        public:
            EGJvD2020FormFactors(const Parameters & p, const Options & o) :
                EGJvD2020FormFactorBase<Process_, typename Process_::Transition, Process_::has_scalar_form_factor>(p, o)
            {
            }

            ~EGJvD2020FormFactors() = default;

            static FormFactors<VacuumToPP> * make(const Parameters & parameters, const Options & options)
            {
                return new EGJvD2020FormFactors(parameters, options);
            }
            //TODO: Check This
            static FormFactors<PToP> * make(const Parameters & parameters, const Options & options)
            {
                return new EGJvD2020FormFactors(parameters, options);
            }
    };
}

#endif
