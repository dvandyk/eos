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
#include <eos/utils/options-impl.hh>
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
        static constexpr const double epsilon_isospin = 0.05;
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
        static constexpr const double epsilon_isospin = 0.05;
    };

    template <typename Process_, typename Transition_, bool has_scalar_form_factor_> class EGJvD2020UnitarityBoundsBase;

    template <typename Process_> class EGJvD2020UnitarityBoundsBase<Process_, VacuumToPP, false> :
        public virtual ParameterUser
    {
        private:
            // parameters for form factors f_+ and f_T
            std::array<UsedParameter, 50u> _a_fp, _a_ft;

            std::string _par_name(const std::string & ff, const std::string & index) const
            {
                return stringify(Process_::label) + "::" + "a_" + ff + "^" + index + "@EGJvD2020";
            }
        public:
            EGJvD2020UnitarityBoundsBase(const Parameters & p, const Options & o) :
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
                    UsedParameter(p[_par_name("+", "9")], *this),

                    UsedParameter(p[_par_name("+", "10")], *this),
                    UsedParameter(p[_par_name("+", "11")], *this),
                    UsedParameter(p[_par_name("+", "12")], *this),
                    UsedParameter(p[_par_name("+", "13")], *this),
                    UsedParameter(p[_par_name("+", "14")], *this),

                    UsedParameter(p[_par_name("+", "15")], *this),
                    UsedParameter(p[_par_name("+", "16")], *this),
                    UsedParameter(p[_par_name("+", "17")], *this),
                    UsedParameter(p[_par_name("+", "18")], *this),
                    UsedParameter(p[_par_name("+", "19")], *this),

                    UsedParameter(p[_par_name("+", "20")], *this),
                    UsedParameter(p[_par_name("+", "21")], *this),
                    UsedParameter(p[_par_name("+", "22")], *this),
                    UsedParameter(p[_par_name("+", "23")], *this),
                    UsedParameter(p[_par_name("+", "24")], *this),

                    UsedParameter(p[_par_name("+", "25")], *this),
                    UsedParameter(p[_par_name("+", "26")], *this),
                    UsedParameter(p[_par_name("+", "27")], *this),
                    UsedParameter(p[_par_name("+", "28")], *this),
                    UsedParameter(p[_par_name("+", "29")], *this),

                    UsedParameter(p[_par_name("+", "30")], *this),
                    UsedParameter(p[_par_name("+", "31")], *this),
                    UsedParameter(p[_par_name("+", "32")], *this),
                    UsedParameter(p[_par_name("+", "33")], *this),
                    UsedParameter(p[_par_name("+", "34")], *this),

                    UsedParameter(p[_par_name("+", "35")], *this),
                    UsedParameter(p[_par_name("+", "36")], *this),
                    UsedParameter(p[_par_name("+", "37")], *this),
                    UsedParameter(p[_par_name("+", "38")], *this),
                    UsedParameter(p[_par_name("+", "39")], *this),

                    UsedParameter(p[_par_name("+", "40")], *this),
                    UsedParameter(p[_par_name("+", "41")], *this),
                    UsedParameter(p[_par_name("+", "42")], *this),
                    UsedParameter(p[_par_name("+", "43")], *this),
                    UsedParameter(p[_par_name("+", "44")], *this),

                    UsedParameter(p[_par_name("+", "45")], *this),
                    UsedParameter(p[_par_name("+", "46")], *this),
                    UsedParameter(p[_par_name("+", "47")], *this),
                    UsedParameter(p[_par_name("+", "48")], *this),
                    UsedParameter(p[_par_name("+", "49")], *this)
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
                    UsedParameter(p[_par_name("T", "9")], *this),

                    UsedParameter(p[_par_name("T", "10")], *this),
                    UsedParameter(p[_par_name("T", "11")], *this),
                    UsedParameter(p[_par_name("T", "12")], *this),
                    UsedParameter(p[_par_name("T", "13")], *this),
                    UsedParameter(p[_par_name("T", "14")], *this),

                    UsedParameter(p[_par_name("T", "15")], *this),
                    UsedParameter(p[_par_name("T", "16")], *this),
                    UsedParameter(p[_par_name("T", "17")], *this),
                    UsedParameter(p[_par_name("T", "18")], *this),
                    UsedParameter(p[_par_name("T", "19")], *this),

                    UsedParameter(p[_par_name("T", "20")], *this),
                    UsedParameter(p[_par_name("T", "21")], *this),
                    UsedParameter(p[_par_name("T", "22")], *this),
                    UsedParameter(p[_par_name("T", "23")], *this),
                    UsedParameter(p[_par_name("T", "24")], *this),

                    UsedParameter(p[_par_name("T", "25")], *this),
                    UsedParameter(p[_par_name("T", "26")], *this),
                    UsedParameter(p[_par_name("T", "27")], *this),
                    UsedParameter(p[_par_name("T", "28")], *this),
                    UsedParameter(p[_par_name("T", "29")], *this),

                    UsedParameter(p[_par_name("T", "30")], *this),
                    UsedParameter(p[_par_name("T", "31")], *this),
                    UsedParameter(p[_par_name("T", "32")], *this),
                    UsedParameter(p[_par_name("T", "33")], *this),
                    UsedParameter(p[_par_name("T", "34")], *this),

                    UsedParameter(p[_par_name("T", "35")], *this),
                    UsedParameter(p[_par_name("T", "36")], *this),
                    UsedParameter(p[_par_name("T", "37")], *this),
                    UsedParameter(p[_par_name("T", "38")], *this),
                    UsedParameter(p[_par_name("T", "39")], *this),

                    UsedParameter(p[_par_name("T", "40")], *this),
                    UsedParameter(p[_par_name("T", "41")], *this),
                    UsedParameter(p[_par_name("T", "42")], *this),
                    UsedParameter(p[_par_name("T", "43")], *this),
                    UsedParameter(p[_par_name("T", "44")], *this),

                    UsedParameter(p[_par_name("T", "45")], *this),
                    UsedParameter(p[_par_name("T", "46")], *this),
                    UsedParameter(p[_par_name("T", "47")], *this),
                    UsedParameter(p[_par_name("T", "48")], *this),
                    UsedParameter(p[_par_name("T", "49")], *this)
                }}
            {
            }

            ~EGJvD2020UnitarityBoundsBase() = default;

            double bound_1m_prior() const
            {
                const double value = bound_1m();

                if (value < 0.0)
                {
                    throw InternalError("Contribution to 1^- unitarity bound must be positive; found to be negative!");
                }
                else if ((0.0 <= value) && (value < 1.0))
                {
                    return 0.0;
                }
                else
                {
                    // add an r-fit like penalty
                    static const double sigma = 0.05; // TODO preliminary assuming 5% uncertainty
                    return -pow((value - 1.0) / sigma, 2) / 2.0;
                }
            }

            double bound_1m() const
            {
                double sum = 0.0;
                for (auto & a : _a_fp)
                {
                    sum += a * a;
                }

                return sum;
            }
    };

    template <typename Process_> class EGJvD2020UnitarityBounds:
        public EGJvD2020UnitarityBoundsBase<Process_, typename Process_::Transition, Process_::has_scalar_form_factor>
    {
        public:
            EGJvD2020UnitarityBounds(const Parameters & p, const Options & o) :
                EGJvD2020UnitarityBoundsBase<Process_, typename Process_::Transition, Process_::has_scalar_form_factor>(p, o)
            {
            }

            ~EGJvD2020UnitarityBounds() = default;
    };

    template <typename Process_, typename Transition_, bool has_scalar_form_factor_> class EGJvD2020FormFactorBase;

    template <typename Process_> class EGJvD2020FormFactorBase<Process_, VacuumToPP, false> :
        public FormFactors<VacuumToPP>
    {
        private:
            // parameters for form factors f_+ and f_T
            std::array<UsedParameter, 50u> _a_fp, _a_ft;

            // parameter for zero point of z
            UsedParameter _t_0;


            std::string _par_name(const std::string & ff, const std::string & index) const
            {
                return stringify(Process_::label) + "::" + "a_" + ff + "^" + index + "@EGJvD2020";
            }

            complex<double> _z(const double & q2, const complex<double> & t_s) const
            {
                const complex<double> t_p = complex<double>{ Process_::t_p, 0.0};
                const complex<double> t   = complex<double>{ q2, 0.0};

                return (sqrt(t_p - t) - sqrt(t_p - t_s)) / (sqrt(t_p - t) + sqrt(t_p - t_s));
            }

            // complex<double> _z(const double & q2, const double & t_0) const
            // {
            //     const double t_p = Process_::t_p;
            //     if (q2 > t_p) {
            //         // assumes that Re(q2) > t_p and Im(q2) < 0 such that Im(z) > 0.
            //         const double re = (q2 + t_0 - 2.0 * t_p) / (q2 - t_0);
            //         const double im = 2.0 * sqrt((q2 - t_p) * (t_p - t_0)) / (q2 - t_0);
            //         return complex<double>{ re, im };
            //     } else {
            //         const double a = sqrt(t_p - t_0);
            //         const double re = (sqrt(t_p - q2) - a) / (sqrt(t_p - q2) + a);
            //         const double im = 0.0;
            //         return complex<double>{ re, im };
            //     }
            // }

            complex<double> _z(const double & q2) const
            {
                const complex<double> t_s = complex<double>(this->_t_0, 0.0);
                return _z(q2, t_s);
            }

        public:
            EGJvD2020FormFactorBase(const Parameters & p, const Options & o) :
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
                    UsedParameter(p[_par_name("+", "9")], *this),

                    UsedParameter(p[_par_name("+", "10")], *this),
                    UsedParameter(p[_par_name("+", "11")], *this),
                    UsedParameter(p[_par_name("+", "12")], *this),
                    UsedParameter(p[_par_name("+", "13")], *this),
                    UsedParameter(p[_par_name("+", "14")], *this),

                    UsedParameter(p[_par_name("+", "15")], *this),
                    UsedParameter(p[_par_name("+", "16")], *this),
                    UsedParameter(p[_par_name("+", "17")], *this),
                    UsedParameter(p[_par_name("+", "18")], *this),
                    UsedParameter(p[_par_name("+", "19")], *this),

                    UsedParameter(p[_par_name("+", "20")], *this),
                    UsedParameter(p[_par_name("+", "21")], *this),
                    UsedParameter(p[_par_name("+", "22")], *this),
                    UsedParameter(p[_par_name("+", "23")], *this),
                    UsedParameter(p[_par_name("+", "24")], *this),

                    UsedParameter(p[_par_name("+", "25")], *this),
                    UsedParameter(p[_par_name("+", "26")], *this),
                    UsedParameter(p[_par_name("+", "27")], *this),
                    UsedParameter(p[_par_name("+", "28")], *this),
                    UsedParameter(p[_par_name("+", "29")], *this),

                    UsedParameter(p[_par_name("+", "30")], *this),
                    UsedParameter(p[_par_name("+", "31")], *this),
                    UsedParameter(p[_par_name("+", "32")], *this),
                    UsedParameter(p[_par_name("+", "33")], *this),
                    UsedParameter(p[_par_name("+", "34")], *this),

                    UsedParameter(p[_par_name("+", "35")], *this),
                    UsedParameter(p[_par_name("+", "36")], *this),
                    UsedParameter(p[_par_name("+", "37")], *this),
                    UsedParameter(p[_par_name("+", "38")], *this),
                    UsedParameter(p[_par_name("+", "39")], *this),

                    UsedParameter(p[_par_name("+", "40")], *this),
                    UsedParameter(p[_par_name("+", "41")], *this),
                    UsedParameter(p[_par_name("+", "42")], *this),
                    UsedParameter(p[_par_name("+", "43")], *this),
                    UsedParameter(p[_par_name("+", "44")], *this),

                    UsedParameter(p[_par_name("+", "45")], *this),
                    UsedParameter(p[_par_name("+", "46")], *this),
                    UsedParameter(p[_par_name("+", "47")], *this),
                    UsedParameter(p[_par_name("+", "48")], *this),
                    UsedParameter(p[_par_name("+", "49")], *this)
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
                    UsedParameter(p[_par_name("T", "9")], *this),

                    UsedParameter(p[_par_name("T", "10")], *this),
                    UsedParameter(p[_par_name("T", "11")], *this),
                    UsedParameter(p[_par_name("T", "12")], *this),
                    UsedParameter(p[_par_name("T", "13")], *this),
                    UsedParameter(p[_par_name("T", "14")], *this),

                    UsedParameter(p[_par_name("T", "15")], *this),
                    UsedParameter(p[_par_name("T", "16")], *this),
                    UsedParameter(p[_par_name("T", "17")], *this),
                    UsedParameter(p[_par_name("T", "18")], *this),
                    UsedParameter(p[_par_name("T", "19")], *this),

                    UsedParameter(p[_par_name("T", "20")], *this),
                    UsedParameter(p[_par_name("T", "21")], *this),
                    UsedParameter(p[_par_name("T", "22")], *this),
                    UsedParameter(p[_par_name("T", "23")], *this),
                    UsedParameter(p[_par_name("T", "24")], *this),

                    UsedParameter(p[_par_name("T", "25")], *this),
                    UsedParameter(p[_par_name("T", "26")], *this),
                    UsedParameter(p[_par_name("T", "27")], *this),
                    UsedParameter(p[_par_name("T", "28")], *this),
                    UsedParameter(p[_par_name("T", "29")], *this),

                    UsedParameter(p[_par_name("T", "30")], *this),
                    UsedParameter(p[_par_name("T", "31")], *this),
                    UsedParameter(p[_par_name("T", "32")], *this),
                    UsedParameter(p[_par_name("T", "33")], *this),
                    UsedParameter(p[_par_name("T", "34")], *this),

                    UsedParameter(p[_par_name("T", "35")], *this),
                    UsedParameter(p[_par_name("T", "36")], *this),
                    UsedParameter(p[_par_name("T", "37")], *this),
                    UsedParameter(p[_par_name("T", "38")], *this),
                    UsedParameter(p[_par_name("T", "39")], *this),

                    UsedParameter(p[_par_name("T", "40")], *this),
                    UsedParameter(p[_par_name("T", "41")], *this),
                    UsedParameter(p[_par_name("T", "42")], *this),
                    UsedParameter(p[_par_name("T", "43")], *this),
                    UsedParameter(p[_par_name("T", "44")], *this),

                    UsedParameter(p[_par_name("T", "45")], *this),
                    UsedParameter(p[_par_name("T", "46")], *this),
                    UsedParameter(p[_par_name("T", "47")], *this),
                    UsedParameter(p[_par_name("T", "48")], *this),
                    UsedParameter(p[_par_name("T", "49")], *this)
                }},
                _t_0(p[stringify(Process_::label) + "::t_0@EGJvD2020"], *this)
            {
                //empty body
            }
            /* f_+ */

            complex<double> phi_p(const complex<double> & z) const //Divided by asymptotic behaviour to improve numeric behaviour of Formfactor
            {
                const double t_p = Process_::t_p;
                const double t_0 = this->_t_0;
                const double tfactor = 1.0 - t_0 / t_p;
                const double chi = 0.00405;                 //GeV^-2 as a makeshift value
                const double Q2 = Process_::Q2;             // = 2.0

                const double part0 = 1.0 / sqrt(12.0 * M_PI * t_p * chi);
                const complex<double> part1 = /*[Asymptotics:] (1.0 + z) * (1.0 + z) * sqrt(1.0 - z)*/ pow(tfactor, 1.25);
                const complex<double> part2 = pow(sqrt(tfactor) * (1.0 + z) + (1.0 - z), -0.5);
                const complex<double> part3 = pow(sqrt(1.0 + Q2 / t_p) * (1.0 - z) + sqrt(tfactor) * (1.0 + z), -3.0);

                return part0 * part1 * part2 * part3;
            }

            complex<double> series_p(const complex<double> & z, std::array<UsedParameter, 50u> c_f) const
            {   
                complex<double> result = 0.0;

                result += c_f[0]() *  1.0 * (1.0 * pow(z, 0));
                result += c_f[1]() *  1.0041580220928046 * (1.0 / 11.0 * pow(z, 0) + 1.0 * pow(z, 1));
                result += c_f[2]() *  1.3915668626887223 * (9.0 / 13.0 * pow(z, 0) + 2.0 / 13.0 * pow(z, 1) + 1.0 * pow(z, 2));
                result += c_f[3]() *  1.394669579723865 * (1.0 / 15.0 * pow(z, 0) + 137.0 / 195.0 * pow(z, 1) + 1.0 / 5.0 * pow(z, 2) + 1.0 * pow(z, 3));
                result += c_f[4]() *  1.6439499152771448 * (9.0 / 17.0 * pow(z, 0) + 44.0 / 255.0 * pow(z, 1) + 274.0 / 255.0 * pow(z, 2) + 4.0 / 17.0 * pow(z, 3) + 1.0 * pow(z, 4));
                result += c_f[5]() *  1.6462315956469282 * (1.0 / 19.0 * pow(z, 0) + 175.0 / 323.0 * pow(z, 1) + 74.0 / 323.0 * pow(z, 2) + 350.0 / 323.0 * pow(z, 3) + 5.0 / 19.0 * pow(z, 4) + 1.0 * pow(z, 5));
                result += c_f[6]() *  1.822044489432169 * (3.0 / 7.0 * pow(z, 0) + 22.0 / 133.0 * pow(z, 1) + 325.0 / 323.0 * pow(z, 2) + 740.0 / 2261.0 * pow(z, 3) + 25.0 / 19.0 * pow(z, 4) + 2.0 / 7.0 * pow(z, 5) + 1.0 * pow(z, 6));
                result += c_f[7]() *  1.8237690941622522 * (1.0 / 23.0 * pow(z, 0) + 71.0 / 161.0 * pow(z, 1) + 681.0 / 3059.0 * pow(z, 2) + 53065.0 / 52003.0 * pow(z, 3) + 1135.0 / 3059.0 * pow(z, 4) + 213.0 / 161.0 * pow(z, 5) + 7.0 / 23.0 * pow(z, 6) + 1.0 * pow(z, 7));
                result += c_f[8]() *  1.9548363704716327 * (9.0 / 25.0 * pow(z, 0) + 88.0 / 575.0 * pow(z, 1) + 3692.0 / 4025.0 * pow(z, 2) + 5448.0 / 15295.0 * pow(z, 3) + 21226.0 / 15295.0 * pow(z, 4) + 1816.0 / 4025.0 * pow(z, 5) + 852.0 / 575.0 * pow(z, 6) + 8.0 / 25.0 * pow(z, 7) + 1.0 * pow(z, 8));
                result += c_f[9]() *  1.9561785171250927 * (1.0 / 27.0 * pow(z, 0) + 251.0 / 675.0 * pow(z, 1) + 1076.0 / 5175.0 * pow(z, 2) + 580.0 / 621.0 * pow(z, 3) + 24046.0 / 58995.0 * pow(z, 4) + 290.0 / 207.0 * pow(z, 5) + 7532.0 / 15525.0 * pow(z, 6) + 1004.0 / 675.0 * pow(z, 7) + 1.0 / 3.0 * pow(z, 8) + 1.0 * pow(z, 9));
                result += c_f[10]() *  2.05778352996703 * (9.0 / 29.0 * pow(z, 0) + 110.0 / 783.0 * pow(z, 1) + 3263.0 / 3915.0 * pow(z, 2) + 2152.0 / 6003.0 * pow(z, 3) + 850.0 / 621.0 * pow(z, 4) + 48092.0 / 90045.0 * pow(z, 5) + 350.0 / 207.0 * pow(z, 6) + 2152.0 / 3915.0 * pow(z, 7) + 1255.0 / 783.0 * pow(z, 8) + 10.0 / 29.0 * pow(z, 9) + 1.0 * pow(z, 10));
                result += c_f[11]() *  2.0588550132627392 * (1.0 / 31.0 * pow(z, 0) + 289.0 / 899.0 * pow(z, 1) + 1555.0 / 8091.0 * pow(z, 2) + 6887.0 / 8091.0 * pow(z, 3) + 76862.0 / 186093.0 * pow(z, 4) + 1289614.0 / 930465.0 * pow(z, 5) + 538034.0 / 930465.0 * pow(z, 6) + 13774.0 / 8091.0 * pow(z, 7) + 1555.0 / 2697.0 * pow(z, 8) + 1445.0 / 899.0 * pow(z, 9) + 11.0 / 31.0 * pow(z, 10) + 1.0 * pow(z, 11));
                result += c_f[12]() *  2.1399786377482064 * (3.0 / 11.0 * pow(z, 0) + 4.0 / 31.0 * pow(z, 1) + 7514.0 / 9889.0 * pow(z, 2) + 31100.0 / 89001.0 * pow(z, 3) + 117079.0 / 89001.0 * pow(z, 4) + 5841512.0 / 10235115.0 * pow(z, 5) + 18054596.0 / 10235115.0 * pow(z, 6) + 307448.0 / 445005.0 * pow(z, 7) + 172175.0 / 89001.0 * pow(z, 8) + 6220.0 / 9889.0 * pow(z, 9) + 578.0 / 341.0 * pow(z, 10) + 4.0 / 11.0 * pow(z, 11) + 1.0 * pow(z, 12));
                result += c_f[13]() *  2.140852633552563 * (1.0 / 35.0 * pow(z, 0) + 109.0 / 385.0 * pow(z, 1) + 2118.0 / 11935.0 * pow(z, 2) + 53842.0 / 69223.0 * pow(z, 3) + 28015.0 / 69223.0 * pow(z, 4) + 2310697.0 / 1730575.0 * pow(z, 5) + 963236.0 / 1550775.0 * pow(z, 6) + 9242788.0 / 5191725.0 * pow(z, 7) + 50427.0 / 69223.0 * pow(z, 8) + 134605.0 / 69223.0 * pow(z, 9) + 706.0 / 1085.0 * pow(z, 10) + 654.0 / 385.0 * pow(z, 11) + 13.0 / 35.0 * pow(z, 12) + 1.0 * pow(z, 13));
                result += c_f[14]() *  2.207143478674807 * (9.0 / 37.0 * pow(z, 0) + 22.0 / 185.0 * pow(z, 1) + 1417.0 / 2035.0 * pow(z, 2) + 4236.0 / 12617.0 * pow(z, 3) + 457657.0 / 365893.0 * pow(z, 4) + 212914.0 / 365893.0 * pow(z, 5) + 16174879.0 / 9147325.0 * pow(z, 6) + 1926472.0 / 2494725.0 * pow(z, 7) + 2310697.0 / 1097679.0 * pow(z, 8) + 302562.0 / 365893.0 * pow(z, 9) + 26921.0 / 12617.0 * pow(z, 10) + 1412.0 / 2035.0 * pow(z, 11) + 327.0 / 185.0 * pow(z, 12) + 14.0 / 37.0 * pow(z, 13) + 1.0 * pow(z, 14));
                result += c_f[15]() *  2.207869393339616 * (1.0 / 39.0 * pow(z, 0) + 365.0 / 1443.0 * pow(z, 1) + 79.0 / 481.0 * pow(z, 2) + 11335.0 / 15873.0 * pow(z, 3) + 192125.0 / 492063.0 * pow(z, 4) + 465415.0 / 365893.0 * pow(z, 5) + 27221635.0 / 42809481.0 * pow(z, 6) + 382730407.0 / 214047405.0 * pow(z, 7) + 3888805.0 / 4756609.0 * pow(z, 8) + 2327075.0 / 1097679.0 * pow(z, 9) + 38425.0 / 44733.0 * pow(z, 10) + 11335.0 / 5291.0 * pow(z, 11) + 79.0 / 111.0 * pow(z, 12) + 2555.0 / 1443.0 * pow(z, 13) + 5.0 / 13.0 * pow(z, 14) + 1.0 * pow(z, 15));
                result += c_f[16]() *  2.2630661281731066 * (9.0 / 41.0 * pow(z, 0) + 176.0 / 1599.0 * pow(z, 1) + 2920.0 / 4551.0 * pow(z, 2) + 6320.0 / 19721.0 * pow(z, 3) + 770780.0 / 650793.0 * pow(z, 4) + 11681200.0 / 20174583.0 * pow(z, 5) + 26063240.0 / 15001613.0 * pow(z, 6) + 1431080240.0 / 1755188721.0 * pow(z, 7) + 3827304070.0 / 1755188721.0 * pow(z, 8) + 186662640.0 / 195020969.0 * pow(z, 9) + 3723320.0 / 1551891.0 * pow(z, 10) + 614800.0 / 650793.0 * pow(z, 11) + 45340.0 / 19721.0 * pow(z, 12) + 44240.0 / 59163.0 * pow(z, 13) + 2920.0 / 1599.0 * pow(z, 14) + 16.0 / 41.0 * pow(z, 15) + 1.0 * pow(z, 16));
                result += c_f[17]() *  2.2636783468041286 * (1.0 / 43.0 * pow(z, 0) + 403.0 / 1763.0 * pow(z, 1) + 3496.0 / 22919.0 * pow(z, 2) + 558840.0 / 848003.0 * pow(z, 3) + 317100.0 / 848003.0 * pow(z, 4) + 1022980.0 / 848003.0 * pow(z, 5) + 16687720.0 / 26288093.0 * pow(z, 6) + 1341455800.0 / 762354697.0 * pow(z, 7) + 1980719830.0 / 2287064091.0 * pow(z, 8) + 1676819750.0 / 762354697.0 * pow(z, 9) + 26223560.0 / 26288093.0 * pow(z, 10) + 2045960.0 / 848003.0 * pow(z, 11) + 63420.0 / 65231.0 * pow(z, 12) + 1955940.0 / 848003.0 * pow(z, 13) + 17480.0 / 22919.0 * pow(z, 14) + 3224.0 / 1763.0 * pow(z, 15) + 17.0 / 43.0 * pow(z, 16) + 1.0 * pow(z, 17));
                result += c_f[18]() *  2.310357038107123 * (1.0 / 5.0 * pow(z, 0) + 22.0 / 215.0 * pow(z, 1) + 5239.0 / 8815.0 * pow(z, 2) + 6992.0 / 22919.0 * pow(z, 3) + 950028.0 / 848003.0 * pow(z, 4) + 481992.0 / 848003.0 * pow(z, 5) + 1432172.0 / 848003.0 * pow(z, 6) + 21932432.0 / 26288093.0 * pow(z, 7) + 1676819750.0 / 762354697.0 * pow(z, 8) + 792287932.0 / 762354697.0 * pow(z, 9) + 67072790.0 / 26288093.0 * pow(z, 10) + 953584.0 / 848003.0 * pow(z, 11) + 2250556.0 / 848003.0 * pow(z, 12) + 887880.0 / 848003.0 * pow(z, 13) + 55884.0 / 22919.0 * pow(z, 14) + 6992.0 / 8815.0 * pow(z, 15) + 403.0 / 215.0 * pow(z, 16) + 2.0 / 5.0 * pow(z, 17) + 1.0 * pow(z, 18));
                result += c_f[19]() *  2.310880157560925 * (1.0 / 47.0 * pow(z, 0) + 49.0 / 235.0 * pow(z, 1) + 1437.0 / 10105.0 * pow(z, 2) + 50645.0 / 82861.0 * pow(z, 3) + 384508.0 / 1077193.0 * pow(z, 4) + 45539196.0 / 39856141.0 * pow(z, 5) + 24904180.0 / 39856141.0 * pow(z, 6) + 68265668.0 / 39856141.0 * pow(z, 7) + 1097897094.0 / 1235540371.0 * pow(z, 8) + 1941532102.0 / 873918799.0 * pow(z, 9) + 1341874226.0 / 1235540371.0 * pow(z, 10) + 102398502.0 / 39856141.0 * pow(z, 11) + 3557740.0 / 3065857.0 * pow(z, 12) + 106258124.0 / 39856141.0 * pow(z, 13) + 1153524.0 / 1077193.0 * pow(z, 14) + 202580.0 / 82861.0 * pow(z, 15) + 8143.0 / 10105.0 * pow(z, 16) + 441.0 / 235.0 * pow(z, 17) + 19.0 / 47.0 * pow(z, 18) + 1.0 * pow(z, 19));
                result += c_f[20]() *  2.350874856721882 * (9.0 / 49.0 * pow(z, 0) + 220.0 / 2303.0 * pow(z, 1) + 26.0 / 47.0 * pow(z, 2) + 28740.0 / 99029.0 * pow(z, 3) + 614975.0 / 580027.0 * pow(z, 4) + 29222608.0 / 52782457.0 * pow(z, 5) + 455391960.0 / 278992987.0 * pow(z, 6) + 1636560400.0 / 1952950909.0 * pow(z, 7) + 4266604250.0 / 1952950909.0 * pow(z, 8) + 9410546520.0 / 8648782597.0 * pow(z, 9) + 3883064204.0 / 1476621419.0 * pow(z, 10) + 348538760.0 / 278992987.0 * pow(z, 11) + 5631917610.0 / 1952950909.0 * pow(z, 12) + 355774000.0 / 278992987.0 * pow(z, 13) + 151797320.0 / 52782457.0 * pow(z, 14) + 4614096.0 / 4060189.0 * pow(z, 15) + 36175.0 / 14147.0 * pow(z, 16) + 1916.0 / 2303.0 * pow(z, 17) + 90.0 / 47.0 * pow(z, 18) + 20.0 / 49.0 * pow(z, 19) + 1.0 * pow(z, 20));
                result += c_f[21]() *  2.351326904578521 * (1.0 / 51.0 * pow(z, 0) + 479.0 / 2499.0 * pow(z, 1) + 5210.0 / 39151.0 * pow(z, 2) + 66890.0 / 117453.0 * pow(z, 3) + 1718965.0 / 5050479.0 * pow(z, 4) + 74720057.0 / 69023213.0 * pow(z, 5) + 2577944.0 / 4225911.0 * pow(z, 6) + 1813904920.0 / 1094510949.0 * pow(z, 7) + 2284525590.0 / 2553858881.0 * pow(z, 8) + 995640670.0 / 450680979.0 * pow(z, 9) + 7315531684.0 / 6419158809.0 * pow(z, 10) + 398256268.0 / 150226993.0 * pow(z, 11) + 9899610890.0 / 7661576643.0 * pow(z, 12) + 453476230.0 / 156358707.0 * pow(z, 13) + 12889720.0 / 9860459.0 * pow(z, 14) + 597760456.0 / 207069639.0 * pow(z, 15) + 343793.0 / 297087.0 * pow(z, 16) + 100335.0 / 39151.0 * pow(z, 17) + 98990.0 / 117453.0 * pow(z, 18) + 4790.0 / 2499.0 * pow(z, 19) + 7.0 / 17.0 * pow(z, 20) + 1.0 * pow(z, 21));
                result += c_f[22]() *  2.3859794553105305 * (9.0 / 53.0 * pow(z, 0) + 242.0 / 2703.0 * pow(z, 1) + 68497.0 / 132447.0 * pow(z, 2) + 573100.0 / 2075003.0 * pow(z, 3) + 367895.0 / 366177.0 * pow(z, 4) + 143705474.0 / 267675387.0 * pow(z, 5) + 821920627.0 / 522604327.0 * pow(z, 6) + 1304439664.0 / 1567812981.0 * pow(z, 7) + 124705963250.0 / 58009080297.0 * pow(z, 8) + 150778688940.0 / 135354520693.0 * pow(z, 9) + 63521874746.0 / 23886091887.0 * pow(z, 10) + 14631063368.0 / 10974690867.0 * pow(z, 11) + 24094504214.0 / 7962030629.0 * pow(z, 12) + 83765938300.0 / 58009080297.0 * pow(z, 13) + 4988238530.0 / 1567812981.0 * pow(z, 14) + 737291984.0 / 522604327.0 * pow(z, 15) + 821920627.0 / 267675387.0 * pow(z, 16) + 7563446.0 / 6225009.0 * pow(z, 17) + 5518425.0 / 2075003.0 * pow(z, 18) + 114620.0 / 132447.0 * pow(z, 19) + 5269.0 / 2703.0 * pow(z, 20) + 22.0 / 53.0 * pow(z, 21) + 1.0 * pow(z, 22));
                result += c_f[23]() *  2.3863739298863123 * (1.0 / 55.0 * pow(z, 0) + 47.0 / 265.0 * pow(z, 1) + 563.0 / 4505.0 * pow(z, 2) + 3361.0 / 6307.0 * pow(z, 3) + 96205.0 / 296429.0 * pow(z, 4) + 1521841.0 / 1482145.0 * pow(z, 5) + 5396241.0 / 9104605.0 * pow(z, 6) + 4176629679.0 / 2613021635.0 * pow(z, 7) + 465044970.0 / 522604327.0 * pow(z, 8) + 42076326770.0 / 19336360099.0 * pow(z, 9) + 113018629394.0 / 96681800495.0 * pow(z, 10) + 2854014392006.0 / 1063499805445.0 * pow(z, 11) + 133567471102.0 / 96681800495.0 * pow(z, 12) + 8415265354.0 / 2762337157.0 * pow(z, 13) + 775074950.0 / 522604327.0 * pow(z, 14) + 8353259358.0 / 2613021635.0 * pow(z, 15) + 5396241.0 / 3748955.0 * pow(z, 16) + 4565523.0 / 1482145.0 * pow(z, 17) + 365579.0 / 296429.0 * pow(z, 18) + 16805.0 / 6307.0 * pow(z, 19) + 3941.0 / 4505.0 * pow(z, 20) + 517.0 / 265.0 * pow(z, 21) + 23.0 / 55.0 * pow(z, 22) + 1.0 * pow(z, 23));
                result += c_f[24]() *  2.416688998884546 * (3.0 / 19.0 * pow(z, 0) + 8.0 / 95.0 * pow(z, 1) + 2444.0 / 5035.0 * pow(z, 2) + 4504.0 / 17119.0 * pow(z, 3) + 6722.0 / 7049.0 * pow(z, 4) + 153928.0 / 296429.0 * pow(z, 5) + 6087364.0 / 4022965.0 * pow(z, 6) + 992908344.0 / 1210912465.0 * pow(z, 7) + 20883148395.0 / 9929482213.0 * pow(z, 8) + 11161079280.0 / 9929482213.0 * pow(z, 9) + 976170781064.0 / 367390841881.0 * pow(z, 10) + 2548056371792.0 / 1836954209405.0 * pow(z, 11) + 5708028784012.0 / 1836954209405.0 * pow(z, 12) + 82195366832.0 / 52484405983.0 * pow(z, 13) + 33661061416.0 / 9929482213.0 * pow(z, 14) + 16121558960.0 / 9929482213.0 * pow(z, 15) + 4176629679.0 / 1210912465.0 * pow(z, 16) + 43169928.0 / 28160755.0 * pow(z, 17) + 18262092.0 / 5632151.0 * pow(z, 18) + 153928.0 / 119833.0 * pow(z, 19) + 47054.0 / 17119.0 * pow(z, 20) + 4504.0 / 5035.0 * pow(z, 21) + 188.0 / 95.0 * pow(z, 22) + 8.0 / 19.0 * pow(z, 23) + 1.0 * pow(z, 24));
                result += c_f[25]() *  2.4170361993746483 * (1.0 / 59.0 * pow(z, 0) + 185.0 / 1121.0 * pow(z, 1) + 132.0 / 1121.0 * pow(z, 2) + 29740.0 / 59413.0 * pow(z, 3) + 312790.0 / 1010021.0 * pow(z, 4) + 6896094.0 / 7070147.0 * pow(z, 5) + 27259340.0 / 47470987.0 * pow(z, 6) + 511450252.0 / 332296909.0 * pow(z, 7) + 12551644395.0 / 14288767087.0 * pow(z, 8) + 65696174435.0 / 30833655293.0 * pow(z, 9) + 692164738936.0 / 585839450567.0 * pow(z, 10) + 58169443650600.0 / 21676059670979.0 * pow(z, 11) + 31208670943948.0 / 21676059670979.0 * pow(z, 12) + 9694907275100.0 / 3096579952997.0 * pow(z, 13) + 943861007640.0 / 585839450567.0 * pow(z, 14) + 105113879096.0 / 30833655293.0 * pow(z, 15) + 1394627155.0 / 840515711.0 * pow(z, 16) + 1150763067.0 / 332296909.0 * pow(z, 17) + 27259340.0 / 17489311.0 * pow(z, 18) + 22986980.0 / 7070147.0 * pow(z, 19) + 1313718.0 / 1010021.0 * pow(z, 20) + 163570.0 / 59413.0 * pow(z, 21) + 1012.0 / 1121.0 * pow(z, 22) + 2220.0 / 1121.0 * pow(z, 23) + 25.0 / 59.0 * pow(z, 24) + 1.0 * pow(z, 25));
                result += c_f[26]() *  2.443781079080322 * (9.0 / 61.0 * pow(z, 0) + 286.0 / 3599.0 * pow(z, 1) + 31265.0 / 68381.0 * pow(z, 2) + 17160.0 / 68381.0 * pow(z, 3) + 3286270.0 / 3624193.0 * pow(z, 4) + 1626508.0 / 3242699.0 * pow(z, 5) + 89649222.0 / 61611281.0 * pow(z, 6) + 16301085320.0 / 20270111449.0 * pow(z, 7) + 41555332975.0 / 20270111449.0 * pow(z, 8) + 979028262810.0 / 871614792307.0 * pow(z, 9) + 4953491552399.0 / 1880852972873.0 * pow(z, 10) + 50716798143856.0 / 35736206484587.0 * pow(z, 11) + 4159115221017900.0 / 1322239639929719.0 * pow(z, 12) + 312086709439480.0 / 188891377132817.0 * pow(z, 13) + 126033794576300.0 / 35736206484587.0 * pow(z, 14) + 63805004116464.0 / 35736206484587.0 * pow(z, 15) + 170810053531.0 / 45874462753.0 * pow(z, 16) + 36260306030.0 / 20270111449.0 * pow(z, 17) + 74799599355.0 / 20270111449.0 * pow(z, 18) + 708742840.0 / 431278967.0 * pow(z, 19) + 209181518.0 / 61611281.0 * pow(z, 20) + 4879524.0 / 3624193.0 * pow(z, 21) + 193310.0 / 68381.0 * pow(z, 22) + 62920.0 / 68381.0 * pow(z, 23) + 7215.0 / 3599.0 * pow(z, 24) + 26.0 / 61.0 * pow(z, 25) + 1.0 * pow(z, 26));
                result += c_f[27]() *  2.4440889958054246 * (1.0 / 63.0 * pow(z, 0) + 593.0 / 3843.0 * pow(z, 1) + 8411.0 / 75579.0 * pow(z, 2) + 2032615.0 / 4308003.0 * pow(z, 3) + 1274390.0 / 4308003.0 * pow(z, 4) + 70638178.0 / 76108053.0 * pow(z, 5) + 308015942.0 / 554501529.0 * pow(z, 6) + 682102538.0 / 460518219.0 * pow(z, 7) + 122418663835.0 / 141890780143.0 * pow(z, 8) + 43512234155.0 / 20934705267.0 * pow(z, 9) + 1224984369323.0 / 1036070413497.0 * pow(z, 10) + 1997711464112689.0 / 750460336176327.0 * pow(z, 11) + 3321192077639228.0 / 2251381008528981.0 * pow(z, 12) + 37744123698600580.0 / 11900156759367471.0 * pow(z, 13) + 1277381568322780.0 / 750460336176327.0 * pow(z, 14) + 7990845856450756.0 / 2251381008528981.0 * pow(z, 15) + 111362215393.0 / 60945318441.0 * pow(z, 16) + 8702446831.0 / 2326078363.0 * pow(z, 17) + 122418663835.0 / 67211422173.0 * pow(z, 18) + 1705256345.0 / 460518219.0 * pow(z, 19) + 308015942.0 / 184833843.0 * pow(z, 20) + 777019958.0 / 228324159.0 * pow(z, 21) + 5862194.0 / 4308003.0 * pow(z, 22) + 4065230.0 / 1436001.0 * pow(z, 23) + 210275.0 / 226737.0 * pow(z, 24) + 7709.0 / 3843.0 * pow(z, 25) + 3.0 / 7.0 * pow(z, 26) + 1.0 * pow(z, 27));
                result += c_f[28]() *  2.467859887040803 * (9.0 / 65.0 * pow(z, 0) + 44.0 / 585.0 * pow(z, 1) + 1186.0 / 2745.0 * pow(z, 2) + 2588.0 / 10797.0 * pow(z, 3) + 531607.0 / 615429.0 * pow(z, 4) + 78424.0 / 161955.0 * pow(z, 5) + 76071884.0 / 54362895.0 * pow(z, 6) + 2179805128.0 / 2772507645.0 * pow(z, 7) + 131173565.0 / 65788317.0 * pow(z, 8) + 22600368708.0 / 20270111449.0 * pow(z, 9) + 38826301246.0 / 14953360905.0 * pow(z, 10) + 1062224208364.0 / 740050295355.0 * pow(z, 11) + 1690371238864583.0 / 536043097268805.0 * pow(z, 12) + 1021905254658224.0 / 597305165528097.0 * pow(z, 13) + 15097649479440232.0 / 4181136158696679.0 * pow(z, 14) + 1021905254658224.0 / 536043097268805.0 * pow(z, 15) + 153670112624053.0 / 39222665653815.0 * pow(z, 16) + 34265297044.0 / 17210471985.0 * pow(z, 17) + 1338837974.0 / 332296909.0 * pow(z, 18) + 7533456236.0 / 3881510703.0 * pow(z, 19) + 183642991.0 / 46991655.0 * pow(z, 20) + 94774136.0 / 54362895.0 * pow(z, 21) + 10867412.0 / 3077145.0 * pow(z, 22) + 862664.0 / 615429.0 * pow(z, 23) + 31271.0 / 10797.0 * pow(z, 24) + 2588.0 / 2745.0 * pow(z, 25) + 1186.0 / 585.0 * pow(z, 26) + 28.0 / 65.0 * pow(z, 27) + 1.0 * pow(z, 28));
                result += c_f[29]() *  2.468134811554221 * (1.0 / 67.0 * pow(z, 0) + 631.0 / 4355.0 * pow(z, 1) + 106.0 / 1005.0 * pow(z, 2) + 5470.0 / 12261.0 * pow(z, 3) + 204667.0 / 723399.0 * pow(z, 4) + 12160111.0 / 13744581.0 * pow(z, 5) + 36900388.0 / 68722905.0 * pow(z, 6) + 5191590364.0 / 3642313965.0 * pow(z, 7) + 10793387.0 / 12780049.0 * pow(z, 8) + 3438317377.0 / 1699746517.0 * pow(z, 9) + 93876107050.0 / 79888086299.0 * pow(z, 10) + 3147042234034.0 / 1198321294485.0 * pow(z, 11) + 76973103461267.0 / 51527815662855.0 * pow(z, 12) + 192061812199303.0 / 60361155490773.0 * pow(z, 13) + 9693553214002888.0 / 5492865149660343.0 * pow(z, 14) + 1536494497594424.0 / 422528088435411.0 * pow(z, 15) + 100657135295503.0 / 51527815662855.0 * pow(z, 16) + 1573521117017.0 / 399440431495.0 * pow(z, 17) + 8534191550.0 / 4204636121.0 * pow(z, 18) + 6876634754.0 / 1699746517.0 * pow(z, 19) + 75553709.0 / 38340147.0 * pow(z, 20) + 14276873501.0 / 3642313965.0 * pow(z, 21) + 121244132.0 / 68722905.0 * pow(z, 22) + 48640444.0 / 13744581.0 * pow(z, 23) + 1023335.0 / 723399.0 * pow(z, 24) + 35555.0 / 12261.0 * pow(z, 25) + 318.0 / 335.0 * pow(z, 26) + 8834.0 / 4355.0 * pow(z, 27) + 29.0 / 67.0 * pow(z, 28) + 1.0 * pow(z, 29));
                result += c_f[30]() *  2.4894020435851436 * (3.0 / 23.0 * pow(z, 0) + 110.0 / 1541.0 * pow(z, 1) + 631.0 / 1541.0 * pow(z, 2) + 1060.0 / 4623.0 * pow(z, 3) + 232475.0 / 282003.0 * pow(z, 4) + 7777346.0 / 16638177.0 * pow(z, 5) + 425603885.0 / 316125363.0 * pow(z, 6) + 10542968.0 / 13744581.0 * pow(z, 7) + 32447439775.0 / 16754644239.0 * pow(z, 8) + 323801610.0 / 293941127.0 * pow(z, 9) + 99711203933.0 / 39094169891.0 * pow(z, 10) + 2645599380500.0 / 1837425984877.0 * pow(z, 11) + 17308732287187.0 / 5512277954631.0 * pow(z, 12) + 2574351286330.0 / 1472223304653.0 * pow(z, 13) + 35531435256871055.0 / 9718146034014453.0 * pow(z, 14) + 19387106428005776.0 / 9718146034014453.0 * pow(z, 15) + 960309060996515.0 / 237027952049133.0 * pow(z, 16) + 514870257266.0 / 239664258897.0 * pow(z, 17) + 7867605585085.0 / 1837425984877.0 * pow(z, 18) + 85341915500.0 / 39094169891.0 * pow(z, 19) + 24068221639.0 / 5584881413.0 * pow(z, 20) + 1834875790.0 / 881823381.0 * pow(z, 21) + 1297897591.0 / 316125363.0 * pow(z, 22) + 579863240.0 / 316125363.0 * pow(z, 23) + 60800555.0 / 16638177.0 * pow(z, 24) + 409334.0 / 282003.0 * pow(z, 25) + 13675.0 / 4623.0 * pow(z, 26) + 1484.0 / 1541.0 * pow(z, 27) + 3155.0 / 1541.0 * pow(z, 28) + 10.0 / 23.0 * pow(z, 29) + 1.0 * pow(z, 30));
                result += c_f[31]() *  2.489648995824587 * (1.0 / 71.0 * pow(z, 0) + 223.0 / 1633.0 * pow(z, 1) + 10965.0 / 109411.0 * pow(z, 2) + 46285.0 / 109411.0 * pow(z, 3) + 29645.0 / 109411.0 * pow(z, 4) + 5638353.0 / 6674071.0 * pow(z, 5) + 204330707.0 / 393770189.0 * pow(z, 6) + 10265913025.0 / 7481633591.0 * pow(z, 7) + 6171521445.0 / 7481633591.0 * pow(z, 8) + 779543621345.0 / 396526580323.0 * pow(z, 9) + 460876593529.0 / 396526580323.0 * pow(z, 10) + 1023548199249.0 / 396526580323.0 * pow(z, 11) + 27957880228655.0 / 18636749275181.0 * pow(z, 12) + 59083905157495.0 / 18636749275181.0 * pow(z, 13) + 10129239538344915.0 / 5609661531829481.0 * pow(z, 14) + 847373003221950227.0 / 229996122805008721.0 * pow(z, 15) + 11479804810124237.0 / 5609661531829481.0 * pow(z, 16) + 531755146417455.0 / 130457244926267.0 * pow(z, 17) + 2150606171435.0 / 980881540799.0 * pow(z, 18) + 1705913665415.0 / 396526580323.0 * pow(z, 19) + 879855314919.0 / 396526580323.0 * pow(z, 20) + 1714995966959.0 / 396526580323.0 * pow(z, 21) + 685724605.0 / 325288417.0 * pow(z, 22) + 30797739075.0 / 7481633591.0 * pow(z, 23) + 729752525.0 / 393770189.0 * pow(z, 24) + 24432863.0 / 6674071.0 * pow(z, 25) + 160083.0 / 109411.0 * pow(z, 26) + 323995.0 / 109411.0 * pow(z, 27) + 105995.0 / 109411.0 * pow(z, 28) + 3345.0 / 1633.0 * pow(z, 29) + 31.0 / 71.0 * pow(z, 30) + 1.0 * pow(z, 31));
                result += c_f[32]() *  2.508788609246518 * (9.0 / 73.0 * pow(z, 0) + 352.0 / 5183.0 * pow(z, 1) + 46384.0 / 119209.0 * pow(z, 2) + 1754400.0 / 7987003.0 * pow(z, 3) + 6294760.0 / 7987003.0 * pow(z, 4) + 3604832.0 / 7987003.0 * pow(z, 5) + 631495536.0 / 487207183.0 * pow(z, 6) + 934083232.0 / 1249792339.0 * pow(z, 7) + 1026591302500.0 / 546159252143.0 * pow(z, 8) + 592466058720.0 / 546159252143.0 * pow(z, 9) + 72341648060816.0 / 28946440363579.0 * pow(z, 10) + 41562689161888.0 / 28946440363579.0 * pow(z, 11) + 90072241533912.0 / 28946440363579.0 * pow(z, 12) + 2408678912007200.0 / 1360482697088213.0 * pow(z, 13) + 34977671853237040.0 / 9523378879617491.0 * pow(z, 14) + 842752729590296928.0 / 409505291823552113.0 * pow(z, 15) + 1694746006443900454.0 / 409505291823552113.0 * pow(z, 16) + 21609044348469152.0 / 9523378879617491.0 * pow(z, 17) + 42540411713396400.0 / 9523378879617491.0 * pow(z, 18) + 68819397485920.0 / 28946440363579.0 * pow(z, 19) + 133743631368536.0 / 28946440363579.0 * pow(z, 20) + 68377327330848.0 / 28946440363579.0 * pow(z, 21) + 2494539588304.0 / 546159252143.0 * pow(z, 22) + 1206875304800.0 / 546159252143.0 * pow(z, 23) + 123190956300.0 / 28745223797.0 * pow(z, 24) + 934083232.0 / 487207183.0 * pow(z, 25) + 30071216.0 / 7987003.0 * pow(z, 26) + 11952864.0 / 7987003.0 * pow(z, 27) + 24068200.0 / 7987003.0 * pow(z, 28) + 116960.0 / 119209.0 * pow(z, 29) + 10704.0 / 5183.0 * pow(z, 30) + 32.0 / 73.0 * pow(z, 31) + 1.0 * pow(z, 32));
                result += c_f[33]() *  2.509011642416648 * (1.0 / 75.0 * pow(z, 0) + 707.0 / 5475.0 * pow(z, 1) + 12368.0 / 129575.0 * pow(z, 2) + 719152.0 / 1788135.0 * pow(z, 3) + 6225928.0 / 23961009.0 * pow(z, 4) + 161353288.0 / 199675075.0 * pow(z, 5) + 300433616.0 / 599025225.0 * pow(z, 6) + 48296248432.0 / 36540538725.0 * pow(z, 7) + 23126460420.0 / 28745223797.0 * pow(z, 8) + 164634153668.0 / 86235671391.0 * pow(z, 9) + 34788357296.0 / 30364673025.0 * pow(z, 10) + 96385981261264.0 / 38087421531025.0 * pow(z, 11) + 171102385184744.0 / 114262264593075.0 * pow(z, 12) + 71834079079256.0 / 22852452918615.0 * pow(z, 13) + 917260940713808.0 / 501230467348289.0 * pow(z, 14) + 19886725062716144.0 / 5370326435874525.0 * pow(z, 15) + 3415852669774535266.0 / 1616468257198232025.0 * pow(z, 16) + 7457521898518554.0 / 1790108811958175.0 * pow(z, 17) + 17427957873562352.0 / 7518457010224335.0 * pow(z, 18) + 143668158158512.0 / 31993434086061.0 * pow(z, 19) + 92132053561016.0 / 38087421531025.0 * pow(z, 20) + 530122896936952.0 / 114262264593075.0 * pow(z, 21) + 3162577936.0 / 1320203175.0 * pow(z, 22) + 658536614672.0 / 143726118985.0 * pow(z, 23) + 192720503500.0 / 86235671391.0 * pow(z, 24) + 156962807404.0 / 36540538725.0 * pow(z, 25) + 386271792.0 / 199675075.0 * pow(z, 26) + 2258946032.0 / 599025225.0 * pow(z, 27) + 180551912.0 / 119805045.0 * pow(z, 28) + 359576.0 / 119209.0 * pow(z, 29) + 383408.0 / 388725.0 * pow(z, 30) + 11312.0 / 5475.0 * pow(z, 31) + 11.0 / 25.0 * pow(z, 32) + 1.0 * pow(z, 33));
                result += c_f[34]() *  2.526327908322109 * (9.0 / 77.0 * pow(z, 0) + 34.0 / 525.0 * pow(z, 1) + 22321.0 / 60225.0 * pow(z, 2) + 420512.0 / 1995455.0 * pow(z, 3) + 14845352.0 / 19669485.0 * pow(z, 4) + 4021949488.0 / 9224988465.0 * pow(z, 5) + 2743005896.0 / 2196425825.0 * pow(z, 6) + 1459248992.0 / 2005432275.0 * pow(z, 7) + 205259055836.0 / 112544859273.0 * pow(z, 8) + 336985566120.0 / 316197461767.0 * pow(z, 9) + 81164637758324.0 / 33200733485535.0 * pow(z, 10) + 476193877792.0 / 334011403275.0 * pow(z, 11) + 819280840720744.0 / 266611950717175.0 * pow(z, 12) + 447498545867792.0 / 251376982104765.0 * pow(z, 13) + 45183635740852024.0 / 12317472123133485.0 * pow(z, 14) + 405429335795503136.0 / 192973729929091265.0 * pow(z, 15) + 1732630921089144046.0 / 413515135562338425.0 * pow(z, 16) + 6831705339549070532.0 / 2894605948936368975.0 * pow(z, 17) + 126777872274815418.0 / 27567675704155895.0 * pow(z, 18) + 31186871984269472.0 / 12317472123133485.0 * pow(z, 19) + 1221179344347352.0 / 251376982104765.0 * pow(z, 20) + 7607475279752464.0 / 2932731457888925.0 * pow(z, 21) + 819280840720744.0 / 166003667427675.0 * pow(z, 22) + 15361092832.0 / 6072934605.0 * pow(z, 23) + 53176831634764.0 / 11066911161845.0 * pow(z, 24) + 37442840680.0 / 16077837039.0 * pow(z, 25) + 205259055836.0 / 46124942325.0 * pow(z, 26) + 4377746976.0 / 2196425825.0 * pow(z, 27) + 35659076648.0 / 9224988465.0 * pow(z, 28) + 211681552.0 / 137686395.0 * pow(z, 29) + 873256.0 / 285065.0 * pow(z, 30) + 420512.0 / 421575.0 * pow(z, 31) + 1717.0 / 825.0 * pow(z, 32) + 34.0 / 77.0 * pow(z, 33) + 1.0 * pow(z, 34));
                result += c_f[35]() *  2.5265303303334368 * (1.0 / 79.0 * pow(z, 0) + 745.0 / 6083.0 * pow(z, 1) + 2771.0 / 30415.0 * pow(z, 2) + 170187.0 / 444059.0 * pow(z, 3) + 7866648.0 / 31528189.0 * pow(z, 4) + 2807053736.0 / 3625741735.0 * pow(z, 5) + 305968040.0 / 630973237.0 * pow(z, 6) + 44215042552.0 / 34703528035.0 * pow(z, 7) + 38089526324.0 / 48584939249.0 * pow(z, 8) + 5492521765268.0 / 2963681294189.0 * pow(z, 9) + 89540174723324.0 / 79480543798705.0 * pow(z, 10) + 18828943412812.0 / 7602486798137.0 * pow(z, 11) + 1301072070964536.0 / 874285981785755.0 * pow(z, 12) + 4111785862091880.0 / 1323918772418429.0 * pow(z, 13) + 2438237631193528.0 / 1323918772418429.0 * pow(z, 14) + 1200231365170526456.0 / 324360099242515105.0 * pow(z, 15) + 6583272526753691134.0 / 3048984932879641987.0 * pow(z, 16) + 5848343058799065018.0 / 1385902242218019085.0 * pow(z, 17) + 7357775176960007738.0 / 3048984932879641987.0 * pow(z, 18) + 300057841292631614.0 / 64872019848503021.0 * pow(z, 19) + 17067663418354696.0 / 6619593862092145.0 * pow(z, 20) + 587397980298840.0 / 120356252038039.0 * pow(z, 21) + 100082466997272.0 / 38012433990685.0 * pow(z, 22) + 37657886825624.0 / 7602486798137.0 * pow(z, 23) + 447700873616620.0 / 174857196357151.0 * pow(z, 24) + 71402782948484.0 / 14818406470945.0 * pow(z, 25) + 114268578972.0 / 48584939249.0 * pow(z, 26) + 154752648932.0 / 34703528035.0 * pow(z, 27) + 1267581880.0 / 630973237.0 * pow(z, 28) + 2807053736.0 / 725148347.0 * pow(z, 29) + 243866088.0 / 157640945.0 * pow(z, 30) + 1361496.0 / 444059.0 * pow(z, 31) + 2771.0 / 2765.0 * pow(z, 32) + 12665.0 / 6083.0 * pow(z, 33) + 35.0 / 79.0 * pow(z, 34) + 1.0 * pow(z, 35));
                result += c_f[36]() *  2.5422721046282533 * (1.0 / 9.0 * pow(z, 0) + 44.0 / 711.0 * pow(z, 1) + 19370.0 / 54747.0 * pow(z, 2) + 11084.0 / 54747.0 * pow(z, 3) + 964393.0 / 1332177.0 * pow(z, 4) + 66429472.0 / 157640945.0 * pow(z, 5) + 5614107472.0 / 4661667945.0 * pow(z, 6) + 174838880.0 / 246902571.0 * pow(z, 7) + 110537606380.0 / 62466350463.0 * pow(z, 8) + 152358105296.0 / 145754817747.0 * pow(z, 9) + 318566262385544.0 / 133365658238505.0 * pow(z, 10) + 11102981665692176.0 / 7868573836071795.0 * pow(z, 11) + 18828943412812.0 / 6220216471203.0 * pow(z, 12) + 44481096443232.0 / 24979599479593.0 * pow(z, 13) + 14489150180704720.0 / 3971756317255287.0 * pow(z, 14) + 126788356822063456.0 / 59576344758829305.0 * pow(z, 15) + 12302371492997896174.0 / 2919240893182635945.0 * pow(z, 16) + 66607227917743227944.0 / 27440864395916777883.0 * pow(z, 17) + 3898895372532710012.0 / 831541345330811451.0 * pow(z, 18) + 1549005300412633208.0 / 583848178636527189.0 * pow(z, 19) + 300057841292631614.0 / 59576344758829305.0 * pow(z, 20) + 165800158921159904.0 / 59576344758829305.0 * pow(z, 21) + 391598653532560.0 / 74938798438779.0 * pow(z, 22) + 44481096443232.0 / 15896108759741.0 * pow(z, 23) + 357749924843428.0 / 68422381183233.0 * pow(z, 24) + 358160698893296.0 / 133365658238505.0 * pow(z, 25) + 10985043530536.0 / 2186322266205.0 * pow(z, 26) + 152358105296.0 / 62466350463.0 * pow(z, 27) + 287397776588.0 / 62466350463.0 * pow(z, 28) + 174838880.0 / 84757599.0 * pow(z, 29) + 5614107472.0 / 1418768505.0 * pow(z, 30) + 3496288.0 / 2220295.0 * pow(z, 31) + 56729.0 / 18249.0 * pow(z, 32) + 55420.0 / 54747.0 * pow(z, 33) + 1490.0 / 711.0 * pow(z, 34) + 4.0 / 9.0 * pow(z, 35) + 1.0 * pow(z, 36));
                result += c_f[37]() *  2.5424566414923317 * (1.0 / 83.0 * pow(z, 0) + 29.0 / 249.0 * pow(z, 1) + 1714.0 / 19671.0 * pow(z, 2) + 237590.0 / 649143.0 * pow(z, 3) + 155737.0 / 649143.0 * pow(z, 4) + 58673137.0 / 78979065.0 * pow(z, 5) + 7890989008.0 / 16822540845.0 * pow(z, 6) + 475587058576.0 / 386918439435.0 * pow(z, 7) + 439871046412.0 / 576078565381.0 * pow(z, 8) + 9326979434836.0 / 5184707088429.0 * pow(z, 9) + 2606115384728.0 / 2356685040195.0 * pow(z, 10) + 17973950688728.0 / 7424111089065.0 * pow(z, 11) + 12502450904018764.0 / 8481709459661805.0 * pow(z, 12) + 5192111834721028.0 / 1696341891932361.0 * pow(z, 13) + 11467391667897328.0 / 6219920270418657.0 * pow(z, 14) + 249380714162993648.0 / 67737487876477155.0 * pow(z, 15) + 10823491457523898462.0 / 4944836614982832315.0 * pow(z, 16) + 16330823181283945214.0 / 3845984033875536245.0 * pow(z, 17) + 808151923495181049964.0 / 325370249265870366327.0 * pow(z, 18) + 32661646362567890428.0 / 6922771260975965241.0 * pow(z, 19) + 4456731776627487602.0 / 1648278871660944105.0 * pow(z, 20) + 31172589270374206.0 / 6157953443316105.0 * pow(z, 21) + 11467391667897328.0 / 4056469741577385.0 * pow(z, 22) + 2966921048412016.0 / 565447297310787.0 * pow(z, 23) + 4808634963084140.0 / 1696341891932361.0 * pow(z, 24) + 116830679476732.0 / 22272333267195.0 * pow(z, 25) + 7818346154184.0 / 2880392826905.0 * pow(z, 26) + 130577712087704.0 / 25923535442145.0 * pow(z, 27) + 12756260345948.0 / 5184707088429.0 * pow(z, 28) + 118896764644.0 / 25794562629.0 * pow(z, 29) + 34945808464.0 / 16822540845.0 * pow(z, 30) + 938770192.0 / 236937195.0 * pow(z, 31) + 155737.0 / 98355.0 * pow(z, 32) + 2019515.0 / 649143.0 * pow(z, 33) + 59990.0 / 59013.0 * pow(z, 34) + 174.0 / 83.0 * pow(z, 35) + 37.0 / 83.0 * pow(z, 36) + 1.0 * pow(z, 37));
                result += c_f[38]() *  2.5568294389699044 * (9.0 / 85.0 * pow(z, 0) + 418.0 / 7055.0 * pow(z, 1) + 7163.0 / 21165.0 * pow(z, 2) + 65132.0 / 334407.0 * pow(z, 3) + 451421.0 / 649143.0 * pow(z, 4) + 6614242.0 / 16228575.0 * pow(z, 5) + 459031013.0 / 394895325.0 * pow(z, 6) + 57955667168.0 / 84112704225.0 * pow(z, 7) + 132884619308.0 / 77383687887.0 * pow(z, 8) + 2949723487704.0 / 2880392826905.0 * pow(z, 9) + 302303862858508.0 / 129617677210725.0 * pow(z, 10) + 180588466071152.0 / 129617677210725.0 * pow(z, 11) + 10044266561348.0 / 3374595949575.0 * pow(z, 12) + 15048198825651544.0 / 8481709459661805.0 * pow(z, 13) + 30672727897553636.0 / 8481709459661805.0 * pow(z, 14) + 333228910820075296.0 / 155498006760466425.0 * pow(z, 15) + 1428438061271853322.0 / 338687439382385775.0 * pow(z, 16) + 1040328531858473534524.0 / 420311112273540746775.0 * pow(z, 17) + 310285640444394959066.0 / 65381728575884116165.0 * pow(z, 18) + 1616303846990362099928.0 / 588435557182957045485.0 * pow(z, 19) + 2171999483110764713462.0 / 420311112273540746775.0 * pow(z, 20) + 24193686787406361268.0 / 8241394358304720525.0 * pow(z, 21) + 34839952713947642.0 / 6390329044950675.0 * pow(z, 22) + 25632993140005792.0 / 8481709459661805.0 * pow(z, 23) + 15750860271716732.0 / 2827236486553935.0 * pow(z, 24) + 2149742689378792.0 / 718788937259475.0 * pow(z, 25) + 10044266561348.0 / 1825601087475.0 * pow(z, 26) + 40778040725744.0 / 14401964134525.0 * pow(z, 27) + 135515524729676.0 / 25923535442145.0 * pow(z, 28) + 983241162568.0 / 386918439435.0 * pow(z, 29) + 132884619308.0 / 28037568075.0 * pow(z, 30) + 2519811616.0 / 1184685975.0 * pow(z, 31) + 65575859.0 / 16228575.0 * pow(z, 32) + 348118.0 / 216381.0 * pow(z, 33) + 3159947.0 / 1003221.0 * pow(z, 34) + 65132.0 / 63495.0 * pow(z, 35) + 14877.0 / 7055.0 * pow(z, 36) + 38.0 / 85.0 * pow(z, 37) + 1.0 * pow(z, 38));
                result += c_f[39]() *  2.5569983571109383 * (1.0 / 87.0 * pow(z, 0) + 821.0 / 7395.0 * pow(z, 1) + 589.0 / 7055.0 * pow(z, 2) + 386935.0 / 1104813.0 * pow(z, 3) + 1185847.0 / 5134131.0 * pow(z, 4) + 1221757.0 / 1711377.0 * pow(z, 5) + 58274083.0 / 128353275.0 * pow(z, 6) + 11120627819.0 / 9369789075.0 * pow(z, 7) + 379017548.0 / 509773965.0 * pow(z, 8) + 5344368232868.0 / 3060173111895.0 * pow(z, 9) + 222286274145628.0 / 205031598496965.0 * pow(z, 10) + 808104195006724.0 / 341719330828275.0 * pow(z, 11) + 1493121770367812.0 / 1025157992484825.0 * pow(z, 12) + 37656009267543716.0 / 12506927508314865.0 * pow(z, 13) + 452148626215611508.0 / 245969574330192345.0 * pow(z, 14) + 2694160320227172124.0 / 737908722990577035.0 * pow(z, 15) + 3552083737635722.0 / 1610451163226925.0 * pow(z, 16) + 92369607954543025574.0 / 21727312399166990475.0 * pow(z, 17) + 337020297290174408062.0 / 132971151882901981707.0 * pow(z, 18) + 22233632755904692277146.0 / 4653990315901569359745.0 * pow(z, 19) + 124165372685853729286.0 / 44323717294300660569.0 * pow(z, 20) + 1016065687499973281314.0 / 195545811592502914275.0 * pow(z, 21) + 208946102213866.0 / 70019615792475.0 * pow(z, 22) + 1347080160113586062.0 / 245969574330192345.0 * pow(z, 23) + 452148626215611508.0 / 147581744598115407.0 * pow(z, 24) + 69932588639724044.0 / 12506927508314865.0 * pow(z, 25) + 114855520797524.0 / 37968814536475.0 * pow(z, 26) + 5656729365047068.0 / 1025157992484825.0 * pow(z, 27) + 20207843104148.0 / 7070055120585.0 * pow(z, 28) + 5344368232868.0 / 1020057703965.0 * pow(z, 29) + 11749543988.0 / 4587965685.0 * pow(z, 30) + 44482511276.0 / 9369789075.0 * pow(z, 31) + 91573559.0 / 42784425.0 * pow(z, 32) + 20769869.0 / 5134131.0 * pow(z, 33) + 8300929.0 / 5134131.0 * pow(z, 34) + 386935.0 / 122757.0 * pow(z, 35) + 21793.0 / 21165.0 * pow(z, 36) + 15599.0 / 7395.0 * pow(z, 37) + 13.0 / 29.0 * pow(z, 38) + 1.0 * pow(z, 39));
                result += c_f[40]() *  2.570173398494093 * (9.0 / 89.0 * pow(z, 0) + 440.0 / 7743.0 * pow(z, 1) + 42692.0 / 131631.0 * pow(z, 2) + 23560.0 / 125579.0 * pow(z, 3) + 3869350.0 / 5784021.0 * pow(z, 4) + 180248744.0 / 456937659.0 * pow(z, 5) + 171045980.0 / 152312553.0 * pow(z, 6) + 1531775896.0 / 2284688295.0 * pow(z, 7) + 55603139095.0 / 33356449107.0 * pow(z, 8) + 3032140384.0 / 3024658859.0 * pow(z, 9) + 21377472931472.0 / 9391565757195.0 * pow(z, 10) + 5011545089828704.0 / 3649562453245977.0 * pow(z, 11) + 17778292290147928.0 / 6082604088743295.0 * pow(z, 12) + 6431909164661344.0 / 3649562453245977.0 * pow(z, 13) + 796155624513781424.0 / 222623309648004597.0 * pow(z, 14) + 47023457126423596832.0 / 21891292115387118705.0 * pow(z, 15) + 55230286564657028542.0 / 13134775269232271223.0 * pow(z, 16) + 71877459161569904.0 / 28666030705439265.0 * pow(z, 17) + 369478431818172102296.0 / 77349232141034486091.0 * pow(z, 18) + 374688444060040146160.0 / 132971151882901981707.0 * pow(z, 19) + 311270858582665691880044.0 / 59172162587891381859615.0 * pow(z, 20) + 7972094554468939280.0 / 2607277487900038857.0 * pow(z, 21) + 369478431818172102296.0 / 65673876346161356115.0 * pow(z, 22) + 18387256994820208.0 / 5733206141087853.0 * pow(z, 23) + 25594523042158135178.0 / 4378258423077423741.0 * pow(z, 24) + 3617189009724892064.0 / 1113116548240022985.0 * pow(z, 25) + 21517719581453552.0 / 3649562453245977.0 * pow(z, 26) + 6431909164661344.0 / 2027534696247765.0 * pow(z, 27) + 21010709070174824.0 / 3649562453245977.0 * pow(z, 28) + 161662744833184.0 / 54471081391731.0 * pow(z, 29) + 21377472931472.0 / 3947179810995.0 * pow(z, 30) + 3032140384.0 / 1150222383.0 * pow(z, 31) + 11120627819.0 / 2284688295.0 * pow(z, 32) + 332994760.0 / 152312553.0 * pow(z, 33) + 1881505780.0 / 456937659.0 * pow(z, 34) + 9486776.0 / 5784021.0 * pow(z, 35) + 11608050.0 / 3641791.0 * pow(z, 36) + 4712.0 / 4539.0 * pow(z, 37) + 16420.0 / 7743.0 * pow(z, 38) + 40.0 / 89.0 * pow(z, 39) + 1.0 * pow(z, 40));
                result += c_f[41]() *  2.570328597515912 * (1.0 / 91.0 * pow(z, 0) + 859.0 / 8099.0 * pow(z, 1) + 18820.0 / 234871.0 * pow(z, 2) + 1340540.0 / 3992807.0 * pow(z, 3) + 4340170.0 / 19494293.0 * pow(z, 4) + 3090578.0 / 4498683.0 * pow(z, 5) + 290224468.0 / 660021063.0 * pow(z, 6) + 757056140.0 / 660021063.0 * pow(z, 7) + 477816617.0 / 660021063.0 * pow(z, 8) + 18856475197.0 / 11118816369.0 * pow(z, 9) + 612294228304.0 / 576554354415.0 * pow(z, 10) + 118475082381616.0 / 51313337542935.0 * pow(z, 11) + 75977275401272.0 / 52892209467333.0 * pow(z, 12) + 156437680448872.0 / 52892209467333.0 * pow(z, 13) + 8794513819791824.0 / 4813191061527303.0 * pow(z, 14) + 1767478251997562512.0 / 489341091255275805.0 * pow(z, 15) + 63871119039352724458.0 / 28871124384061272495.0 * pow(z, 16) + 24483484459704940550.0 / 5774224876812254499.0 * pow(z, 17) + 3423536899728300392.0 / 1332513433110520269.0 * pow(z, 18) + 4416396882682987805336.0 / 918101755413148465341.0 * pow(z, 19) + 32058549484045530870868.0 / 11148378458588231364855.0 * pow(z, 20) + 24290182854756432929348.0 / 4590508777065742326705.0 * pow(z, 21) + 4144281510197416264.0 / 1332513433110520269.0 * pow(z, 22) + 97933937838819762200.0 / 17322674630436763497.0 * pow(z, 23) + 18785623246868448370.0 / 5774224876812254499.0 * pow(z, 24) + 220934781499695314.0 / 37641622404251985.0 * pow(z, 25) + 8794513819791824.0 / 2673995034181835.0 * pow(z, 26) + 312875360897744.0 / 52892209467333.0 * pow(z, 27) + 75977275401272.0 / 23710300795701.0 * pow(z, 28) + 59237541190808.0 / 10262667508587.0 * pow(z, 29) + 1725556461584.0 / 576554354415.0 * pow(z, 30) + 301703603152.0 / 55594081845.0 * pow(z, 31) + 5255982787.0 / 1980063189.0 * pow(z, 32) + 3217488595.0 / 660021063.0 * pow(z, 33) + 1451122340.0 / 660021063.0 * pow(z, 34) + 6181156.0 / 1499561.0 * pow(z, 35) + 32117258.0 / 19494293.0 * pow(z, 36) + 12735130.0 / 3992807.0 * pow(z, 37) + 18820.0 / 18067.0 * pow(z, 38) + 17180.0 / 8099.0 * pow(z, 39) + 41.0 / 91.0 * pow(z, 40) + 1.0 * pow(z, 41));
                result += c_f[42]() *  2.582449679876761 * (3.0 / 31.0 * pow(z, 0) + 22.0 / 403.0 * pow(z, 1) + 859.0 / 2759.0 * pow(z, 2) + 188200.0 / 1040143.0 * pow(z, 3) + 670270.0 / 1040143.0 * pow(z, 4) + 32985292.0 / 86331869.0 * pow(z, 5) + 21634046.0 / 19922739.0 * pow(z, 6) + 13350325528.0 / 20460652953.0 * pow(z, 7) + 33121206125.0 / 20460652953.0 * pow(z, 8) + 6689432638.0 / 6820217651.0 * pow(z, 9) + 131995326379.0 / 59428156455.0 * pow(z, 10) + 779283563296.0 / 576554354415.0 * pow(z, 11) + 4561290671692216.0 / 1590713463830985.0 * pow(z, 12) + 37228864946623280.0 / 21315560415335199.0 * pow(z, 13) + 5788194176608264.0 / 1639658493487323.0 * pow(z, 14) + 17589027639583648.0 / 8198292467436615.0 * pow(z, 15) + 9058326041487507874.0 / 2167081975559078565.0 * pow(z, 16) + 323112719846137311964.0 / 127857836557985635335.0 * pow(z, 17) + 122417422298524702750.0 / 25571567311597127067.0 * pow(z, 18) + 16937498346024222992.0 / 5901130918060875477.0 * pow(z, 19) + 108201723625733201230732.0 / 20329396012719716018265.0 * pow(z, 20) + 64117098968091061741736.0 / 20329396012719716018265.0 * pow(z, 21) + 2208198441341493902668.0 / 383573509673956906005.0 * pow(z, 22) + 19820476787900686480.0 / 5901130918060875477.0 * pow(z, 23) + 465186204734393870450.0 / 76714701934791381201.0 * pow(z, 24) + 7514249298747379348.0 / 2167081975559078565.0 * pow(z, 25) + 220934781499695314.0 / 35525934025558665.0 * pow(z, 26) + 123123193477085536.0 / 35525934025558665.0 * pow(z, 27) + 10168449229176680.0 / 1639658493487323.0 * pow(z, 28) + 1063681855617808.0 / 318142692766197.0 * pow(z, 29) + 9537244131720088.0 / 1590713463830985.0 * pow(z, 30) + 779283563296.0 / 251734999815.0 * pow(z, 31) + 131995326379.0 / 23608445715.0 * pow(z, 32) + 167235815950.0 / 61381958859.0 * pow(z, 33) + 102013314865.0 / 20460652953.0 * pow(z, 34) + 580448936.0 / 258995607.0 * pow(z, 35) + 27815202.0 / 6640913.0 * pow(z, 36) + 1736068.0 / 1040143.0 * pow(z, 37) + 3351350.0 / 1040143.0 * pow(z, 38) + 37640.0 / 35867.0 * pow(z, 39) + 859.0 / 403.0 * pow(z, 40) + 14.0 / 31.0 * pow(z, 41) + 1.0 * pow(z, 42));
                result += c_f[43]() *  2.5825927637719244 * (1.0 / 95.0 * pow(z, 0) + 299.0 / 2945.0 * pow(z, 1) + 2949.0 / 38285.0 * pow(z, 2) + 219701.0 / 681473.0 * pow(z, 3) + 136970.0 / 637507.0 * pow(z, 4) + 3442722.0 / 5200715.0 * pow(z, 5) + 183957914.0 / 431659345.0 * pow(z, 6) + 478920978.0 / 431659345.0 * pow(z, 7) + 4808049965.0 / 6820217651.0 * pow(z, 8) + 101123995385.0 / 61381958859.0 * pow(z, 9) + 1595434935283.0 / 1534548971475.0 * pow(z, 10) + 84153839498197.0 / 37340691639225.0 * pow(z, 11) + 865552074855944.0 / 611812870704225.0 * pow(z, 12) + 4617273927251048.0 / 1590713463830985.0 * pow(z, 13) + 12873445648801784.0 / 7105186805111733.0 * pow(z, 14) + 1900603611630909832.0 / 532889010383379975.0 * pow(z, 15) + 1178171235757099538.0 / 532889010383379975.0 * pow(z, 16) + 15229039004211589954.0 / 3611803292598464275.0 * pow(z, 17) + 993821643998116876442.0 / 383573509673956906005.0 * pow(z, 18) + 7029326310664449841598.0 / 1457579336761036242819.0 * pow(z, 19) + 428904621397755304732.0 / 146343306903718498275.0 * pow(z, 20) + 10343280843412745178661276.0 / 1931292621208373021735175.0 * pow(z, 21) + 1409258041735481715548.0 / 439029920711155494825.0 * pow(z, 22) + 14058652621328899683196.0 / 2429298894601727071365.0 * pow(z, 23) + 261532011578451809590.0 / 76714701934791381201.0 * pow(z, 24) + 15229039004211589954.0 / 2500479202568167575.0 * pow(z, 25) + 207912571015958742.0 / 59209890042597775.0 * pow(z, 26) + 3326056320354092206.0 / 532889010383379975.0 * pow(z, 27) + 12873445648801784.0 / 3675096623333655.0 * pow(z, 28) + 659610561035864.0 / 106047564255399.0 * pow(z, 29) + 865552074855944.0 / 256566687714675.0 * pow(z, 30) + 673230715985576.0 / 112022074917675.0 * pow(z, 31) + 1595434935283.0 / 511516323825.0 * pow(z, 32) + 343821584309.0 / 61381958859.0 * pow(z, 33) + 168281748775.0 / 61381958859.0 * pow(z, 34) + 2155144401.0 / 431659345.0 * pow(z, 35) + 972348974.0 / 431659345.0 * pow(z, 36) + 21803906.0 / 5200715.0 * pow(z, 37) + 82182.0 / 49039.0 * pow(z, 38) + 2197010.0 / 681473.0 * pow(z, 39) + 40303.0 / 38285.0 * pow(z, 40) + 6279.0 / 2945.0 * pow(z, 41) + 43.0 / 95.0 * pow(z, 42) + 1.0 * pow(z, 43));
                result += c_f[44]() *  2.593781542064307 * (9.0 / 97.0 * pow(z, 0) + 484.0 / 9215.0 * pow(z, 1) + 85514.0 / 285665.0 * pow(z, 2) + 129756.0 / 742729.0 * pow(z, 3) + 41084087.0 / 66102881.0 * pow(z, 4) + 1205336.0 / 3254641.0 * pow(z, 5) + 530179188.0 / 504469355.0 * pow(z, 6) + 26595058424.0 / 41870956465.0 * pow(z, 7) + 13170326895.0 / 8374191293.0 * pow(z, 8) + 634662595380.0 / 661561112147.0 * pow(z, 9) + 444945579694.0 / 205312069287.0 * pow(z, 10) + 6381739741132.0 / 4801653233325.0 * pow(z, 11) + 10182614579281837.0 / 3622047089004825.0 * pow(z, 12) + 266590039055630752.0 / 154299205991605545.0 * pow(z, 13) + 536922996683193296.0 / 154299205991605545.0 * pow(z, 14) + 566431608547278496.0 / 265078123113783885.0 * pow(z, 15) + 214293057211385083558.0 / 51690234007187857575.0 * pow(z, 16) + 131123528120731313288.0 / 51690234007187857575.0 * pow(z, 17) + 335038858092654978988.0 / 70068983876410206935.0 * pow(z, 18) + 108169639988847668446424.0 / 37206630438373819882485.0 * pow(z, 19) + 3788806881448138464621322.0 / 706925978329102577767215.0 * pow(z, 20) + 45831522400788709705648.0 / 14195300769660694332675.0 * pow(z, 21) + 20686561686825490357322552.0 / 3534629891645512888836075.0 * pow(z, 22) + 29655690965216223927184.0 / 8517180461796416599605.0 * pow(z, 23) + 77322589417308948257578.0 / 12402210146124606627495.0 * pow(z, 24) + 2301481701890375924392.0 / 630620854887691862415.0 * pow(z, 25) + 335038858092654978988.0 / 51690234007187857575.0 * pow(z, 26) + 21345690624305097512.0 / 5743359334131984175.0 * pow(z, 27) + 5226659931985002038.0 / 795234369341351655.0 * pow(z, 28) + 566431608547278496.0 / 154299205991605545.0 * pow(z, 29) + 333762943884147184.0 / 51433068663868515.0 * pow(z, 30) + 38084291293661536.0 / 10866141267014475.0 * pow(z, 31) + 925692234480167.0 / 148851250233075.0 * pow(z, 32) + 6381739741132.0 / 1984683336441.0 * pow(z, 33) + 34260809636438.0 / 5954050009323.0 * pow(z, 34) + 211554198460.0 / 75367721637.0 * pow(z, 35) + 213359295699.0 / 41870956465.0 * pow(z, 36) + 1156306888.0 / 504469355.0 * pow(z, 37) + 429192676.0 / 100893871.0 * pow(z, 38) + 3616008.0 / 2132351.0 * pow(z, 39) + 2416711.0 / 742729.0 * pow(z, 40) + 302764.0 / 285665.0 * pow(z, 41) + 19734.0 / 9215.0 * pow(z, 42) + 44.0 / 97.0 * pow(z, 43) + 1.0 * pow(z, 44));
                result += c_f[45]() *  2.5939138744815287 * (1.0 / 99.0 * pow(z, 0) + 85.0 / 873.0 * pow(z, 1) + 410.0 / 5529.0 * pow(z, 2) + 159430.0 / 514197.0 * pow(z, 3) + 1387505.0 / 6684561.0 * pow(z, 4) + 68719.0 / 107601.0 * pow(z, 5) + 375306260.0 / 908044839.0 * pow(z, 6) + 975346300.0 / 908044839.0 * pow(z, 7) + 5750040565.0 / 8374191293.0 * pow(z, 8) + 1086028678355.0 / 678309494733.0 * pow(z, 9) + 54522289283638.0 / 53586450083907.0 * pow(z, 10) + 432194659508290.0 / 196483650307659.0 * pow(z, 11) + 59148587701.0 / 42495202287.0 * pow(z, 12) + 11135712077601335.0 / 3911810856125211.0 * pow(z, 13) + 166022440594908400.0 / 92579523594963327.0 * pow(z, 14) + 75135396196060720.0 / 21364505444991537.0 * pow(z, 15) + 3153760866736849370.0 / 1431421864814432979.0 * pow(z, 16) + 279011063923699826.0 / 66697076138306913.0 * pow(z, 17) + 583407147515665060.0 / 224198605332381069.0 * pow(z, 18) + 5469474622953562949140.0 / 1135117538797845352347.0 * pow(z, 19) + 66307649255441491036174.0 / 22323978263024291929491.0 * pow(z, 20) + 6864605520245066446972130.0 / 1272466760992384639980987.0 * pow(z, 21) + 46019176815546718932450056.0 / 13997134370916231039790857.0 * pow(z, 22) + 2496220189180024162535320.0 / 424155586997461546660329.0 * pow(z, 23) + 236813033055148182272050.0 / 66971934789072875788473.0 * pow(z, 24) + 546947462295356294914.0 / 87316733753680411719.0 * pow(z, 25) + 30705639342929740.0 / 8303652049347447.0 * pow(z, 26) + 3906154894931797564.0 / 600273685244762217.0 * pow(z, 27) + 185515345102167610.0 / 49359374648773551.0 * pow(z, 28) + 46959622622537950.0 / 7121501814997179.0 * pow(z, 29) + 33204488118981680.0 / 8959308734996451.0 * pow(z, 30) + 25453056177374480.0 / 3911810856125211.0 * pow(z, 31) + 650634464711.0 / 184145876577.0 * pow(z, 32) + 3673654605820465.0 / 589450950922977.0 * pow(z, 33) + 173480011357030.0 / 53586450083907.0 * pow(z, 34) + 434411471342.0 / 75367721637.0 * pow(z, 35) + 212751500905.0 / 75367721637.0 * pow(z, 36) + 4632894925.0 / 908044839.0 * pow(z, 37) + 53615180.0 / 23283201.0 * pow(z, 38) + 1374380.0 / 322803.0 * pow(z, 39) + 11377541.0 / 6684561.0 * pow(z, 40) + 558005.0 / 171399.0 * pow(z, 41) + 17630.0 / 16587.0 * pow(z, 42) + 1870.0 / 873.0 * pow(z, 43) + 5.0 / 11.0 * pow(z, 44) + 1.0 * pow(z, 45));
                result += c_f[46]() *  2.6042739671849557 * (9.0 / 101.0 * pow(z, 0) + 46.0 / 909.0 * pow(z, 1) + 25415.0 / 88173.0 * pow(z, 2) + 94300.0 / 558429.0 * pow(z, 3) + 31168565.0 / 51933897.0 * pow(z, 4) + 12765046.0 / 35533719.0 * pow(z, 5) + 11063759.0 / 10867701.0 * pow(z, 6) + 56724860440.0 / 91712528739.0 * pow(z, 7) + 140206030625.0 / 91712528739.0 * pow(z, 8) + 793505597970.0 / 845793320593.0 * pow(z, 9) + 4995731920433.0 / 2362388240277.0 * pow(z, 10) + 228002300640668.0 / 174588111563697.0 * pow(z, 11) + 4970238584345335.0 / 1804077152824869.0 * pow(z, 12) + 95229226198610.0 / 55796200602831.0 * pow(z, 13) + 1353784425434105155.0 / 395092896468646311.0 * pow(z, 14) + 1527406453473157280.0 / 719271683314715079.0 * pow(z, 15) + 8856584826610657370.0 / 2157815049944145237.0 * pow(z, 16) + 366949352612087532580.0 / 144573608346257730879.0 * pow(z, 17) + 32086272351225479990.0 / 6736404689968998213.0 * pow(z, 18) + 66385592259414097880.0 / 22644059138570487969.0 * pow(z, 19) + 616409790006866544368078.0 / 114646871418582380587047.0 * pow(z, 20) + 7407511673965035141469724.0 / 2254721804565453484878591.0 * pow(z, 21) + 760723102652612363532638770.0 / 128519142860230848638079687.0 * pow(z, 22) + 460191768155467189324500560.0 / 128519142860230848638079687.0 * pow(z, 23) + 14353266087785138934578090.0 / 2254721804565453484878591.0 * pow(z, 24) + 435735980821472655380572.0 / 114646871418582380587047.0 * pow(z, 25) + 12579791632793194783022.0 / 1879456908501350501427.0 * pow(z, 26) + 9887215868423376280.0 / 2516006570952276441.0 * pow(z, 27) + 32086272351225479990.0 / 4663664785363152609.0 * pow(z, 28) + 8533705874699710060.0 / 2157815049944145237.0 * pow(z, 29) + 4968328073464515110.0 / 719271683314715079.0 * pow(z, 30) + 1527406453473157280.0 / 395092896468646311.0 * pow(z, 31) + 36588768254975815.0 / 5412231458474607.0 * pow(z, 32) + 68020875856150.0 / 18598733534277.0 * pow(z, 33) + 34791670090417345.0 / 5412231458474607.0 * pow(z, 34) + 228002300640668.0 / 68509258968033.0 * pow(z, 35) + 4995731920433.0 / 845793320593.0 * pow(z, 36) + 264501865990.0 / 91712528739.0 * pow(z, 37) + 476700504125.0 / 91712528739.0 * pow(z, 38) + 2466298280.0 / 1054166997.0 * pow(z, 39) + 1580537.0 / 366327.0 * pow(z, 40) + 89355322.0 / 51933897.0 * pow(z, 41) + 1833445.0 / 558429.0 * pow(z, 42) + 94300.0 / 88173.0 * pow(z, 43) + 1955.0 / 909.0 * pow(z, 44) + 46.0 / 101.0 * pow(z, 45) + 1.0 * pow(z, 46));
                result += c_f[47]() *  2.604396714762013 * (1.0 / 103.0 * pow(z, 0) + 973.0 / 10403.0 * pow(z, 1) + 2231.0 / 31209.0 * pow(z, 2) + 904015.0 / 3027273.0 * pow(z, 3) + 11546345.0 / 57518187.0 * pow(z, 4) + 57889781.0 / 93845463.0 * pow(z, 5) + 37643893.0 / 93845463.0 * pow(z, 6) + 8692639573.0 / 8352246207.0 * pow(z, 7) + 54011633585.0 / 80738380001.0 * pow(z, 8) + 125689940345.0 / 80738380001.0 * pow(z, 9) + 6671292962411.0 / 6701285540083.0 * pow(z, 10) + 388468039666001.0 / 180934709582241.0 * pow(z, 11) + 19559051693232571.0 / 14293842056997039.0 * pow(z, 12) + 39887123011201235.0 / 14293842056997039.0 * pow(z, 13) + 25333947877058675.0 / 14293842056997039.0 * pow(z, 14) + 3614543648030409955.0 / 1043450470160783847.0 * pow(z, 15) + 162291192781199714950.0 / 74084983381415653137.0 * pow(z, 16) + 1401464144620391810.0 / 338287595348929923.0 * pow(z, 17) + 4310050862437000637270.0 / 1654564628851616253393.0 * pow(z, 18) + 78653546484799769870.0 / 16381828008431844093.0 * pow(z, 19) + 153346608474365402474.0 / 51172101923245863507.0 * pow(z, 20) + 1639126778244326326289554.0 / 302785327079845774370919.0 * pow(z, 21) + 59794382192783366038919974.0 / 17864334297710900687884221.0 * pow(z, 22) + 69685827888041149985133770.0 / 11704219022638176312751731.0 * pow(z, 23) + 64993893687808006564043450.0 / 17864334297710900687884221.0 * pow(z, 24) + 1937149828834203840160382.0 / 302785327079845774370919.0 * pow(z, 25) + 21906658353480771782.0 / 5685789102582873723.0 * pow(z, 26) + 110114965078719677818.0 / 16381828008431844093.0 * pow(z, 27) + 226844782233526349330.0 / 57053952719021250117.0 * pow(z, 28) + 7007320723101959050.0 / 1014862786046789769.0 * pow(z, 29) + 9546540751835277350.0 / 2389838173594053327.0 * pow(z, 30) + 7229087296060819910.0 / 1043450470160783847.0 * pow(z, 31) + 55734685329529085.0 / 14293842056997039.0 * pow(z, 32) + 96868727312917285.0 / 14293842056997039.0 * pow(z, 33) + 52658985327933845.0 / 14293842056997039.0 * pow(z, 34) + 388468039666001.0 / 60311569860747.0 * pow(z, 35) + 22439803600837.0 / 6701285540083.0 * pow(z, 36) + 477621773311.0 / 80738380001.0 * pow(z, 37) + 702151236605.0 / 242215140003.0 * pow(z, 38) + 43463197865.0 / 8352246207.0 * pow(z, 39) + 220485659.0 / 93845463.0 * pow(z, 40) + 405228467.0 / 93845463.0 * pow(z, 41) + 99298567.0 / 57518187.0 * pow(z, 42) + 9944165.0 / 3027273.0 * pow(z, 43) + 11155.0 / 10403.0 * pow(z, 44) + 22379.0 / 10403.0 * pow(z, 45) + 47.0 / 103.0 * pow(z, 46) + 1.0 * pow(z, 47));
                result += c_f[48]() *  2.614016928281532 * (3.0 / 35.0 * pow(z, 0) + 176.0 / 3605.0 * pow(z, 1) + 14456.0 / 52015.0 * pow(z, 2) + 35696.0 / 218463.0 * pow(z, 3) + 1756372.0 / 3027273.0 * pow(z, 4) + 36948304.0 / 105954555.0 * pow(z, 5) + 463118248.0 / 469227315.0 * pow(z, 6) + 1978993232.0 / 3284591205.0 * pow(z, 7) + 86926395730.0 / 58465723449.0 * pow(z, 8) + 74073097488.0 / 80738380001.0 * pow(z, 9) + 201103904552.0 / 97442872415.0 * pow(z, 10) + 1386242693488.0 / 1080852506465.0 * pow(z, 11) + 17092593745304044.0 / 6332714835378435.0 * pow(z, 12) + 24072679007055472.0 / 14293842056997039.0 * pow(z, 13) + 337331097466159016.0 / 100056894398979273.0 * pow(z, 14) + 210778446337128176.0 / 100056894398979273.0 * pow(z, 15) + 4234179701978480233.0 / 1043450470160783847.0 * pow(z, 16) + 1313604007452534163360.0 / 518594883669909571959.0 * pow(z, 17) + 1601673308137590640.0 / 338287595348929923.0 * pow(z, 18) + 34117455247922362939232.0 / 11581952401961313773751.0 * pow(z, 19) + 440459860314878711272.0 / 81909140042159220465.0 * pow(z, 20) + 5958611072146769924704.0 / 1791023567313605222745.0 * pow(z, 21) + 63180886725054032940615536.0 / 10597486447794602102982165.0 * pow(z, 22) + 65365287366024052315837984.0 / 17864334297710900687884221.0 * pow(z, 23) + 27874331155216459994053508.0 / 4312080692550907062592743.0 * pow(z, 24) + 5942298851456732028712544.0 / 1513926635399228871854595.0 * pow(z, 25) + 1192092202359510055483312.0 / 173729286029419706606265.0 * pow(z, 26) + 350506533655692348512.0 / 85286836538743105845.0 * pow(z, 27) + 817996883441917606648.0 / 114672796059022908651.0 * pow(z, 28) + 725903303147284317856.0 / 172864961223303190653.0 * pow(z, 29) + 7367697217432916944.0 / 1014862786046789769.0 * pow(z, 30) + 30548930405872887520.0 / 7304153291125486929.0 * pow(z, 31) + 103272675658011713.0 / 14293842056997039.0 * pow(z, 32) + 405343166032938800.0 / 100056894398979273.0 * pow(z, 33) + 100287623571020248.0 / 14293842056997039.0 * pow(z, 34) + 24072679007055472.0 / 6332714835378435.0 * pow(z, 35) + 1553872158664004.0 / 234544993902905.0 * pow(z, 36) + 1386242693488.0 / 403691900005.0 * pow(z, 37) + 3418766377384.0 / 565168660007.0 * pow(z, 38) + 24691032496.0 / 8352246207.0 * pow(z, 39) + 17385279146.0 / 3284591205.0 * pow(z, 40) + 1118561392.0 / 469227315.0 * pow(z, 41) + 463118248.0 / 105954555.0 * pow(z, 42) + 36948304.0 / 21190911.0 * pow(z, 43) + 103316.0 / 31209.0 * pow(z, 44) + 392656.0 / 364105.0 * pow(z, 45) + 1112.0 / 515.0 * pow(z, 46) + 16.0 / 35.0 * pow(z, 47) + 1.0 * pow(z, 48));
                result += c_f[49]() *  2.614131094941481 * (1.0 / 107.0 * pow(z, 0) + 337.0 / 3745.0 * pow(z, 1) + 26616.0 / 385735.0 * pow(z, 2) + 2244040.0 / 7791847.0 * pow(z, 3) + 1514228.0 / 7791847.0 * pow(z, 4) + 450823644.0 / 755809159.0 * pow(z, 5) + 210313656.0 / 539863685.0 * pow(z, 6) + 16890737976.0 / 16735774235.0 * pow(z, 7) + 15275836998.0 / 23430083929.0 * pow(z, 8) + 3157987190194.0 / 2085277469681.0 * pow(z, 9) + 58899516395896.0 / 60473046620749.0 * pow(z, 10) + 633729114679272.0 / 302365233103745.0 * pow(z, 11) + 33741041258761876.0 / 25096314347610835.0 * pow(z, 12) + 1960825618787924.0 / 717037552788881.0 * pow(z, 13) + 99113491752813176.0 / 56645966670321599.0 * pow(z, 14) + 4055530066101328168.0 / 1189565300076753579.0 * pow(z, 15) + 2586244720853199647.0 / 1189565300076753579.0 * pow(z, 16) + 1332480036282877479.0 / 325236954702633001.0 * pow(z, 17) + 16035613563947708448656.0 / 6165516950297813799957.0 * pow(z, 18) + 29433531590728775913712.0 / 6165516950297813799957.0 * pow(z, 19) + 59256911059608357345544.0 / 19670935031902548790339.0 * pow(z, 20) + 1598016281395832130174056.0 / 295064025478538231855085.0 * pow(z, 21) + 98648695649719616751152.0 / 29090819413377008774445.0 * pow(z, 22) + 50384822011414656658124816.0 / 8399489258622388333474753.0 * pow(z, 23) + 5529661760072588058265577972.0 / 1486709598776162735025031281.0 * pow(z, 24) + 163750671537097634138905652.0 / 25198467775867165000424259.0 * pow(z, 25) + 38601663515107676120016.0 / 9696939804459002924815.0 * pow(z, 26) + 2033838903594695438403344.0 / 295064025478538231855085.0 * pow(z, 27) + 8465273008515479620792.0 / 2034924313645091254173.0 * pow(z, 28) + 14716765795364387956856.0 / 2055172316765937933319.0 * pow(z, 29) + 843979661260405707824.0 / 198887643557993993547.0 * pow(z, 30) + 7106560193508679888.0 / 975710864107899003.0 * pow(z, 31) + 1673452466434423301.0 / 396521766692251193.0 * pow(z, 32) + 8618001390465322357.0 / 1189565300076753579.0 * pow(z, 33) + 693794442269692232.0 / 169937900010964797.0 * pow(z, 34) + 5042123019740376.0 / 717037552788881.0 * pow(z, 35) + 96032194351860724.0 / 25096314347610835.0 * pow(z, 36) + 2006808863151028.0 / 302365233103745.0 * pow(z, 37) + 208825558130904.0 / 60473046620749.0 * pow(z, 38) + 12631948760776.0 / 2085277469681.0 * pow(z, 39) + 69589924102.0 / 23430083929.0 * pow(z, 40) + 88676374374.0 / 16735774235.0 * pow(z, 41) + 1291926744.0 / 539863685.0 * pow(z, 42) + 3306040056.0 / 755809159.0 * pow(z, 43) + 13628052.0 / 7791847.0 * pow(z, 44) + 25806460.0 / 7791847.0 * pow(z, 45) + 416984.0 / 385735.0 * pow(z, 46) + 8088.0 / 3745.0 * pow(z, 47) + 49.0 / 107.0 * pow(z, 48) + 1.0 * pow(z, 49));

                return result;
            }

            virtual complex<double> f_p(const double & q2) const
            {
                const auto z = _z(q2);
                //return this->phi_p(z);
                const auto phi      = this->phi_p(z);

                const double epsilon = Process_::epsilon_isospin;
                const auto series_total   = this->series_p(z, _a_fp) + (epsilon / 6) * this->series_p(z, _a_ft);

                complex<double> correction_to_asymptotics = pow(1.0 - z, 5.0 / 2.0 - 1.0 / 2.0);
                complex<double> K = (315 * M_PI) / 8192;

                return 1.0 / (phi) * correction_to_asymptotics * sqrt(K) * series_total;
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

    template <typename Process_> class EGJvD2020FormFactorBase<Process_, PToP, false> :
        public FormFactors<PToP>
    {
        private:
            // parameters for form factors f_+ and f_T
            std::array<UsedParameter, 50u> _a_fp, _a_ft;

            // parameter for zero point of z
            UsedParameter _t_0;

            std::string _par_name(const std::string & ff, const std::string & index) const
            {
                return stringify(Process_::label) + "::" + "a_" + ff + "^" + index + "@EGJvD2020";
            }

            double _z(const double & q2, const complex<double> & t_s) const
            {
                const complex<double> t_p   = complex<double>{ Process_::t_p, 0.0};
                const complex<double> t     = complex<double>{ q2, 0.0};
                const complex<double> result= (sqrt(t_p - t) - sqrt(t_p - t_s)) / (sqrt(t_p - t) + sqrt(t_p - t_s));

                return result.real();
            }

            // double _z(const double & q2, const double & t_0) const
            // {
            //     const double t_p = Process_::t_p;
            //     const double a = sqrt(t_p - t_0);
            //     const double z = (sqrt(t_p - q2) - a) / (sqrt(t_p - q2) + a);

            //     return z;
            // }

            double _z(const double & q2) const
            {
                const complex<double> t_s = complex<double>(this->_t_0, 0.0);
                return _z(q2, t_s);
            }

        public:
            EGJvD2020FormFactorBase(const Parameters & p, const Options & o) :
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
                    UsedParameter(p[_par_name("+", "9")], *this),

                    UsedParameter(p[_par_name("+", "10")], *this),
                    UsedParameter(p[_par_name("+", "11")], *this),
                    UsedParameter(p[_par_name("+", "12")], *this),
                    UsedParameter(p[_par_name("+", "13")], *this),
                    UsedParameter(p[_par_name("+", "14")], *this),

                    UsedParameter(p[_par_name("+", "15")], *this),
                    UsedParameter(p[_par_name("+", "16")], *this),
                    UsedParameter(p[_par_name("+", "17")], *this),
                    UsedParameter(p[_par_name("+", "18")], *this),
                    UsedParameter(p[_par_name("+", "19")], *this),

                    UsedParameter(p[_par_name("+", "20")], *this),
                    UsedParameter(p[_par_name("+", "21")], *this),
                    UsedParameter(p[_par_name("+", "22")], *this),
                    UsedParameter(p[_par_name("+", "23")], *this),
                    UsedParameter(p[_par_name("+", "24")], *this),

                    UsedParameter(p[_par_name("+", "25")], *this),
                    UsedParameter(p[_par_name("+", "26")], *this),
                    UsedParameter(p[_par_name("+", "27")], *this),
                    UsedParameter(p[_par_name("+", "28")], *this),
                    UsedParameter(p[_par_name("+", "29")], *this),

                    UsedParameter(p[_par_name("+", "30")], *this),
                    UsedParameter(p[_par_name("+", "31")], *this),
                    UsedParameter(p[_par_name("+", "32")], *this),
                    UsedParameter(p[_par_name("+", "33")], *this),
                    UsedParameter(p[_par_name("+", "34")], *this),

                    UsedParameter(p[_par_name("+", "35")], *this),
                    UsedParameter(p[_par_name("+", "36")], *this),
                    UsedParameter(p[_par_name("+", "37")], *this),
                    UsedParameter(p[_par_name("+", "38")], *this),
                    UsedParameter(p[_par_name("+", "39")], *this),

                    UsedParameter(p[_par_name("+", "40")], *this),
                    UsedParameter(p[_par_name("+", "41")], *this),
                    UsedParameter(p[_par_name("+", "42")], *this),
                    UsedParameter(p[_par_name("+", "43")], *this),
                    UsedParameter(p[_par_name("+", "44")], *this),

                    UsedParameter(p[_par_name("+", "45")], *this),
                    UsedParameter(p[_par_name("+", "46")], *this),
                    UsedParameter(p[_par_name("+", "47")], *this),
                    UsedParameter(p[_par_name("+", "48")], *this),
                    UsedParameter(p[_par_name("+", "49")], *this)
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
                    UsedParameter(p[_par_name("T", "9")], *this),

                    UsedParameter(p[_par_name("T", "10")], *this),
                    UsedParameter(p[_par_name("T", "11")], *this),
                    UsedParameter(p[_par_name("T", "12")], *this),
                    UsedParameter(p[_par_name("T", "13")], *this),
                    UsedParameter(p[_par_name("T", "14")], *this),

                    UsedParameter(p[_par_name("T", "15")], *this),
                    UsedParameter(p[_par_name("T", "16")], *this),
                    UsedParameter(p[_par_name("T", "17")], *this),
                    UsedParameter(p[_par_name("T", "18")], *this),
                    UsedParameter(p[_par_name("T", "19")], *this),

                    UsedParameter(p[_par_name("T", "20")], *this),
                    UsedParameter(p[_par_name("T", "21")], *this),
                    UsedParameter(p[_par_name("T", "22")], *this),
                    UsedParameter(p[_par_name("T", "23")], *this),
                    UsedParameter(p[_par_name("T", "24")], *this),

                    UsedParameter(p[_par_name("T", "25")], *this),
                    UsedParameter(p[_par_name("T", "26")], *this),
                    UsedParameter(p[_par_name("T", "27")], *this),
                    UsedParameter(p[_par_name("T", "28")], *this),
                    UsedParameter(p[_par_name("T", "29")], *this),

                    UsedParameter(p[_par_name("T", "30")], *this),
                    UsedParameter(p[_par_name("T", "31")], *this),
                    UsedParameter(p[_par_name("T", "32")], *this),
                    UsedParameter(p[_par_name("T", "33")], *this),
                    UsedParameter(p[_par_name("T", "34")], *this),

                    UsedParameter(p[_par_name("T", "35")], *this),
                    UsedParameter(p[_par_name("T", "36")], *this),
                    UsedParameter(p[_par_name("T", "37")], *this),
                    UsedParameter(p[_par_name("T", "38")], *this),
                    UsedParameter(p[_par_name("T", "39")], *this),

                    UsedParameter(p[_par_name("T", "40")], *this),
                    UsedParameter(p[_par_name("T", "41")], *this),
                    UsedParameter(p[_par_name("T", "42")], *this),
                    UsedParameter(p[_par_name("T", "43")], *this),
                    UsedParameter(p[_par_name("T", "44")], *this),

                    UsedParameter(p[_par_name("T", "45")], *this),
                    UsedParameter(p[_par_name("T", "46")], *this),
                    UsedParameter(p[_par_name("T", "47")], *this),
                    UsedParameter(p[_par_name("T", "48")], *this),
                    UsedParameter(p[_par_name("T", "49")], *this)
                }},
                _t_0(p[stringify(Process_::label) + "::t_0@EGJvD2020"], *this)
            {
                // TODO: empty body?
            }
            /* f_+ */

            double phi_p(const double & z) const //Divided by asymptotic behaviour to improve numeric behaviour of Formfactor
            {
                const double t_p = Process_::t_p;
                const double t_0 = this->_t_0;
                const double tfactor = 1.0 - t_0 / t_p;
                const double chi = 0.00405; //GeV^-2 as a makeshift value
                const double Q2 = Process_::Q2;

                const double part0 = 1.0 / sqrt(12.0 * M_PI * t_p * chi);
                const double part1 = pow(tfactor, 1.25); /*[Asymptotics:] (1.0 + z) * (1.0 + z) * sqrt(1.0 - z) */
                const double part2 = pow(sqrt(tfactor) * (1.0 + z) + (1.0 - z), -0.5);
                const double part3 = pow(sqrt(1.0 + Q2 / t_p) * (1.0 - z) + sqrt(tfactor) * (1.0 + z), -3.0);

                return part0 * part1 * part2 * part3;
            }

            double series_p(const double & z, std::array<UsedParameter, 50u> c_f) const
            {
                double result = 0.0;

                result += c_f[0]() *  1.0 * (1.0 * pow(z, 0));
                result += c_f[1]() *  1.0041580220928046 * (1.0 / 11.0 * pow(z, 0) + 1.0 * pow(z, 1));
                result += c_f[2]() *  1.3915668626887223 * (9.0 / 13.0 * pow(z, 0) + 2.0 / 13.0 * pow(z, 1) + 1.0 * pow(z, 2));
                result += c_f[3]() *  1.394669579723865 * (1.0 / 15.0 * pow(z, 0) + 137.0 / 195.0 * pow(z, 1) + 1.0 / 5.0 * pow(z, 2) + 1.0 * pow(z, 3));
                result += c_f[4]() *  1.6439499152771448 * (9.0 / 17.0 * pow(z, 0) + 44.0 / 255.0 * pow(z, 1) + 274.0 / 255.0 * pow(z, 2) + 4.0 / 17.0 * pow(z, 3) + 1.0 * pow(z, 4));
                result += c_f[5]() *  1.6462315956469282 * (1.0 / 19.0 * pow(z, 0) + 175.0 / 323.0 * pow(z, 1) + 74.0 / 323.0 * pow(z, 2) + 350.0 / 323.0 * pow(z, 3) + 5.0 / 19.0 * pow(z, 4) + 1.0 * pow(z, 5));
                result += c_f[6]() *  1.822044489432169 * (3.0 / 7.0 * pow(z, 0) + 22.0 / 133.0 * pow(z, 1) + 325.0 / 323.0 * pow(z, 2) + 740.0 / 2261.0 * pow(z, 3) + 25.0 / 19.0 * pow(z, 4) + 2.0 / 7.0 * pow(z, 5) + 1.0 * pow(z, 6));
                result += c_f[7]() *  1.8237690941622522 * (1.0 / 23.0 * pow(z, 0) + 71.0 / 161.0 * pow(z, 1) + 681.0 / 3059.0 * pow(z, 2) + 53065.0 / 52003.0 * pow(z, 3) + 1135.0 / 3059.0 * pow(z, 4) + 213.0 / 161.0 * pow(z, 5) + 7.0 / 23.0 * pow(z, 6) + 1.0 * pow(z, 7));
                result += c_f[8]() *  1.9548363704716327 * (9.0 / 25.0 * pow(z, 0) + 88.0 / 575.0 * pow(z, 1) + 3692.0 / 4025.0 * pow(z, 2) + 5448.0 / 15295.0 * pow(z, 3) + 21226.0 / 15295.0 * pow(z, 4) + 1816.0 / 4025.0 * pow(z, 5) + 852.0 / 575.0 * pow(z, 6) + 8.0 / 25.0 * pow(z, 7) + 1.0 * pow(z, 8));
                result += c_f[9]() *  1.9561785171250927 * (1.0 / 27.0 * pow(z, 0) + 251.0 / 675.0 * pow(z, 1) + 1076.0 / 5175.0 * pow(z, 2) + 580.0 / 621.0 * pow(z, 3) + 24046.0 / 58995.0 * pow(z, 4) + 290.0 / 207.0 * pow(z, 5) + 7532.0 / 15525.0 * pow(z, 6) + 1004.0 / 675.0 * pow(z, 7) + 1.0 / 3.0 * pow(z, 8) + 1.0 * pow(z, 9));
                result += c_f[10]() *  2.05778352996703 * (9.0 / 29.0 * pow(z, 0) + 110.0 / 783.0 * pow(z, 1) + 3263.0 / 3915.0 * pow(z, 2) + 2152.0 / 6003.0 * pow(z, 3) + 850.0 / 621.0 * pow(z, 4) + 48092.0 / 90045.0 * pow(z, 5) + 350.0 / 207.0 * pow(z, 6) + 2152.0 / 3915.0 * pow(z, 7) + 1255.0 / 783.0 * pow(z, 8) + 10.0 / 29.0 * pow(z, 9) + 1.0 * pow(z, 10));
                result += c_f[11]() *  2.0588550132627392 * (1.0 / 31.0 * pow(z, 0) + 289.0 / 899.0 * pow(z, 1) + 1555.0 / 8091.0 * pow(z, 2) + 6887.0 / 8091.0 * pow(z, 3) + 76862.0 / 186093.0 * pow(z, 4) + 1289614.0 / 930465.0 * pow(z, 5) + 538034.0 / 930465.0 * pow(z, 6) + 13774.0 / 8091.0 * pow(z, 7) + 1555.0 / 2697.0 * pow(z, 8) + 1445.0 / 899.0 * pow(z, 9) + 11.0 / 31.0 * pow(z, 10) + 1.0 * pow(z, 11));
                result += c_f[12]() *  2.1399786377482064 * (3.0 / 11.0 * pow(z, 0) + 4.0 / 31.0 * pow(z, 1) + 7514.0 / 9889.0 * pow(z, 2) + 31100.0 / 89001.0 * pow(z, 3) + 117079.0 / 89001.0 * pow(z, 4) + 5841512.0 / 10235115.0 * pow(z, 5) + 18054596.0 / 10235115.0 * pow(z, 6) + 307448.0 / 445005.0 * pow(z, 7) + 172175.0 / 89001.0 * pow(z, 8) + 6220.0 / 9889.0 * pow(z, 9) + 578.0 / 341.0 * pow(z, 10) + 4.0 / 11.0 * pow(z, 11) + 1.0 * pow(z, 12));
                result += c_f[13]() *  2.140852633552563 * (1.0 / 35.0 * pow(z, 0) + 109.0 / 385.0 * pow(z, 1) + 2118.0 / 11935.0 * pow(z, 2) + 53842.0 / 69223.0 * pow(z, 3) + 28015.0 / 69223.0 * pow(z, 4) + 2310697.0 / 1730575.0 * pow(z, 5) + 963236.0 / 1550775.0 * pow(z, 6) + 9242788.0 / 5191725.0 * pow(z, 7) + 50427.0 / 69223.0 * pow(z, 8) + 134605.0 / 69223.0 * pow(z, 9) + 706.0 / 1085.0 * pow(z, 10) + 654.0 / 385.0 * pow(z, 11) + 13.0 / 35.0 * pow(z, 12) + 1.0 * pow(z, 13));
                result += c_f[14]() *  2.207143478674807 * (9.0 / 37.0 * pow(z, 0) + 22.0 / 185.0 * pow(z, 1) + 1417.0 / 2035.0 * pow(z, 2) + 4236.0 / 12617.0 * pow(z, 3) + 457657.0 / 365893.0 * pow(z, 4) + 212914.0 / 365893.0 * pow(z, 5) + 16174879.0 / 9147325.0 * pow(z, 6) + 1926472.0 / 2494725.0 * pow(z, 7) + 2310697.0 / 1097679.0 * pow(z, 8) + 302562.0 / 365893.0 * pow(z, 9) + 26921.0 / 12617.0 * pow(z, 10) + 1412.0 / 2035.0 * pow(z, 11) + 327.0 / 185.0 * pow(z, 12) + 14.0 / 37.0 * pow(z, 13) + 1.0 * pow(z, 14));
                result += c_f[15]() *  2.207869393339616 * (1.0 / 39.0 * pow(z, 0) + 365.0 / 1443.0 * pow(z, 1) + 79.0 / 481.0 * pow(z, 2) + 11335.0 / 15873.0 * pow(z, 3) + 192125.0 / 492063.0 * pow(z, 4) + 465415.0 / 365893.0 * pow(z, 5) + 27221635.0 / 42809481.0 * pow(z, 6) + 382730407.0 / 214047405.0 * pow(z, 7) + 3888805.0 / 4756609.0 * pow(z, 8) + 2327075.0 / 1097679.0 * pow(z, 9) + 38425.0 / 44733.0 * pow(z, 10) + 11335.0 / 5291.0 * pow(z, 11) + 79.0 / 111.0 * pow(z, 12) + 2555.0 / 1443.0 * pow(z, 13) + 5.0 / 13.0 * pow(z, 14) + 1.0 * pow(z, 15));
                result += c_f[16]() *  2.2630661281731066 * (9.0 / 41.0 * pow(z, 0) + 176.0 / 1599.0 * pow(z, 1) + 2920.0 / 4551.0 * pow(z, 2) + 6320.0 / 19721.0 * pow(z, 3) + 770780.0 / 650793.0 * pow(z, 4) + 11681200.0 / 20174583.0 * pow(z, 5) + 26063240.0 / 15001613.0 * pow(z, 6) + 1431080240.0 / 1755188721.0 * pow(z, 7) + 3827304070.0 / 1755188721.0 * pow(z, 8) + 186662640.0 / 195020969.0 * pow(z, 9) + 3723320.0 / 1551891.0 * pow(z, 10) + 614800.0 / 650793.0 * pow(z, 11) + 45340.0 / 19721.0 * pow(z, 12) + 44240.0 / 59163.0 * pow(z, 13) + 2920.0 / 1599.0 * pow(z, 14) + 16.0 / 41.0 * pow(z, 15) + 1.0 * pow(z, 16));
                result += c_f[17]() *  2.2636783468041286 * (1.0 / 43.0 * pow(z, 0) + 403.0 / 1763.0 * pow(z, 1) + 3496.0 / 22919.0 * pow(z, 2) + 558840.0 / 848003.0 * pow(z, 3) + 317100.0 / 848003.0 * pow(z, 4) + 1022980.0 / 848003.0 * pow(z, 5) + 16687720.0 / 26288093.0 * pow(z, 6) + 1341455800.0 / 762354697.0 * pow(z, 7) + 1980719830.0 / 2287064091.0 * pow(z, 8) + 1676819750.0 / 762354697.0 * pow(z, 9) + 26223560.0 / 26288093.0 * pow(z, 10) + 2045960.0 / 848003.0 * pow(z, 11) + 63420.0 / 65231.0 * pow(z, 12) + 1955940.0 / 848003.0 * pow(z, 13) + 17480.0 / 22919.0 * pow(z, 14) + 3224.0 / 1763.0 * pow(z, 15) + 17.0 / 43.0 * pow(z, 16) + 1.0 * pow(z, 17));
                result += c_f[18]() *  2.310357038107123 * (1.0 / 5.0 * pow(z, 0) + 22.0 / 215.0 * pow(z, 1) + 5239.0 / 8815.0 * pow(z, 2) + 6992.0 / 22919.0 * pow(z, 3) + 950028.0 / 848003.0 * pow(z, 4) + 481992.0 / 848003.0 * pow(z, 5) + 1432172.0 / 848003.0 * pow(z, 6) + 21932432.0 / 26288093.0 * pow(z, 7) + 1676819750.0 / 762354697.0 * pow(z, 8) + 792287932.0 / 762354697.0 * pow(z, 9) + 67072790.0 / 26288093.0 * pow(z, 10) + 953584.0 / 848003.0 * pow(z, 11) + 2250556.0 / 848003.0 * pow(z, 12) + 887880.0 / 848003.0 * pow(z, 13) + 55884.0 / 22919.0 * pow(z, 14) + 6992.0 / 8815.0 * pow(z, 15) + 403.0 / 215.0 * pow(z, 16) + 2.0 / 5.0 * pow(z, 17) + 1.0 * pow(z, 18));
                result += c_f[19]() *  2.310880157560925 * (1.0 / 47.0 * pow(z, 0) + 49.0 / 235.0 * pow(z, 1) + 1437.0 / 10105.0 * pow(z, 2) + 50645.0 / 82861.0 * pow(z, 3) + 384508.0 / 1077193.0 * pow(z, 4) + 45539196.0 / 39856141.0 * pow(z, 5) + 24904180.0 / 39856141.0 * pow(z, 6) + 68265668.0 / 39856141.0 * pow(z, 7) + 1097897094.0 / 1235540371.0 * pow(z, 8) + 1941532102.0 / 873918799.0 * pow(z, 9) + 1341874226.0 / 1235540371.0 * pow(z, 10) + 102398502.0 / 39856141.0 * pow(z, 11) + 3557740.0 / 3065857.0 * pow(z, 12) + 106258124.0 / 39856141.0 * pow(z, 13) + 1153524.0 / 1077193.0 * pow(z, 14) + 202580.0 / 82861.0 * pow(z, 15) + 8143.0 / 10105.0 * pow(z, 16) + 441.0 / 235.0 * pow(z, 17) + 19.0 / 47.0 * pow(z, 18) + 1.0 * pow(z, 19));
                result += c_f[20]() *  2.350874856721882 * (9.0 / 49.0 * pow(z, 0) + 220.0 / 2303.0 * pow(z, 1) + 26.0 / 47.0 * pow(z, 2) + 28740.0 / 99029.0 * pow(z, 3) + 614975.0 / 580027.0 * pow(z, 4) + 29222608.0 / 52782457.0 * pow(z, 5) + 455391960.0 / 278992987.0 * pow(z, 6) + 1636560400.0 / 1952950909.0 * pow(z, 7) + 4266604250.0 / 1952950909.0 * pow(z, 8) + 9410546520.0 / 8648782597.0 * pow(z, 9) + 3883064204.0 / 1476621419.0 * pow(z, 10) + 348538760.0 / 278992987.0 * pow(z, 11) + 5631917610.0 / 1952950909.0 * pow(z, 12) + 355774000.0 / 278992987.0 * pow(z, 13) + 151797320.0 / 52782457.0 * pow(z, 14) + 4614096.0 / 4060189.0 * pow(z, 15) + 36175.0 / 14147.0 * pow(z, 16) + 1916.0 / 2303.0 * pow(z, 17) + 90.0 / 47.0 * pow(z, 18) + 20.0 / 49.0 * pow(z, 19) + 1.0 * pow(z, 20));
                result += c_f[21]() *  2.351326904578521 * (1.0 / 51.0 * pow(z, 0) + 479.0 / 2499.0 * pow(z, 1) + 5210.0 / 39151.0 * pow(z, 2) + 66890.0 / 117453.0 * pow(z, 3) + 1718965.0 / 5050479.0 * pow(z, 4) + 74720057.0 / 69023213.0 * pow(z, 5) + 2577944.0 / 4225911.0 * pow(z, 6) + 1813904920.0 / 1094510949.0 * pow(z, 7) + 2284525590.0 / 2553858881.0 * pow(z, 8) + 995640670.0 / 450680979.0 * pow(z, 9) + 7315531684.0 / 6419158809.0 * pow(z, 10) + 398256268.0 / 150226993.0 * pow(z, 11) + 9899610890.0 / 7661576643.0 * pow(z, 12) + 453476230.0 / 156358707.0 * pow(z, 13) + 12889720.0 / 9860459.0 * pow(z, 14) + 597760456.0 / 207069639.0 * pow(z, 15) + 343793.0 / 297087.0 * pow(z, 16) + 100335.0 / 39151.0 * pow(z, 17) + 98990.0 / 117453.0 * pow(z, 18) + 4790.0 / 2499.0 * pow(z, 19) + 7.0 / 17.0 * pow(z, 20) + 1.0 * pow(z, 21));
                result += c_f[22]() *  2.3859794553105305 * (9.0 / 53.0 * pow(z, 0) + 242.0 / 2703.0 * pow(z, 1) + 68497.0 / 132447.0 * pow(z, 2) + 573100.0 / 2075003.0 * pow(z, 3) + 367895.0 / 366177.0 * pow(z, 4) + 143705474.0 / 267675387.0 * pow(z, 5) + 821920627.0 / 522604327.0 * pow(z, 6) + 1304439664.0 / 1567812981.0 * pow(z, 7) + 124705963250.0 / 58009080297.0 * pow(z, 8) + 150778688940.0 / 135354520693.0 * pow(z, 9) + 63521874746.0 / 23886091887.0 * pow(z, 10) + 14631063368.0 / 10974690867.0 * pow(z, 11) + 24094504214.0 / 7962030629.0 * pow(z, 12) + 83765938300.0 / 58009080297.0 * pow(z, 13) + 4988238530.0 / 1567812981.0 * pow(z, 14) + 737291984.0 / 522604327.0 * pow(z, 15) + 821920627.0 / 267675387.0 * pow(z, 16) + 7563446.0 / 6225009.0 * pow(z, 17) + 5518425.0 / 2075003.0 * pow(z, 18) + 114620.0 / 132447.0 * pow(z, 19) + 5269.0 / 2703.0 * pow(z, 20) + 22.0 / 53.0 * pow(z, 21) + 1.0 * pow(z, 22));
                result += c_f[23]() *  2.3863739298863123 * (1.0 / 55.0 * pow(z, 0) + 47.0 / 265.0 * pow(z, 1) + 563.0 / 4505.0 * pow(z, 2) + 3361.0 / 6307.0 * pow(z, 3) + 96205.0 / 296429.0 * pow(z, 4) + 1521841.0 / 1482145.0 * pow(z, 5) + 5396241.0 / 9104605.0 * pow(z, 6) + 4176629679.0 / 2613021635.0 * pow(z, 7) + 465044970.0 / 522604327.0 * pow(z, 8) + 42076326770.0 / 19336360099.0 * pow(z, 9) + 113018629394.0 / 96681800495.0 * pow(z, 10) + 2854014392006.0 / 1063499805445.0 * pow(z, 11) + 133567471102.0 / 96681800495.0 * pow(z, 12) + 8415265354.0 / 2762337157.0 * pow(z, 13) + 775074950.0 / 522604327.0 * pow(z, 14) + 8353259358.0 / 2613021635.0 * pow(z, 15) + 5396241.0 / 3748955.0 * pow(z, 16) + 4565523.0 / 1482145.0 * pow(z, 17) + 365579.0 / 296429.0 * pow(z, 18) + 16805.0 / 6307.0 * pow(z, 19) + 3941.0 / 4505.0 * pow(z, 20) + 517.0 / 265.0 * pow(z, 21) + 23.0 / 55.0 * pow(z, 22) + 1.0 * pow(z, 23));
                result += c_f[24]() *  2.416688998884546 * (3.0 / 19.0 * pow(z, 0) + 8.0 / 95.0 * pow(z, 1) + 2444.0 / 5035.0 * pow(z, 2) + 4504.0 / 17119.0 * pow(z, 3) + 6722.0 / 7049.0 * pow(z, 4) + 153928.0 / 296429.0 * pow(z, 5) + 6087364.0 / 4022965.0 * pow(z, 6) + 992908344.0 / 1210912465.0 * pow(z, 7) + 20883148395.0 / 9929482213.0 * pow(z, 8) + 11161079280.0 / 9929482213.0 * pow(z, 9) + 976170781064.0 / 367390841881.0 * pow(z, 10) + 2548056371792.0 / 1836954209405.0 * pow(z, 11) + 5708028784012.0 / 1836954209405.0 * pow(z, 12) + 82195366832.0 / 52484405983.0 * pow(z, 13) + 33661061416.0 / 9929482213.0 * pow(z, 14) + 16121558960.0 / 9929482213.0 * pow(z, 15) + 4176629679.0 / 1210912465.0 * pow(z, 16) + 43169928.0 / 28160755.0 * pow(z, 17) + 18262092.0 / 5632151.0 * pow(z, 18) + 153928.0 / 119833.0 * pow(z, 19) + 47054.0 / 17119.0 * pow(z, 20) + 4504.0 / 5035.0 * pow(z, 21) + 188.0 / 95.0 * pow(z, 22) + 8.0 / 19.0 * pow(z, 23) + 1.0 * pow(z, 24));
                result += c_f[25]() *  2.4170361993746483 * (1.0 / 59.0 * pow(z, 0) + 185.0 / 1121.0 * pow(z, 1) + 132.0 / 1121.0 * pow(z, 2) + 29740.0 / 59413.0 * pow(z, 3) + 312790.0 / 1010021.0 * pow(z, 4) + 6896094.0 / 7070147.0 * pow(z, 5) + 27259340.0 / 47470987.0 * pow(z, 6) + 511450252.0 / 332296909.0 * pow(z, 7) + 12551644395.0 / 14288767087.0 * pow(z, 8) + 65696174435.0 / 30833655293.0 * pow(z, 9) + 692164738936.0 / 585839450567.0 * pow(z, 10) + 58169443650600.0 / 21676059670979.0 * pow(z, 11) + 31208670943948.0 / 21676059670979.0 * pow(z, 12) + 9694907275100.0 / 3096579952997.0 * pow(z, 13) + 943861007640.0 / 585839450567.0 * pow(z, 14) + 105113879096.0 / 30833655293.0 * pow(z, 15) + 1394627155.0 / 840515711.0 * pow(z, 16) + 1150763067.0 / 332296909.0 * pow(z, 17) + 27259340.0 / 17489311.0 * pow(z, 18) + 22986980.0 / 7070147.0 * pow(z, 19) + 1313718.0 / 1010021.0 * pow(z, 20) + 163570.0 / 59413.0 * pow(z, 21) + 1012.0 / 1121.0 * pow(z, 22) + 2220.0 / 1121.0 * pow(z, 23) + 25.0 / 59.0 * pow(z, 24) + 1.0 * pow(z, 25));
                result += c_f[26]() *  2.443781079080322 * (9.0 / 61.0 * pow(z, 0) + 286.0 / 3599.0 * pow(z, 1) + 31265.0 / 68381.0 * pow(z, 2) + 17160.0 / 68381.0 * pow(z, 3) + 3286270.0 / 3624193.0 * pow(z, 4) + 1626508.0 / 3242699.0 * pow(z, 5) + 89649222.0 / 61611281.0 * pow(z, 6) + 16301085320.0 / 20270111449.0 * pow(z, 7) + 41555332975.0 / 20270111449.0 * pow(z, 8) + 979028262810.0 / 871614792307.0 * pow(z, 9) + 4953491552399.0 / 1880852972873.0 * pow(z, 10) + 50716798143856.0 / 35736206484587.0 * pow(z, 11) + 4159115221017900.0 / 1322239639929719.0 * pow(z, 12) + 312086709439480.0 / 188891377132817.0 * pow(z, 13) + 126033794576300.0 / 35736206484587.0 * pow(z, 14) + 63805004116464.0 / 35736206484587.0 * pow(z, 15) + 170810053531.0 / 45874462753.0 * pow(z, 16) + 36260306030.0 / 20270111449.0 * pow(z, 17) + 74799599355.0 / 20270111449.0 * pow(z, 18) + 708742840.0 / 431278967.0 * pow(z, 19) + 209181518.0 / 61611281.0 * pow(z, 20) + 4879524.0 / 3624193.0 * pow(z, 21) + 193310.0 / 68381.0 * pow(z, 22) + 62920.0 / 68381.0 * pow(z, 23) + 7215.0 / 3599.0 * pow(z, 24) + 26.0 / 61.0 * pow(z, 25) + 1.0 * pow(z, 26));
                result += c_f[27]() *  2.4440889958054246 * (1.0 / 63.0 * pow(z, 0) + 593.0 / 3843.0 * pow(z, 1) + 8411.0 / 75579.0 * pow(z, 2) + 2032615.0 / 4308003.0 * pow(z, 3) + 1274390.0 / 4308003.0 * pow(z, 4) + 70638178.0 / 76108053.0 * pow(z, 5) + 308015942.0 / 554501529.0 * pow(z, 6) + 682102538.0 / 460518219.0 * pow(z, 7) + 122418663835.0 / 141890780143.0 * pow(z, 8) + 43512234155.0 / 20934705267.0 * pow(z, 9) + 1224984369323.0 / 1036070413497.0 * pow(z, 10) + 1997711464112689.0 / 750460336176327.0 * pow(z, 11) + 3321192077639228.0 / 2251381008528981.0 * pow(z, 12) + 37744123698600580.0 / 11900156759367471.0 * pow(z, 13) + 1277381568322780.0 / 750460336176327.0 * pow(z, 14) + 7990845856450756.0 / 2251381008528981.0 * pow(z, 15) + 111362215393.0 / 60945318441.0 * pow(z, 16) + 8702446831.0 / 2326078363.0 * pow(z, 17) + 122418663835.0 / 67211422173.0 * pow(z, 18) + 1705256345.0 / 460518219.0 * pow(z, 19) + 308015942.0 / 184833843.0 * pow(z, 20) + 777019958.0 / 228324159.0 * pow(z, 21) + 5862194.0 / 4308003.0 * pow(z, 22) + 4065230.0 / 1436001.0 * pow(z, 23) + 210275.0 / 226737.0 * pow(z, 24) + 7709.0 / 3843.0 * pow(z, 25) + 3.0 / 7.0 * pow(z, 26) + 1.0 * pow(z, 27));
                result += c_f[28]() *  2.467859887040803 * (9.0 / 65.0 * pow(z, 0) + 44.0 / 585.0 * pow(z, 1) + 1186.0 / 2745.0 * pow(z, 2) + 2588.0 / 10797.0 * pow(z, 3) + 531607.0 / 615429.0 * pow(z, 4) + 78424.0 / 161955.0 * pow(z, 5) + 76071884.0 / 54362895.0 * pow(z, 6) + 2179805128.0 / 2772507645.0 * pow(z, 7) + 131173565.0 / 65788317.0 * pow(z, 8) + 22600368708.0 / 20270111449.0 * pow(z, 9) + 38826301246.0 / 14953360905.0 * pow(z, 10) + 1062224208364.0 / 740050295355.0 * pow(z, 11) + 1690371238864583.0 / 536043097268805.0 * pow(z, 12) + 1021905254658224.0 / 597305165528097.0 * pow(z, 13) + 15097649479440232.0 / 4181136158696679.0 * pow(z, 14) + 1021905254658224.0 / 536043097268805.0 * pow(z, 15) + 153670112624053.0 / 39222665653815.0 * pow(z, 16) + 34265297044.0 / 17210471985.0 * pow(z, 17) + 1338837974.0 / 332296909.0 * pow(z, 18) + 7533456236.0 / 3881510703.0 * pow(z, 19) + 183642991.0 / 46991655.0 * pow(z, 20) + 94774136.0 / 54362895.0 * pow(z, 21) + 10867412.0 / 3077145.0 * pow(z, 22) + 862664.0 / 615429.0 * pow(z, 23) + 31271.0 / 10797.0 * pow(z, 24) + 2588.0 / 2745.0 * pow(z, 25) + 1186.0 / 585.0 * pow(z, 26) + 28.0 / 65.0 * pow(z, 27) + 1.0 * pow(z, 28));
                result += c_f[29]() *  2.468134811554221 * (1.0 / 67.0 * pow(z, 0) + 631.0 / 4355.0 * pow(z, 1) + 106.0 / 1005.0 * pow(z, 2) + 5470.0 / 12261.0 * pow(z, 3) + 204667.0 / 723399.0 * pow(z, 4) + 12160111.0 / 13744581.0 * pow(z, 5) + 36900388.0 / 68722905.0 * pow(z, 6) + 5191590364.0 / 3642313965.0 * pow(z, 7) + 10793387.0 / 12780049.0 * pow(z, 8) + 3438317377.0 / 1699746517.0 * pow(z, 9) + 93876107050.0 / 79888086299.0 * pow(z, 10) + 3147042234034.0 / 1198321294485.0 * pow(z, 11) + 76973103461267.0 / 51527815662855.0 * pow(z, 12) + 192061812199303.0 / 60361155490773.0 * pow(z, 13) + 9693553214002888.0 / 5492865149660343.0 * pow(z, 14) + 1536494497594424.0 / 422528088435411.0 * pow(z, 15) + 100657135295503.0 / 51527815662855.0 * pow(z, 16) + 1573521117017.0 / 399440431495.0 * pow(z, 17) + 8534191550.0 / 4204636121.0 * pow(z, 18) + 6876634754.0 / 1699746517.0 * pow(z, 19) + 75553709.0 / 38340147.0 * pow(z, 20) + 14276873501.0 / 3642313965.0 * pow(z, 21) + 121244132.0 / 68722905.0 * pow(z, 22) + 48640444.0 / 13744581.0 * pow(z, 23) + 1023335.0 / 723399.0 * pow(z, 24) + 35555.0 / 12261.0 * pow(z, 25) + 318.0 / 335.0 * pow(z, 26) + 8834.0 / 4355.0 * pow(z, 27) + 29.0 / 67.0 * pow(z, 28) + 1.0 * pow(z, 29));
                result += c_f[30]() *  2.4894020435851436 * (3.0 / 23.0 * pow(z, 0) + 110.0 / 1541.0 * pow(z, 1) + 631.0 / 1541.0 * pow(z, 2) + 1060.0 / 4623.0 * pow(z, 3) + 232475.0 / 282003.0 * pow(z, 4) + 7777346.0 / 16638177.0 * pow(z, 5) + 425603885.0 / 316125363.0 * pow(z, 6) + 10542968.0 / 13744581.0 * pow(z, 7) + 32447439775.0 / 16754644239.0 * pow(z, 8) + 323801610.0 / 293941127.0 * pow(z, 9) + 99711203933.0 / 39094169891.0 * pow(z, 10) + 2645599380500.0 / 1837425984877.0 * pow(z, 11) + 17308732287187.0 / 5512277954631.0 * pow(z, 12) + 2574351286330.0 / 1472223304653.0 * pow(z, 13) + 35531435256871055.0 / 9718146034014453.0 * pow(z, 14) + 19387106428005776.0 / 9718146034014453.0 * pow(z, 15) + 960309060996515.0 / 237027952049133.0 * pow(z, 16) + 514870257266.0 / 239664258897.0 * pow(z, 17) + 7867605585085.0 / 1837425984877.0 * pow(z, 18) + 85341915500.0 / 39094169891.0 * pow(z, 19) + 24068221639.0 / 5584881413.0 * pow(z, 20) + 1834875790.0 / 881823381.0 * pow(z, 21) + 1297897591.0 / 316125363.0 * pow(z, 22) + 579863240.0 / 316125363.0 * pow(z, 23) + 60800555.0 / 16638177.0 * pow(z, 24) + 409334.0 / 282003.0 * pow(z, 25) + 13675.0 / 4623.0 * pow(z, 26) + 1484.0 / 1541.0 * pow(z, 27) + 3155.0 / 1541.0 * pow(z, 28) + 10.0 / 23.0 * pow(z, 29) + 1.0 * pow(z, 30));
                result += c_f[31]() *  2.489648995824587 * (1.0 / 71.0 * pow(z, 0) + 223.0 / 1633.0 * pow(z, 1) + 10965.0 / 109411.0 * pow(z, 2) + 46285.0 / 109411.0 * pow(z, 3) + 29645.0 / 109411.0 * pow(z, 4) + 5638353.0 / 6674071.0 * pow(z, 5) + 204330707.0 / 393770189.0 * pow(z, 6) + 10265913025.0 / 7481633591.0 * pow(z, 7) + 6171521445.0 / 7481633591.0 * pow(z, 8) + 779543621345.0 / 396526580323.0 * pow(z, 9) + 460876593529.0 / 396526580323.0 * pow(z, 10) + 1023548199249.0 / 396526580323.0 * pow(z, 11) + 27957880228655.0 / 18636749275181.0 * pow(z, 12) + 59083905157495.0 / 18636749275181.0 * pow(z, 13) + 10129239538344915.0 / 5609661531829481.0 * pow(z, 14) + 847373003221950227.0 / 229996122805008721.0 * pow(z, 15) + 11479804810124237.0 / 5609661531829481.0 * pow(z, 16) + 531755146417455.0 / 130457244926267.0 * pow(z, 17) + 2150606171435.0 / 980881540799.0 * pow(z, 18) + 1705913665415.0 / 396526580323.0 * pow(z, 19) + 879855314919.0 / 396526580323.0 * pow(z, 20) + 1714995966959.0 / 396526580323.0 * pow(z, 21) + 685724605.0 / 325288417.0 * pow(z, 22) + 30797739075.0 / 7481633591.0 * pow(z, 23) + 729752525.0 / 393770189.0 * pow(z, 24) + 24432863.0 / 6674071.0 * pow(z, 25) + 160083.0 / 109411.0 * pow(z, 26) + 323995.0 / 109411.0 * pow(z, 27) + 105995.0 / 109411.0 * pow(z, 28) + 3345.0 / 1633.0 * pow(z, 29) + 31.0 / 71.0 * pow(z, 30) + 1.0 * pow(z, 31));
                result += c_f[32]() *  2.508788609246518 * (9.0 / 73.0 * pow(z, 0) + 352.0 / 5183.0 * pow(z, 1) + 46384.0 / 119209.0 * pow(z, 2) + 1754400.0 / 7987003.0 * pow(z, 3) + 6294760.0 / 7987003.0 * pow(z, 4) + 3604832.0 / 7987003.0 * pow(z, 5) + 631495536.0 / 487207183.0 * pow(z, 6) + 934083232.0 / 1249792339.0 * pow(z, 7) + 1026591302500.0 / 546159252143.0 * pow(z, 8) + 592466058720.0 / 546159252143.0 * pow(z, 9) + 72341648060816.0 / 28946440363579.0 * pow(z, 10) + 41562689161888.0 / 28946440363579.0 * pow(z, 11) + 90072241533912.0 / 28946440363579.0 * pow(z, 12) + 2408678912007200.0 / 1360482697088213.0 * pow(z, 13) + 34977671853237040.0 / 9523378879617491.0 * pow(z, 14) + 842752729590296928.0 / 409505291823552113.0 * pow(z, 15) + 1694746006443900454.0 / 409505291823552113.0 * pow(z, 16) + 21609044348469152.0 / 9523378879617491.0 * pow(z, 17) + 42540411713396400.0 / 9523378879617491.0 * pow(z, 18) + 68819397485920.0 / 28946440363579.0 * pow(z, 19) + 133743631368536.0 / 28946440363579.0 * pow(z, 20) + 68377327330848.0 / 28946440363579.0 * pow(z, 21) + 2494539588304.0 / 546159252143.0 * pow(z, 22) + 1206875304800.0 / 546159252143.0 * pow(z, 23) + 123190956300.0 / 28745223797.0 * pow(z, 24) + 934083232.0 / 487207183.0 * pow(z, 25) + 30071216.0 / 7987003.0 * pow(z, 26) + 11952864.0 / 7987003.0 * pow(z, 27) + 24068200.0 / 7987003.0 * pow(z, 28) + 116960.0 / 119209.0 * pow(z, 29) + 10704.0 / 5183.0 * pow(z, 30) + 32.0 / 73.0 * pow(z, 31) + 1.0 * pow(z, 32));
                result += c_f[33]() *  2.509011642416648 * (1.0 / 75.0 * pow(z, 0) + 707.0 / 5475.0 * pow(z, 1) + 12368.0 / 129575.0 * pow(z, 2) + 719152.0 / 1788135.0 * pow(z, 3) + 6225928.0 / 23961009.0 * pow(z, 4) + 161353288.0 / 199675075.0 * pow(z, 5) + 300433616.0 / 599025225.0 * pow(z, 6) + 48296248432.0 / 36540538725.0 * pow(z, 7) + 23126460420.0 / 28745223797.0 * pow(z, 8) + 164634153668.0 / 86235671391.0 * pow(z, 9) + 34788357296.0 / 30364673025.0 * pow(z, 10) + 96385981261264.0 / 38087421531025.0 * pow(z, 11) + 171102385184744.0 / 114262264593075.0 * pow(z, 12) + 71834079079256.0 / 22852452918615.0 * pow(z, 13) + 917260940713808.0 / 501230467348289.0 * pow(z, 14) + 19886725062716144.0 / 5370326435874525.0 * pow(z, 15) + 3415852669774535266.0 / 1616468257198232025.0 * pow(z, 16) + 7457521898518554.0 / 1790108811958175.0 * pow(z, 17) + 17427957873562352.0 / 7518457010224335.0 * pow(z, 18) + 143668158158512.0 / 31993434086061.0 * pow(z, 19) + 92132053561016.0 / 38087421531025.0 * pow(z, 20) + 530122896936952.0 / 114262264593075.0 * pow(z, 21) + 3162577936.0 / 1320203175.0 * pow(z, 22) + 658536614672.0 / 143726118985.0 * pow(z, 23) + 192720503500.0 / 86235671391.0 * pow(z, 24) + 156962807404.0 / 36540538725.0 * pow(z, 25) + 386271792.0 / 199675075.0 * pow(z, 26) + 2258946032.0 / 599025225.0 * pow(z, 27) + 180551912.0 / 119805045.0 * pow(z, 28) + 359576.0 / 119209.0 * pow(z, 29) + 383408.0 / 388725.0 * pow(z, 30) + 11312.0 / 5475.0 * pow(z, 31) + 11.0 / 25.0 * pow(z, 32) + 1.0 * pow(z, 33));
                result += c_f[34]() *  2.526327908322109 * (9.0 / 77.0 * pow(z, 0) + 34.0 / 525.0 * pow(z, 1) + 22321.0 / 60225.0 * pow(z, 2) + 420512.0 / 1995455.0 * pow(z, 3) + 14845352.0 / 19669485.0 * pow(z, 4) + 4021949488.0 / 9224988465.0 * pow(z, 5) + 2743005896.0 / 2196425825.0 * pow(z, 6) + 1459248992.0 / 2005432275.0 * pow(z, 7) + 205259055836.0 / 112544859273.0 * pow(z, 8) + 336985566120.0 / 316197461767.0 * pow(z, 9) + 81164637758324.0 / 33200733485535.0 * pow(z, 10) + 476193877792.0 / 334011403275.0 * pow(z, 11) + 819280840720744.0 / 266611950717175.0 * pow(z, 12) + 447498545867792.0 / 251376982104765.0 * pow(z, 13) + 45183635740852024.0 / 12317472123133485.0 * pow(z, 14) + 405429335795503136.0 / 192973729929091265.0 * pow(z, 15) + 1732630921089144046.0 / 413515135562338425.0 * pow(z, 16) + 6831705339549070532.0 / 2894605948936368975.0 * pow(z, 17) + 126777872274815418.0 / 27567675704155895.0 * pow(z, 18) + 31186871984269472.0 / 12317472123133485.0 * pow(z, 19) + 1221179344347352.0 / 251376982104765.0 * pow(z, 20) + 7607475279752464.0 / 2932731457888925.0 * pow(z, 21) + 819280840720744.0 / 166003667427675.0 * pow(z, 22) + 15361092832.0 / 6072934605.0 * pow(z, 23) + 53176831634764.0 / 11066911161845.0 * pow(z, 24) + 37442840680.0 / 16077837039.0 * pow(z, 25) + 205259055836.0 / 46124942325.0 * pow(z, 26) + 4377746976.0 / 2196425825.0 * pow(z, 27) + 35659076648.0 / 9224988465.0 * pow(z, 28) + 211681552.0 / 137686395.0 * pow(z, 29) + 873256.0 / 285065.0 * pow(z, 30) + 420512.0 / 421575.0 * pow(z, 31) + 1717.0 / 825.0 * pow(z, 32) + 34.0 / 77.0 * pow(z, 33) + 1.0 * pow(z, 34));
                result += c_f[35]() *  2.5265303303334368 * (1.0 / 79.0 * pow(z, 0) + 745.0 / 6083.0 * pow(z, 1) + 2771.0 / 30415.0 * pow(z, 2) + 170187.0 / 444059.0 * pow(z, 3) + 7866648.0 / 31528189.0 * pow(z, 4) + 2807053736.0 / 3625741735.0 * pow(z, 5) + 305968040.0 / 630973237.0 * pow(z, 6) + 44215042552.0 / 34703528035.0 * pow(z, 7) + 38089526324.0 / 48584939249.0 * pow(z, 8) + 5492521765268.0 / 2963681294189.0 * pow(z, 9) + 89540174723324.0 / 79480543798705.0 * pow(z, 10) + 18828943412812.0 / 7602486798137.0 * pow(z, 11) + 1301072070964536.0 / 874285981785755.0 * pow(z, 12) + 4111785862091880.0 / 1323918772418429.0 * pow(z, 13) + 2438237631193528.0 / 1323918772418429.0 * pow(z, 14) + 1200231365170526456.0 / 324360099242515105.0 * pow(z, 15) + 6583272526753691134.0 / 3048984932879641987.0 * pow(z, 16) + 5848343058799065018.0 / 1385902242218019085.0 * pow(z, 17) + 7357775176960007738.0 / 3048984932879641987.0 * pow(z, 18) + 300057841292631614.0 / 64872019848503021.0 * pow(z, 19) + 17067663418354696.0 / 6619593862092145.0 * pow(z, 20) + 587397980298840.0 / 120356252038039.0 * pow(z, 21) + 100082466997272.0 / 38012433990685.0 * pow(z, 22) + 37657886825624.0 / 7602486798137.0 * pow(z, 23) + 447700873616620.0 / 174857196357151.0 * pow(z, 24) + 71402782948484.0 / 14818406470945.0 * pow(z, 25) + 114268578972.0 / 48584939249.0 * pow(z, 26) + 154752648932.0 / 34703528035.0 * pow(z, 27) + 1267581880.0 / 630973237.0 * pow(z, 28) + 2807053736.0 / 725148347.0 * pow(z, 29) + 243866088.0 / 157640945.0 * pow(z, 30) + 1361496.0 / 444059.0 * pow(z, 31) + 2771.0 / 2765.0 * pow(z, 32) + 12665.0 / 6083.0 * pow(z, 33) + 35.0 / 79.0 * pow(z, 34) + 1.0 * pow(z, 35));
                result += c_f[36]() *  2.5422721046282533 * (1.0 / 9.0 * pow(z, 0) + 44.0 / 711.0 * pow(z, 1) + 19370.0 / 54747.0 * pow(z, 2) + 11084.0 / 54747.0 * pow(z, 3) + 964393.0 / 1332177.0 * pow(z, 4) + 66429472.0 / 157640945.0 * pow(z, 5) + 5614107472.0 / 4661667945.0 * pow(z, 6) + 174838880.0 / 246902571.0 * pow(z, 7) + 110537606380.0 / 62466350463.0 * pow(z, 8) + 152358105296.0 / 145754817747.0 * pow(z, 9) + 318566262385544.0 / 133365658238505.0 * pow(z, 10) + 11102981665692176.0 / 7868573836071795.0 * pow(z, 11) + 18828943412812.0 / 6220216471203.0 * pow(z, 12) + 44481096443232.0 / 24979599479593.0 * pow(z, 13) + 14489150180704720.0 / 3971756317255287.0 * pow(z, 14) + 126788356822063456.0 / 59576344758829305.0 * pow(z, 15) + 12302371492997896174.0 / 2919240893182635945.0 * pow(z, 16) + 66607227917743227944.0 / 27440864395916777883.0 * pow(z, 17) + 3898895372532710012.0 / 831541345330811451.0 * pow(z, 18) + 1549005300412633208.0 / 583848178636527189.0 * pow(z, 19) + 300057841292631614.0 / 59576344758829305.0 * pow(z, 20) + 165800158921159904.0 / 59576344758829305.0 * pow(z, 21) + 391598653532560.0 / 74938798438779.0 * pow(z, 22) + 44481096443232.0 / 15896108759741.0 * pow(z, 23) + 357749924843428.0 / 68422381183233.0 * pow(z, 24) + 358160698893296.0 / 133365658238505.0 * pow(z, 25) + 10985043530536.0 / 2186322266205.0 * pow(z, 26) + 152358105296.0 / 62466350463.0 * pow(z, 27) + 287397776588.0 / 62466350463.0 * pow(z, 28) + 174838880.0 / 84757599.0 * pow(z, 29) + 5614107472.0 / 1418768505.0 * pow(z, 30) + 3496288.0 / 2220295.0 * pow(z, 31) + 56729.0 / 18249.0 * pow(z, 32) + 55420.0 / 54747.0 * pow(z, 33) + 1490.0 / 711.0 * pow(z, 34) + 4.0 / 9.0 * pow(z, 35) + 1.0 * pow(z, 36));
                result += c_f[37]() *  2.5424566414923317 * (1.0 / 83.0 * pow(z, 0) + 29.0 / 249.0 * pow(z, 1) + 1714.0 / 19671.0 * pow(z, 2) + 237590.0 / 649143.0 * pow(z, 3) + 155737.0 / 649143.0 * pow(z, 4) + 58673137.0 / 78979065.0 * pow(z, 5) + 7890989008.0 / 16822540845.0 * pow(z, 6) + 475587058576.0 / 386918439435.0 * pow(z, 7) + 439871046412.0 / 576078565381.0 * pow(z, 8) + 9326979434836.0 / 5184707088429.0 * pow(z, 9) + 2606115384728.0 / 2356685040195.0 * pow(z, 10) + 17973950688728.0 / 7424111089065.0 * pow(z, 11) + 12502450904018764.0 / 8481709459661805.0 * pow(z, 12) + 5192111834721028.0 / 1696341891932361.0 * pow(z, 13) + 11467391667897328.0 / 6219920270418657.0 * pow(z, 14) + 249380714162993648.0 / 67737487876477155.0 * pow(z, 15) + 10823491457523898462.0 / 4944836614982832315.0 * pow(z, 16) + 16330823181283945214.0 / 3845984033875536245.0 * pow(z, 17) + 808151923495181049964.0 / 325370249265870366327.0 * pow(z, 18) + 32661646362567890428.0 / 6922771260975965241.0 * pow(z, 19) + 4456731776627487602.0 / 1648278871660944105.0 * pow(z, 20) + 31172589270374206.0 / 6157953443316105.0 * pow(z, 21) + 11467391667897328.0 / 4056469741577385.0 * pow(z, 22) + 2966921048412016.0 / 565447297310787.0 * pow(z, 23) + 4808634963084140.0 / 1696341891932361.0 * pow(z, 24) + 116830679476732.0 / 22272333267195.0 * pow(z, 25) + 7818346154184.0 / 2880392826905.0 * pow(z, 26) + 130577712087704.0 / 25923535442145.0 * pow(z, 27) + 12756260345948.0 / 5184707088429.0 * pow(z, 28) + 118896764644.0 / 25794562629.0 * pow(z, 29) + 34945808464.0 / 16822540845.0 * pow(z, 30) + 938770192.0 / 236937195.0 * pow(z, 31) + 155737.0 / 98355.0 * pow(z, 32) + 2019515.0 / 649143.0 * pow(z, 33) + 59990.0 / 59013.0 * pow(z, 34) + 174.0 / 83.0 * pow(z, 35) + 37.0 / 83.0 * pow(z, 36) + 1.0 * pow(z, 37));
                result += c_f[38]() *  2.5568294389699044 * (9.0 / 85.0 * pow(z, 0) + 418.0 / 7055.0 * pow(z, 1) + 7163.0 / 21165.0 * pow(z, 2) + 65132.0 / 334407.0 * pow(z, 3) + 451421.0 / 649143.0 * pow(z, 4) + 6614242.0 / 16228575.0 * pow(z, 5) + 459031013.0 / 394895325.0 * pow(z, 6) + 57955667168.0 / 84112704225.0 * pow(z, 7) + 132884619308.0 / 77383687887.0 * pow(z, 8) + 2949723487704.0 / 2880392826905.0 * pow(z, 9) + 302303862858508.0 / 129617677210725.0 * pow(z, 10) + 180588466071152.0 / 129617677210725.0 * pow(z, 11) + 10044266561348.0 / 3374595949575.0 * pow(z, 12) + 15048198825651544.0 / 8481709459661805.0 * pow(z, 13) + 30672727897553636.0 / 8481709459661805.0 * pow(z, 14) + 333228910820075296.0 / 155498006760466425.0 * pow(z, 15) + 1428438061271853322.0 / 338687439382385775.0 * pow(z, 16) + 1040328531858473534524.0 / 420311112273540746775.0 * pow(z, 17) + 310285640444394959066.0 / 65381728575884116165.0 * pow(z, 18) + 1616303846990362099928.0 / 588435557182957045485.0 * pow(z, 19) + 2171999483110764713462.0 / 420311112273540746775.0 * pow(z, 20) + 24193686787406361268.0 / 8241394358304720525.0 * pow(z, 21) + 34839952713947642.0 / 6390329044950675.0 * pow(z, 22) + 25632993140005792.0 / 8481709459661805.0 * pow(z, 23) + 15750860271716732.0 / 2827236486553935.0 * pow(z, 24) + 2149742689378792.0 / 718788937259475.0 * pow(z, 25) + 10044266561348.0 / 1825601087475.0 * pow(z, 26) + 40778040725744.0 / 14401964134525.0 * pow(z, 27) + 135515524729676.0 / 25923535442145.0 * pow(z, 28) + 983241162568.0 / 386918439435.0 * pow(z, 29) + 132884619308.0 / 28037568075.0 * pow(z, 30) + 2519811616.0 / 1184685975.0 * pow(z, 31) + 65575859.0 / 16228575.0 * pow(z, 32) + 348118.0 / 216381.0 * pow(z, 33) + 3159947.0 / 1003221.0 * pow(z, 34) + 65132.0 / 63495.0 * pow(z, 35) + 14877.0 / 7055.0 * pow(z, 36) + 38.0 / 85.0 * pow(z, 37) + 1.0 * pow(z, 38));
                result += c_f[39]() *  2.5569983571109383 * (1.0 / 87.0 * pow(z, 0) + 821.0 / 7395.0 * pow(z, 1) + 589.0 / 7055.0 * pow(z, 2) + 386935.0 / 1104813.0 * pow(z, 3) + 1185847.0 / 5134131.0 * pow(z, 4) + 1221757.0 / 1711377.0 * pow(z, 5) + 58274083.0 / 128353275.0 * pow(z, 6) + 11120627819.0 / 9369789075.0 * pow(z, 7) + 379017548.0 / 509773965.0 * pow(z, 8) + 5344368232868.0 / 3060173111895.0 * pow(z, 9) + 222286274145628.0 / 205031598496965.0 * pow(z, 10) + 808104195006724.0 / 341719330828275.0 * pow(z, 11) + 1493121770367812.0 / 1025157992484825.0 * pow(z, 12) + 37656009267543716.0 / 12506927508314865.0 * pow(z, 13) + 452148626215611508.0 / 245969574330192345.0 * pow(z, 14) + 2694160320227172124.0 / 737908722990577035.0 * pow(z, 15) + 3552083737635722.0 / 1610451163226925.0 * pow(z, 16) + 92369607954543025574.0 / 21727312399166990475.0 * pow(z, 17) + 337020297290174408062.0 / 132971151882901981707.0 * pow(z, 18) + 22233632755904692277146.0 / 4653990315901569359745.0 * pow(z, 19) + 124165372685853729286.0 / 44323717294300660569.0 * pow(z, 20) + 1016065687499973281314.0 / 195545811592502914275.0 * pow(z, 21) + 208946102213866.0 / 70019615792475.0 * pow(z, 22) + 1347080160113586062.0 / 245969574330192345.0 * pow(z, 23) + 452148626215611508.0 / 147581744598115407.0 * pow(z, 24) + 69932588639724044.0 / 12506927508314865.0 * pow(z, 25) + 114855520797524.0 / 37968814536475.0 * pow(z, 26) + 5656729365047068.0 / 1025157992484825.0 * pow(z, 27) + 20207843104148.0 / 7070055120585.0 * pow(z, 28) + 5344368232868.0 / 1020057703965.0 * pow(z, 29) + 11749543988.0 / 4587965685.0 * pow(z, 30) + 44482511276.0 / 9369789075.0 * pow(z, 31) + 91573559.0 / 42784425.0 * pow(z, 32) + 20769869.0 / 5134131.0 * pow(z, 33) + 8300929.0 / 5134131.0 * pow(z, 34) + 386935.0 / 122757.0 * pow(z, 35) + 21793.0 / 21165.0 * pow(z, 36) + 15599.0 / 7395.0 * pow(z, 37) + 13.0 / 29.0 * pow(z, 38) + 1.0 * pow(z, 39));
                result += c_f[40]() *  2.570173398494093 * (9.0 / 89.0 * pow(z, 0) + 440.0 / 7743.0 * pow(z, 1) + 42692.0 / 131631.0 * pow(z, 2) + 23560.0 / 125579.0 * pow(z, 3) + 3869350.0 / 5784021.0 * pow(z, 4) + 180248744.0 / 456937659.0 * pow(z, 5) + 171045980.0 / 152312553.0 * pow(z, 6) + 1531775896.0 / 2284688295.0 * pow(z, 7) + 55603139095.0 / 33356449107.0 * pow(z, 8) + 3032140384.0 / 3024658859.0 * pow(z, 9) + 21377472931472.0 / 9391565757195.0 * pow(z, 10) + 5011545089828704.0 / 3649562453245977.0 * pow(z, 11) + 17778292290147928.0 / 6082604088743295.0 * pow(z, 12) + 6431909164661344.0 / 3649562453245977.0 * pow(z, 13) + 796155624513781424.0 / 222623309648004597.0 * pow(z, 14) + 47023457126423596832.0 / 21891292115387118705.0 * pow(z, 15) + 55230286564657028542.0 / 13134775269232271223.0 * pow(z, 16) + 71877459161569904.0 / 28666030705439265.0 * pow(z, 17) + 369478431818172102296.0 / 77349232141034486091.0 * pow(z, 18) + 374688444060040146160.0 / 132971151882901981707.0 * pow(z, 19) + 311270858582665691880044.0 / 59172162587891381859615.0 * pow(z, 20) + 7972094554468939280.0 / 2607277487900038857.0 * pow(z, 21) + 369478431818172102296.0 / 65673876346161356115.0 * pow(z, 22) + 18387256994820208.0 / 5733206141087853.0 * pow(z, 23) + 25594523042158135178.0 / 4378258423077423741.0 * pow(z, 24) + 3617189009724892064.0 / 1113116548240022985.0 * pow(z, 25) + 21517719581453552.0 / 3649562453245977.0 * pow(z, 26) + 6431909164661344.0 / 2027534696247765.0 * pow(z, 27) + 21010709070174824.0 / 3649562453245977.0 * pow(z, 28) + 161662744833184.0 / 54471081391731.0 * pow(z, 29) + 21377472931472.0 / 3947179810995.0 * pow(z, 30) + 3032140384.0 / 1150222383.0 * pow(z, 31) + 11120627819.0 / 2284688295.0 * pow(z, 32) + 332994760.0 / 152312553.0 * pow(z, 33) + 1881505780.0 / 456937659.0 * pow(z, 34) + 9486776.0 / 5784021.0 * pow(z, 35) + 11608050.0 / 3641791.0 * pow(z, 36) + 4712.0 / 4539.0 * pow(z, 37) + 16420.0 / 7743.0 * pow(z, 38) + 40.0 / 89.0 * pow(z, 39) + 1.0 * pow(z, 40));
                result += c_f[41]() *  2.570328597515912 * (1.0 / 91.0 * pow(z, 0) + 859.0 / 8099.0 * pow(z, 1) + 18820.0 / 234871.0 * pow(z, 2) + 1340540.0 / 3992807.0 * pow(z, 3) + 4340170.0 / 19494293.0 * pow(z, 4) + 3090578.0 / 4498683.0 * pow(z, 5) + 290224468.0 / 660021063.0 * pow(z, 6) + 757056140.0 / 660021063.0 * pow(z, 7) + 477816617.0 / 660021063.0 * pow(z, 8) + 18856475197.0 / 11118816369.0 * pow(z, 9) + 612294228304.0 / 576554354415.0 * pow(z, 10) + 118475082381616.0 / 51313337542935.0 * pow(z, 11) + 75977275401272.0 / 52892209467333.0 * pow(z, 12) + 156437680448872.0 / 52892209467333.0 * pow(z, 13) + 8794513819791824.0 / 4813191061527303.0 * pow(z, 14) + 1767478251997562512.0 / 489341091255275805.0 * pow(z, 15) + 63871119039352724458.0 / 28871124384061272495.0 * pow(z, 16) + 24483484459704940550.0 / 5774224876812254499.0 * pow(z, 17) + 3423536899728300392.0 / 1332513433110520269.0 * pow(z, 18) + 4416396882682987805336.0 / 918101755413148465341.0 * pow(z, 19) + 32058549484045530870868.0 / 11148378458588231364855.0 * pow(z, 20) + 24290182854756432929348.0 / 4590508777065742326705.0 * pow(z, 21) + 4144281510197416264.0 / 1332513433110520269.0 * pow(z, 22) + 97933937838819762200.0 / 17322674630436763497.0 * pow(z, 23) + 18785623246868448370.0 / 5774224876812254499.0 * pow(z, 24) + 220934781499695314.0 / 37641622404251985.0 * pow(z, 25) + 8794513819791824.0 / 2673995034181835.0 * pow(z, 26) + 312875360897744.0 / 52892209467333.0 * pow(z, 27) + 75977275401272.0 / 23710300795701.0 * pow(z, 28) + 59237541190808.0 / 10262667508587.0 * pow(z, 29) + 1725556461584.0 / 576554354415.0 * pow(z, 30) + 301703603152.0 / 55594081845.0 * pow(z, 31) + 5255982787.0 / 1980063189.0 * pow(z, 32) + 3217488595.0 / 660021063.0 * pow(z, 33) + 1451122340.0 / 660021063.0 * pow(z, 34) + 6181156.0 / 1499561.0 * pow(z, 35) + 32117258.0 / 19494293.0 * pow(z, 36) + 12735130.0 / 3992807.0 * pow(z, 37) + 18820.0 / 18067.0 * pow(z, 38) + 17180.0 / 8099.0 * pow(z, 39) + 41.0 / 91.0 * pow(z, 40) + 1.0 * pow(z, 41));
                result += c_f[42]() *  2.582449679876761 * (3.0 / 31.0 * pow(z, 0) + 22.0 / 403.0 * pow(z, 1) + 859.0 / 2759.0 * pow(z, 2) + 188200.0 / 1040143.0 * pow(z, 3) + 670270.0 / 1040143.0 * pow(z, 4) + 32985292.0 / 86331869.0 * pow(z, 5) + 21634046.0 / 19922739.0 * pow(z, 6) + 13350325528.0 / 20460652953.0 * pow(z, 7) + 33121206125.0 / 20460652953.0 * pow(z, 8) + 6689432638.0 / 6820217651.0 * pow(z, 9) + 131995326379.0 / 59428156455.0 * pow(z, 10) + 779283563296.0 / 576554354415.0 * pow(z, 11) + 4561290671692216.0 / 1590713463830985.0 * pow(z, 12) + 37228864946623280.0 / 21315560415335199.0 * pow(z, 13) + 5788194176608264.0 / 1639658493487323.0 * pow(z, 14) + 17589027639583648.0 / 8198292467436615.0 * pow(z, 15) + 9058326041487507874.0 / 2167081975559078565.0 * pow(z, 16) + 323112719846137311964.0 / 127857836557985635335.0 * pow(z, 17) + 122417422298524702750.0 / 25571567311597127067.0 * pow(z, 18) + 16937498346024222992.0 / 5901130918060875477.0 * pow(z, 19) + 108201723625733201230732.0 / 20329396012719716018265.0 * pow(z, 20) + 64117098968091061741736.0 / 20329396012719716018265.0 * pow(z, 21) + 2208198441341493902668.0 / 383573509673956906005.0 * pow(z, 22) + 19820476787900686480.0 / 5901130918060875477.0 * pow(z, 23) + 465186204734393870450.0 / 76714701934791381201.0 * pow(z, 24) + 7514249298747379348.0 / 2167081975559078565.0 * pow(z, 25) + 220934781499695314.0 / 35525934025558665.0 * pow(z, 26) + 123123193477085536.0 / 35525934025558665.0 * pow(z, 27) + 10168449229176680.0 / 1639658493487323.0 * pow(z, 28) + 1063681855617808.0 / 318142692766197.0 * pow(z, 29) + 9537244131720088.0 / 1590713463830985.0 * pow(z, 30) + 779283563296.0 / 251734999815.0 * pow(z, 31) + 131995326379.0 / 23608445715.0 * pow(z, 32) + 167235815950.0 / 61381958859.0 * pow(z, 33) + 102013314865.0 / 20460652953.0 * pow(z, 34) + 580448936.0 / 258995607.0 * pow(z, 35) + 27815202.0 / 6640913.0 * pow(z, 36) + 1736068.0 / 1040143.0 * pow(z, 37) + 3351350.0 / 1040143.0 * pow(z, 38) + 37640.0 / 35867.0 * pow(z, 39) + 859.0 / 403.0 * pow(z, 40) + 14.0 / 31.0 * pow(z, 41) + 1.0 * pow(z, 42));
                result += c_f[43]() *  2.5825927637719244 * (1.0 / 95.0 * pow(z, 0) + 299.0 / 2945.0 * pow(z, 1) + 2949.0 / 38285.0 * pow(z, 2) + 219701.0 / 681473.0 * pow(z, 3) + 136970.0 / 637507.0 * pow(z, 4) + 3442722.0 / 5200715.0 * pow(z, 5) + 183957914.0 / 431659345.0 * pow(z, 6) + 478920978.0 / 431659345.0 * pow(z, 7) + 4808049965.0 / 6820217651.0 * pow(z, 8) + 101123995385.0 / 61381958859.0 * pow(z, 9) + 1595434935283.0 / 1534548971475.0 * pow(z, 10) + 84153839498197.0 / 37340691639225.0 * pow(z, 11) + 865552074855944.0 / 611812870704225.0 * pow(z, 12) + 4617273927251048.0 / 1590713463830985.0 * pow(z, 13) + 12873445648801784.0 / 7105186805111733.0 * pow(z, 14) + 1900603611630909832.0 / 532889010383379975.0 * pow(z, 15) + 1178171235757099538.0 / 532889010383379975.0 * pow(z, 16) + 15229039004211589954.0 / 3611803292598464275.0 * pow(z, 17) + 993821643998116876442.0 / 383573509673956906005.0 * pow(z, 18) + 7029326310664449841598.0 / 1457579336761036242819.0 * pow(z, 19) + 428904621397755304732.0 / 146343306903718498275.0 * pow(z, 20) + 10343280843412745178661276.0 / 1931292621208373021735175.0 * pow(z, 21) + 1409258041735481715548.0 / 439029920711155494825.0 * pow(z, 22) + 14058652621328899683196.0 / 2429298894601727071365.0 * pow(z, 23) + 261532011578451809590.0 / 76714701934791381201.0 * pow(z, 24) + 15229039004211589954.0 / 2500479202568167575.0 * pow(z, 25) + 207912571015958742.0 / 59209890042597775.0 * pow(z, 26) + 3326056320354092206.0 / 532889010383379975.0 * pow(z, 27) + 12873445648801784.0 / 3675096623333655.0 * pow(z, 28) + 659610561035864.0 / 106047564255399.0 * pow(z, 29) + 865552074855944.0 / 256566687714675.0 * pow(z, 30) + 673230715985576.0 / 112022074917675.0 * pow(z, 31) + 1595434935283.0 / 511516323825.0 * pow(z, 32) + 343821584309.0 / 61381958859.0 * pow(z, 33) + 168281748775.0 / 61381958859.0 * pow(z, 34) + 2155144401.0 / 431659345.0 * pow(z, 35) + 972348974.0 / 431659345.0 * pow(z, 36) + 21803906.0 / 5200715.0 * pow(z, 37) + 82182.0 / 49039.0 * pow(z, 38) + 2197010.0 / 681473.0 * pow(z, 39) + 40303.0 / 38285.0 * pow(z, 40) + 6279.0 / 2945.0 * pow(z, 41) + 43.0 / 95.0 * pow(z, 42) + 1.0 * pow(z, 43));
                result += c_f[44]() *  2.593781542064307 * (9.0 / 97.0 * pow(z, 0) + 484.0 / 9215.0 * pow(z, 1) + 85514.0 / 285665.0 * pow(z, 2) + 129756.0 / 742729.0 * pow(z, 3) + 41084087.0 / 66102881.0 * pow(z, 4) + 1205336.0 / 3254641.0 * pow(z, 5) + 530179188.0 / 504469355.0 * pow(z, 6) + 26595058424.0 / 41870956465.0 * pow(z, 7) + 13170326895.0 / 8374191293.0 * pow(z, 8) + 634662595380.0 / 661561112147.0 * pow(z, 9) + 444945579694.0 / 205312069287.0 * pow(z, 10) + 6381739741132.0 / 4801653233325.0 * pow(z, 11) + 10182614579281837.0 / 3622047089004825.0 * pow(z, 12) + 266590039055630752.0 / 154299205991605545.0 * pow(z, 13) + 536922996683193296.0 / 154299205991605545.0 * pow(z, 14) + 566431608547278496.0 / 265078123113783885.0 * pow(z, 15) + 214293057211385083558.0 / 51690234007187857575.0 * pow(z, 16) + 131123528120731313288.0 / 51690234007187857575.0 * pow(z, 17) + 335038858092654978988.0 / 70068983876410206935.0 * pow(z, 18) + 108169639988847668446424.0 / 37206630438373819882485.0 * pow(z, 19) + 3788806881448138464621322.0 / 706925978329102577767215.0 * pow(z, 20) + 45831522400788709705648.0 / 14195300769660694332675.0 * pow(z, 21) + 20686561686825490357322552.0 / 3534629891645512888836075.0 * pow(z, 22) + 29655690965216223927184.0 / 8517180461796416599605.0 * pow(z, 23) + 77322589417308948257578.0 / 12402210146124606627495.0 * pow(z, 24) + 2301481701890375924392.0 / 630620854887691862415.0 * pow(z, 25) + 335038858092654978988.0 / 51690234007187857575.0 * pow(z, 26) + 21345690624305097512.0 / 5743359334131984175.0 * pow(z, 27) + 5226659931985002038.0 / 795234369341351655.0 * pow(z, 28) + 566431608547278496.0 / 154299205991605545.0 * pow(z, 29) + 333762943884147184.0 / 51433068663868515.0 * pow(z, 30) + 38084291293661536.0 / 10866141267014475.0 * pow(z, 31) + 925692234480167.0 / 148851250233075.0 * pow(z, 32) + 6381739741132.0 / 1984683336441.0 * pow(z, 33) + 34260809636438.0 / 5954050009323.0 * pow(z, 34) + 211554198460.0 / 75367721637.0 * pow(z, 35) + 213359295699.0 / 41870956465.0 * pow(z, 36) + 1156306888.0 / 504469355.0 * pow(z, 37) + 429192676.0 / 100893871.0 * pow(z, 38) + 3616008.0 / 2132351.0 * pow(z, 39) + 2416711.0 / 742729.0 * pow(z, 40) + 302764.0 / 285665.0 * pow(z, 41) + 19734.0 / 9215.0 * pow(z, 42) + 44.0 / 97.0 * pow(z, 43) + 1.0 * pow(z, 44));
                result += c_f[45]() *  2.5939138744815287 * (1.0 / 99.0 * pow(z, 0) + 85.0 / 873.0 * pow(z, 1) + 410.0 / 5529.0 * pow(z, 2) + 159430.0 / 514197.0 * pow(z, 3) + 1387505.0 / 6684561.0 * pow(z, 4) + 68719.0 / 107601.0 * pow(z, 5) + 375306260.0 / 908044839.0 * pow(z, 6) + 975346300.0 / 908044839.0 * pow(z, 7) + 5750040565.0 / 8374191293.0 * pow(z, 8) + 1086028678355.0 / 678309494733.0 * pow(z, 9) + 54522289283638.0 / 53586450083907.0 * pow(z, 10) + 432194659508290.0 / 196483650307659.0 * pow(z, 11) + 59148587701.0 / 42495202287.0 * pow(z, 12) + 11135712077601335.0 / 3911810856125211.0 * pow(z, 13) + 166022440594908400.0 / 92579523594963327.0 * pow(z, 14) + 75135396196060720.0 / 21364505444991537.0 * pow(z, 15) + 3153760866736849370.0 / 1431421864814432979.0 * pow(z, 16) + 279011063923699826.0 / 66697076138306913.0 * pow(z, 17) + 583407147515665060.0 / 224198605332381069.0 * pow(z, 18) + 5469474622953562949140.0 / 1135117538797845352347.0 * pow(z, 19) + 66307649255441491036174.0 / 22323978263024291929491.0 * pow(z, 20) + 6864605520245066446972130.0 / 1272466760992384639980987.0 * pow(z, 21) + 46019176815546718932450056.0 / 13997134370916231039790857.0 * pow(z, 22) + 2496220189180024162535320.0 / 424155586997461546660329.0 * pow(z, 23) + 236813033055148182272050.0 / 66971934789072875788473.0 * pow(z, 24) + 546947462295356294914.0 / 87316733753680411719.0 * pow(z, 25) + 30705639342929740.0 / 8303652049347447.0 * pow(z, 26) + 3906154894931797564.0 / 600273685244762217.0 * pow(z, 27) + 185515345102167610.0 / 49359374648773551.0 * pow(z, 28) + 46959622622537950.0 / 7121501814997179.0 * pow(z, 29) + 33204488118981680.0 / 8959308734996451.0 * pow(z, 30) + 25453056177374480.0 / 3911810856125211.0 * pow(z, 31) + 650634464711.0 / 184145876577.0 * pow(z, 32) + 3673654605820465.0 / 589450950922977.0 * pow(z, 33) + 173480011357030.0 / 53586450083907.0 * pow(z, 34) + 434411471342.0 / 75367721637.0 * pow(z, 35) + 212751500905.0 / 75367721637.0 * pow(z, 36) + 4632894925.0 / 908044839.0 * pow(z, 37) + 53615180.0 / 23283201.0 * pow(z, 38) + 1374380.0 / 322803.0 * pow(z, 39) + 11377541.0 / 6684561.0 * pow(z, 40) + 558005.0 / 171399.0 * pow(z, 41) + 17630.0 / 16587.0 * pow(z, 42) + 1870.0 / 873.0 * pow(z, 43) + 5.0 / 11.0 * pow(z, 44) + 1.0 * pow(z, 45));
                result += c_f[46]() *  2.6042739671849557 * (9.0 / 101.0 * pow(z, 0) + 46.0 / 909.0 * pow(z, 1) + 25415.0 / 88173.0 * pow(z, 2) + 94300.0 / 558429.0 * pow(z, 3) + 31168565.0 / 51933897.0 * pow(z, 4) + 12765046.0 / 35533719.0 * pow(z, 5) + 11063759.0 / 10867701.0 * pow(z, 6) + 56724860440.0 / 91712528739.0 * pow(z, 7) + 140206030625.0 / 91712528739.0 * pow(z, 8) + 793505597970.0 / 845793320593.0 * pow(z, 9) + 4995731920433.0 / 2362388240277.0 * pow(z, 10) + 228002300640668.0 / 174588111563697.0 * pow(z, 11) + 4970238584345335.0 / 1804077152824869.0 * pow(z, 12) + 95229226198610.0 / 55796200602831.0 * pow(z, 13) + 1353784425434105155.0 / 395092896468646311.0 * pow(z, 14) + 1527406453473157280.0 / 719271683314715079.0 * pow(z, 15) + 8856584826610657370.0 / 2157815049944145237.0 * pow(z, 16) + 366949352612087532580.0 / 144573608346257730879.0 * pow(z, 17) + 32086272351225479990.0 / 6736404689968998213.0 * pow(z, 18) + 66385592259414097880.0 / 22644059138570487969.0 * pow(z, 19) + 616409790006866544368078.0 / 114646871418582380587047.0 * pow(z, 20) + 7407511673965035141469724.0 / 2254721804565453484878591.0 * pow(z, 21) + 760723102652612363532638770.0 / 128519142860230848638079687.0 * pow(z, 22) + 460191768155467189324500560.0 / 128519142860230848638079687.0 * pow(z, 23) + 14353266087785138934578090.0 / 2254721804565453484878591.0 * pow(z, 24) + 435735980821472655380572.0 / 114646871418582380587047.0 * pow(z, 25) + 12579791632793194783022.0 / 1879456908501350501427.0 * pow(z, 26) + 9887215868423376280.0 / 2516006570952276441.0 * pow(z, 27) + 32086272351225479990.0 / 4663664785363152609.0 * pow(z, 28) + 8533705874699710060.0 / 2157815049944145237.0 * pow(z, 29) + 4968328073464515110.0 / 719271683314715079.0 * pow(z, 30) + 1527406453473157280.0 / 395092896468646311.0 * pow(z, 31) + 36588768254975815.0 / 5412231458474607.0 * pow(z, 32) + 68020875856150.0 / 18598733534277.0 * pow(z, 33) + 34791670090417345.0 / 5412231458474607.0 * pow(z, 34) + 228002300640668.0 / 68509258968033.0 * pow(z, 35) + 4995731920433.0 / 845793320593.0 * pow(z, 36) + 264501865990.0 / 91712528739.0 * pow(z, 37) + 476700504125.0 / 91712528739.0 * pow(z, 38) + 2466298280.0 / 1054166997.0 * pow(z, 39) + 1580537.0 / 366327.0 * pow(z, 40) + 89355322.0 / 51933897.0 * pow(z, 41) + 1833445.0 / 558429.0 * pow(z, 42) + 94300.0 / 88173.0 * pow(z, 43) + 1955.0 / 909.0 * pow(z, 44) + 46.0 / 101.0 * pow(z, 45) + 1.0 * pow(z, 46));
                result += c_f[47]() *  2.604396714762013 * (1.0 / 103.0 * pow(z, 0) + 973.0 / 10403.0 * pow(z, 1) + 2231.0 / 31209.0 * pow(z, 2) + 904015.0 / 3027273.0 * pow(z, 3) + 11546345.0 / 57518187.0 * pow(z, 4) + 57889781.0 / 93845463.0 * pow(z, 5) + 37643893.0 / 93845463.0 * pow(z, 6) + 8692639573.0 / 8352246207.0 * pow(z, 7) + 54011633585.0 / 80738380001.0 * pow(z, 8) + 125689940345.0 / 80738380001.0 * pow(z, 9) + 6671292962411.0 / 6701285540083.0 * pow(z, 10) + 388468039666001.0 / 180934709582241.0 * pow(z, 11) + 19559051693232571.0 / 14293842056997039.0 * pow(z, 12) + 39887123011201235.0 / 14293842056997039.0 * pow(z, 13) + 25333947877058675.0 / 14293842056997039.0 * pow(z, 14) + 3614543648030409955.0 / 1043450470160783847.0 * pow(z, 15) + 162291192781199714950.0 / 74084983381415653137.0 * pow(z, 16) + 1401464144620391810.0 / 338287595348929923.0 * pow(z, 17) + 4310050862437000637270.0 / 1654564628851616253393.0 * pow(z, 18) + 78653546484799769870.0 / 16381828008431844093.0 * pow(z, 19) + 153346608474365402474.0 / 51172101923245863507.0 * pow(z, 20) + 1639126778244326326289554.0 / 302785327079845774370919.0 * pow(z, 21) + 59794382192783366038919974.0 / 17864334297710900687884221.0 * pow(z, 22) + 69685827888041149985133770.0 / 11704219022638176312751731.0 * pow(z, 23) + 64993893687808006564043450.0 / 17864334297710900687884221.0 * pow(z, 24) + 1937149828834203840160382.0 / 302785327079845774370919.0 * pow(z, 25) + 21906658353480771782.0 / 5685789102582873723.0 * pow(z, 26) + 110114965078719677818.0 / 16381828008431844093.0 * pow(z, 27) + 226844782233526349330.0 / 57053952719021250117.0 * pow(z, 28) + 7007320723101959050.0 / 1014862786046789769.0 * pow(z, 29) + 9546540751835277350.0 / 2389838173594053327.0 * pow(z, 30) + 7229087296060819910.0 / 1043450470160783847.0 * pow(z, 31) + 55734685329529085.0 / 14293842056997039.0 * pow(z, 32) + 96868727312917285.0 / 14293842056997039.0 * pow(z, 33) + 52658985327933845.0 / 14293842056997039.0 * pow(z, 34) + 388468039666001.0 / 60311569860747.0 * pow(z, 35) + 22439803600837.0 / 6701285540083.0 * pow(z, 36) + 477621773311.0 / 80738380001.0 * pow(z, 37) + 702151236605.0 / 242215140003.0 * pow(z, 38) + 43463197865.0 / 8352246207.0 * pow(z, 39) + 220485659.0 / 93845463.0 * pow(z, 40) + 405228467.0 / 93845463.0 * pow(z, 41) + 99298567.0 / 57518187.0 * pow(z, 42) + 9944165.0 / 3027273.0 * pow(z, 43) + 11155.0 / 10403.0 * pow(z, 44) + 22379.0 / 10403.0 * pow(z, 45) + 47.0 / 103.0 * pow(z, 46) + 1.0 * pow(z, 47));
                result += c_f[48]() *  2.614016928281532 * (3.0 / 35.0 * pow(z, 0) + 176.0 / 3605.0 * pow(z, 1) + 14456.0 / 52015.0 * pow(z, 2) + 35696.0 / 218463.0 * pow(z, 3) + 1756372.0 / 3027273.0 * pow(z, 4) + 36948304.0 / 105954555.0 * pow(z, 5) + 463118248.0 / 469227315.0 * pow(z, 6) + 1978993232.0 / 3284591205.0 * pow(z, 7) + 86926395730.0 / 58465723449.0 * pow(z, 8) + 74073097488.0 / 80738380001.0 * pow(z, 9) + 201103904552.0 / 97442872415.0 * pow(z, 10) + 1386242693488.0 / 1080852506465.0 * pow(z, 11) + 17092593745304044.0 / 6332714835378435.0 * pow(z, 12) + 24072679007055472.0 / 14293842056997039.0 * pow(z, 13) + 337331097466159016.0 / 100056894398979273.0 * pow(z, 14) + 210778446337128176.0 / 100056894398979273.0 * pow(z, 15) + 4234179701978480233.0 / 1043450470160783847.0 * pow(z, 16) + 1313604007452534163360.0 / 518594883669909571959.0 * pow(z, 17) + 1601673308137590640.0 / 338287595348929923.0 * pow(z, 18) + 34117455247922362939232.0 / 11581952401961313773751.0 * pow(z, 19) + 440459860314878711272.0 / 81909140042159220465.0 * pow(z, 20) + 5958611072146769924704.0 / 1791023567313605222745.0 * pow(z, 21) + 63180886725054032940615536.0 / 10597486447794602102982165.0 * pow(z, 22) + 65365287366024052315837984.0 / 17864334297710900687884221.0 * pow(z, 23) + 27874331155216459994053508.0 / 4312080692550907062592743.0 * pow(z, 24) + 5942298851456732028712544.0 / 1513926635399228871854595.0 * pow(z, 25) + 1192092202359510055483312.0 / 173729286029419706606265.0 * pow(z, 26) + 350506533655692348512.0 / 85286836538743105845.0 * pow(z, 27) + 817996883441917606648.0 / 114672796059022908651.0 * pow(z, 28) + 725903303147284317856.0 / 172864961223303190653.0 * pow(z, 29) + 7367697217432916944.0 / 1014862786046789769.0 * pow(z, 30) + 30548930405872887520.0 / 7304153291125486929.0 * pow(z, 31) + 103272675658011713.0 / 14293842056997039.0 * pow(z, 32) + 405343166032938800.0 / 100056894398979273.0 * pow(z, 33) + 100287623571020248.0 / 14293842056997039.0 * pow(z, 34) + 24072679007055472.0 / 6332714835378435.0 * pow(z, 35) + 1553872158664004.0 / 234544993902905.0 * pow(z, 36) + 1386242693488.0 / 403691900005.0 * pow(z, 37) + 3418766377384.0 / 565168660007.0 * pow(z, 38) + 24691032496.0 / 8352246207.0 * pow(z, 39) + 17385279146.0 / 3284591205.0 * pow(z, 40) + 1118561392.0 / 469227315.0 * pow(z, 41) + 463118248.0 / 105954555.0 * pow(z, 42) + 36948304.0 / 21190911.0 * pow(z, 43) + 103316.0 / 31209.0 * pow(z, 44) + 392656.0 / 364105.0 * pow(z, 45) + 1112.0 / 515.0 * pow(z, 46) + 16.0 / 35.0 * pow(z, 47) + 1.0 * pow(z, 48));
                result += c_f[49]() *  2.614131094941481 * (1.0 / 107.0 * pow(z, 0) + 337.0 / 3745.0 * pow(z, 1) + 26616.0 / 385735.0 * pow(z, 2) + 2244040.0 / 7791847.0 * pow(z, 3) + 1514228.0 / 7791847.0 * pow(z, 4) + 450823644.0 / 755809159.0 * pow(z, 5) + 210313656.0 / 539863685.0 * pow(z, 6) + 16890737976.0 / 16735774235.0 * pow(z, 7) + 15275836998.0 / 23430083929.0 * pow(z, 8) + 3157987190194.0 / 2085277469681.0 * pow(z, 9) + 58899516395896.0 / 60473046620749.0 * pow(z, 10) + 633729114679272.0 / 302365233103745.0 * pow(z, 11) + 33741041258761876.0 / 25096314347610835.0 * pow(z, 12) + 1960825618787924.0 / 717037552788881.0 * pow(z, 13) + 99113491752813176.0 / 56645966670321599.0 * pow(z, 14) + 4055530066101328168.0 / 1189565300076753579.0 * pow(z, 15) + 2586244720853199647.0 / 1189565300076753579.0 * pow(z, 16) + 1332480036282877479.0 / 325236954702633001.0 * pow(z, 17) + 16035613563947708448656.0 / 6165516950297813799957.0 * pow(z, 18) + 29433531590728775913712.0 / 6165516950297813799957.0 * pow(z, 19) + 59256911059608357345544.0 / 19670935031902548790339.0 * pow(z, 20) + 1598016281395832130174056.0 / 295064025478538231855085.0 * pow(z, 21) + 98648695649719616751152.0 / 29090819413377008774445.0 * pow(z, 22) + 50384822011414656658124816.0 / 8399489258622388333474753.0 * pow(z, 23) + 5529661760072588058265577972.0 / 1486709598776162735025031281.0 * pow(z, 24) + 163750671537097634138905652.0 / 25198467775867165000424259.0 * pow(z, 25) + 38601663515107676120016.0 / 9696939804459002924815.0 * pow(z, 26) + 2033838903594695438403344.0 / 295064025478538231855085.0 * pow(z, 27) + 8465273008515479620792.0 / 2034924313645091254173.0 * pow(z, 28) + 14716765795364387956856.0 / 2055172316765937933319.0 * pow(z, 29) + 843979661260405707824.0 / 198887643557993993547.0 * pow(z, 30) + 7106560193508679888.0 / 975710864107899003.0 * pow(z, 31) + 1673452466434423301.0 / 396521766692251193.0 * pow(z, 32) + 8618001390465322357.0 / 1189565300076753579.0 * pow(z, 33) + 693794442269692232.0 / 169937900010964797.0 * pow(z, 34) + 5042123019740376.0 / 717037552788881.0 * pow(z, 35) + 96032194351860724.0 / 25096314347610835.0 * pow(z, 36) + 2006808863151028.0 / 302365233103745.0 * pow(z, 37) + 208825558130904.0 / 60473046620749.0 * pow(z, 38) + 12631948760776.0 / 2085277469681.0 * pow(z, 39) + 69589924102.0 / 23430083929.0 * pow(z, 40) + 88676374374.0 / 16735774235.0 * pow(z, 41) + 1291926744.0 / 539863685.0 * pow(z, 42) + 3306040056.0 / 755809159.0 * pow(z, 43) + 13628052.0 / 7791847.0 * pow(z, 44) + 25806460.0 / 7791847.0 * pow(z, 45) + 416984.0 / 385735.0 * pow(z, 46) + 8088.0 / 3745.0 * pow(z, 47) + 49.0 / 107.0 * pow(z, 48) + 1.0 * pow(z, 49));

                return result;
            }

            virtual double f_p(const double & q2) const
            {
                const auto z = _z(q2);

                const auto phi = this->phi_p(z);

                const double epsilon = Process_::epsilon_isospin;
                const auto series_total   = this->series_p(z, _a_fp) + (epsilon / 6) * this->series_p(z, _a_ft);

                double correction_to_asymptotics= pow(1.0 - z, 5.0 / 2.0 - 1.0 / 2.0);
                double K = (315 * M_PI) / 8192;

                return 1.0 / (phi) * correction_to_asymptotics * sqrt(K) * series_total;
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