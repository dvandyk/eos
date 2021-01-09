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
        static constexpr const auto asymptotic_case_switch = 1; // makeshift implementation of switching the different asymptotic behaviours
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
        static constexpr const auto asymptotic_case_switch = 1; // makeshift implementation of switching the different asymptotic behaviours
    };

    template <typename Process_, typename Transition_, bool has_scalar_form_factor_> class EGJvD2020UnitarityBoundsBase;

    template <typename Process_> class EGJvD2020UnitarityBoundsBase<Process_, VacuumToPP, false> :
        public virtual ParameterUser
    {
        private:
            // parameters for form factors f_+ and f_T
            std::array<UsedParameter, 20u> _a_fp, _a_ft;

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
                    UsedParameter(p[_par_name("+", "19")], *this)
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
                    UsedParameter(p[_par_name("T", "19")], *this)
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
            std::array<UsedParameter, 20u> _a_fp, _a_ft;

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
                    UsedParameter(p[_par_name("+", "19")], *this)
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
                    UsedParameter(p[_par_name("T", "19")], *this)
                }},
                _t_0(p[stringify(Process_::label) + "::t_0@EGJvD2020"], *this)
            {
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

            complex<double> blaschke_p(const double & q2) const
            {
                constexpr const complex<double> t_s_p = complex<double>{0.77 * 0.77, 0.77 * 0.1};
                constexpr const complex<double> t_s_m = complex<double>{0.77 * 0.77, - 0.77 * 0.1};
                return 1.0;

                return 1.0 / (_z(q2, t_s_p) * _z(q2, t_s_m));
            }

            complex<double> series_p(const complex<double> & z) const
            {
                complex<double> result = 0.0;

                //TODO: properly implement the options for the different polynomial choices

               switch (Process_::asymptotic_case_switch) //intermediate solution: DONT FORGET TO APPLY THIS OPTION BELOW TOO
                {
                    case 0:
                        result += _a_fp[0]() *  1.0 * (1.0 * pow(z, 0));
                        result += _a_fp[1]() *  1.1067971810589328 * (-3.0 / 7.0 * pow(z, 0) + 1.0 * pow(z, 1));
                        result += _a_fp[2]() *  1.3311179511974134 * (5.0 / 9.0 * pow(z, 0) + -2.0 / 3.0 * pow(z, 1) + 1.0 * pow(z, 2));
                        result += _a_fp[3]() *  1.3835670610779949 * (-3.0 / 11.0 * pow(z, 0) + 73.0 / 99.0 * pow(z, 1) + -9.0 / 11.0 * pow(z, 2) + 1.0 * pow(z, 3));
                        result += _a_fp[4]() *  1.4988643161678277 * (5.0 / 13.0 * pow(z, 0) + -84.0 / 143.0 * pow(z, 1) + 146.0 / 143.0 * pow(z, 2) + -12.0 / 13.0 * pow(z, 3) + 1.0 * pow(z, 4));
                        result += _a_fp[5]() *  1.5297719867820068 * (-1.0 / 5.0 * pow(z, 0) + 37.0 / 65.0 * pow(z, 1) + -566.0 / 715.0 * pow(z, 2) + 74.0 / 65.0 * pow(z, 3) + -1.0 / 1.0 * pow(z, 4) + 1.0 * pow(z, 5));
                        result += _a_fp[6]() *  1.600566275047013 * (5.0 / 17.0 * pow(z, 0) + -42.0 / 85.0 * pow(z, 1) + 999.0 / 1105.0 * pow(z, 2) + -1132.0 / 1105.0 * pow(z, 3) + 111.0 / 85.0 * pow(z, 4) + -18.0 / 17.0 * pow(z, 5) + 1.0 * pow(z, 6));
                        result += _a_fp[7]() *  1.6208989129696074 * (-3.0 / 19.0 * pow(z, 0) + 149.0 / 323.0 * pow(z, 1) + -1131.0 / 1615.0 * pow(z, 2) + 22377.0 / 20995.0 * pow(z, 3) + -377.0 / 323.0 * pow(z, 4) + 447.0 / 323.0 * pow(z, 5) + -21.0 / 19.0 * pow(z, 6) + 1.0 * pow(z, 7));
                        result += _a_fp[8]() *  1.6688932588717365 * (5.0 / 21.0 * pow(z, 0) + -8.0 / 19.0 * pow(z, 1) + 1788.0 / 2261.0 * pow(z, 2) + -33176.0 / 33915.0 * pow(z, 3) + 14918.0 / 11305.0 * pow(z, 4) + -3016.0 / 2261.0 * pow(z, 5) + 596.0 / 399.0 * pow(z, 6) + -8.0 / 7.0 * pow(z, 7) + 1.0 * pow(z, 8));
                        result += _a_fp[9]() *  1.683273648470033 * (-3.0 / 23.0 * pow(z, 0) + 187.0 / 483.0 * pow(z, 1) + -1884.0 / 3059.0 * pow(z, 2) + 50172.0 / 52003.0 * pow(z, 3) + -179462.0 / 156009.0 * pow(z, 4) + 75258.0 / 52003.0 * pow(z, 5) + -628.0 / 437.0 * pow(z, 6) + 748.0 / 483.0 * pow(z, 7) + -27.0 / 23.0 * pow(z, 8) + 1.0 * pow(z, 9));
                        result += _a_fp[10]() *  1.7179839734269011 * (1.0 / 5.0 * pow(z, 0) + -42.0 / 115.0 * pow(z, 1) + 561.0 / 805.0 * pow(z, 2) + -13816.0 / 15295.0 * pow(z, 3) + 326118.0 / 260015.0 * pow(z, 4) + -358924.0 / 260015.0 * pow(z, 5) + 25086.0 / 15295.0 * pow(z, 6) + -1256.0 / 805.0 * pow(z, 7) + 187.0 / 115.0 * pow(z, 8) + -6.0 / 5.0 * pow(z, 9) + 1.0 * pow(z, 10));
                        result += _a_fp[11]() *  1.7286880269769846 * (-1.0 / 9.0 * pow(z, 0) + 1.0 / 3.0 * pow(z, 1) + -113.0 / 207.0 * pow(z, 2) + 1261.0 / 1449.0 * pow(z, 3) + -9962.0 / 9177.0 * pow(z, 4) + 3293986.0 / 2340135.0 * pow(z, 5) + -9962.0 / 6555.0 * pow(z, 6) + 2522.0 / 1449.0 * pow(z, 7) + -113.0 / 69.0 * pow(z, 8) + 5.0 / 3.0 * pow(z, 9) + -11.0 / 9.0 * pow(z, 10) + 1.0 * pow(z, 11));
                        result += _a_fp[12]() *  1.7549693763633214 * (5.0 / 29.0 * pow(z, 0) + -28.0 / 87.0 * pow(z, 1) + 18.0 / 29.0 * pow(z, 2) + -4972.0 / 6003.0 * pow(z, 3) + 16393.0 / 14007.0 * pow(z, 4) + -119544.0 / 88711.0 * pow(z, 5) + 6587972.0 / 3991995.0 * pow(z, 6) + -39848.0 / 23345.0 * pow(z, 7) + 1261.0 / 667.0 * pow(z, 8) + -452.0 / 261.0 * pow(z, 9) + 50.0 / 29.0 * pow(z, 10) + -36.0 / 29.0 * pow(z, 11) + 1.0 * pow(z, 12));
                        result += _a_fp[13]() *  1.763245410942023 * (-3.0 / 31.0 * pow(z, 0) + 263.0 / 899.0 * pow(z, 1) + -1318.0 / 2697.0 * pow(z, 2) + 2126.0 / 2697.0 * pow(z, 3) + -188179.0 / 186093.0 * pow(z, 4) + 414221.0 / 310155.0 * pow(z, 5) + -8882276.0 / 5892945.0 * pow(z, 6) + 1656884.0 / 930465.0 * pow(z, 7) + -188179.0 / 103385.0 * pow(z, 8) + 5315.0 / 2697.0 * pow(z, 9) + -14498.0 / 8091.0 * pow(z, 10) + 1578.0 / 899.0 * pow(z, 11) + -39.0 / 31.0 * pow(z, 12) + 1.0 * pow(z, 13));
                        result += _a_fp[14]() *  1.7838399589145673 * (5.0 / 33.0 * pow(z, 0) + -98.0 / 341.0 * pow(z, 1) + 5523.0 / 9889.0 * pow(z, 2) + -18452.0 / 24273.0 * pow(z, 3) + 96733.0 / 89001.0 * pow(z, 4) + -2634506.0 / 2047023.0 * pow(z, 5) + 49292299.0 / 30705345.0 * pow(z, 6) + -17764552.0 / 10235115.0 * pow(z, 7) + 20296829.0 / 10235115.0 * pow(z, 8) + -2634506.0 / 1335015.0 * pow(z, 9) + 186025.0 / 89001.0 * pow(z, 10) + -18452.0 / 9889.0 * pow(z, 11) + 1841.0 / 1023.0 * pow(z, 12) + -14.0 / 11.0 * pow(z, 13) + 1.0 * pow(z, 14));
                        result += _a_fp[15]() *  1.7904291706088293 * (-3.0 / 35.0 * pow(z, 0) + 43.0 / 165.0 * pow(z, 1) + -753.0 / 1705.0 * pow(z, 2) + 35523.0 / 49445.0 * pow(z, 3) + -250807.0 / 267003.0 * pow(z, 4) + 931561.0 / 741675.0 * pow(z, 5) + -2571079.0 / 1764675.0 * pow(z, 6) + 1885111433.0 / 1074687075.0 * pow(z, 7) + -367297.0 / 196075.0 * pow(z, 8) + 931561.0 / 445005.0 * pow(z, 9) + -250807.0 / 121365.0 * pow(z, 10) + 106569.0 / 49445.0 * pow(z, 11) + -3263.0 / 1705.0 * pow(z, 12) + 301.0 / 165.0 * pow(z, 13) + -9.0 / 7.0 * pow(z, 14) + 1.0 * pow(z, 15));
                        result += _a_fp[16]() *  1.8070045025509684 * (5.0 / 37.0 * pow(z, 0) + -48.0 / 185.0 * pow(z, 1) + 1032.0 / 2035.0 * pow(z, 2) + -4016.0 / 5735.0 * pow(z, 3) + 1847196.0 / 1829465.0 * pow(z, 4) + -4012912.0 / 3293037.0 * pow(z, 5) + 126692296.0 / 82325925.0 * pow(z, 6) + -111658288.0 / 65292975.0 * pow(z, 7) + 3770222866.0 / 1893496275.0 * pow(z, 8) + -5876752.0 / 2838825.0 * pow(z, 9) + 7452488.0 / 3293037.0 * pow(z, 10) + -4012912.0 / 1829465.0 * pow(z, 11) + 142092.0 / 63085.0 * pow(z, 12) + -4016.0 / 2035.0 * pow(z, 13) + 344.0 / 185.0 * pow(z, 14) + -48.0 / 37.0 * pow(z, 15) + 1.0 * pow(z, 16));
                        result += _a_fp[17]() *  1.8123745129112863 * (-1.0 / 13.0 * pow(z, 0) + 113.0 / 481.0 * pow(z, 1) + -968.0 / 2405.0 * pow(z, 2) + 17432.0 / 26455.0 * pow(z, 3) + -143276.0 / 164021.0 * pow(z, 4) + 509572.0 / 432419.0 * pow(z, 5) + -19873448.0 / 14269827.0 * pow(z, 6) + 201936184.0 / 118915225.0 * pow(z, 7) + -15288465814.0 / 8205150525.0 * pow(z, 8) + 50484046.0 / 23783045.0 * pow(z, 9) + -2839064.0 / 1297257.0 * pow(z, 10) + 1019144.0 / 432419.0 * pow(z, 11) + -143276.0 / 63085.0 * pow(z, 12) + 61012.0 / 26455.0 * pow(z, 13) + -968.0 / 481.0 * pow(z, 14) + 904.0 / 481.0 * pow(z, 15) + -17.0 / 13.0 * pow(z, 16) + 1.0 * pow(z, 17));
                        result += _a_fp[18]() *  1.826003630370946 * (5.0 / 41.0 * pow(z, 0) + -126.0 / 533.0 * pow(z, 1) + 9153.0 / 19721.0 * pow(z, 2) + -63888.0 / 98605.0 * pow(z, 3) + 78444.0 / 83435.0 * pow(z, 4) + -7736904.0 / 6724861.0 * pow(z, 5) + 25988172.0 / 17729179.0 * pow(z, 6) + -323653296.0 / 195020969.0 * pow(z, 7) + 9541484694.0 / 4875524225.0 * pow(z, 8) + -30576931628.0 / 14626572675.0 * pow(z, 9) + 454356414.0 / 195020969.0 * pow(z, 10) + -459928368.0 / 195020969.0 * pow(z, 11) + 1528716.0 / 611351.0 * pow(z, 12) + -2578968.0 / 1084655.0 * pow(z, 13) + 235332.0 / 98605.0 * pow(z, 14) + -40656.0 / 19721.0 * pow(z, 15) + 1017.0 / 533.0 * pow(z, 16) + -54.0 / 41.0 * pow(z, 17) + 1.0 * pow(z, 18));
                        result += _a_fp[19]() *  1.8304639525740638 * (-3.0 / 43.0 * pow(z, 0) + 377.0 / 1763.0 * pow(z, 1) + -8469.0 / 22919.0 * pow(z, 2) + 515547.0 / 848003.0 * pow(z, 3) + -690636.0 / 848003.0 * pow(z, 4) + 239940.0 / 216931.0 * pow(z, 5) + -383134500.0 / 289169023.0 * pow(z, 6) + 333468060.0 / 204534187.0 * pow(z, 7) + -15280160970.0 / 8385901667.0 * pow(z, 8) + 6782473438.0 / 3225346795.0 * pow(z, 9) + -5093386990.0 / 2287064091.0 * pow(z, 10) + 500202090.0 / 204534187.0 * pow(z, 11) + -54733500.0 / 22243771.0 * pow(z, 12) + 559860.0 / 216931.0 * pow(z, 13) + -2071908.0 / 848003.0 * pow(z, 14) + 2062188.0 / 848003.0 * pow(z, 15) + -47991.0 / 22919.0 * pow(z, 16) + 3393.0 / 1763.0 * pow(z, 17) + -57.0 / 43.0 * pow(z, 18) + 1.0 * pow(z, 19));

                        break;

                    case 1:
                        result += _a_fp[0]() *  1.0 * (1.0 * pow(z, 0));
                        result += _a_fp[1]() *  1.0041580220928046 * (1.0 / 11.0 * pow(z, 0) + 1.0 * pow(z, 1));
                        result += _a_fp[2]() *  1.3915668626887223 * (9.0 / 13.0 * pow(z, 0) + 2.0 / 13.0 * pow(z, 1) + 1.0 * pow(z, 2));
                        result += _a_fp[3]() *  1.394669579723865 * (1.0 / 15.0 * pow(z, 0) + 137.0 / 195.0 * pow(z, 1) + 1.0 / 5.0 * pow(z, 2) + 1.0 * pow(z, 3));
                        result += _a_fp[4]() *  1.6439499152771448 * (9.0 / 17.0 * pow(z, 0) + 44.0 / 255.0 * pow(z, 1) + 274.0 / 255.0 * pow(z, 2) + 4.0 / 17.0 * pow(z, 3) + 1.0 * pow(z, 4));
                        result += _a_fp[5]() *  1.6462315956469282 * (1.0 / 19.0 * pow(z, 0) + 175.0 / 323.0 * pow(z, 1) + 74.0 / 323.0 * pow(z, 2) + 350.0 / 323.0 * pow(z, 3) + 5.0 / 19.0 * pow(z, 4) + 1.0 * pow(z, 5));
                        result += _a_fp[6]() *  1.822044489432169 * (3.0 / 7.0 * pow(z, 0) + 22.0 / 133.0 * pow(z, 1) + 325.0 / 323.0 * pow(z, 2) + 740.0 / 2261.0 * pow(z, 3) + 25.0 / 19.0 * pow(z, 4) + 2.0 / 7.0 * pow(z, 5) + 1.0 * pow(z, 6));
                        result += _a_fp[7]() *  1.8237690941622522 * (1.0 / 23.0 * pow(z, 0) + 71.0 / 161.0 * pow(z, 1) + 681.0 / 3059.0 * pow(z, 2) + 53065.0 / 52003.0 * pow(z, 3) + 1135.0 / 3059.0 * pow(z, 4) + 213.0 / 161.0 * pow(z, 5) + 7.0 / 23.0 * pow(z, 6) + 1.0 * pow(z, 7));
                        result += _a_fp[8]() *  1.9548363704716327 * (9.0 / 25.0 * pow(z, 0) + 88.0 / 575.0 * pow(z, 1) + 3692.0 / 4025.0 * pow(z, 2) + 5448.0 / 15295.0 * pow(z, 3) + 21226.0 / 15295.0 * pow(z, 4) + 1816.0 / 4025.0 * pow(z, 5) + 852.0 / 575.0 * pow(z, 6) + 8.0 / 25.0 * pow(z, 7) + 1.0 * pow(z, 8));
                        result += _a_fp[9]() *  1.9561785171250927 * (1.0 / 27.0 * pow(z, 0) + 251.0 / 675.0 * pow(z, 1) + 1076.0 / 5175.0 * pow(z, 2) + 580.0 / 621.0 * pow(z, 3) + 24046.0 / 58995.0 * pow(z, 4) + 290.0 / 207.0 * pow(z, 5) + 7532.0 / 15525.0 * pow(z, 6) + 1004.0 / 675.0 * pow(z, 7) + 1.0 / 3.0 * pow(z, 8) + 1.0 * pow(z, 9));
                        result += _a_fp[10]() *  2.05778352996703 * (9.0 / 29.0 * pow(z, 0) + 110.0 / 783.0 * pow(z, 1) + 3263.0 / 3915.0 * pow(z, 2) + 2152.0 / 6003.0 * pow(z, 3) + 850.0 / 621.0 * pow(z, 4) + 48092.0 / 90045.0 * pow(z, 5) + 350.0 / 207.0 * pow(z, 6) + 2152.0 / 3915.0 * pow(z, 7) + 1255.0 / 783.0 * pow(z, 8) + 10.0 / 29.0 * pow(z, 9) + 1.0 * pow(z, 10));
                        result += _a_fp[11]() *  2.0588550132627392 * (1.0 / 31.0 * pow(z, 0) + 289.0 / 899.0 * pow(z, 1) + 1555.0 / 8091.0 * pow(z, 2) + 6887.0 / 8091.0 * pow(z, 3) + 76862.0 / 186093.0 * pow(z, 4) + 1289614.0 / 930465.0 * pow(z, 5) + 538034.0 / 930465.0 * pow(z, 6) + 13774.0 / 8091.0 * pow(z, 7) + 1555.0 / 2697.0 * pow(z, 8) + 1445.0 / 899.0 * pow(z, 9) + 11.0 / 31.0 * pow(z, 10) + 1.0 * pow(z, 11));
                        result += _a_fp[12]() *  2.1399786377482064 * (3.0 / 11.0 * pow(z, 0) + 4.0 / 31.0 * pow(z, 1) + 7514.0 / 9889.0 * pow(z, 2) + 31100.0 / 89001.0 * pow(z, 3) + 117079.0 / 89001.0 * pow(z, 4) + 5841512.0 / 10235115.0 * pow(z, 5) + 18054596.0 / 10235115.0 * pow(z, 6) + 307448.0 / 445005.0 * pow(z, 7) + 172175.0 / 89001.0 * pow(z, 8) + 6220.0 / 9889.0 * pow(z, 9) + 578.0 / 341.0 * pow(z, 10) + 4.0 / 11.0 * pow(z, 11) + 1.0 * pow(z, 12));
                        result += _a_fp[13]() *  2.140852633552563 * (1.0 / 35.0 * pow(z, 0) + 109.0 / 385.0 * pow(z, 1) + 2118.0 / 11935.0 * pow(z, 2) + 53842.0 / 69223.0 * pow(z, 3) + 28015.0 / 69223.0 * pow(z, 4) + 2310697.0 / 1730575.0 * pow(z, 5) + 963236.0 / 1550775.0 * pow(z, 6) + 9242788.0 / 5191725.0 * pow(z, 7) + 50427.0 / 69223.0 * pow(z, 8) + 134605.0 / 69223.0 * pow(z, 9) + 706.0 / 1085.0 * pow(z, 10) + 654.0 / 385.0 * pow(z, 11) + 13.0 / 35.0 * pow(z, 12) + 1.0 * pow(z, 13));
                        result += _a_fp[14]() *  2.207143478674807 * (9.0 / 37.0 * pow(z, 0) + 22.0 / 185.0 * pow(z, 1) + 1417.0 / 2035.0 * pow(z, 2) + 4236.0 / 12617.0 * pow(z, 3) + 457657.0 / 365893.0 * pow(z, 4) + 212914.0 / 365893.0 * pow(z, 5) + 16174879.0 / 9147325.0 * pow(z, 6) + 1926472.0 / 2494725.0 * pow(z, 7) + 2310697.0 / 1097679.0 * pow(z, 8) + 302562.0 / 365893.0 * pow(z, 9) + 26921.0 / 12617.0 * pow(z, 10) + 1412.0 / 2035.0 * pow(z, 11) + 327.0 / 185.0 * pow(z, 12) + 14.0 / 37.0 * pow(z, 13) + 1.0 * pow(z, 14));
                        result += _a_fp[15]() *  2.207869393339616 * (1.0 / 39.0 * pow(z, 0) + 365.0 / 1443.0 * pow(z, 1) + 79.0 / 481.0 * pow(z, 2) + 11335.0 / 15873.0 * pow(z, 3) + 192125.0 / 492063.0 * pow(z, 4) + 465415.0 / 365893.0 * pow(z, 5) + 27221635.0 / 42809481.0 * pow(z, 6) + 382730407.0 / 214047405.0 * pow(z, 7) + 3888805.0 / 4756609.0 * pow(z, 8) + 2327075.0 / 1097679.0 * pow(z, 9) + 38425.0 / 44733.0 * pow(z, 10) + 11335.0 / 5291.0 * pow(z, 11) + 79.0 / 111.0 * pow(z, 12) + 2555.0 / 1443.0 * pow(z, 13) + 5.0 / 13.0 * pow(z, 14) + 1.0 * pow(z, 15));
                        result += _a_fp[16]() *  2.2630661281731066 * (9.0 / 41.0 * pow(z, 0) + 176.0 / 1599.0 * pow(z, 1) + 2920.0 / 4551.0 * pow(z, 2) + 6320.0 / 19721.0 * pow(z, 3) + 770780.0 / 650793.0 * pow(z, 4) + 11681200.0 / 20174583.0 * pow(z, 5) + 26063240.0 / 15001613.0 * pow(z, 6) + 1431080240.0 / 1755188721.0 * pow(z, 7) + 3827304070.0 / 1755188721.0 * pow(z, 8) + 186662640.0 / 195020969.0 * pow(z, 9) + 3723320.0 / 1551891.0 * pow(z, 10) + 614800.0 / 650793.0 * pow(z, 11) + 45340.0 / 19721.0 * pow(z, 12) + 44240.0 / 59163.0 * pow(z, 13) + 2920.0 / 1599.0 * pow(z, 14) + 16.0 / 41.0 * pow(z, 15) + 1.0 * pow(z, 16));
                        result += _a_fp[17]() *  2.2636783468041286 * (1.0 / 43.0 * pow(z, 0) + 403.0 / 1763.0 * pow(z, 1) + 3496.0 / 22919.0 * pow(z, 2) + 558840.0 / 848003.0 * pow(z, 3) + 317100.0 / 848003.0 * pow(z, 4) + 1022980.0 / 848003.0 * pow(z, 5) + 16687720.0 / 26288093.0 * pow(z, 6) + 1341455800.0 / 762354697.0 * pow(z, 7) + 1980719830.0 / 2287064091.0 * pow(z, 8) + 1676819750.0 / 762354697.0 * pow(z, 9) + 26223560.0 / 26288093.0 * pow(z, 10) + 2045960.0 / 848003.0 * pow(z, 11) + 63420.0 / 65231.0 * pow(z, 12) + 1955940.0 / 848003.0 * pow(z, 13) + 17480.0 / 22919.0 * pow(z, 14) + 3224.0 / 1763.0 * pow(z, 15) + 17.0 / 43.0 * pow(z, 16) + 1.0 * pow(z, 17));
                        result += _a_fp[18]() *  2.310357038107123 * (1.0 / 5.0 * pow(z, 0) + 22.0 / 215.0 * pow(z, 1) + 5239.0 / 8815.0 * pow(z, 2) + 6992.0 / 22919.0 * pow(z, 3) + 950028.0 / 848003.0 * pow(z, 4) + 481992.0 / 848003.0 * pow(z, 5) + 1432172.0 / 848003.0 * pow(z, 6) + 21932432.0 / 26288093.0 * pow(z, 7) + 1676819750.0 / 762354697.0 * pow(z, 8) + 792287932.0 / 762354697.0 * pow(z, 9) + 67072790.0 / 26288093.0 * pow(z, 10) + 953584.0 / 848003.0 * pow(z, 11) + 2250556.0 / 848003.0 * pow(z, 12) + 887880.0 / 848003.0 * pow(z, 13) + 55884.0 / 22919.0 * pow(z, 14) + 6992.0 / 8815.0 * pow(z, 15) + 403.0 / 215.0 * pow(z, 16) + 2.0 / 5.0 * pow(z, 17) + 1.0 * pow(z, 18));
                        result += _a_fp[19]() *  2.310880157560925 * (1.0 / 47.0 * pow(z, 0) + 49.0 / 235.0 * pow(z, 1) + 1437.0 / 10105.0 * pow(z, 2) + 50645.0 / 82861.0 * pow(z, 3) + 384508.0 / 1077193.0 * pow(z, 4) + 45539196.0 / 39856141.0 * pow(z, 5) + 24904180.0 / 39856141.0 * pow(z, 6) + 68265668.0 / 39856141.0 * pow(z, 7) + 1097897094.0 / 1235540371.0 * pow(z, 8) + 1941532102.0 / 873918799.0 * pow(z, 9) + 1341874226.0 / 1235540371.0 * pow(z, 10) + 102398502.0 / 39856141.0 * pow(z, 11) + 3557740.0 / 3065857.0 * pow(z, 12) + 106258124.0 / 39856141.0 * pow(z, 13) + 1153524.0 / 1077193.0 * pow(z, 14) + 202580.0 / 82861.0 * pow(z, 15) + 8143.0 / 10105.0 * pow(z, 16) + 441.0 / 235.0 * pow(z, 17) + 19.0 / 47.0 * pow(z, 18) + 1.0 * pow(z, 19));

                        break;
                }

                return result;
            }

            virtual complex<double> f_p(const double & q2) const
            {
                const auto z = _z(q2);
                //return this->phi_p(z);
                const auto phi      = this->phi_p(z);
                const auto blaschke = this->blaschke_p(q2);
                const auto series   = this->series_p(z);
                complex<double> correction_to_asymptotics;
                complex<double> K;
                switch (Process_::asymptotic_case_switch)
                {
                    case 0: // The asymptotics used by Buck and Lebed
                        correction_to_asymptotics = 1.0;
                        K = (5 * M_PI) / 64;
                        break;
                    case 1: // correct asyptotic behaviour
                        correction_to_asymptotics = pow(1.0 - z, 5.0 / 2.0 - 1.0 / 2.0);
                        K = (315 * M_PI) / 8192;
                        break;
                    default:
                        correction_to_asymptotics = 1.0;
                        K = (5 * M_PI) / 64;
                        break;
                }

                return 1.0 / (phi * blaschke) * correction_to_asymptotics * sqrt(K) * series; //todo combine asymptotics and phi into one
                // return 1.0 / (phi * blaschke) /* asymptotics*/ * sqrt(K) * series; //todo combine asymptotics and phi into one
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
            std::array<UsedParameter, 20u> _a_fp, _a_ft;

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
                    UsedParameter(p[_par_name("+", "19")], *this)
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
                    UsedParameter(p[_par_name("T", "19")], *this)
                }},
                _t_0(p[stringify(Process_::label) + "::t_0@EGJvD2020"], *this)
            {
            }

            /* f_+ */

            double phi_p(const double & z) const //Divided by asymptotic behaviour to improve numeric behaviour of Formfactor
            {
                const double t_p = Process_::t_p;
                const double t_0 = this->_t_0;
                const double tfactor = 1.0 - t_0 / t_p;
                const double chi = 0.00405; //GeV^-2 as a makeshift value
                const double Q2 = Process_::Q2;

                const double part0 = 1.0 / sqrt(12.0 * M_PI * t_p * chi); // = 
                const double part1 = /*[Asymptotics:] (1.0 + z) * (1.0 + z) * sqrt(1.0 - z) */ pow(tfactor, 1.25);
                const double part2 = pow(sqrt(tfactor) * (1.0 + z) + (1.0 - z), -0.5);
                const double part3 = pow(sqrt(1.0 + Q2 / t_p) * (1.0 - z) + sqrt(tfactor) * (1.0 + z), -3.0);

                return part0 * part1 * part2 * part3;
            }

            double blaschke_p(const double & q2) const
            {
                constexpr const complex<double> t_s_p = complex<double>{0.77 * 0.77, 0.77 * 0.1};
                constexpr const complex<double> t_s_m = complex<double>{0.77 * 0.77, - 0.77 * 0.1};
                return 1.0;

                return 1.0 / (_z(q2, t_s_p) * _z(q2, t_s_m));
            }

            double series_p(const double & z) const
            {
                double result = 0.0;

                switch (Process_::asymptotic_case_switch) //intermediate solution: DONT FORGET TO APPLY THIS OPTION ABOVE TOO
                {
                    case 0:
                        result += _a_fp[0]() *  1.0 * (1.0 * pow(z, 0));
                        result += _a_fp[1]() *  1.1067971810589328 * (-3.0 / 7.0 * pow(z, 0) + 1.0 * pow(z, 1));
                        result += _a_fp[2]() *  1.3311179511974134 * (5.0 / 9.0 * pow(z, 0) + -2.0 / 3.0 * pow(z, 1) + 1.0 * pow(z, 2));
                        result += _a_fp[3]() *  1.3835670610779949 * (-3.0 / 11.0 * pow(z, 0) + 73.0 / 99.0 * pow(z, 1) + -9.0 / 11.0 * pow(z, 2) + 1.0 * pow(z, 3));
                        result += _a_fp[4]() *  1.4988643161678277 * (5.0 / 13.0 * pow(z, 0) + -84.0 / 143.0 * pow(z, 1) + 146.0 / 143.0 * pow(z, 2) + -12.0 / 13.0 * pow(z, 3) + 1.0 * pow(z, 4));
                        result += _a_fp[5]() *  1.5297719867820068 * (-1.0 / 5.0 * pow(z, 0) + 37.0 / 65.0 * pow(z, 1) + -566.0 / 715.0 * pow(z, 2) + 74.0 / 65.0 * pow(z, 3) + -1.0 / 1.0 * pow(z, 4) + 1.0 * pow(z, 5));
                        result += _a_fp[6]() *  1.600566275047013 * (5.0 / 17.0 * pow(z, 0) + -42.0 / 85.0 * pow(z, 1) + 999.0 / 1105.0 * pow(z, 2) + -1132.0 / 1105.0 * pow(z, 3) + 111.0 / 85.0 * pow(z, 4) + -18.0 / 17.0 * pow(z, 5) + 1.0 * pow(z, 6));
                        result += _a_fp[7]() *  1.6208989129696074 * (-3.0 / 19.0 * pow(z, 0) + 149.0 / 323.0 * pow(z, 1) + -1131.0 / 1615.0 * pow(z, 2) + 22377.0 / 20995.0 * pow(z, 3) + -377.0 / 323.0 * pow(z, 4) + 447.0 / 323.0 * pow(z, 5) + -21.0 / 19.0 * pow(z, 6) + 1.0 * pow(z, 7));
                        result += _a_fp[8]() *  1.6688932588717365 * (5.0 / 21.0 * pow(z, 0) + -8.0 / 19.0 * pow(z, 1) + 1788.0 / 2261.0 * pow(z, 2) + -33176.0 / 33915.0 * pow(z, 3) + 14918.0 / 11305.0 * pow(z, 4) + -3016.0 / 2261.0 * pow(z, 5) + 596.0 / 399.0 * pow(z, 6) + -8.0 / 7.0 * pow(z, 7) + 1.0 * pow(z, 8));
                        result += _a_fp[9]() *  1.683273648470033 * (-3.0 / 23.0 * pow(z, 0) + 187.0 / 483.0 * pow(z, 1) + -1884.0 / 3059.0 * pow(z, 2) + 50172.0 / 52003.0 * pow(z, 3) + -179462.0 / 156009.0 * pow(z, 4) + 75258.0 / 52003.0 * pow(z, 5) + -628.0 / 437.0 * pow(z, 6) + 748.0 / 483.0 * pow(z, 7) + -27.0 / 23.0 * pow(z, 8) + 1.0 * pow(z, 9));
                        result += _a_fp[10]() *  1.7179839734269011 * (1.0 / 5.0 * pow(z, 0) + -42.0 / 115.0 * pow(z, 1) + 561.0 / 805.0 * pow(z, 2) + -13816.0 / 15295.0 * pow(z, 3) + 326118.0 / 260015.0 * pow(z, 4) + -358924.0 / 260015.0 * pow(z, 5) + 25086.0 / 15295.0 * pow(z, 6) + -1256.0 / 805.0 * pow(z, 7) + 187.0 / 115.0 * pow(z, 8) + -6.0 / 5.0 * pow(z, 9) + 1.0 * pow(z, 10));
                        result += _a_fp[11]() *  1.7286880269769846 * (-1.0 / 9.0 * pow(z, 0) + 1.0 / 3.0 * pow(z, 1) + -113.0 / 207.0 * pow(z, 2) + 1261.0 / 1449.0 * pow(z, 3) + -9962.0 / 9177.0 * pow(z, 4) + 3293986.0 / 2340135.0 * pow(z, 5) + -9962.0 / 6555.0 * pow(z, 6) + 2522.0 / 1449.0 * pow(z, 7) + -113.0 / 69.0 * pow(z, 8) + 5.0 / 3.0 * pow(z, 9) + -11.0 / 9.0 * pow(z, 10) + 1.0 * pow(z, 11));
                        result += _a_fp[12]() *  1.7549693763633214 * (5.0 / 29.0 * pow(z, 0) + -28.0 / 87.0 * pow(z, 1) + 18.0 / 29.0 * pow(z, 2) + -4972.0 / 6003.0 * pow(z, 3) + 16393.0 / 14007.0 * pow(z, 4) + -119544.0 / 88711.0 * pow(z, 5) + 6587972.0 / 3991995.0 * pow(z, 6) + -39848.0 / 23345.0 * pow(z, 7) + 1261.0 / 667.0 * pow(z, 8) + -452.0 / 261.0 * pow(z, 9) + 50.0 / 29.0 * pow(z, 10) + -36.0 / 29.0 * pow(z, 11) + 1.0 * pow(z, 12));
                        result += _a_fp[13]() *  1.763245410942023 * (-3.0 / 31.0 * pow(z, 0) + 263.0 / 899.0 * pow(z, 1) + -1318.0 / 2697.0 * pow(z, 2) + 2126.0 / 2697.0 * pow(z, 3) + -188179.0 / 186093.0 * pow(z, 4) + 414221.0 / 310155.0 * pow(z, 5) + -8882276.0 / 5892945.0 * pow(z, 6) + 1656884.0 / 930465.0 * pow(z, 7) + -188179.0 / 103385.0 * pow(z, 8) + 5315.0 / 2697.0 * pow(z, 9) + -14498.0 / 8091.0 * pow(z, 10) + 1578.0 / 899.0 * pow(z, 11) + -39.0 / 31.0 * pow(z, 12) + 1.0 * pow(z, 13));
                        result += _a_fp[14]() *  1.7838399589145673 * (5.0 / 33.0 * pow(z, 0) + -98.0 / 341.0 * pow(z, 1) + 5523.0 / 9889.0 * pow(z, 2) + -18452.0 / 24273.0 * pow(z, 3) + 96733.0 / 89001.0 * pow(z, 4) + -2634506.0 / 2047023.0 * pow(z, 5) + 49292299.0 / 30705345.0 * pow(z, 6) + -17764552.0 / 10235115.0 * pow(z, 7) + 20296829.0 / 10235115.0 * pow(z, 8) + -2634506.0 / 1335015.0 * pow(z, 9) + 186025.0 / 89001.0 * pow(z, 10) + -18452.0 / 9889.0 * pow(z, 11) + 1841.0 / 1023.0 * pow(z, 12) + -14.0 / 11.0 * pow(z, 13) + 1.0 * pow(z, 14));
                        result += _a_fp[15]() *  1.7904291706088293 * (-3.0 / 35.0 * pow(z, 0) + 43.0 / 165.0 * pow(z, 1) + -753.0 / 1705.0 * pow(z, 2) + 35523.0 / 49445.0 * pow(z, 3) + -250807.0 / 267003.0 * pow(z, 4) + 931561.0 / 741675.0 * pow(z, 5) + -2571079.0 / 1764675.0 * pow(z, 6) + 1885111433.0 / 1074687075.0 * pow(z, 7) + -367297.0 / 196075.0 * pow(z, 8) + 931561.0 / 445005.0 * pow(z, 9) + -250807.0 / 121365.0 * pow(z, 10) + 106569.0 / 49445.0 * pow(z, 11) + -3263.0 / 1705.0 * pow(z, 12) + 301.0 / 165.0 * pow(z, 13) + -9.0 / 7.0 * pow(z, 14) + 1.0 * pow(z, 15));
                        result += _a_fp[16]() *  1.8070045025509684 * (5.0 / 37.0 * pow(z, 0) + -48.0 / 185.0 * pow(z, 1) + 1032.0 / 2035.0 * pow(z, 2) + -4016.0 / 5735.0 * pow(z, 3) + 1847196.0 / 1829465.0 * pow(z, 4) + -4012912.0 / 3293037.0 * pow(z, 5) + 126692296.0 / 82325925.0 * pow(z, 6) + -111658288.0 / 65292975.0 * pow(z, 7) + 3770222866.0 / 1893496275.0 * pow(z, 8) + -5876752.0 / 2838825.0 * pow(z, 9) + 7452488.0 / 3293037.0 * pow(z, 10) + -4012912.0 / 1829465.0 * pow(z, 11) + 142092.0 / 63085.0 * pow(z, 12) + -4016.0 / 2035.0 * pow(z, 13) + 344.0 / 185.0 * pow(z, 14) + -48.0 / 37.0 * pow(z, 15) + 1.0 * pow(z, 16));
                        result += _a_fp[17]() *  1.8123745129112863 * (-1.0 / 13.0 * pow(z, 0) + 113.0 / 481.0 * pow(z, 1) + -968.0 / 2405.0 * pow(z, 2) + 17432.0 / 26455.0 * pow(z, 3) + -143276.0 / 164021.0 * pow(z, 4) + 509572.0 / 432419.0 * pow(z, 5) + -19873448.0 / 14269827.0 * pow(z, 6) + 201936184.0 / 118915225.0 * pow(z, 7) + -15288465814.0 / 8205150525.0 * pow(z, 8) + 50484046.0 / 23783045.0 * pow(z, 9) + -2839064.0 / 1297257.0 * pow(z, 10) + 1019144.0 / 432419.0 * pow(z, 11) + -143276.0 / 63085.0 * pow(z, 12) + 61012.0 / 26455.0 * pow(z, 13) + -968.0 / 481.0 * pow(z, 14) + 904.0 / 481.0 * pow(z, 15) + -17.0 / 13.0 * pow(z, 16) + 1.0 * pow(z, 17));
                        result += _a_fp[18]() *  1.826003630370946 * (5.0 / 41.0 * pow(z, 0) + -126.0 / 533.0 * pow(z, 1) + 9153.0 / 19721.0 * pow(z, 2) + -63888.0 / 98605.0 * pow(z, 3) + 78444.0 / 83435.0 * pow(z, 4) + -7736904.0 / 6724861.0 * pow(z, 5) + 25988172.0 / 17729179.0 * pow(z, 6) + -323653296.0 / 195020969.0 * pow(z, 7) + 9541484694.0 / 4875524225.0 * pow(z, 8) + -30576931628.0 / 14626572675.0 * pow(z, 9) + 454356414.0 / 195020969.0 * pow(z, 10) + -459928368.0 / 195020969.0 * pow(z, 11) + 1528716.0 / 611351.0 * pow(z, 12) + -2578968.0 / 1084655.0 * pow(z, 13) + 235332.0 / 98605.0 * pow(z, 14) + -40656.0 / 19721.0 * pow(z, 15) + 1017.0 / 533.0 * pow(z, 16) + -54.0 / 41.0 * pow(z, 17) + 1.0 * pow(z, 18));
                        result += _a_fp[19]() *  1.8304639525740638 * (-3.0 / 43.0 * pow(z, 0) + 377.0 / 1763.0 * pow(z, 1) + -8469.0 / 22919.0 * pow(z, 2) + 515547.0 / 848003.0 * pow(z, 3) + -690636.0 / 848003.0 * pow(z, 4) + 239940.0 / 216931.0 * pow(z, 5) + -383134500.0 / 289169023.0 * pow(z, 6) + 333468060.0 / 204534187.0 * pow(z, 7) + -15280160970.0 / 8385901667.0 * pow(z, 8) + 6782473438.0 / 3225346795.0 * pow(z, 9) + -5093386990.0 / 2287064091.0 * pow(z, 10) + 500202090.0 / 204534187.0 * pow(z, 11) + -54733500.0 / 22243771.0 * pow(z, 12) + 559860.0 / 216931.0 * pow(z, 13) + -2071908.0 / 848003.0 * pow(z, 14) + 2062188.0 / 848003.0 * pow(z, 15) + -47991.0 / 22919.0 * pow(z, 16) + 3393.0 / 1763.0 * pow(z, 17) + -57.0 / 43.0 * pow(z, 18) + 1.0 * pow(z, 19));

                        break;

                    case 1:
                        result += _a_fp[0]() *  1.0 * (1.0 * pow(z, 0));
                        result += _a_fp[1]() *  1.0041580220928046 * (1.0 / 11.0 * pow(z, 0) + 1.0 * pow(z, 1));
                        result += _a_fp[2]() *  1.3915668626887223 * (9.0 / 13.0 * pow(z, 0) + 2.0 / 13.0 * pow(z, 1) + 1.0 * pow(z, 2));
                        result += _a_fp[3]() *  1.394669579723865 * (1.0 / 15.0 * pow(z, 0) + 137.0 / 195.0 * pow(z, 1) + 1.0 / 5.0 * pow(z, 2) + 1.0 * pow(z, 3));
                        result += _a_fp[4]() *  1.6439499152771448 * (9.0 / 17.0 * pow(z, 0) + 44.0 / 255.0 * pow(z, 1) + 274.0 / 255.0 * pow(z, 2) + 4.0 / 17.0 * pow(z, 3) + 1.0 * pow(z, 4));
                        result += _a_fp[5]() *  1.6462315956469282 * (1.0 / 19.0 * pow(z, 0) + 175.0 / 323.0 * pow(z, 1) + 74.0 / 323.0 * pow(z, 2) + 350.0 / 323.0 * pow(z, 3) + 5.0 / 19.0 * pow(z, 4) + 1.0 * pow(z, 5));
                        result += _a_fp[6]() *  1.822044489432169 * (3.0 / 7.0 * pow(z, 0) + 22.0 / 133.0 * pow(z, 1) + 325.0 / 323.0 * pow(z, 2) + 740.0 / 2261.0 * pow(z, 3) + 25.0 / 19.0 * pow(z, 4) + 2.0 / 7.0 * pow(z, 5) + 1.0 * pow(z, 6));
                        result += _a_fp[7]() *  1.8237690941622522 * (1.0 / 23.0 * pow(z, 0) + 71.0 / 161.0 * pow(z, 1) + 681.0 / 3059.0 * pow(z, 2) + 53065.0 / 52003.0 * pow(z, 3) + 1135.0 / 3059.0 * pow(z, 4) + 213.0 / 161.0 * pow(z, 5) + 7.0 / 23.0 * pow(z, 6) + 1.0 * pow(z, 7));
                        result += _a_fp[8]() *  1.9548363704716327 * (9.0 / 25.0 * pow(z, 0) + 88.0 / 575.0 * pow(z, 1) + 3692.0 / 4025.0 * pow(z, 2) + 5448.0 / 15295.0 * pow(z, 3) + 21226.0 / 15295.0 * pow(z, 4) + 1816.0 / 4025.0 * pow(z, 5) + 852.0 / 575.0 * pow(z, 6) + 8.0 / 25.0 * pow(z, 7) + 1.0 * pow(z, 8));
                        result += _a_fp[9]() *  1.9561785171250927 * (1.0 / 27.0 * pow(z, 0) + 251.0 / 675.0 * pow(z, 1) + 1076.0 / 5175.0 * pow(z, 2) + 580.0 / 621.0 * pow(z, 3) + 24046.0 / 58995.0 * pow(z, 4) + 290.0 / 207.0 * pow(z, 5) + 7532.0 / 15525.0 * pow(z, 6) + 1004.0 / 675.0 * pow(z, 7) + 1.0 / 3.0 * pow(z, 8) + 1.0 * pow(z, 9));
                        result += _a_fp[10]() *  2.05778352996703 * (9.0 / 29.0 * pow(z, 0) + 110.0 / 783.0 * pow(z, 1) + 3263.0 / 3915.0 * pow(z, 2) + 2152.0 / 6003.0 * pow(z, 3) + 850.0 / 621.0 * pow(z, 4) + 48092.0 / 90045.0 * pow(z, 5) + 350.0 / 207.0 * pow(z, 6) + 2152.0 / 3915.0 * pow(z, 7) + 1255.0 / 783.0 * pow(z, 8) + 10.0 / 29.0 * pow(z, 9) + 1.0 * pow(z, 10));
                        result += _a_fp[11]() *  2.0588550132627392 * (1.0 / 31.0 * pow(z, 0) + 289.0 / 899.0 * pow(z, 1) + 1555.0 / 8091.0 * pow(z, 2) + 6887.0 / 8091.0 * pow(z, 3) + 76862.0 / 186093.0 * pow(z, 4) + 1289614.0 / 930465.0 * pow(z, 5) + 538034.0 / 930465.0 * pow(z, 6) + 13774.0 / 8091.0 * pow(z, 7) + 1555.0 / 2697.0 * pow(z, 8) + 1445.0 / 899.0 * pow(z, 9) + 11.0 / 31.0 * pow(z, 10) + 1.0 * pow(z, 11));
                        result += _a_fp[12]() *  2.1399786377482064 * (3.0 / 11.0 * pow(z, 0) + 4.0 / 31.0 * pow(z, 1) + 7514.0 / 9889.0 * pow(z, 2) + 31100.0 / 89001.0 * pow(z, 3) + 117079.0 / 89001.0 * pow(z, 4) + 5841512.0 / 10235115.0 * pow(z, 5) + 18054596.0 / 10235115.0 * pow(z, 6) + 307448.0 / 445005.0 * pow(z, 7) + 172175.0 / 89001.0 * pow(z, 8) + 6220.0 / 9889.0 * pow(z, 9) + 578.0 / 341.0 * pow(z, 10) + 4.0 / 11.0 * pow(z, 11) + 1.0 * pow(z, 12));
                        result += _a_fp[13]() *  2.140852633552563 * (1.0 / 35.0 * pow(z, 0) + 109.0 / 385.0 * pow(z, 1) + 2118.0 / 11935.0 * pow(z, 2) + 53842.0 / 69223.0 * pow(z, 3) + 28015.0 / 69223.0 * pow(z, 4) + 2310697.0 / 1730575.0 * pow(z, 5) + 963236.0 / 1550775.0 * pow(z, 6) + 9242788.0 / 5191725.0 * pow(z, 7) + 50427.0 / 69223.0 * pow(z, 8) + 134605.0 / 69223.0 * pow(z, 9) + 706.0 / 1085.0 * pow(z, 10) + 654.0 / 385.0 * pow(z, 11) + 13.0 / 35.0 * pow(z, 12) + 1.0 * pow(z, 13));
                        result += _a_fp[14]() *  2.207143478674807 * (9.0 / 37.0 * pow(z, 0) + 22.0 / 185.0 * pow(z, 1) + 1417.0 / 2035.0 * pow(z, 2) + 4236.0 / 12617.0 * pow(z, 3) + 457657.0 / 365893.0 * pow(z, 4) + 212914.0 / 365893.0 * pow(z, 5) + 16174879.0 / 9147325.0 * pow(z, 6) + 1926472.0 / 2494725.0 * pow(z, 7) + 2310697.0 / 1097679.0 * pow(z, 8) + 302562.0 / 365893.0 * pow(z, 9) + 26921.0 / 12617.0 * pow(z, 10) + 1412.0 / 2035.0 * pow(z, 11) + 327.0 / 185.0 * pow(z, 12) + 14.0 / 37.0 * pow(z, 13) + 1.0 * pow(z, 14));
                        result += _a_fp[15]() *  2.207869393339616 * (1.0 / 39.0 * pow(z, 0) + 365.0 / 1443.0 * pow(z, 1) + 79.0 / 481.0 * pow(z, 2) + 11335.0 / 15873.0 * pow(z, 3) + 192125.0 / 492063.0 * pow(z, 4) + 465415.0 / 365893.0 * pow(z, 5) + 27221635.0 / 42809481.0 * pow(z, 6) + 382730407.0 / 214047405.0 * pow(z, 7) + 3888805.0 / 4756609.0 * pow(z, 8) + 2327075.0 / 1097679.0 * pow(z, 9) + 38425.0 / 44733.0 * pow(z, 10) + 11335.0 / 5291.0 * pow(z, 11) + 79.0 / 111.0 * pow(z, 12) + 2555.0 / 1443.0 * pow(z, 13) + 5.0 / 13.0 * pow(z, 14) + 1.0 * pow(z, 15));
                        result += _a_fp[16]() *  2.2630661281731066 * (9.0 / 41.0 * pow(z, 0) + 176.0 / 1599.0 * pow(z, 1) + 2920.0 / 4551.0 * pow(z, 2) + 6320.0 / 19721.0 * pow(z, 3) + 770780.0 / 650793.0 * pow(z, 4) + 11681200.0 / 20174583.0 * pow(z, 5) + 26063240.0 / 15001613.0 * pow(z, 6) + 1431080240.0 / 1755188721.0 * pow(z, 7) + 3827304070.0 / 1755188721.0 * pow(z, 8) + 186662640.0 / 195020969.0 * pow(z, 9) + 3723320.0 / 1551891.0 * pow(z, 10) + 614800.0 / 650793.0 * pow(z, 11) + 45340.0 / 19721.0 * pow(z, 12) + 44240.0 / 59163.0 * pow(z, 13) + 2920.0 / 1599.0 * pow(z, 14) + 16.0 / 41.0 * pow(z, 15) + 1.0 * pow(z, 16));
                        result += _a_fp[17]() *  2.2636783468041286 * (1.0 / 43.0 * pow(z, 0) + 403.0 / 1763.0 * pow(z, 1) + 3496.0 / 22919.0 * pow(z, 2) + 558840.0 / 848003.0 * pow(z, 3) + 317100.0 / 848003.0 * pow(z, 4) + 1022980.0 / 848003.0 * pow(z, 5) + 16687720.0 / 26288093.0 * pow(z, 6) + 1341455800.0 / 762354697.0 * pow(z, 7) + 1980719830.0 / 2287064091.0 * pow(z, 8) + 1676819750.0 / 762354697.0 * pow(z, 9) + 26223560.0 / 26288093.0 * pow(z, 10) + 2045960.0 / 848003.0 * pow(z, 11) + 63420.0 / 65231.0 * pow(z, 12) + 1955940.0 / 848003.0 * pow(z, 13) + 17480.0 / 22919.0 * pow(z, 14) + 3224.0 / 1763.0 * pow(z, 15) + 17.0 / 43.0 * pow(z, 16) + 1.0 * pow(z, 17));
                        result += _a_fp[18]() *  2.310357038107123 * (1.0 / 5.0 * pow(z, 0) + 22.0 / 215.0 * pow(z, 1) + 5239.0 / 8815.0 * pow(z, 2) + 6992.0 / 22919.0 * pow(z, 3) + 950028.0 / 848003.0 * pow(z, 4) + 481992.0 / 848003.0 * pow(z, 5) + 1432172.0 / 848003.0 * pow(z, 6) + 21932432.0 / 26288093.0 * pow(z, 7) + 1676819750.0 / 762354697.0 * pow(z, 8) + 792287932.0 / 762354697.0 * pow(z, 9) + 67072790.0 / 26288093.0 * pow(z, 10) + 953584.0 / 848003.0 * pow(z, 11) + 2250556.0 / 848003.0 * pow(z, 12) + 887880.0 / 848003.0 * pow(z, 13) + 55884.0 / 22919.0 * pow(z, 14) + 6992.0 / 8815.0 * pow(z, 15) + 403.0 / 215.0 * pow(z, 16) + 2.0 / 5.0 * pow(z, 17) + 1.0 * pow(z, 18));
                        result += _a_fp[19]() *  2.310880157560925 * (1.0 / 47.0 * pow(z, 0) + 49.0 / 235.0 * pow(z, 1) + 1437.0 / 10105.0 * pow(z, 2) + 50645.0 / 82861.0 * pow(z, 3) + 384508.0 / 1077193.0 * pow(z, 4) + 45539196.0 / 39856141.0 * pow(z, 5) + 24904180.0 / 39856141.0 * pow(z, 6) + 68265668.0 / 39856141.0 * pow(z, 7) + 1097897094.0 / 1235540371.0 * pow(z, 8) + 1941532102.0 / 873918799.0 * pow(z, 9) + 1341874226.0 / 1235540371.0 * pow(z, 10) + 102398502.0 / 39856141.0 * pow(z, 11) + 3557740.0 / 3065857.0 * pow(z, 12) + 106258124.0 / 39856141.0 * pow(z, 13) + 1153524.0 / 1077193.0 * pow(z, 14) + 202580.0 / 82861.0 * pow(z, 15) + 8143.0 / 10105.0 * pow(z, 16) + 441.0 / 235.0 * pow(z, 17) + 19.0 / 47.0 * pow(z, 18) + 1.0 * pow(z, 19));

                        break;
                }

                return result;
            }

            virtual double f_p(const double & q2) const
            {
                const auto z = _z(q2);

                const auto phi      = this->phi_p(z);
                const auto blaschke = this->blaschke_p(q2);
                const auto series   = this->series_p(z);
                double correction_to_asymptotics;
                double K;
                switch (Process_::asymptotic_case_switch)
                {   
                    case 0: // The asymptotics used by Buck and Lebed
                        correction_to_asymptotics = 1.0;
                        K = (5 * M_PI) / 64;
                        break;
                    case 1: // correct asyptotic behaviour
                        correction_to_asymptotics = pow(1.0 - z, 5.0 / 2.0 - 1.0 / 2.0);
                        K = (315 * M_PI) / 8192;
                        break;
                    default:
                        correction_to_asymptotics = 1.0;
                        K = (5 * M_PI) / 64;
                        break;
                }

                return 1.0 / (phi * blaschke) * correction_to_asymptotics * sqrt(K) * series;
                // return 1.0 / (phi * blaschke) /* asymptotics*/ * sqrt(K) * series;
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