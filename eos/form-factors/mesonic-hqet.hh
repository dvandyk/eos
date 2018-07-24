/* vim: set sw=4 sts=4 tw=120 et foldmethod=syntax : */

#ifndef MASTER_GUARD_EOS_FORM_FACTORS_MESONIC_HQET_HH
#define MASTER_GUARD_EOS_FORM_FACTORS_MESONIC_HQET_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/utils/derivative.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/polylog.hh>
#include <eos/utils/power_of.hh>

#include <limits>

#include <iostream>

namespace eos
{
    /* HQET Form Factors, based on [BLPR2017] and [JS2018] */
    template <typename Process_, typename Transition_> class HQETFormFactors;

    class HQETFormFactorBase :
        public ParameterUser
    {
        protected:
            std::shared_ptr<Model> _model;

            // parameters for the leading Isgur-Wise function xi
            UsedParameter _rho2, _c;

            // parameters for the subleading Isgur-Wise function chi_2
            UsedParameter _chi2one, _chi2pone;

            // parameters for the subleading Isgur-Wise function chi_3
            UsedParameter _chi3pone;

            // parameters for the subleading Isgur-Wise function eta
            UsedParameter _etaone, _etapone;

            // parameters for subsubleading 1/m_c corrections in h_A1 (B->D^*)
            UsedParameter _delta_a1;

            // parameters for subsubleading 1/m_c corrections in h_+ (B->D)
            UsedParameter _delta_fp;

        public:
            HQETFormFactorBase(const Parameters & p, const Options & o) :
                _model(Model::make("SM", p, o)),
                _rho2(p["B(*)->D(*)::rho^2@HQET"], *this),
                _c(p["B(*)->D(*)::c@HQET"], *this),
                _chi2one(p["B(*)->D(*)::chi_2(1)@HQET"], *this),
                _chi2pone(p["B(*)->D(*)::chi_2'(1)@HQET"], *this),
                _chi3pone(p["B(*)->D(*)::chi_3'(1)@HQET"], *this),
                _etaone(p["B(*)->D(*)::eta(1)@HQET"], *this),
                _etapone(p["B(*)->D(*)::eta'(1)@HQET"], *this),
                _delta_a1(p["B(*)->D(*)::delta_a1@HQET"], *this),
                _delta_fp(p["B(*)->D(*)::delta_f+@HQET"], *this)
            {
            }

            ~HQETFormFactorBase() = default;

        protected:
            /*
             * HQET parameters following [BLPR2017]
             */
            inline double _mu() const { return 2.31; } // mu^2 = m_b * m_c
            inline double _alpha_s() const { return 0.26; }
            inline double _m_b_1S() const { return 4.71; }
            inline double _m_c_1S() const { return 1.38; }
            inline double _m_b_pole() const { return _m_b_1S() * (1 + 2.0 / 9.0 * power_of<2>(_alpha_s())); }
            inline double _lambda_1() const { return -0.30; }
            inline double _LambdaBar() const { return 5.313 - _m_b_pole() + _lambda_1() / (2.0 * _m_b_1S()); }

            /*
             * Interface to Process_-specific kinematics.
             */
            virtual double _w(const double & q2) const = 0;
            virtual double _q2(const double & w) const = 0;
            virtual double _z(const double & q2) const = 0;

            /*
             * Isgur-Wise functions
             */
            double _xi(const double & q2) const
            {
                const double z = _z(q2);

                return 1.0 - 8.0 * _rho2 * z + (64.0 * _c - 16.0 * _rho2) * z * z;
            }

            double _chi2(const double & q2) const
            {
                const double w = _w(q2);

                return _chi2one + _chi2pone * (w - 1.0);
            }

            double _chi3(const double & q2) const
            {
                const double w = _w(q2);

                return 0.0 + _chi3pone * (w - 1.0);
            }

            double _eta(const double & q2) const
            {
                const double w = _w(q2);

                return _etaone + _etapone * (w - 1.0);
            }

            /*
             * Auxilliary functions for the HQET Wilson coefficients
             *
             * We use a fixed scale mu = sqrt(m_b * m_c), with m_b = 4.2 and m_c = 1.27,
             * which yields mu = 2.31 GeV.
             */

            inline double _wz(const double & z) const
            {
                return 0.5 * (z + 1.0 / z);
            }

            inline double _wp(const double & w) const { return w + std::sqrt(w * w - 1.0); }
            inline double _wm(const double & w) const { return w - std::sqrt(w * w - 1.0); }

            double _r(const double & w) const
            {
                if (w < 1.0)
                    return std::numeric_limits<double>::quiet_NaN();

                if (w - 1.0 < 1.0e-5)
                    return 1.0 - (w - 1.0) / 3.0;

                return std::log(_wp(w)) / std::sqrt(w * w - 1.0);
            }

            inline double _Omega(const double & w, const double & z) const
            {
                if (w < 1.0)
                    return std::numeric_limits<double>::quiet_NaN();

                const double lnz = std::log(z);

                if (w - 1.0 < 1.0e-5)
                    return -1.0 - (1.0 + z) / (1.0 - z) * lnz;

                const double wm = _wm(w);
                const double wp = _wp(w);

                const complex<double> li2wmz = dilog(1.0 - wm * z);
                const complex<double> li2wpz = dilog(1.0 - wp * z);
                const complex<double> li2wm2 = dilog(1.0 - wm * wm);
                const complex<double> li2wp2 = dilog(1.0 - wp * wp);

                return w * real(2.0 * (li2wmz - li2wpz) + li2wp2 - li2wm2) / (2.0 * std::sqrt(w * w - 1.0))
                    - w * _r(w) * lnz + 1.0;
            }

            /* Wilson Coefficients */

            inline double _CV1(const double & w, const double & z) const
            {
                const double z2  = z * z;
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = 2.0 * (w + 1.0) * ((3.0 * w - 1.0) * z - z2 - 1.0) * _r(w);
                result += (12.0 * z * (wz - w) - (z2 - 1.0) * lnz);
                result += 4.0 * z * (w - wz) * _Omega(w, z);

                return result / (6.0 * z * (w - wz));
            }

            inline double _CV2(const double & w, const double & z) const
            {
                const double z2  = z * z, z3 = z2 * z;
                const double w2  = w * w;
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = ((4.0 * w2 + 2.0 * w) * z2 - (2.0 * w2 + 5.0 * w - 1.0) * z - (1.0 + w) * z3 + 2.0) * _r(w);
                result += z * (2.0 * (z - 1.0) * (wz - w) + (z2 - (4.0 * w - 2.0) * z + (-2.0 * w + 3)) * lnz);

                return -1.0 * result / (6.0 * z2 * power_of<2>(w - wz));
            }

            inline double _CV3(const double & w, const double & z) const
            {
                const double z2  = z * z, z3 = z2 * z;
                const double w2  = w * w;
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = (-2.0 * z3 + (2.0 * w2 + 5.0 * w - 1.0) * z2 - (4.0 * w2 + 2.0 * w) * z + w + 1.0) * _r(w);
                result += 2.0 * z * (z - 1.0) * (wz - w) + ((-2.0 * w + 3.0) * z2 + (-4.0 * w + 2.0) * z + 1.0) * lnz;

                return +1.0 * result / (6.0 * z * power_of<2>(w - wz));
            }

            inline double _CA1(const double & w, const double & z) const
            {
                const double z2  = z * z;
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = 2.0 * (w - 1.0) * ((3.0 * w + 1.0) * z - z2 - 1.0) * _r(w);
                result += (12.0 * z * (wz - w) - (z2 - 1.0) * lnz);
                result += 4.0 * z * (w - wz) * _Omega(w, z);

                return result / (6.0 * z * (w - wz));
            }

            inline double _CA2(const double & w, const double & z) const
            {
                const double z2  = z * z, z3 = z2 * z;
                const double w2  = w * w;
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = ((4.0 * w2 - 2.0 * w) * z2 + (2.0 * w2 - 5.0 * w - 1.0) * z + (1.0 - w) * z3 + 2.0) * _r(w);
                result += z * (2.0 * (z + 1.0) * (wz - w) + (z2 - (4.0 * w + 2.0) * z + (2.0 * w + 3)) * lnz);

                return -1.0 * result / (6.0 * z2 * power_of<2>(w - wz));
            }

            inline double _CA3(const double & w, const double & z) const
            {
                const double z2  = z * z, z3 = z2 * z;
                const double w2  = w * w;
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = (2.0 * z3 + (2.0 * w2 - 5.0 * w - 1.0) * z2 + (4.0 * w2 - 2.0 * w) * z - w + 1.0) * _r(w);
                result += 2.0 * z * (z + 1.0) * (wz - w) - ((2.0 * w + 3.0) * z2 - (4.0 * w + 2.0) * z + 1.0) * lnz;

                return +1.0 * result / (6.0 * z * power_of<2>(w - wz));
            }
    };

    template <typename Process_> class HQETFormFactors<Process_, PToP> :
        public HQETFormFactorBase,
        public FormFactors<PToP>
    {
        private:
            /*
             * Kinematics
             */
            virtual double _w(const double & q2) const override
            {
                static constexpr double mB = Process_::m_B, mB2 = power_of<2>(Process_::m_B);
                static constexpr double mP = Process_::m_P, mP2 = power_of<2>(Process_::m_P);

                return (mB2 + mP2 - q2) / (2.0 * mB * mP);
            }

            virtual double _q2(const double & w) const override
            {
                static constexpr double mB = Process_::m_B, mB2 = power_of<2>(Process_::m_B);
                static constexpr double mP = Process_::m_P, mP2 = power_of<2>(Process_::m_P);

                return mB2 + mP2 - 2.0 * mB * mP * w;
            }

            virtual double _z(const double & q2) const override
            {
                const double w = _w(q2);

                return (std::sqrt(w + 1.0) - std::sqrt(2.0)) / (std::sqrt(w + 1.0) + std::sqrt(2.0));
            }

            /* HQET form factors h_i */

            double _h_p(const double & q2) const
            {
                const double m_b_1S = _m_b_1S();
                const double m_c_1S = _m_c_1S();

                const double w = this->_w(q2);
                const double z = m_c_1S / m_b_1S;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_1S);
                const double eps_c = _LambdaBar() / (2.0 * m_c_1S);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;

                double result = (1.0 + as * (_CV1(w, z) + (w + 1.0) / 2.0 * (_CV2(w, z) + _CV3(w, z))));
                result += eps_c * (L1);
                result += eps_b * (L1);
                result += eps_c * eps_c * _delta_fp;

                return result * xi;
            }

            double _h_m(const double & q2) const
            {
                const double m_b_1S = _m_b_1S();
                const double m_c_1S = _m_c_1S();

                const double w = this->_w(q2);
                const double z = m_c_1S / m_b_1S;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_1S);
                const double eps_c = _LambdaBar() / (2.0 * m_c_1S);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L4 = 2.0 * eta - 1.0;

                double result = (0.0 + as * (w + 1.0) / 2.0 * (_CV2(w, z) - _CV3(w, z)));
                result += eps_c * L4;
                result -= eps_b * L4;

                return result * xi;
            }

        public:
            HQETFormFactors(const Parameters & p, const Options & o) :
                HQETFormFactorBase(p, o)
            {
            }

            ~HQETFormFactors() = default;

            static FormFactors<PToP> * make(const Parameters & parameters, unsigned)
            {
                return new HQETFormFactors(parameters, Options());
            }

            virtual double f_p(const double & q2) const
            {
                constexpr double r = Process_::m_P / Process_::m_B;

                // cf. [FKKM2008], eq. (22)
                return 1.0 / (2.0 * sqrt(r)) * ((1.0 + r) * _h_p(q2) - (1.0 - r) * _h_m(q2));
            }

            double f_m(const double & q2) const
            {
                constexpr double r = Process_::m_P / Process_::m_B;

                // cf. [FKKM2008], eq. (22)
                return 1.0 / (2.0 * sqrt(r)) * ((1.0 + r) * _h_m(q2) - (1.0 - r) * _h_p(q2));
            }

            virtual double f_0(const double & q2) const
            {
                return -1.0;
            }

            virtual double f_t(const double & q2) const
            {
                return -1.0;
            }

            Diagnostics diagnostics() const
            {
                Diagnostics results;

                // Inputs
                {
                    const double m_b = _model->m_b_msbar(this->_mu());
                    const double m_c = _model->m_c_msbar(this->_mu());
                    const double z   = m_c / m_b;
                    const double wz  = _wz(z);

                    results.add(Diagnostics::Entry{ z,  "z = m_c / m_b" });
                    results.add(Diagnostics::Entry{ wz, "w_z"           });
                }

                // xi
                {
                    results.add(Diagnostics::Entry{ _xi(_q2(1.10)), "xi(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.05)), "xi(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.00)), "xi(w = 1.00)" });
                }

                // chi2
                {
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.10)), "chi2(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.05)), "chi2(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.00)), "chi2(w = 1.00)" });
                }

                // chi3
                {
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.10)), "chi3(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.05)), "chi3(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.00)), "chi3(w = 1.00)" });
                }

                // eta
                {
                    results.add(Diagnostics::Entry{ _eta(_q2(1.10)), "eta(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.05)), "eta(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.00)), "eta(w = 1.00)" });
                }

                // r(w)
                {
                    results.add(Diagnostics::Entry{ _r(1.1),     "r(w = 1.1)"     });
                    results.add(Diagnostics::Entry{ _r(1.0007),  "r(w = 1.0007)"  });
                    results.add(Diagnostics::Entry{ _r(1.0001),  "r(w = 1.0001)"  });
                    results.add(Diagnostics::Entry{ _r(1.00005), "r(w = 1.00005)" });
                    results.add(Diagnostics::Entry{ _r(1.0),     "r(w = 1.0)"     });
                }

                // Omega(w, z = 0.25)
                {
                    results.add(Diagnostics::Entry{ _Omega(1.1,     0.25), "Omega(w = 1.1,     z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0007,  0.25), "Omega(w = 1.0007,  z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0001,  0.25), "Omega(w = 1.0001,  z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.00005, 0.25), "Omega(w = 1.00005, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0,     0.25), "Omega(w = 1.0,     z = 0.25)" });
                }

                // Omega(w, z = 0.20)
                {
                    results.add(Diagnostics::Entry{ _Omega(1.1,     0.20), "Omega(w = 1.1,     z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0007,  0.20), "Omega(w = 1.0007,  z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0001,  0.20), "Omega(w = 1.0001,  z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.00005, 0.20), "Omega(w = 1.00005, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0,     0.20), "Omega(w = 1.0,     z = 0.20)" });
                }

                // WCs at w = 1.2, z = 0.20
                {
                    results.add(Diagnostics::Entry{ _CV1(1.2, 0.20), "C_{V_1}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CV2(1.2, 0.20), "C_{V_2}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CV3(1.2, 0.20), "C_{V_3}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA1(1.2, 0.20), "C_{A_1}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA2(1.2, 0.20), "C_{A_2}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA3(1.2, 0.20), "C_{A_3}(w = 1.2, z = 0.20)" });
                }

                // WCs at w = 1.0, z = 0.25
                {
                    results.add(Diagnostics::Entry{ _CV1(1.0, 0.25), "C_{V_1}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CV2(1.0, 0.25), "C_{V_2}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CV3(1.0, 0.25), "C_{V_3}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA1(1.0, 0.25), "C_{A_1}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA2(1.0, 0.25), "C_{A_2}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA3(1.0, 0.25), "C_{A_3}(w = 1.0, z = 0.25)" });
                }

                // HQET definition of the form factors
                {
                    results.add(Diagnostics::Entry{ _h_p(_q2(1.4)), "h_p(w = 1.4)" });
                    results.add(Diagnostics::Entry{ _h_m(_q2(1.4)), "h_m(w = 1.4)" });

                    results.add(Diagnostics::Entry{ _h_p(_q2(1.2)), "h_p(w = 1.2)" });
                    results.add(Diagnostics::Entry{ _h_m(_q2(1.2)), "h_m(w = 1.2)" });

                    results.add(Diagnostics::Entry{ _h_p(_q2(1.0)), "h_p(w = 1.0)" });
                    results.add(Diagnostics::Entry{ _h_m(_q2(1.0)), "h_m(w = 1.0)" });
                }

                return results;
            }
    };

    template <typename Process_> class HQETFormFactors<Process_, PToV> :
        public HQETFormFactorBase,
        public FormFactors<PToV>
    {
        private:
            /*
             * Kinematics
             */
            virtual double _w(const double & q2) const override
            {
                static constexpr double mB = Process_::mB, mB2 = power_of<2>(Process_::mB);
                static constexpr double mV = Process_::mV, mV2 = power_of<2>(Process_::mV);

                return (mB2 + mV2 - q2) / (2.0 * mB * mV);
            }

            virtual double _q2(const double & w) const override
            {
                static constexpr double mB = Process_::mB, mB2 = power_of<2>(Process_::mB);
                static constexpr double mV = Process_::mV, mV2 = power_of<2>(Process_::mV);

                return mB2 + mV2 - 2.0 * mB * mV * w;
            }

            virtual double _z(const double & q2) const override
            {
                const double w = _w(q2);

                return (std::sqrt(w + 1.0) - std::sqrt(2.0)) / (std::sqrt(w + 1.0) + std::sqrt(2.0));
            }

            /* HQET form factors h_i */

            double _h_a1(const double & q2) const
            {
                const double m_b_1S = _m_b_1S();
                const double m_c_1S = _m_c_1S();

                const double w = this->_w(q2);
                const double z = m_c_1S / m_b_1S;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_1S);
                const double eps_c = _LambdaBar() / (2.0 * m_c_1S);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
                const double L2 = -4.0 * chi3;
                const double L4 = 2.0 * eta - 1.0;
                const double L5 = -1.0;

                double result = (1.0 + as * _CA1(w, z));
                result += eps_c * (L2 - L5 * (w - 1.0) / (w + 1.0));
                result += eps_b * (L1 - L4 * (w - 1.0) / (w + 1.0));
                result += eps_c * eps_c * _delta_a1;

                return result * xi;
            }

            double _h_a2(const double & q2) const
            {
                const double m_b_1S = _m_b_1S();
                const double m_c_1S = _m_c_1S();

                const double w = this->_w(q2);
                const double z = m_c_1S / m_b_1S;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_1S);
                const double eps_c = _LambdaBar() / (2.0 * m_c_1S);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L3 = 4.0 * chi2;
                const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

                double result = (0.0 + as * _CA2(w, z));
                result += eps_c * (L3 + L6);

                return result * xi;
            }

            double _h_a3(const double & q2) const
            {
                const double m_b_1S = _m_b_1S();
                const double m_c_1S = _m_c_1S();

                const double w = this->_w(q2);
                const double z = m_c_1S / m_b_1S;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_1S);
                const double eps_c = _LambdaBar() / (2.0 * m_c_1S);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
                const double L2 = -4.0 * chi3;
                const double L3 = 4.0 * chi2;
                const double L4 = 2.0 * eta - 1.0;
                const double L5 = -1.0;
                const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

                double result = (1.0 + as * (_CA1(w, z) +_CA3(w, z)));
                result += eps_c * (L2 - L3 + L6 - L5);
                result += eps_b * (L1 - L4);

                return result * xi;
            }

            double _h_v(const double & q2) const
            {
                auto print = [](const std::string & l, const double & x) -> const double &
                {
                    std::cout << l << ": " << x << std::endl;
                };

                const double m_b_1S = _m_b_1S();
                const double m_c_1S = _m_c_1S();

                const double w = this->_w(q2);
                const double z = m_c_1S / m_b_1S;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_1S);
                const double eps_c = _LambdaBar() / (2.0 * m_c_1S);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
                const double L2 = -4.0 * chi3;
                const double L4 = 2.0 * eta - 1.0;
                const double L5 = -1.0;

                double result = (1.0 + as * _CV1(w, z));
                result += eps_c * (L2 - L5);
                result += eps_b * (L1 - L4);

                return result * xi;
            }

        public:
            HQETFormFactors(const Parameters & p, const Options & o) :
                HQETFormFactorBase(p, o)
            {
            }

            ~HQETFormFactors() = default;

            static FormFactors<PToV> * make(const Parameters & parameters, unsigned)
            {
                return new HQETFormFactors(parameters, Options());
            }

            virtual double v(const double & q2) const
            {
                constexpr double r = Process_::mV / Process_::mB;

                // cf. [FKKM2008], eq. (22)
                return (1.0 + r) / 2.0 * sqrt(r) * _h_v(q2);
            }

            virtual double a_0(const double & q2) const
            {
                constexpr double r = Process_::mV / Process_::mB;
                const     double w = _w(q2);

                // cf. [FKKM2008], eq. (22)
                const double a_30 = (1.0 + r * r - 2.0 * r * w) / (4.0 * r * sqrt(r)) * (r * _h_a2(q2) - _h_a3(q2));
                return a_3(q2) - a_30;
            }

            virtual double a_1(const double & q2) const
            {
                constexpr double r = Process_::mV / Process_::mB;
                const     double w = _w(q2);

                // cf. [FKKM2008], eq. (22)
                return sqrt(r) * (1.0 + w) / (1 + r) * _h_a1(q2);
            }

            virtual double a_2(const double & q2) const
            {
                constexpr double r = Process_::mV / Process_::mB;
                const     double w = _w(q2);

                // cf. [FKKM2008], eq. (22)
                return (1.0 + r) / (2.0 * sqrt(r)) * (r * _h_a2(q2) + _h_a3(q2));
            }

            double a_3(const double & q2) const
            {
                constexpr double r = Process_::mV / Process_::mB;

                // cf. [FKKM2008], below eq. (6)
                return ((1.0 + r) * a_1(q2) - (1.0 - r) * a_2(q2)) / (2.0 * r);
            }

            virtual double a_12(const double & q2) const
            {
                constexpr double mB = Process_::mB, mB2 = mB * mB;
                constexpr double mV = Process_::mV, mV2 = mV * mV;

                double result = (mB + mV) * (mB + mV) * (mB2 - mV2 - q2) * a_1(q2) - eos::lambda(mB2, mV2, q2) * a_2(q2);
                result /= 8.0 * mB * mV2 * (mB - mV);

                return result;
            }

            virtual double t_1(const double & s) const
            {
                return 0.0;
            }

            virtual double t_2(const double & s) const
            {
                return 0.0;
            }

            virtual double t_3(const double & s) const
            {
                return 0.0;
            }

            virtual double t_23(const double & s) const
            {
                return 0.0;
            }

            Diagnostics diagnostics() const
            {
                Diagnostics results;

                // Inputs
                {
                    const double m_b = _model->m_b_msbar(this->_mu());
                    const double m_c = _model->m_c_msbar(this->_mu());
                    const double z   = m_c / m_b;
                    const double wz  = _wz(z);

                    results.add(Diagnostics::Entry{ z,  "z = m_c / m_b" });
                    results.add(Diagnostics::Entry{ wz, "w_z"           });
                }

                // xi
                {
                    results.add(Diagnostics::Entry{ _xi(_q2(1.40)), "xi(w = 1.40)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.20)), "xi(w = 1.20)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.10)), "xi(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.05)), "xi(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.00)), "xi(w = 1.00)" });
                }

                // chi2
                {
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.10)), "chi2(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.05)), "chi2(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.00)), "chi2(w = 1.00)" });
                }

                // chi3
                {
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.10)), "chi3(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.05)), "chi3(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.00)), "chi3(w = 1.00)" });
                }

                // eta
                {
                    results.add(Diagnostics::Entry{ _eta(_q2(1.10)), "eta(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.05)), "eta(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.00)), "eta(w = 1.00)" });
                }

                // r(w)
                {
                    results.add(Diagnostics::Entry{ _r(1.1),     "r(w = 1.1)"     });
                    results.add(Diagnostics::Entry{ _r(1.0007),  "r(w = 1.0007)"  });
                    results.add(Diagnostics::Entry{ _r(1.0001),  "r(w = 1.0001)"  });
                    results.add(Diagnostics::Entry{ _r(1.00005), "r(w = 1.00005)" });
                    results.add(Diagnostics::Entry{ _r(1.0),     "r(w = 1.0)"     });
                }

                // Omega(w, z = 0.25)
                {
                    results.add(Diagnostics::Entry{ _Omega(1.1,     0.25), "Omega(w = 1.1,     z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0007,  0.25), "Omega(w = 1.0007,  z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0001,  0.25), "Omega(w = 1.0001,  z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.00005, 0.25), "Omega(w = 1.00005, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0,     0.25), "Omega(w = 1.0,     z = 0.25)" });
                }

                // Omega(w, z = 0.20)
                {
                    results.add(Diagnostics::Entry{ _Omega(1.1,     0.20), "Omega(w = 1.1,     z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0007,  0.20), "Omega(w = 1.0007,  z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0001,  0.20), "Omega(w = 1.0001,  z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.00005, 0.20), "Omega(w = 1.00005, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0,     0.20), "Omega(w = 1.0,     z = 0.20)" });
                }

                // WCs at w = 1.2, z = 0.20
                {
                    results.add(Diagnostics::Entry{ _CV1(1.2, 0.20), "C_{V_1}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CV2(1.2, 0.20), "C_{V_2}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CV3(1.2, 0.20), "C_{V_3}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA1(1.2, 0.20), "C_{A_1}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA2(1.2, 0.20), "C_{A_2}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA3(1.2, 0.20), "C_{A_3}(w = 1.2, z = 0.20)" });
                }

                // WCs at w = 1.0, z = 0.25
                {
                    results.add(Diagnostics::Entry{ _CV1(1.0, 0.25), "C_{V_1}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CV2(1.0, 0.25), "C_{V_2}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CV3(1.0, 0.25), "C_{V_3}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA1(1.0, 0.25), "C_{A_1}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA2(1.0, 0.25), "C_{A_2}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA3(1.0, 0.25), "C_{A_3}(w = 1.0, z = 0.25)" });
                }

                // HQET definition of the form factors
                {
                    results.add(Diagnostics::Entry{ _h_a1(_q2(1.4)), "h_A1(w = 1.4)" });
                    results.add(Diagnostics::Entry{ _h_a2(_q2(1.4)), "h_A2(w = 1.4)" });
                    results.add(Diagnostics::Entry{ _h_a3(_q2(1.4)), "h_A3(w = 1.4)" });
                    results.add(Diagnostics::Entry{ _h_v (_q2(1.4)), "h_V (w = 1.4)" });

                    results.add(Diagnostics::Entry{ _h_a1(_q2(1.2)), "h_A1(w = 1.2)" });
                    results.add(Diagnostics::Entry{ _h_a2(_q2(1.2)), "h_A2(w = 1.2)" });
                    results.add(Diagnostics::Entry{ _h_a3(_q2(1.2)), "h_A3(w = 1.2)" });
                    results.add(Diagnostics::Entry{ _h_v (_q2(1.2)), "h_V (w = 1.2)" });

                    results.add(Diagnostics::Entry{ _h_a1(_q2(1.0)), "h_A1(w = 1.0)" });
                    results.add(Diagnostics::Entry{ _h_a2(_q2(1.0)), "h_A2(w = 1.0)" });
                    results.add(Diagnostics::Entry{ _h_a3(_q2(1.0)), "h_A3(w = 1.0)" });
                    results.add(Diagnostics::Entry{ _h_v (_q2(1.0)), "h_V (w = 1.0)" });
                }

                return results;
            }
    };
}

#endif
