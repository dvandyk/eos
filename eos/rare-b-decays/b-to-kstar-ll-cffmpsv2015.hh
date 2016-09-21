/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_CFFMPSV2015_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_CFFMPSV2015_HH 1

#include <eos/rare-b-decays/b-to-kstar-ll-base.hh>
#include <eos/rare-b-decays/b-to-kstar-ll-cffmpsv2015.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>

namespace eos
{
    struct Common_Pararameters
    {
        BToKstarDilepton::Amplitudes amp;
        complex<double> Y_top_c ,
            Y_top_b ,
            Y_top_0 ,
            Y_top_ ,
            Y_top,
            Y_up;
        complex<double>    c7eff,
            c9eff ,
            c910_mi_r ,
            c910_mi_l ,
            c910_pl_r ,
            c910_pl_l ,
            // here b-qu
            c7_mi ,
            c7_pl ;
    };

    /* store intermediate results for potential reuse
     * -> mostly filled in 'amp_semileptonic', which is always called
     */
    struct Commons_buffer
    {
        double
            shat, sqrt_s,
            mKhat2, m_K2,
            m_B2,
            m_sum, m_diff, m2_diff,
            norm_s,
            lam_s, sqrt_lam;

        double
            ff_V,
            ff_A0, ff_A1, ff_A2,
            ff_T1, ff_T2, ff_T3;

        double
            pre_long, pre_perp, pre_par,
            kin_long_1, kin_long_2;

        complex<double>
            c910_mi_r, c910_mi_l,
            c910_pl_r, c910_pl_l,
            // here b-qu
            c7_mi, c7_pl;
    };

    template <>
    class BToKstarDileptonAmplitudes<tag::CFFMPSV2015> :
        public BToKstarDilepton::AmplitudeGenerator
    {
        public:
            UsedParameter hbar;

            UsedParameter m_b_MSbar;
            UsedParameter m_c;
            UsedParameter m_s_MSbar;

            UsedParameter mu;
            UsedParameter alpha_e;
            UsedParameter g_fermi;
            UsedParameter tau;

            UsedParameter f_B;
            UsedParameter f_Kstar_par;
            UsedParameter f_Kstar_perp;
            UsedParameter lambda_B_p;
            UsedParameter a_1_par;
            UsedParameter a_2_par;
            UsedParameter a_1_perp;
            UsedParameter a_2_perp;

            UsedParameter uncertainty_xi_perp;
            UsedParameter uncertainty_xi_par;

            UsedParameter re_h_00;
            UsedParameter im_h_00;
            UsedParameter re_h_01;
            UsedParameter im_h_01;
            UsedParameter re_h_02;
            UsedParameter im_h_02;
            UsedParameter re_h_p0;
            UsedParameter im_h_p0;
            UsedParameter re_h_p1;
            UsedParameter im_h_p1;
            UsedParameter re_h_p2;
            UsedParameter im_h_p2;
            UsedParameter re_h_m0;
            UsedParameter im_h_m0;
            UsedParameter re_h_m1;
            UsedParameter im_h_m1;
            UsedParameter re_h_m2;
            UsedParameter im_h_m2;

            UsedParameter abs_h_00;
            UsedParameter arg_h_00;
            UsedParameter abs_h_01;
            UsedParameter arg_h_01;
            UsedParameter abs_h_02;
            UsedParameter arg_h_02;
            UsedParameter abs_h_p0;
            UsedParameter arg_h_p0;
            UsedParameter abs_h_p1;
            UsedParameter arg_h_p1;
            UsedParameter abs_h_p2;
            UsedParameter arg_h_p2;
            UsedParameter abs_h_m0;
            UsedParameter arg_h_m0;
            UsedParameter abs_h_m1;
            UsedParameter arg_h_m1;
            UsedParameter abs_h_m2;
            UsedParameter arg_h_m2;

            double e_q;

            char q;

            std::string coordinates;

            std::function<complex<double> (const double &)> h_0;
            std::function<complex<double> (const double &)> h_p;
            std::function<complex<double> (const double &)> h_m;

            std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                    const double &, const double &, const double &, const double &,
                    const double &, const double &)> qcdf_dilepton_massless_case;
            std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                    const double &, const double &, const double &, const double &,
                    const double &, const double &, const double &)> qcdf_dilepton_charm_case;
            std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                    const double &, const double &, const double &, const double &,
                    const double &, const double &, const double &)> qcdf_dilepton_bottom_case;

            //mutable Commons_buffer cbuf;

            BToKstarDileptonAmplitudes(const Parameters & p, const Options & o);
            ~BToKstarDileptonAmplitudes();

            virtual BToKstarDilepton::Amplitudes amplitudes(const double & q2) const;

            BToKstarDilepton::Amplitudes amp_semileptonic(const double & s, const WilsonCoefficients<BToS> & wc) const;
            WilsonCoefficients<BToS> wilson_coefficients() const;
            double m_b_PS() const;
            double mu_f() const;
            BToKstarDilepton::DipoleFormFactors dipole_form_factors(const double & q2, const WilsonCoefficients<BToS> & wc) const;
            double norm(const double & q2) const;
            double s_hat(const double & q2) const;
            double xi_perp(const double & q2) const;
            double xi_par(const double & q2) const;

            complex<double> h_0_cartesian(const double & q2) const;
            complex<double> h_p_cartesian(const double & q2) const;
            complex<double> h_m_cartesian(const double & q2) const;

            complex<double> h_0_polar(const double & q2) const;
            complex<double> h_p_polar(const double & q2) const;
            complex<double> h_m_polar(const double & q2) const;
    };

}

#endif
