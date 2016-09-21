/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_BASE_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_BASE_HH 1

#include <eos/utils/model.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton.hh>

namespace eos
{
    class BToKstarDilepton::AmplitudeGenerator :
        public ParameterUser
    {
        public:
            std::shared_ptr<Model> model;
            std::shared_ptr<FormFactors<PToV>> form_factors;

            UsedParameter mu;
            UsedParameter alpha_e;
            UsedParameter g_fermi;
            UsedParameter tau;

            UsedParameter m_B;
            UsedParameter m_Kstar;
            UsedParameter m_l;

            bool cp_conjugate;
            std::string lepton_flavour;

            AmplitudeGenerator(const Parameters &, const Options &);

            double s_hat(const double & q2) const;
            double beta_l(const double & q2) const;
            double energy(const double & q2) const;
            double lambda(const double & q2) const;

            virtual ~AmplitudeGenerator();
            virtual BToKstarDilepton::Amplitudes amplitudes(const double & q2) const = 0;
    };

    struct BToKstarDilepton::DipoleFormFactors
    {
        complex<double> calT_perp_left;
        complex<double> calT_perp_right;
        complex<double> calT_parallel;
    };

    template <typename Tag_> class BToKstarDileptonAmplitudes;

    namespace tag
    {
        struct ABBBSW2008;
        struct CFFMPSV2015;
        struct BFS2004;
        struct GP2004;
    }
}

#endif
