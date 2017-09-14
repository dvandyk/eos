/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef BTOGAMMA_FF_GUARD_EOS_B_DECAYS_B_TO_GAMMA_L_NU_HH
#define BTOGAMMA_FF_GUARD_EOS_B_DECAYS_B_TO_GAMMA_L_NU_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    /*
     * Decay: B -> gamma l nu
     */
    class BToGammaLeptonNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<BToGammaLeptonNeutrino>
    {
        public:
            BToGammaLeptonNeutrino(const Parameters & parameters, const Options & options);
            ~BToGammaLeptonNeutrino();

            // Differential Observables
            double differential_branching_ratio(const double & Egamma) const;

            // Integrated Observables
            double integrated_branching_ratio(const double & Egamma_min, const double & Egamma_max) const;
            double integrated_forward_backward_asymmetry(const double & Egamma_min, const double & Egamma_max) const;
            double integrated_photon_energy_moment_1(const double & Egamma_min, const double & Egamma_max) const;
            double integrated_photon_energy_moment_2(const double & Egamma_min, const double & Egamma_max) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_Egamma;
    };
}


#endif
