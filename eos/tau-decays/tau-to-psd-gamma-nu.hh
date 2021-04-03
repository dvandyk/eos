/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_EOS_TAU_DECAYS_TAU_TO_PSD_GAMMA_NU_HH
#define EOS_GUARD_EOS_TAU_DECAYS_TAU_TO_PSD_GAMMA_NU_HH 1

#include <eos/decays.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    class TauToPseudoscalarGammaNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<TauToPseudoscalarGammaNeutrino>
        
    {
        public:
            TauToPseudoscalarGammaNeutrino(const Parameters & parameters, const Options & options);
            ~TauToPseudoscalarGammaNeutrino();

            // Observables
            double differential_ratio(const double & E_gamma, const double & E_pi) const;
          /*  double branching_ratio(const double & q2) const;*/
            double dummy(const double & ) const;
              
            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;
            static const std::string description;
            static const std::string kinematics_description_E_gamma;
            static const std::string kinematics_description_E_pi;
    };
}

#endif
