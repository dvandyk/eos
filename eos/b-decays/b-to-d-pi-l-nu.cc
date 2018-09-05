/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <eos/b-decays/b-to-d-pi-l-nu.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    template <>
    struct Implementation<BToDPiLeptonNeutrino>
    {
        // angular observables of the P-wave D-pi resonance (D*)
        UsedParameter a_fb, a_t2, a_im, f_l;

        // angular observables of the S-wave D-pi background, and the
        // interference with the P-wave resonance
        UsedParameter a_s, f_s;

        Implementation(const Parameters & p, const Options &, ParameterUser & u) :
            a_fb(p["B->Dpilnu::A_FB"], u),
            a_t2(p["B->Dpilnu::A_T^2"], u),
            a_im(p["B->Dpilnu::A_im"], u),
            f_l(p["B->Dpilnu::F_L"], u),
            a_s(p["B->Dpilnu::A_S"], u),
            f_s(p["B->Dpilnu::F_S"], u)
        {
        }

        double pdf_d(const double & c_d) const
        {
            const double c_d2 = c_d * c_d;

            return f_s / 2.0 + a_s * c_d
                + (1.0 - f_s) * (
                    3.0 / 2.0 * c_d2 * f_l
                    - 3.0 / 4.0 * (c_d2 - 1.0) * (1.0 - f_l)
                );
        }

        double pdf_l(const double & c_l) const
        {
            const double c_l2 = c_l * c_l;

            return -3.0 / 4.0 * (c_l2 - 1.0) * f_s
                + (1.0 - f_s) * (
                    a_fb * c_l
                    - 3.0 / 4.0 * (c_l2 - 1.0) * f_l
                    + 3.0 / 8.0 * (c_l2 + 1.0) * (1.0 - f_l)
                );
        }

        double pdf_phi(const double & phi) const
        {
            const double c2 = std::cos(2.0 * phi), s2 = std::sin(2.0 * phi);

            double result = f_s;
            result += (1.0 - f_s) * (a_t2 * (1.0 - f_l) * c2 / 2.0 + 1.0 + a_im * s2);
            result /= 2.0 * M_PI;

            return result;
        }
    };

    BToDPiLeptonNeutrino::BToDPiLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToDPiLeptonNeutrino>(new Implementation<BToDPiLeptonNeutrino>(parameters, options, *this))
    {
    }

    BToDPiLeptonNeutrino::~BToDPiLeptonNeutrino()
    {
    }

    double
    BToDPiLeptonNeutrino::differential_pdf_d(const double & c_d) const
    {
        return _imp->pdf_d(c_d);
    }

    double
    BToDPiLeptonNeutrino::differential_pdf_l(const double & c_l) const
    {
        return _imp->pdf_l(c_l);
    }

    double
    BToDPiLeptonNeutrino::differential_pdf_phi(const double & phi) const
    {
        return _imp->pdf_phi(phi);
    }

    double
    BToDPiLeptonNeutrino::integrated_pdf_d(const double & c_d_min, const double & c_d_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_d, _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, c_d_min, c_d_max);
    }

    double
    BToDPiLeptonNeutrino::integrated_pdf_l(const double & c_l_min, const double & c_l_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_l, _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, c_l_min, c_l_max);
    }

    double
    BToDPiLeptonNeutrino::integrated_pdf_phi(const double & phi_min, const double & phi_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_phi, _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, phi_min, phi_max);
    }

    const std::string
    BToDPiLeptonNeutrino::description = "\
The decay B->D pi l nu, where l is a massless lepton.";

    const std::string
    BToDPiLeptonNeutrino::kinematics_description_c_d = "\
The cosine of the helicity angles theta_D in the D-pi rest frame.";

    const std::string
    BToDPiLeptonNeutrino::kinematics_description_c_l = "\
The cosine of the helicity angles theta_L in the l-nu rest frame.";

    const std::string
    BToDPiLeptonNeutrino::kinematics_description_phi = "\
The azimuthal angles between the decay planes.";

}
