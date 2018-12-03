/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <eos/b-decays/b-to-d-pi-l-nu.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    using std::conj;
    using std::norm;
    using std::real;
    using std::sqrt;

    template <>
    struct Implementation<BToDPiLeptonNeutrino>
    {
        // model
        SwitchOption opt_model;

        // meson masses
        UsedParameter m_B, m_Dstar;

        // lepton mass
        SwitchOption opt_l;
        UsedParameter m_l;

        // form factors
        std::shared_ptr<FormFactors<PToV>> ff;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_model(o, "model", { "SM", "CKMScan" }, "SM"),
            m_B(p["mass::B_d"], u),
            m_Dstar(p["mass::D^*_d"], u),
            opt_l(o, "l", { "e", "mu", "tau" }, "mu"),
            m_l(p["mass::" + opt_l.value()], u),
            ff(FormFactorFactory<PToV>::create("B->D^*::" + o.get("form-factors", "HQET"), p, o))
        {
            if (! ff.get())
                throw InternalError("Form factors not found!");

            u.uses(*ff);
        }

        inline double pdf_normalization(const double & q2) const
        {
            // we only require q2-dependent terms of the normalization, since
            // the constants parts drop out in the PDF.

            const double m_B     = this->m_B(),     m_B2 = m_B * m_B;
            const double m_Dstar = this->m_Dstar(), m_Dstar2 = m_Dstar * m_Dstar;

            const double pDstar = sqrt(eos::lambda(m_B2, m_Dstar2, q2)) / (2.0 * m_B);
            const double beta   = 1.0 - m_l() * m_l() / q2;

            return pDstar * q2 * beta * beta;
        }

        inline double a_long(const double & q2) const
        {
            const double m_B     = this->m_B(),     m_B2 = m_B * m_B;
            const double m_Dstar = this->m_Dstar(), m_Dstar2 = m_Dstar * m_Dstar;
            const double lambda  = eos::lambda(m_B2, m_Dstar2, q2);

            double result = (m_B + m_Dstar) * (m_B2 - m_Dstar2 - q2) * ff->a_1(q2);
            result -= lambda / (m_B + m_Dstar) * ff->a_2(q2);
            result /= (2.0 * m_Dstar * sqrt(q2));

            return result;
        }

        inline double a_para(const double & q2) const
        {
            const double m_B     = this->m_B();
            const double m_Dstar = this->m_Dstar();

            return sqrt(2.0) * (m_B + m_Dstar) * ff->a_1(q2);
        }

        inline double a_perp(const double & q2) const
        {
            const double m_B     = this->m_B(),     m_B2 = m_B * m_B;
            const double m_Dstar = this->m_Dstar(), m_Dstar2 = m_Dstar * m_Dstar;
            const double lambda  = eos::lambda(m_B2, m_Dstar2, q2);

            return -sqrt(2.0) * sqrt(lambda) / (m_B + m_Dstar) * ff->v(q2);
        }

        inline double a_time(const double & q2) const
        {
            const double m_B     = this->m_B(),     m_B2 = m_B * m_B;
            const double m_Dstar = this->m_Dstar(), m_Dstar2 = m_Dstar * m_Dstar;
            const double lambda  = eos::lambda(m_B2, m_Dstar2, q2);

            // cf. [DSD2014], eq. (22), p. 17
            return sqrt(lambda / q2) * ff->a_0(q2);
        }

        double pdf_q2(const double & q2) const
        {
            const double nf = pdf_normalization(q2);

            const double m_l2    = m_l() * m_l();
            const double a_long2 = norm(a_long(q2));
            const double a_para2 = norm(a_para(q2));
            const double a_perp2 = norm(a_perp(q2));
            const double a_time2 = norm(a_time(q2));

            const double a = 2.0 * (a_long2 + a_para2 + a_perp2) * (1.0 + m_l2 / (2.0 * q2))
                           + 3.0 * a_time2 * m_l2 / q2;

            return nf * a;
        }

        double pdf_q2d(const double & q2, const double & c_d) const
        {
            const double nf = pdf_normalization(q2);

            const double m_l2    = m_l() * m_l();

            const double a_long2 = norm(a_long(q2));
            const double a_para2 = norm(a_para(q2));
            const double a_perp2 = norm(a_perp(q2));
            const double a_time2 = norm(a_time(q2));

            const double a = (a_para2 + a_perp2) * (1.0 + m_l2 / (2.0 * q2));
            const double b = (2.0 * a_long2 - a_para2 - a_perp2) * (1.0 + m_l2 / (2.0 * q2))
                           + 3.0 * m_l2 / q2 * a_time2;

            return 3.0 / 2.0 * nf * (a + b * c_d * c_d);
        }

        double pdf_d(const double & c_d) const
        {
            const double q2_min = power_of<2>(m_l());
            const double q2_max = 10.68;

            std::function<double (const double &)> num   = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_q2d, this, std::placeholders::_1, c_d);
            std::function<double (const double &)> denom = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_q2,  this, std::placeholders::_1);

            return integrate<GSL::QAGS>(num, q2_min, q2_max) / integrate<GSL::QAGS>(denom, q2_min, q2_max);
        }

        double pdf_q2l(const double & q2, const double & c_l) const
        {
            const double nf = pdf_normalization(q2);

            const double m_l2         = m_l() * m_l();

            const double a_long       = this->a_long(q2); 
            const double a_para       = this->a_para(q2); 
            const double a_perp       = this->a_perp(q2); 
            const double a_time       = this->a_time(q2); 

            const double a_long2      = norm(a_long);
            const double a_para2      = norm(a_para);
            const double a_perp2      = norm(a_perp);
            const double a_time2      = norm(a_time);

            const double re_para_perp = real(a_para * conj(a_perp));
            const double re_time_long = real(a_time * conj(a_long));

            const double a = 2.0 * a_long2 + (a_para2 + a_perp2) * (1.0 + m_l2 / q2)
                           + 2.0 * a_time2 * m_l2 / q2;
            const double b = -4.0 * (re_para_perp + re_time_long * m_l2 / q2);
            const double c = -(2.0 * a_long2 - a_para2 - a_perp2) * (1.0 - m_l2 / q2);

            return 3.0 / 4.0 * nf * (a + b * c_l + c * c_l * c_l);
        }

        double pdf_l(const double & c_l) const
        {
            const double q2_min = power_of<2>(m_l());
            const double q2_max = 10.68;

            std::function<double (const double &)> num   = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_q2l, this, std::placeholders::_1, c_l);
            std::function<double (const double &)> denom = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_q2,  this, std::placeholders::_1);

            return integrate<GSL::QAGS>(num, q2_min, q2_max) / integrate<GSL::QAGS>(denom, q2_min, q2_max);
        }

        double pdf_q2chi(const double & q2, const double & c_chi) const
        {
            const double c_chi2 = c_chi * c_chi;

            const double nf = pdf_normalization(q2);

            const double m_l2         = m_l() * m_l();

            const double a_long       = this->a_long(q2); 
            const double a_para       = this->a_para(q2); 
            const double a_perp       = this->a_perp(q2); 
            const double a_time       = this->a_time(q2); 

            const double a_long2      = norm(a_long);
            const double a_para2      = norm(a_para);
            const double a_perp2      = norm(a_perp);
            const double a_time2      = norm(a_time);

            const double re_para_long = real(a_para * conj(a_long));
            const double re_para_time = real(a_para * conj(a_time));
            const double re_perp_long = real(a_perp * conj(a_long));

            const double a = 2.0 * a_long2 + 3.0 * a_para2 + a_perp2
                           + m_l2 / q2 * (a_long2 + 2.0 * a_perp2 + 3.0 * a_time2);
            const double b = 3.0 * M_PI / 10.0 * (re_perp_long - m_l2 / q2 * (re_para_time));
            const double c = -2.0 * (a_para2 - a_perp2) * (1.0 - m_l2 / q2);

            return nf / (2.0 * M_PI) * (a + b * c_chi + c * c_chi2);
        }

        double pdf_chi(const double & chi) const
        {
            const double q2_min = power_of<2>(m_l());
            const double q2_max = 10.68;

            std::function<double (const double &)> num   = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_q2chi, this, std::placeholders::_1, cos(chi));
            std::function<double (const double &)> denom = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_q2,    this, std::placeholders::_1);

            return integrate<GSL::QAGS>(num, q2_min, q2_max) / integrate<GSL::QAGS>(denom, q2_min, q2_max);
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
    BToDPiLeptonNeutrino::differential_pdf_chi(const double & chi) const
    {
        return _imp->pdf_chi(chi);
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
    BToDPiLeptonNeutrino::integrated_pdf_chi(const double & chi_min, const double & chi_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToDPiLeptonNeutrino>::pdf_chi, _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, chi_min, chi_max);
    }

    const std::string
    BToDPiLeptonNeutrino::description = "\
The decay B->D pi l nu, where l is a massless lepton.";

    const std::string
    BToDPiLeptonNeutrino::kinematics_description_c_d = "\
The cosine of the helicity angle theta_D in the D-pi rest frame.";

    const std::string
    BToDPiLeptonNeutrino::kinematics_description_c_l = "\
The cosine of the helicity angle theta_L in the l-nu rest frame.";

    const std::string
    BToDPiLeptonNeutrino::kinematics_description_chi = "\
The azimuthal angle between the decay planes.";

}
