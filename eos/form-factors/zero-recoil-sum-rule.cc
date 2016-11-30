/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <eos/form-factors/zero-recoil-sum-rule.hh>
#include <eos/utils/model.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/power_of.hh>

namespace eos
{
    template <>
    struct Implementation<ZeroRecoilSumRule<LambdaBToLambdaC>>
    {
        std::shared_ptr<Model> model;

        // renormalization scale (kinetic scheme)
        UsedParameter mu;

        // excitation energy cut off
        UsedParameter wM;

        // matrix elements of dim 5
        UsedParameter mu2_pi;

        // matrix elements of dim 6
        UsedParameter rho3_D;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            mu(p["Lambda_b->Lambda_c::mu@ZRSR"], u),
            wM(p["Lambda_b->Lambda_c::wM@ZRSR"], u),
            mu2_pi(p["Lambda_b->Lambda_b::mu_pi^2@1GeV"], u),
            rho3_D(p["Lambda_b->Lambda_b::rho_D^3@1GeV"], u)
        {
            u.uses(*model);
        }

        inline double xiA() const
        {
            const double mb = model->m_b_kin(mu), mb2 = mb * mb;
            const double mc = model->m_c_kin(mu), mc2 = mc * mc;

            const double a_s = model->alpha_s(sqrt(mb * mc));
            const double mu2 = mu() * mu();
            const double mu_p = mu() - mc + std::sqrt(mc2 + mu2), mu_p2 = mu_p * mu_p;

            const double etaApert = 1.0 + a_s / M_PI * ((mb + mc) / (mb - mc) * log(mb / mc) - 8.0 / 3.0);
            const double etaAsoft = - a_s * mu2 / (3.0 * M_PI) * (1.0 / mc2 + 2.0 / (3.0 * mb * mc) + 1.0 / mb2);
            const double etaAspec = 4.0 * a_s / (3.0 * M_PI) * (
                        (wM - mu_p) * (wM + 2.0 * mc + mu_p) / (24.0 * mb2 * power_of<2>(wM + mc)) *
                        (2.0 * wM * (wM + 2.0 * mc) + mc2 * (mc2 - 3.0 * mb2 - 2.0 * mb * mc + 4.0 * mc * mu_p + 2.0 * mu_p2) / power_of<2>(mc + mu_p))
                        - (3.0 * mb - mc) * (mb + mc) / (12.0 * mb2) * log((mc + mu_p) / (mc + wM))
                    );

            // cf. [GMU2012], eq. (10), p. 7. Note that our etaSpec
            // correspond to eq. (15), and therefore enters with a
            // relative sign compared to the original eq. (10).
            return power_of<2>(etaApert) - 2.0 * etaAsoft + etaAspec;
        }

        inline double deltaA() const
        {
            const double mb = model->m_b_kin(mu), mb2 = mb * mb;
            const double mc = model->m_c_kin(mu), mc2 = mc * mc, mc3 = mc2 * mc;

            const double deltaA2 = mu2_pi / 4.0 * (1.0 / mc2 + 1.0 / mb2 + 2.0 / (3.0 * mb * mc));
            const double deltaA3 = rho3_D / (4.0 * mc3) + rho3_D / (12.0 * mb) * (1.0 / mc2 + 1.0 / (mc * mb) + 3.0 / mb2);

            return deltaA2 + deltaA3;
        }

        double xiV() const
        {
            // cf. [U2003], eq. (27)

            const double mb = model->m_b_kin(mu), mb2 = mb * mb;
            const double mc = model->m_c_kin(mu), mc2 = mc * mc;

            const double a_s = model->alpha_s(sqrt(mb * mc));
            const double mu = this->mu(), mu2 = mu * mu;
            const double wb = std::sqrt(mu2 + mb2), wb3 = wb * wb * wb;
            const double wc = std::sqrt(mu2 + mc2), wc3 = wc * wc * wc;

            const double xiVnlo1 =
                (3.0 * mb2 + 2.0 * mc * mb + 3.0 * mc2) / (2.0 * (mb2 - mc2)) * log((mu + wb) / (mu + wc))
                - 2.0;
            const double xiVnlo2 =
                4.0 / (3.0 * mu2) * (mc * wb - mb * wc) / (mb - mc)
                + 2.0 / 3.0 * (mc / wb - mb / wc) / (mb - mc)
                - 1.0 / 3.0 * (wb / mb - wc / mc) / (mb - mc)
                + 2.0 * mc * mb / (mb2 - mc2) * (1.0 / wb - 1.0 / wc)
                + 1.0 / (6.0 * (mc + mb)) * (wc / mc * (3.0 - mb / mc) + wb / mb * (3.0 - mc / mb))
                + 4.0 * mc * mb / (3.0 * (mc + mb)) * (mb / wb3 + mc / wc3)
                + mu / 6.0 * power_of<2>(1.0 / mc - 1.0 / mb)
                - 2.0 * mu2 / (3.0 * mc * mb) * (mb2 / wc + mc2 / wb) / (mb2 - mc2);

            return 1.0 + 2.0 * a_s / (3.0 * M_PI) * (xiVnlo1 - mu * xiVnlo2);
        }

        inline double deltaV() const
        {
            const double mb = model->m_b_kin(mu);
            const double mc = model->m_c_kin(mu);

            const double deltaV2 = mu2_pi / 4.0 * power_of<2>(1.0 / mc - 1.0 / mb);
            const double deltaV3 = rho3_D / 4.0 * power_of<2>(1.0 / mc - 1.0 / mb) * (1.0 / mc + 1.0 / mb);

            return deltaV2 + deltaV3;
        }
    };

    ZeroRecoilSumRule<LambdaBToLambdaC>::ZeroRecoilSumRule(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<ZeroRecoilSumRule<LambdaBToLambdaC>>(new Implementation<ZeroRecoilSumRule<LambdaBToLambdaC>>(parameters, options, *this))
    {
    }

    ZeroRecoilSumRule<LambdaBToLambdaC>::~ZeroRecoilSumRule()
    {
    }

    double
    ZeroRecoilSumRule<LambdaBToLambdaC>::axialvector_current() const
    {
        return _imp->xiA() - _imp->deltaA();
    }

    double
    ZeroRecoilSumRule<LambdaBToLambdaC>::vector_current() const
    {
        return _imp->xiV() - _imp->deltaV();
    }

    template <>
    struct Implementation<ZeroRecoilSumRule<BToD>>
    {
        std::shared_ptr<Model> model;

        // renormalization scale (kinetic scheme)
        UsedParameter mu;

        // excitation energy cut off
        UsedParameter wM;

        // matrix elements of dim 5
        UsedParameter mu2_pi;
        UsedParameter mu2_G;

        // matrix elements of dim 6
        UsedParameter rho3_D;
        UsedParameter rho3_LS;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            mu(p["B->D::mu@ZRSR"], u),
            wM(p["B->D::wM@ZRSR"], u),
            mu2_pi(p["B->B::mu_pi^2@1GeV"], u),
            mu2_G(p["B->B::mu_G^2@1GeV"], u),
            rho3_D(p["B->B::rho_D^3@1GeV"], u),
            rho3_LS(p["B->B::rho_LS^3@1GeV"], u)
        {
            u.uses(*model);
        }

        double xiV() const
        {
            // cf. [U2003], eq. (27)

            const double mb = model->m_b_kin(mu), mb2 = mb * mb;
            const double mc = model->m_c_kin(mu), mc2 = mc * mc;

            const double a_s = model->alpha_s(sqrt(mb * mc));
            const double mu = this->mu(), mu2 = mu * mu;
            const double wb = std::sqrt(mu2 + mb2), wb3 = wb * wb * wb;
            const double wc = std::sqrt(mu2 + mc2), wc3 = wc * wc * wc;

            const double xiVnlo1 =
                (3.0 * mb2 + 2.0 * mc * mb + 3.0 * mc2) / (2.0 * (mb2 - mc2)) * log((mu + wb) / (mu + wc))
                - 2.0;
            const double xiVnlo2 =
                4.0 / (3.0 * mu2) * (mc * wb - mb * wc) / (mb - mc)
                + 2.0 / 3.0 * (mc / wb - mb / wc) / (mb - mc)
                - 1.0 / 3.0 * (wb / mb - wc / mc) / (mb - mc)
                + 2.0 * mc * mb / (mb2 - mc2) * (1.0 / wb - 1.0 / wc)
                + 1.0 / (6.0 * (mc + mb)) * (wc / mc * (3.0 - mb / mc) + wb / mb * (3.0 - mc / mb))
                + 4.0 * mc * mb / (3.0 * (mc + mb)) * (mb / wb3 + mc / wc3)
                + mu / 6.0 * power_of<2>(1.0 / mc - 1.0 / mb)
                - 2.0 * mu2 / (3.0 * mc * mb) * (mb2 / wc + mc2 / wb) / (mb2 - mc2);

            return 1.0 + 2.0 * a_s / (3.0 * M_PI) * (xiVnlo1 - mu * xiVnlo2);
        }

        inline double deltaV() const
        {
            const double mb = model->m_b_kin(mu);
            const double mc = model->m_c_kin(mu);

            const double deltaV2 = (mu2_pi - mu2_G) / 4.0 * power_of<2>(1.0 / mc - 1.0 / mb);
            const double deltaV3 = (rho3_D + rho3_LS) / 4.0 * power_of<2>(1.0 / mc - 1.0 / mb) * (1.0 / mc + 1.0 / mb);

            return deltaV2 + deltaV3;
        }
    };

    ZeroRecoilSumRule<BToD>::ZeroRecoilSumRule(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<ZeroRecoilSumRule<BToD>>(new Implementation<ZeroRecoilSumRule<BToD>>(parameters, options, *this))
    {
    }

    ZeroRecoilSumRule<BToD>::~ZeroRecoilSumRule()
    {
    }

    double
    ZeroRecoilSumRule<BToD>::vector_current() const
    {
        return _imp->xiV() - _imp->deltaV();
    }

    template <>
    struct Implementation<ZeroRecoilSumRule<BToDstar>>
    {
        std::shared_ptr<Model> model;

        // renormalization scale (kinetic scheme)
        UsedParameter mu;

        // excitation energy cut off
        UsedParameter wM;

        // matrix elements of dim 5
        UsedParameter mu2_pi;
        UsedParameter mu2_G;

        // matrix elements of dim 6
        UsedParameter rho3_D;
        UsedParameter rho3_LS;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            mu(p["B->D::mu@ZRSR"], u),
            wM(p["B->D::wM@ZRSR"], u),
            mu2_pi(p["B->B::mu_pi^2@1GeV"], u),
            mu2_G(p["B->B::mu_G^2@1GeV"], u),
            rho3_D(p["B->B::rho_D^3@1GeV"], u),
            rho3_LS(p["B->B::rho_LS^3@1GeV"], u)
        {
            u.uses(*model);
        }

        inline double xiA() const
        {
            const double mb = model->m_b_kin(mu), mb2 = mb * mb;
            const double mc = model->m_c_kin(mu), mc2 = mc * mc;

            const double a_s = model->alpha_s(sqrt(mb * mc));
            const double mu2 = mu() * mu();
            const double mu_p = mu() - mc + std::sqrt(mc2 + mu2), mu_p2 = mu_p * mu_p;

            const double etaApert = 1.0 + a_s / M_PI * ((mb + mc) / (mb - mc) * log(mb / mc) - 8.0 / 3.0);
            const double etaAsoft = - a_s * mu2 / (3.0 * M_PI) * (1.0 / mc2 + 2.0 / (3.0 * mb * mc) + 1.0 / mb2);
            const double etaAspec = 4.0 * a_s / (3.0 * M_PI) * (
                        (wM - mu_p) * (wM + 2.0 * mc + mu_p) / (24.0 * mb2 * power_of<2>(wM + mc)) *
                        (2.0 * wM * (wM + 2.0 * mc) + mc2 * (mc2 - 3.0 * mb2 - 2.0 * mb * mc + 4.0 * mc * mu_p + 2.0 * mu_p2) / power_of<2>(mc + mu_p))
                        - (3.0 * mb - mc) * (mb + mc) / (12.0 * mb2) * log((mc + mu_p) / (mc + wM))
                    );

            // cf. [GMU2012], eq. (10), p. 7. Note that our etaSpec
            // correspond to eq. (15), and therefore enters with a
            // relative sign compared to the original eq. (10).
            return power_of<2>(etaApert) - 2 * etaAsoft + etaAspec;
        }

        inline double deltaA() const
        {
            const double mb = model->m_b_kin(mu), mb2 = mb * mb;
            const double mc = model->m_c_kin(mu), mc2 = mc * mc, mc3 = mc2 * mc;

            const double deltaA2 = mu2_G / (3.0 * mc2) + (mu2_pi - mu2_G) / 4.0 * (1.0 / mc2 + 1.0 / mb2 + 2.0 / (3.0 * mb * mc));
            const double deltaA3 = (3.0 * rho3_D - rho3_LS) / (12.0 * mc3)
                    + (rho3_D + rho3_LS) / (12.0 * mb) * (1.0 / mc2 + 1.0 / (mc * mb) + 3.0 / mb2);

            return deltaA2 + deltaA3;
        }
    };

    ZeroRecoilSumRule<BToDstar>::ZeroRecoilSumRule(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<ZeroRecoilSumRule<BToDstar>>(new Implementation<ZeroRecoilSumRule<BToDstar>>(parameters, options, *this))
    {
    }

    ZeroRecoilSumRule<BToDstar>::~ZeroRecoilSumRule()
    {
    }

    double
    ZeroRecoilSumRule<BToDstar>::axialvector_current() const
    {
        return _imp->xiA() - _imp->deltaA();
    }
}
