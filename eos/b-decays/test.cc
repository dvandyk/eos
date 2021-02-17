/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <eos/b-decays/test.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    using std::norm;

    /*
     * Decay: B_q -> l nubar, cf. [FMvD2015]
     */
    template <>
    struct Implementation<TestCacheableObservableProvider>
    {
        std::shared_ptr<Model> model;

        UsedParameter hbar;

        UsedParameter g_fermi;

        UsedParameter m_B;

        UsedParameter f_B;

        UsedParameter tau_B;

        SwitchOption opt_l;

        UsedParameter m_l;

        using IntermediateResult = TestCacheableObservableProvider::IntermediateResult;

        IntermediateResult _intermediate_result;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            hbar(p["hbar"], u),
            g_fermi(p["G_Fermi"], u),
            m_B(p["mass::B_u"], u),
            f_B(p["decay-constant::B_u"], u),
            tau_B(p["life_time::B_u"], u),
            opt_l(o, "l", {"e", "mu", "tau"}, "tau"),
            m_l(p["mass::" + opt_l.value()], u)
        {
            u.uses(*model);
        }

        const IntermediateResult * prepare(const double & q2)
        {
            _intermediate_result.a = 2.0;
            _intermediate_result.b = 13.0;

            _intermediate_result.q2 = q2;

            return &_intermediate_result;
        }

        double evaluate(const IntermediateResult * intermediate_result)
        {
            return intermediate_result->b - intermediate_result->a * intermediate_result->q2;
        }
    };

    TestCacheableObservableProvider::TestCacheableObservableProvider(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<TestCacheableObservableProvider>(new Implementation<TestCacheableObservableProvider>(parameters, options, *this))
    {
    }

    TestCacheableObservableProvider::~TestCacheableObservableProvider()
    {
    }

    const TestCacheableObservableProvider::IntermediateResult *
    TestCacheableObservableProvider::prepare(const double & q2) const
    {
        return _imp->prepare(q2);
    }

    double
    TestCacheableObservableProvider::evaluate(const TestCacheableObservableProvider::IntermediateResult * ir) const
    {
        return _imp->evaluate(ir);
    }

    const std::set<ReferenceName>
    TestCacheableObservableProvider::references
    {
    };
}
