/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef BSTOKSTARLNU_GUARD_EOS_B_DECAYS_TEST_HH
#define BSTOKSTARLNU_GUARD_EOS_B_DECAYS_TEST_HH 1

#include <eos/observable.hh>
#include <eos/decays.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    class TestCacheableObservableProvider :
        public ParameterUser,
        public PrivateImplementationPattern<TestCacheableObservableProvider>
    {
        public:
            struct IntermediateResult :
                public CacheableObservable::IntermediateResult
            {
                // e.g. amplitudes
                double a;
                double b;

                // e.g. 
                double q2;
            };

            TestCacheableObservableProvider(const Parameters & parameters, const Options & options);
            ~TestCacheableObservableProvider();

            // Observables
            const IntermediateResult * prepare(const double & q2) const;

            double evaluate(const IntermediateResult *) const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;
    };
}

#endif

