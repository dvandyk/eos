/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Danny van Dyk
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef EOS_GUARD_SRC_UTILS_CONCRETE_COMPOSITE_OBSERVABLE_HH
#define EOS_GUARD_SRC_UTILS_CONCRETE_COMPOSITE_OBSERVABLE_HH 1

#include <eos/observable-impl.hh>
#include <eos/utils/apply.hh>
#include <eos/utils/join.hh>
#include <eos/utils/tuple-maker.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <array>
#include <functional>
#include <string>

namespace eos
{
    template <typename Decay_, typename ... Args_>
    class ConcreteCacheableObservable;

    template <typename Decay_, typename ... Args_>
    class ConcreteCachedObservable :
        public Observable
    {
        private:
            QualifiedName _name;

            Parameters _parameters;

            Kinematics _kinematics;

            Options _options;

            std::shared_ptr<Decay_> _decay;

            const typename Decay_::IntermediateResult * _intermediate_result;

            std::function<const typename Decay_::IntermediateResult * (const Decay_ *, const Args_ & ...)> _prepare_fn;

            std::function<double (const Decay_ *, const typename Decay_::IntermediateResult *)> _evaluate_fn;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names;

        public:
            ConcreteCachedObservable(const QualifiedName & name,
                    const Parameters & parameters,
                    const Kinematics & kinematics,
                    const Options & options,
                    const std::shared_ptr<Decay_> & decay,
                    const typename Decay_::IntermediateResult * intermediate_result,
                    const std::function<const typename Decay_::IntermediateResult * (const Decay_ *, const Args_ & ...)> & prepare_fn,
                    const std::function<double (const Decay_ *, const typename Decay_::IntermediateResult *)> & evaluate_fn,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names) :
                _name(name),
                _parameters(parameters),
                _kinematics(kinematics),
                _options(options),
                _decay(decay),
                _intermediate_result(intermediate_result),
                _prepare_fn(prepare_fn),
                _evaluate_fn(evaluate_fn),
                _kinematics_names(kinematics_names)
            {
                uses(*_decay);
            }

            ~ConcreteCachedObservable() = default;

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual double evaluate() const
            {
                return _evaluate_fn(_decay.get(), _intermediate_result);
            };

            virtual Parameters parameters()
            {
                return _parameters;
            };

            virtual Kinematics kinematics()
            {
                return _kinematics;
            };

            virtual Options options()
            {
                return _options;
            }

            virtual ObservablePtr clone() const
            {
                return ObservablePtr(new ConcreteCacheableObservable<Decay_, Args_ ...>(_name, _parameters.clone(), _kinematics.clone(), _options, _prepare_fn, _evaluate_fn, _kinematics_names));
            }

            virtual ObservablePtr clone(const Parameters & parameters) const
            {
                return ObservablePtr(new ConcreteCacheableObservable<Decay_, Args_ ...>(_name, parameters, _kinematics.clone(), _options, _prepare_fn, _evaluate_fn, _kinematics_names));
            }
    };

    template <typename Decay_, typename ... Args_>
    class ConcreteCacheableObservable :
        public CacheableObservable
    {
        private:
            QualifiedName _name;

            Parameters _parameters;

            Kinematics _kinematics;

            Options _options;

            std::shared_ptr<Decay_> _decay;

            std::function<const typename Decay_::IntermediateResult * (const Decay_ *, const Args_ & ...)> _prepare_fn;

            std::function<double (const Decay_ *, const typename Decay_::IntermediateResult *)> _evaluate_fn;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names;

            std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, KinematicVariable>::Type ...> _argument_tuple;

        public:
            ConcreteCacheableObservable(const QualifiedName & name,
                    const Parameters & parameters,
                    const Kinematics & kinematics,
                    const Options & options,
                    const std::function<const typename Decay_::IntermediateResult * (const Decay_ *, const Args_ & ...)> & prepare_fn,
                    const std::function<double (const Decay_ *, const typename Decay_::IntermediateResult *)> & evaluate_fn,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names) :
                _name(name),
                _parameters(parameters),
                _kinematics(kinematics),
                _options(options),
                _decay(new Decay_(parameters, options)),
                _prepare_fn(prepare_fn),
                _evaluate_fn(evaluate_fn),
                _kinematics_names(kinematics_names),
                _argument_tuple(impl::TupleMaker<sizeof...(Args_)>::make(_kinematics, _kinematics_names, _decay.get()))
            {
                uses(*_decay);
            }

            ~ConcreteCacheableObservable() = default;

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual double evaluate() const
            {
                std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, double>::Type ...> values = _argument_tuple;

                const typename Decay_::IntermediateResult * intermediate_result = apply(_prepare_fn, values);

                return _evaluate_fn(_decay.get(), intermediate_result);
            };

            virtual const CacheableObservable::IntermediateResult * prepare() const
            {
                std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, double>::Type ...> values = _argument_tuple;

                return apply(_prepare_fn, values);
            }

            virtual double evaluate(const CacheableObservable::IntermediateResult * intermediate_result) const
            {
                return _evaluate_fn(_decay.get(), static_cast<const typename Decay_::IntermediateResult *>(intermediate_result));
            }

            virtual Parameters parameters()
            {
                return _parameters;
            };

            virtual Kinematics kinematics()
            {
                return _kinematics;
            };

            virtual Options options()
            {
                return _options;
            }

            virtual ObservablePtr make_cached_observable(const CacheableObservable * _other) const
            {
                auto other = dynamic_cast<decltype(this)>(_other);
                if (nullptr == other)
                    return { nullptr };

                if (other->_parameters != this->_parameters)
                    return { nullptr };

                if (other->_kinematics != this->_kinematics)
                    return { nullptr };

                if (other->_options != this->_options)
                    return { nullptr };

                std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, double>::Type ...> values = _argument_tuple;

                return ObservablePtr(new ConcreteCachedObservable<Decay_, Args_ ...>(_name, _parameters, _kinematics, _options, other->_decay, apply(other->_prepare_fn, values), _prepare_fn, _evaluate_fn, _kinematics_names));
            }

            virtual ObservablePtr clone() const
            {
                return ObservablePtr(new ConcreteCacheableObservable(_name, _parameters.clone(), _kinematics.clone(), _options, _prepare_fn, _evaluate_fn, _kinematics_names));
            }

            virtual ObservablePtr clone(const Parameters & parameters) const
            {
                return ObservablePtr(new ConcreteCacheableObservable(_name, parameters, _kinematics.clone(), _options, _prepare_fn, _evaluate_fn, _kinematics_names));
            }
    };

    template <typename Decay_, typename ... Args_>
    class ConcreteCacheableObservableEntry :
        public ObservableEntry
    {
        private:
            QualifiedName _name;

            std::string _latex;

            std::function<const typename Decay_::IntermediateResult * (const Decay_ *, const Args_ & ...)> _prepare_fn;

            std::function<double (const Decay_ *, const typename Decay_::IntermediateResult *)> _evaluate_fn;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names;

            std::array<const std::string, sizeof...(Args_)> _kinematics_names_array;

            Options _forced_options;

        public:
            ConcreteCacheableObservableEntry(const QualifiedName & name, const std::string & latex,
                    const std::function<const typename Decay_::IntermediateResult * (const Decay_ *, const Args_ & ...)> & prepare_fn,
                    const std::function<double (const Decay_ *, const typename Decay_::IntermediateResult *)> & evaluate_fn,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names,
                    const Options & forced_options) :
                _name(name),
                _latex(latex),
                _prepare_fn(prepare_fn),
                _evaluate_fn(evaluate_fn),
                _kinematics_names(kinematics_names),
                _kinematics_names_array(impl::make_array<const std::string>(kinematics_names)),
                _forced_options(forced_options)
            {
            }

            ~ConcreteCacheableObservableEntry()
            {
            }

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual const std::string & latex() const
            {
                return _latex;
            }

            virtual ObservableEntry::KinematicVariableIterator begin_kinematic_variables() const
            {
                return _kinematics_names_array.begin();
            }

            virtual ObservableEntry::KinematicVariableIterator end_kinematic_variables() const
            {
                return _kinematics_names_array.end();
            }

            virtual ObservablePtr make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
            {
                return ObservablePtr(new ConcreteCacheableObservable<Decay_, Args_ ...>(_name, parameters, kinematics, options + _forced_options, _prepare_fn, _evaluate_fn, _kinematics_names));
            }

            virtual std::ostream & insert(std::ostream & os) const
            {
                os << "    type: composite observable" << std::endl;

                if (sizeof...(Args_) > 0)
                {
                    os << "    kinematic variables: " << join(std::begin(_kinematics_names_array), std::end(_kinematics_names_array)) << std::endl;
                }

                return os;
            }
    };

    template <typename Decay_, typename Tuple_, typename ... Args_>
    ObservableEntryPtr make_concrete_composite_observable_entry(const QualifiedName & name, const std::string & latex,
            const typename Decay_::IntermediateResult * (Decay_::* prepare_fn)(const Args_ & ...) const,
            double (Decay_::* evaluate_fn)(const typename Decay_::IntermediateResult *) const,
            const Tuple_ & kinematics_names,
            const Options & forced_options)
    {
        static_assert(sizeof...(Args_) == impl::TupleSize<Tuple_>::size, "Need as many function arguments as kinematics names!");

        return std::make_shared<ConcreteCacheableObservableEntry<Decay_, Args_ ...>>(name, latex,
                std::function<const typename Decay_::IntermediateResult * (const Decay_ *, const Args_ & ...)>(std::mem_fn(prepare_fn)),
                std::function<double (const Decay_ *, const typename Decay_::IntermediateResult *)>(std::mem_fn(evaluate_fn)),
                kinematics_names, forced_options);
    }
}


#endif
