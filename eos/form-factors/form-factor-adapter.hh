/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013, 2016, 2017 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_FORM_FACTOR_ADAPTER_HH
#define EOS_GUARD_EOS_FORM_FACTORS_FORM_FACTOR_ADAPTER_HH 1

#include <eos/observable.hh>
#include <eos/form-factors/form-factors.hh>

namespace eos
{
    /* Form factor adapter class for interfacing Observable */
    template <typename Transition_>
    class FormFactorAdapter :
        public Observable
    {
        private:
            QualifiedName _name;

            qnp::Prefix _process;

            std::string _kinematic_variable;

            Parameters _parameters;

            Kinematics _kinematics;

            KinematicVariable _s;

            Options _options;

            std::shared_ptr<FormFactors<Transition_>> _form_factors;

            std::function<double (const FormFactors<Transition_> *, const double &)> _form_factor_function;

        public:
            FormFactorAdapter(const QualifiedName & name,
                    const qnp::Prefix & process,
                    const std::string & kinematic_variable,
                    const Parameters & parameters,
                    const Kinematics & kinematics,
                    const Options & options,
                    const std::function<double (const FormFactors<Transition_> *, const double &)> & form_factor_function) :
                _name(name),
                _process(process),
                _kinematic_variable(kinematic_variable),
                _parameters(parameters),
                _kinematics(kinematics),
                _s(kinematics[_kinematic_variable]),
                _options(options),
                _form_factor_function(form_factor_function)
            {
                if (! _options.has("form-factors"))
                    throw UnknownOptionError("form-factors");

                _form_factors = FormFactorFactory<Transition_>::create(process.str() + '@' + options["form-factors"], _parameters);
                if (! _form_factors)
                    throw NoSuchFormFactorError(process.str(), options["form-factors"]);
            }

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual double evaluate() const
            {
                return _form_factor_function(_form_factors.get(), _s());
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
                return ObservablePtr(new FormFactorAdapter(_name, _process, _kinematic_variable, _parameters.clone(), _kinematics.clone(), _options, _form_factor_function));
            }

            virtual ObservablePtr clone(const Parameters & parameters) const
            {
                return ObservablePtr(new FormFactorAdapter(_name, _process, _kinematic_variable, parameters, _kinematics.clone(), _options, _form_factor_function));
            }
    };

    template <typename Transition_>
    class FormFactorAdapterEntry :
        public ObservableEntry
    {
        private:
            QualifiedName _name;

            qnp::Prefix _process;

            std::string _kinematic_variable;

            std::function<double (const FormFactors<Transition_> *, const double &)> _form_factor_function;

        public:
            FormFactorAdapterEntry(const QualifiedName & name,
                    const qnp::Prefix & process,
                    const std::string & kinematic_variable,
                    const std::function<double (const FormFactors<Transition_> *, const double &)> & form_factor_function) :
                _name(name),
                _process(process),
                _kinematic_variable(kinematic_variable),
                _form_factor_function(form_factor_function)
            {
            }

            ~FormFactorAdapterEntry()
            {
            }

            virtual ObservablePtr make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
            {
                return ObservablePtr(new FormFactorAdapter<Transition_>(_name, _process, _kinematic_variable, parameters, kinematics, options, _form_factor_function));
            }

            virtual std::ostream & insert(std::ostream & os) const
            {
                os << "    type: adapter for one form factor" << std::endl;

                return os;
            }
    };

    /* Form factor ratio adapter class for interfacing Observable */
    template <typename Transition_>
    class FormFactorRatioAdapter :
        public Observable
    {
        private:
            QualifiedName _name;

            qnp::Prefix _process;

            Parameters _parameters;

            Kinematics _kinematics;

            KinematicVariable _s;

            Options _options;

            std::shared_ptr<FormFactors<Transition_>> _form_factors;

            std::function<double (const FormFactors<Transition_> *, const double &)> _form_factor_numerator;

            std::function<double (const FormFactors<Transition_> *, const double &)> _form_factor_denominator;

        public:
            FormFactorRatioAdapter(const QualifiedName & name,
                    const qnp::Prefix & process,
                    const Parameters & parameters,
                    const Kinematics & kinematics,
                    const Options & options,
                    const std::function<double (const FormFactors<Transition_> *, const double &)> & form_factor_numerator,
                    const std::function<double (const FormFactors<Transition_> *, const double &)> & form_factor_denominator) :
                _name(name),
                _process(process),
                _parameters(parameters),
                _kinematics(kinematics),
                _s(kinematics["s"]),
                _options(options),
                _form_factor_numerator(form_factor_numerator),
                _form_factor_denominator(form_factor_denominator)
            {
                if (! _options.has("form-factors"))
                    throw UnknownOptionError("form-factors");

                _form_factors = FormFactorFactory<Transition_>::create(process.str() + '@' + options["form-factors"], _parameters);
            }

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual double evaluate() const
            {
                double s = _s();

                return _form_factor_numerator(_form_factors.get(), s) / _form_factor_denominator(_form_factors.get(), s);
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
                return ObservablePtr(new FormFactorRatioAdapter(_name, _process, _parameters.clone(), _kinematics.clone(), _options, _form_factor_numerator, _form_factor_denominator));
            }

            virtual ObservablePtr clone(const Parameters & parameters) const
            {
                return ObservablePtr(new FormFactorRatioAdapter(_name, _process, parameters, _kinematics.clone(), _options, _form_factor_numerator, _form_factor_denominator));
            }
    };

    template <typename Transition_>
    class FormFactorRatioAdapterEntry :
        public ObservableEntry
    {
        private:
            QualifiedName _name;

            qnp::Prefix _process;

            std::function<double (const FormFactors<Transition_> *, const double &)> _form_factor_numerator;

            std::function<double (const FormFactors<Transition_> *, const double &)> _form_factor_denominator;

        public:
            FormFactorRatioAdapterEntry(const QualifiedName & name,
                    const qnp::Prefix & process,
                    const std::function<double (const FormFactors<Transition_> *, const double &)> & form_factor_numerator,
                    const std::function<double (const FormFactors<Transition_> *, const double &)> & form_factor_denominator) :
                _name(name),
                _process(process),
                _form_factor_numerator(form_factor_numerator),
                _form_factor_denominator(form_factor_denominator)
            {
            }

            ~FormFactorRatioAdapterEntry()
            {
            }

            virtual ObservablePtr make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
            {
                return ObservablePtr(new FormFactorRatioAdapter<Transition_>(_name, _process, parameters, kinematics, options, _form_factor_numerator, _form_factor_denominator));
            }

            virtual std::ostream & insert(std::ostream & os) const
            {
                os << "    type: adapter for a ratio of two form factors" << std::endl;

                return os;
            }
    };
}

#endif
