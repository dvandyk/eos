/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017, 2018 Danny van Dyk
 * Copyright (c) 2018 Ahmet Kokulu
 * Copyright (c) 2018 Christoph Bobeth
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

#ifndef EOS_GUARD_EOS_UTILS_OPTIONS_IMPL_HH
#define EOS_GUARD_EOS_UTILS_OPTIONS_IMPL_HH 1

#include <eos/utils/join.hh>
#include <eos/utils/options.hh>
#include <eos/utils/qualified-name.hh>

#include <algorithm>
#include <map>

namespace eos
{
    class NameOption
    {
        public:
            std::unique_ptr<qnp::Name> _value;

        public:
            NameOption(const Options & options, const std::string & key) :
                _value(nullptr)
            {
                if (! options.has(key))
                    throw UnspecifiedOptionError(key);

                std::string raw_value(options[key]);
                try
                {
                    _value = std::make_unique<qnp::Name>(raw_value);
                }
                catch (QualifiedNameSyntaxError & e)
                {
                    throw InvalidOptionValueError(key, raw_value);
                }
            }

            NameOption(const Options & options, const std::string & key, const qnp::Name & default_value) :
                _value(std::make_unique<qnp::Name>(default_value))
            {
                if (options.has(key))
                {
                    std::string raw_value(options[key]);
                    try
                    {
                        _value = std::make_unique<qnp::Name>(raw_value);
                    }
                    catch (QualifiedNameSyntaxError & e)
                    {
                        throw InvalidOptionValueError(key, raw_value);
                    }
                }
            }

            const qnp::Name & value() const { return *_value; };
    };

    class SwitchOption
    {
        public:
            std::string _value;

        public:
            SwitchOption(const Options & options, const std::string & key,
                    const std::initializer_list<std::string> & allowed_values)
            {
                if (allowed_values.begin() == allowed_values.end())
                {
                    throw InternalError("SwitchOption: The list of allowed_values is empty");
                }

                if (! options.has(key))
                {
                    throw UnspecifiedOptionError(key, join(allowed_values.begin(), allowed_values.end()));
                }

                _value = options[key];

                if (std::find(allowed_values.begin(), allowed_values.end(), _value) == allowed_values.end())
                {
                    throw InvalidOptionValueError(key, _value, join(allowed_values.begin(), allowed_values.end()));
                }
            }

            SwitchOption(const Options & options, const std::string & key,
                    const std::initializer_list<std::string> & allowed_values,
                    const std::string & default_value) :
                _value(options.get(key, default_value))
            {
                if (allowed_values.begin() == allowed_values.end())
                {
                    throw InternalError("SwitchOption: The list of allowed_values is empty");
                }

                if (std::find(allowed_values.begin(), allowed_values.end(), default_value) == allowed_values.end())
                {
                    throw InternalError("SwitchOption: The default value '" + default_value + "'is not in the list of allowed values: '" + join(allowed_values.begin(), allowed_values.end()) + "'");
                }

                if (std::find(allowed_values.begin(), allowed_values.end(), _value) == allowed_values.end())
                {
                    throw InvalidOptionValueError(key, _value, join(allowed_values.begin(), allowed_values.end()));
                }
            }

            ~SwitchOption() = default;

            const std::string & value() const { return _value; };
    };

    template <typename T_>
    class MapOption
    {
        public:
            T_ _value;

        public:
            MapOption(const Options & options, const std::string & key,
                    const std::initializer_list<typename std::map<std::string, T_>::value_type> & _map)
            {
                std::map<std::string, T_> map{_map};

                if (map.begin() == map.end())
                {
                    throw InternalError("MappedOption: The map is empty");
                }

                if (! options.has(key))
                {
                    throw UnspecifiedOptionError(key);
                }

                std::string option_value = options[key];

                const auto entry = map.find(option_value);
                if (entry ==  map.end())
                {
                    std::vector<std::string> allowed_values;
                    for (const auto & e : map)
                    {
                        allowed_values.push_back(e.first);
                    }

                    throw InvalidOptionValueError(key, option_value, join(allowed_values.begin(), allowed_values.end()));
                }

                _value = entry->second;
            }

            ~MapOption() = default;

            const T_ & value() const { return _value; };
    };
}

#endif
