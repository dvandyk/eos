/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013, 2015, 2016 Danny van Dyk
 * Copyright (c) 2014 Frederik Beaujean
 * Copyright (c) 2014 Christoph Bobeth
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

#include <eos/utils/complex.hh>
#include <eos/utils/log.hh>
#include <eos/utils/matrix.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/wilson_scan_model.hh>

#include <cmath>

namespace eos
{
    using std::sqrt;

    namespace wcimplementation
    {
        complex<double> polar(const Parameter & abs, const Parameter & arg) { return std::polar(abs(), arg()); }
        complex<double> cartesian(const Parameter & re, const Parameter & im) { return complex<double>(re(), im()); }
        complex<double> polar_negative(const Parameter & abs, const Parameter & arg) { return std::polar(abs(), arg() + M_PI); }
        complex<double> cartesian_negative(const Parameter & re, const Parameter & im) { return complex<double>(-re(), -im()); }
        complex<double> zero() { return complex<double>(0.0, 0.0); }
    }

    /* b->d Wilson coefficients */
    WilsonScanComponent<components::DeltaBD1>::WilsonScanComponent(const Parameters & p, const Options &, ParameterUser & u) :
        _alpha_s_Z__deltabd1(p["QCD::alpha_s(MZ)"], u),
        _mu_b__deltabd1(p["QCD::mu_b"], u),
        _m_Z__deltabd1(p["mass::Z"], u),
        _mu__deltabd1(p["mu"], u),
        /* b->d */
        _bd_c1(p["b->d::c1"], u),
        _bd_c2(p["b->d::c2"], u),
        _bd_c3(p["b->d::c3"], u),
        _bd_c4(p["b->d::c4"], u),
        _bd_c5(p["b->d::c5"], u),
        _bd_c6(p["b->d::c6"], u),
        _bd_re_c7(p["b->d::Re{c7}"], u),
        _bd_im_c7(p["b->d::Im{c7}"], u),
        _bd_re_c7prime(p["b->d::Re{c7'}"], u),
        _bd_im_c7prime(p["b->d::Im{c7'}"], u),
        _bd_c8(p["b->d::c8"], u),
        _bd_c8prime(p["b->d::c8'"], u),
        /* b->dee */
        _bdee_re_c9(p["b->dee::Re{c9}"], u),
        _bdee_im_c9(p["b->dee::Im{c9}"], u),
        _bdee_re_c10(p["b->dee::Re{c10}"], u),
        _bdee_im_c10(p["b->dee::Im{c10}"], u),
        _bdee_re_c9prime(p["b->dee::Re{c9'}"], u),
        _bdee_im_c9prime(p["b->dee::Im{c9'}"], u),
        _bdee_re_c10prime(p["b->dee::Re{c10'}"], u),
        _bdee_im_c10prime(p["b->dee::Im{c10'}"], u),
        _bdee_re_cS(p["b->dee::Re{cS}"], u),
        _bdee_im_cS(p["b->dee::Im{cS}"], u),
        _bdee_re_cSprime(p["b->dee::Re{cS'}"], u),
        _bdee_im_cSprime(p["b->dee::Im{cS'}"], u),
        _bdee_re_cP(p["b->dee::Re{cP}"], u),
        _bdee_im_cP(p["b->dee::Im{cP}"], u),
        _bdee_re_cPprime(p["b->dee::Re{cP'}"], u),
        _bdee_im_cPprime(p["b->dee::Im{cP'}"], u),
        _bdee_re_cT(p["b->dee::Re{cT}"], u),
        _bdee_im_cT(p["b->dee::Im{cT}"], u),
        _bdee_re_cT5(p["b->dee::Re{cT5}"], u),
        _bdee_im_cT5(p["b->dee::Im{cT5}"], u),
        /* b->dmumu */
        _bdmumu_re_c9(p["b->dmumu::Re{c9}"], u),
        _bdmumu_im_c9(p["b->dmumu::Im{c9}"], u),
        _bdmumu_re_c10(p["b->dmumu::Re{c10}"], u),
        _bdmumu_im_c10(p["b->dmumu::Im{c10}"], u),
        _bdmumu_re_c9prime(p["b->dmumu::Re{c9'}"], u),
        _bdmumu_im_c9prime(p["b->dmumu::Im{c9'}"], u),
        _bdmumu_re_c10prime(p["b->dmumu::Re{c10'}"], u),
        _bdmumu_im_c10prime(p["b->dmumu::Im{c10'}"], u),
        _bdmumu_re_cS(p["b->dmumu::Re{cS}"], u),
        _bdmumu_im_cS(p["b->dmumu::Im{cS}"], u),
        _bdmumu_re_cSprime(p["b->dmumu::Re{cS'}"], u),
        _bdmumu_im_cSprime(p["b->dmumu::Im{cS'}"], u),
        _bdmumu_re_cP(p["b->dmumu::Re{cP}"], u),
        _bdmumu_im_cP(p["b->dmumu::Im{cP}"], u),
        _bdmumu_re_cPprime(p["b->dmumu::Re{cP'}"], u),
        _bdmumu_im_cPprime(p["b->dmumu::Im{cP'}"], u),
        _bdmumu_re_cT(p["b->dmumu::Re{cT}"], u),
        _bdmumu_im_cT(p["b->dmumu::Im{cT}"], u),
        _bdmumu_re_cT5(p["b->dmumu::Re{cT5}"], u),
        _bdmumu_im_cT5(p["b->dmumu::Im{cT5}"], u),


        /* functions for b->dgamma */
        _bd_c7(std::bind(&wcimplementation::cartesian,           _bd_re_c7,           _bd_im_c7)),
        _bd_c7prime(std::bind(&wcimplementation::cartesian,      _bd_re_c7prime,      _bd_im_c7prime)),

        /* functions for b->dee */
        _bdee_c9(std::bind(&wcimplementation::cartesian,         _bdee_re_c9,         _bdee_im_c9)),
        _bdee_c10(std::bind(&wcimplementation::cartesian,        _bdee_re_c10,        _bdee_im_c10)),
        _bdee_c9prime(std::bind(&wcimplementation::cartesian,    _bdee_re_c9prime,    _bdee_im_c9prime)),
        _bdee_c10prime(std::bind(&wcimplementation::cartesian,   _bdee_re_c10prime,   _bdee_im_c10prime)),
        _bdee_cS(std::bind(&wcimplementation::cartesian,         _bdee_re_cS,         _bdee_im_cS)),
        _bdee_cSprime(std::bind(&wcimplementation::cartesian,    _bdee_re_cSprime,    _bdee_im_cSprime)),
        _bdee_cP(std::bind(&wcimplementation::cartesian,         _bdee_re_cP,         _bdee_im_cP)),
        _bdee_cPprime(std::bind(&wcimplementation::cartesian,    _bdee_re_cPprime,    _bdee_im_cPprime)),
        _bdee_cT(std::bind(&wcimplementation::cartesian,         _bdee_re_cT,         _bdee_im_cT)),
        _bdee_cT5(std::bind(&wcimplementation::cartesian,        _bdee_re_cT5,        _bdee_im_cT5)),

        /* functions for b->dmumu */
        _bdmumu_c9(std::bind(&wcimplementation::cartesian,       _bdmumu_re_c9,       _bdmumu_im_c9)),
        _bdmumu_c10(std::bind(&wcimplementation::cartesian,      _bdmumu_re_c10,      _bdmumu_im_c10)),
        _bdmumu_c9prime(std::bind(&wcimplementation::cartesian,  _bdmumu_re_c9prime,  _bdmumu_im_c9prime)),
        _bdmumu_c10prime(std::bind(&wcimplementation::cartesian, _bdmumu_re_c10prime, _bdmumu_im_c10prime)),
        _bdmumu_cS(std::bind(&wcimplementation::cartesian,       _bdmumu_re_cS,       _bdmumu_im_cS)),
        _bdmumu_cSprime(std::bind(&wcimplementation::cartesian,  _bdmumu_re_cSprime,  _bdmumu_im_cSprime)),
        _bdmumu_cP(std::bind(&wcimplementation::cartesian,       _bdmumu_re_cP,       _bdmumu_im_cP)),
        _bdmumu_cPprime(std::bind(&wcimplementation::cartesian,  _bdmumu_re_cPprime,  _bdmumu_im_cPprime)),
        _bdmumu_cT(std::bind(&wcimplementation::cartesian,       _bdmumu_re_cT,       _bdmumu_im_cT)),
        _bdmumu_cT5(std::bind(&wcimplementation::cartesian,      _bdmumu_re_cT5,      _bdmumu_im_cT5))
    {
    }

    WilsonCoefficients<BToD>
    WilsonScanComponent<components::DeltaBD1>::wilson_coefficients_b_to_d(const std::string & lepton_flavour, const bool & cp_conjugate) const
    {
        std::function<complex<double> ()> bd_c9,  bd_c9prime;
        std::function<complex<double> ()> bd_c10, bd_c10prime;
        std::function<complex<double> ()> bd_cS,  bd_cSprime;
        std::function<complex<double> ()> bd_cP,  bd_cPprime;
        std::function<complex<double> ()> bd_cT,  bd_cT5;

        if ("e" == lepton_flavour)
        {
            bd_c9  = _bdee_c9;    bd_c9prime  = _bdee_c9prime;
            bd_c10 = _bdee_c10;   bd_c10prime = _bdee_c10prime;
            bd_cS  = _bdee_cS;    bd_cSprime  = _bdee_cSprime;
            bd_cP  = _bdee_cP;    bd_cPprime  = _bdee_cPprime;
            bd_cT  = _bdee_cT;    bd_cT5      = _bdee_cT5;
        }
        else if ("mu" == lepton_flavour)
        {
            bd_c9  = _bdmumu_c9;    bd_c9prime  = _bdmumu_c9prime;
            bd_c10 = _bdmumu_c10;   bd_c10prime = _bdmumu_c10prime;
            bd_cS  = _bdmumu_cS;    bd_cSprime  = _bdmumu_cSprime;
            bd_cP  = _bdmumu_cP;    bd_cPprime  = _bdmumu_cPprime;
            bd_cT  = _bdmumu_cT;    bd_cT5      = _bdmumu_cT5;
        }
        else
        {
            throw InternalError("WilsonScan presently only implements 'e' and 'mu' lepton flavours");
        }

        double alpha_s = 0.0;
        if (_mu__deltabd1 < _mu_b__deltabd1)
        {
            alpha_s = QCD::alpha_s(_mu_b__deltabd1, _alpha_s_Z__deltabd1, _m_Z__deltabd1, QCD::beta_function_nf_5);
            alpha_s = QCD::alpha_s(_mu__deltabd1, alpha_s, _mu_b__deltabd1, QCD::beta_function_nf_4);
        }
        else
        {
            alpha_s = QCD::alpha_s(_mu__deltabd1, _alpha_s_Z__deltabd1, _m_Z__deltabd1, QCD::beta_function_nf_5);
        }

        complex<double> a_s = alpha_s / 4.0 / M_PI;
        WilsonCoefficients<BToD> result;
        result._sm_like_coefficients = std::array<std::complex<double>, 15>
        {
            _bd_c1(), _bd_c2(), _bd_c3(), _bd_c4(), _bd_c5(), _bd_c6(),
            0.0, 0.0, 0.0, 0.0, 0.0,
            a_s * _bd_c7(), a_s * _bd_c8(), a_s * bd_c9(), a_s * bd_c10()
        };
        result._primed_coefficients = std::array<std::complex<double>, 15>
        {
            /* we only consider c7', c8', c9' and c10' */
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            a_s * _bd_c7prime(), a_s * _bd_c8prime(), a_s * bd_c9prime(), a_s * bd_c10prime()
        };
        result._scalar_tensor_coefficients = std::array<std::complex<double>, 6>
        {
            bd_cS(), bd_cSprime(), bd_cP(), bd_cPprime(), bd_cT(), bd_cT5()
        };
        result._alpha_s = alpha_s;

        if (cp_conjugate)
        {
            for (auto & c : result._sm_like_coefficients)
            {
                c = conj(c);
            }

            for (auto & c : result._primed_coefficients)
            {
                c = conj(c);
            }

            for (auto & c : result._scalar_tensor_coefficients)
            {
                c = conj(c);
            }
        }

        return result;
    }

    /* b->s Wilson coefficients */
    WilsonScanComponent<components::DeltaBS1>::WilsonScanComponent(const Parameters & p, const Options &, ParameterUser & u) :
        _alpha_s_Z__deltabs1(p["QCD::alpha_s(MZ)"], u),
        _mu_b__deltabs1(p["QCD::mu_b"], u),
        _m_Z__deltabs1(p["mass::Z"], u),
        _mu__deltabs1(p["mu"], u),
        /* b->s */
        _c1(p["b->s::c1"], u),
        _c2(p["b->s::c2"], u),
        _c3(p["b->s::c3"], u),
        _c4(p["b->s::c4"], u),
        _c5(p["b->s::c5"], u),
        _c6(p["b->s::c6"], u),
        _re_c7(p["b->s::Re{c7}"], u),
        _im_c7(p["b->s::Im{c7}"], u),
        _re_c7prime(p["b->s::Re{c7'}"], u),
        _im_c7prime(p["b->s::Im{c7'}"], u),
        _c8(p["b->s::c8"], u),
        _c8prime(p["b->s::c8'"], u),
        /* b->see */
        _e_re_c9(p["b->see::Re{c9}"], u),
        _e_im_c9(p["b->see::Im{c9}"], u),
        _e_re_c10(p["b->see::Re{c10}"], u),
        _e_im_c10(p["b->see::Im{c10}"], u),
        _e_re_c9prime(p["b->see::Re{c9'}"], u),
        _e_im_c9prime(p["b->see::Im{c9'}"], u),
        _e_re_c10prime(p["b->see::Re{c10'}"], u),
        _e_im_c10prime(p["b->see::Im{c10'}"], u),
        _e_re_cS(p["b->see::Re{cS}"], u),
        _e_im_cS(p["b->see::Im{cS}"], u),
        _e_re_cSprime(p["b->see::Re{cS'}"], u),
        _e_im_cSprime(p["b->see::Im{cS'}"], u),
        _e_re_cP(p["b->see::Re{cP}"], u),
        _e_im_cP(p["b->see::Im{cP}"], u),
        _e_re_cPprime(p["b->see::Re{cP'}"], u),
        _e_im_cPprime(p["b->see::Im{cP'}"], u),
        _e_re_cT(p["b->see::Re{cT}"], u),
        _e_im_cT(p["b->see::Im{cT}"], u),
        _e_re_cT5(p["b->see::Re{cT5}"], u),
        _e_im_cT5(p["b->see::Im{cT5}"], u),
        /* b->smumu */
        _mu_re_c9(p["b->smumu::Re{c9}"], u),
        _mu_im_c9(p["b->smumu::Im{c9}"], u),
        _mu_re_c10(p["b->smumu::Re{c10}"], u),
        _mu_im_c10(p["b->smumu::Im{c10}"], u),
        _mu_re_c9prime(p["b->smumu::Re{c9'}"], u),
        _mu_im_c9prime(p["b->smumu::Im{c9'}"], u),
        _mu_re_c10prime(p["b->smumu::Re{c10'}"], u),
        _mu_im_c10prime(p["b->smumu::Im{c10'}"], u),
        _mu_re_cS(p["b->smumu::Re{cS}"], u),
        _mu_im_cS(p["b->smumu::Im{cS}"], u),
        _mu_re_cSprime(p["b->smumu::Re{cS'}"], u),
        _mu_im_cSprime(p["b->smumu::Im{cS'}"], u),
        _mu_re_cP(p["b->smumu::Re{cP}"], u),
        _mu_im_cP(p["b->smumu::Im{cP}"], u),
        _mu_re_cPprime(p["b->smumu::Re{cP'}"], u),
        _mu_im_cPprime(p["b->smumu::Im{cP'}"], u),
        _mu_re_cT(p["b->smumu::Re{cT}"], u),
        _mu_im_cT(p["b->smumu::Im{cT}"], u),
        _mu_re_cT5(p["b->smumu::Re{cT5}"], u),
        _mu_im_cT5(p["b->smumu::Im{cT5}"], u),


        /* functions for b->sgamma */
        _c7(std::bind(&wcimplementation::cartesian,          _re_c7,          _im_c7)),
        _c7prime(std::bind(&wcimplementation::cartesian,     _re_c7prime,     _im_c7prime)),

        /* functions for b->see */
        _e_c9(std::bind(&wcimplementation::cartesian,        _e_re_c9,        _e_im_c9)),
        _e_c10(std::bind(&wcimplementation::cartesian,       _e_re_c10,       _e_im_c10)),
        _e_c9prime(std::bind(&wcimplementation::cartesian,   _e_re_c9prime,   _e_im_c9prime)),
        _e_c10prime(std::bind(&wcimplementation::cartesian,  _e_re_c10prime,  _e_im_c10prime)),
        _e_cS(std::bind(&wcimplementation::cartesian,        _e_re_cS,        _e_im_cS)),
        _e_cSprime(std::bind(&wcimplementation::cartesian,   _e_re_cSprime,   _e_im_cSprime)),
        _e_cP(std::bind(&wcimplementation::cartesian,        _e_re_cP,        _e_im_cP)),
        _e_cPprime(std::bind(&wcimplementation::cartesian,   _e_re_cPprime,   _e_im_cPprime)),
        _e_cT(std::bind(&wcimplementation::cartesian,        _e_re_cT,        _e_im_cT)),
        _e_cT5(std::bind(&wcimplementation::cartesian,       _e_re_cT5,       _e_im_cT5)),

        /* functions for b->smumu */
        _mu_c9(std::bind(&wcimplementation::cartesian,       _mu_re_c9,       _mu_im_c9)),
        _mu_c10(std::bind(&wcimplementation::cartesian,      _mu_re_c10,      _mu_im_c10)),
        _mu_c9prime(std::bind(&wcimplementation::cartesian,  _mu_re_c9prime,  _mu_im_c9prime)),
        _mu_c10prime(std::bind(&wcimplementation::cartesian, _mu_re_c10prime, _mu_im_c10prime)),
        _mu_cS(std::bind(&wcimplementation::cartesian,       _mu_re_cS,       _mu_im_cS)),
        _mu_cSprime(std::bind(&wcimplementation::cartesian,  _mu_re_cSprime,  _mu_im_cSprime)),
        _mu_cP(std::bind(&wcimplementation::cartesian,       _mu_re_cP,       _mu_im_cP)),
        _mu_cPprime(std::bind(&wcimplementation::cartesian,  _mu_re_cPprime,  _mu_im_cPprime)),
        _mu_cT(std::bind(&wcimplementation::cartesian,       _mu_re_cT,       _mu_im_cT)),
        _mu_cT5(std::bind(&wcimplementation::cartesian,      _mu_re_cT5,      _mu_im_cT5))
    {
    }

    WilsonCoefficients<BToS>
    WilsonScanComponent<components::DeltaBS1>::wilson_coefficients_b_to_s(const std::string & lepton_flavour, const bool & cp_conjugate) const
    {
        std::function<complex<double> ()> c9,  c9prime;
        std::function<complex<double> ()> c10, c10prime;
        std::function<complex<double> ()> cS,  cSprime;
        std::function<complex<double> ()> cP,  cPprime;
        std::function<complex<double> ()> cT,  cT5;

        if ("e" == lepton_flavour)
        {
            c9 = _e_c9;     c9prime = _e_c9prime;
            c10 = _e_c10;   c10prime = _e_c10prime;
            cS = _e_cS;     cSprime = _e_cSprime;
            cP = _e_cP;     cPprime = _e_cPprime;
            cT = _e_cT;     cT5 = _e_cT5;
        }
        else if ("mu" == lepton_flavour)
        {
            c9 = _mu_c9;    c9prime = _mu_c9prime;
            c10 = _mu_c10;  c10prime = _mu_c10prime;
            cS = _mu_cS;    cSprime = _mu_cSprime;
            cP = _mu_cP;    cPprime = _mu_cPprime;
            cT = _mu_cT;    cT5 = _mu_cT5;
        }
        else
        {
            throw InternalError("WilsonScan presently only implements 'e' and 'mu' lepton flavours");
        }

        double alpha_s = 0.0;
        if (_mu__deltabs1 < _mu_b__deltabs1)
        {
            alpha_s = QCD::alpha_s(_mu_b__deltabs1, _alpha_s_Z__deltabs1, _m_Z__deltabs1, QCD::beta_function_nf_5);
            alpha_s = QCD::alpha_s(_mu__deltabs1, alpha_s, _mu_b__deltabs1, QCD::beta_function_nf_4);
        }
        else
        {
            alpha_s = QCD::alpha_s(_mu__deltabs1, _alpha_s_Z__deltabs1, _m_Z__deltabs1, QCD::beta_function_nf_5);
        }

        complex<double> a_s = alpha_s / 4.0 / M_PI;
        WilsonCoefficients<BToS> result;
        result._sm_like_coefficients = std::array<std::complex<double>, 15>
        {
            _c1(), _c2(), _c3(), _c4(), _c5(), _c6(),
            0.0, 0.0, 0.0, 0.0, 0.0,
            a_s * _c7(), a_s * _c8(), a_s * c9(), a_s * c10()
        };
        result._primed_coefficients = std::array<std::complex<double>, 15>
        {
            /* we only consider c7', c8', c9' and c10' */
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            a_s * _c7prime(), a_s * _c8prime(), a_s * c9prime(), a_s * c10prime()
        };
        result._scalar_tensor_coefficients = std::array<std::complex<double>, 6>
        {
            cS(), cSprime(), cP(), cPprime(), cT(), cT5()
        };
        result._alpha_s = alpha_s;

        if (cp_conjugate)
        {
            for (auto & c : result._sm_like_coefficients)
            {
                c = conj(c);
            }

            for (auto & c : result._primed_coefficients)
            {
                c = conj(c);
            }

            for (auto & c : result._scalar_tensor_coefficients)
            {
                c = conj(c);
            }
        }

        return result;
    }

    /* b->u Wilson coefficients */
    WilsonScanComponent<components::DeltaBU1>::WilsonScanComponent(const Parameters & p, const Options &, ParameterUser & u) :
        _re_csl(p["b->ulnu::Re{cSL}"], u),
        _im_csl(p["b->ulnu::Im{cSL}"], u),
        _re_csr(p["b->ulnu::Re{cSR}"], u),
        _im_csr(p["b->ulnu::Im{cSR}"], u),
        _re_cvl(p["b->ulnu::Re{cVL}"], u),
        _im_cvl(p["b->ulnu::Im{cVL}"], u),
        _re_cvr(p["b->ulnu::Re{cVR}"], u),
        _im_cvr(p["b->ulnu::Im{cVR}"], u),
        _re_ct(p["b->ulnu::Re{cT}"], u),
        _im_ct(p["b->ulnu::Im{cT}"], u),
        _csl(std::bind(&wcimplementation::cartesian, _re_csl, _im_csl)),
        _csr(std::bind(&wcimplementation::cartesian, _re_csr, _im_csr)),
        _cvl(std::bind(&wcimplementation::cartesian, _re_cvl, _im_cvl)),
        _cvr(std::bind(&wcimplementation::cartesian, _re_cvr, _im_cvr)),
        _ct(std::bind(&wcimplementation::cartesian, _re_ct, _im_ct))
    {
    }

    WilsonCoefficients<BToU>
    WilsonScanComponent<components::DeltaBU1>::wilson_coefficients_b_to_u(const bool & cp_conjugate) const
    {
        WilsonCoefficients<BToU> result
        {
            {{
                _cvl(), _cvr(), _csl(), _csr(), _ct()
            }},
        };

        if (cp_conjugate)
        {
            for (auto c = result._coefficients.begin(), c_end = result._coefficients.end() ; c != c_end ; ++c)
            {
                *c = conj(*c);
            }
        }

        return result;
    }

    ConstrainedWilsonScanComponent::ConstrainedWilsonScanComponent(const Parameters & p, const Options & o, ParameterUser & u) :
        WilsonScanComponent<components::DeltaBD1>(p, o, u),
        WilsonScanComponent<components::DeltaBS1>(p, o, u)
    {
        /* b->smee */
        _e_cT = std::bind(&wcimplementation::zero);
        _e_cT5 = std::bind(&wcimplementation::zero);
        _e_cP = std::bind(&wcimplementation::cartesian_negative, _e_re_cS,        _e_im_cS);
        _e_cPprime = std::bind(&wcimplementation::cartesian,     _e_re_cSprime,   _e_im_cSprime);

        u.drop(_e_re_cP.id());       u.drop(_e_im_cP.id());
        u.drop(_e_re_cPprime.id());  u.drop(_e_im_cPprime.id());
        u.drop(_e_re_cT.id());       u.drop(_e_im_cT.id());
        u.drop(_e_re_cT5.id());      u.drop(_e_im_cT5.id());

        /* b->smumu */
        _mu_cT = std::bind(&wcimplementation::zero);
        _mu_cT5 = std::bind(&wcimplementation::zero);
        _mu_cP = std::bind(&wcimplementation::cartesian_negative, _mu_re_cS,      _mu_im_cS);
        _mu_cPprime = std::bind(&wcimplementation::cartesian,     _mu_re_cSprime, _mu_im_cSprime);

        u.drop(_mu_re_cP.id());      u.drop(_mu_im_cP.id());
        u.drop(_mu_re_cPprime.id()); u.drop(_mu_im_cPprime.id());
        u.drop(_mu_re_cT.id());      u.drop(_mu_im_cT.id());
        u.drop(_mu_re_cT5.id());     u.drop(_mu_im_cT5.id());
    }

    WilsonScanModel::WilsonScanModel(const Parameters & parameters, const Options & options) :
        SMComponent<components::CKM>(parameters, *this),
        SMComponent<components::QCD>(parameters, *this),
        WilsonScanComponent<components::DeltaBD1>(parameters, options, *this),
        WilsonScanComponent<components::DeltaBS1>(parameters, options, *this),
        WilsonScanComponent<components::DeltaBU1>(parameters, options, *this)
    {
    }

    WilsonScanModel::~WilsonScanModel()
    {
    }

    std::shared_ptr<Model>
    WilsonScanModel::make(const Parameters & parameters, const Options & options)
    {
        return std::shared_ptr<Model>(new WilsonScanModel(parameters, options));
    }

    ConstrainedWilsonScanModel::ConstrainedWilsonScanModel(const Parameters & parameters, const Options & options) :
        SMComponent<components::CKM>(parameters, *this),
        SMComponent<components::QCD>(parameters, *this),
        ConstrainedWilsonScanComponent(parameters, options, *this),
        WilsonScanComponent<components::DeltaBU1>(parameters, options, *this)
    {
    }

    ConstrainedWilsonScanModel::~ConstrainedWilsonScanModel()
    {
    }

    std::shared_ptr<Model>
    ConstrainedWilsonScanModel::make(const Parameters & parameters, const Options & options)
    {
        return std::shared_ptr<Model>(new ConstrainedWilsonScanModel(parameters, options));
    }
}
