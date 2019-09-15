/* vim: set sw=4 sts=4 tw=140 et foldmethod=marker : */

/*
 * Copyright (c) 2019 Danny van Dyk
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

#include <eos/form-factors/unitarity-bounds.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <cmath>

#include <iostream>

namespace eos
{
    template <> struct Implementation<BGLCoefficients>
    {
        // parameters for the leading Isgur-Wise function xi
        UsedParameter xipone, xippone, xipppone;

        // parameters for the subleading Isgur-Wise function chi2
        UsedParameter chi2one, chi2pone, chi2ppone;

        // parameters for the subleading Isgur-Wise function chi3
        UsedParameter chi3pone, chi3ppone;

        // parameters for the subleading Isgur-Wise function eta
        UsedParameter etaone, etapone, etappone;

        // parameters for subsubleading 1/mc corrections in h+ (B->D), equal to delta{h+}
        UsedParameter l1one, l1pone, l1ppone;

        // parameters for subsubleading 1/mc corrections in hA1 (B->D^*), equal to delta{A1}
        UsedParameter l2one, l2pone, l2ppone;

        // parameters for subsubleading 1/m_c corrections
        UsedParameter l3one, l3pone, l3ppone;
        UsedParameter l4one, l4pone, l4ppone;
        UsedParameter l5one, l5pone, l5ppone;
        UsedParameter l6one, l6pone, l6ppone;

        Implementation(const Parameters & p, const Options & /*o*/, ParameterUser & u) :
            xipone(p["B(*)->D(*)::xi'(1)@HQET"], u),
            xippone(p["B(*)->D(*)::xi''(1)@HQET"], u),
            xipppone(p["B(*)->D(*)::xi'''(1)@HQET"], u),
            chi2one(p["B(*)->D(*)::chi_2(1)@HQET"], u),
            chi2pone(p["B(*)->D(*)::chi_2'(1)@HQET"], u),
            chi2ppone(p["B(*)->D(*)::chi_2''(1)@HQET"], u),
            chi3pone(p["B(*)->D(*)::chi_3'(1)@HQET"], u),
            chi3ppone(p["B(*)->D(*)::chi_3''(1)@HQET"], u),
            etaone(p["B(*)->D(*)::eta(1)@HQET"], u),
            etapone(p["B(*)->D(*)::eta'(1)@HQET"], u),
            etappone(p["B(*)->D(*)::eta''(1)@HQET"], u),
            l1one(p["B(*)->D(*)::l_1(1)@HQET"], u),
            l1pone(p["B(*)->D(*)::l_1'(1)@HQET"], u),
            l1ppone(p["B(*)->D(*)::l_1''(1)@HQET"], u),
            l2one(p["B(*)->D(*)::l_2(1)@HQET"], u),
            l2pone(p["B(*)->D(*)::l_2'(1)@HQET"], u),
            l2ppone(p["B(*)->D(*)::l_2''(1)@HQET"], u),
            l3one(p["B(*)->D(*)::l_3(1)@HQET"], u),
            l3pone(p["B(*)->D(*)::l_3'(1)@HQET"], u),
            l3ppone(p["B(*)->D(*)::l_3''(1)@HQET"], u),
            l4one(p["B(*)->D(*)::l_4(1)@HQET"], u),
            l4pone(p["B(*)->D(*)::l_4'(1)@HQET"], u),
            l4ppone(p["B(*)->D(*)::l_4''(1)@HQET"], u),
            l5one(p["B(*)->D(*)::l_5(1)@HQET"], u),
            l5pone(p["B(*)->D(*)::l_5'(1)@HQET"], u),
            l5ppone(p["B(*)->D(*)::l_5''(1)@HQET"], u),
            l6one(p["B(*)->D(*)::l_6(1)@HQET"], u),
            l6pone(p["B(*)->D(*)::l_6'(1)@HQET"], u),
            l6ppone(p["B(*)->D(*)::l_6''(1)@HQET"], u)
        {
        }

        ~Implementation() = default;

        /*
         * HQET parameters following [BLPR2017]
         */
        inline double _mu() const { return 2.31; } // mu^2 = m_b * m_c
        inline double _alpha_s() const { return 0.26; }
        inline double _m_b_1S() const { return 4.71; }
        inline double _m_b_pole() const { return _m_b_1S() * (1 + 2.0 / 9.0 * power_of<2>(_alpha_s())); }
        inline double _m_c_pole() const { return _m_b_pole() - 3.40; }
        inline double _lambda_1() const { return -0.30; }
        inline double _LambdaBar() const { return 5.313 - _m_b_pole() + _lambda_1() / (2.0 * _m_b_1S()); }
        inline double _eps_b() const { return _LambdaBar() / (2.0 * _m_b_pole()); }
        inline double _eps_c() const { return _LambdaBar() / (2.0 * _m_c_pole()); }

        // B -> D form factors
        // {{{
        double V1_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.008703312831206976 + 0.003212569766727432*as + -0.004155423451479522*epsb + 0.004155423451479522*epsc +
            0.008310846902959045*epsb*etaone + -0.008310846902959045*epsc*etaone + 0.008703312831206976*l1one*pow(epsc,2) +
            -0.004155423451479522*l4one*pow(epsc,2);
        }

        double V1_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.060227992698537264 + 0.029238681829861973*as + -0.02875604015951836*epsb + -0.27850601059862323*chi2one*epsb +
            0.8355180317958697*chi3pone*epsb + 0.02875604015951836*epsc + -0.27850601059862323*chi2one*epsc +
            0.8355180317958697*chi3pone*epsc + 0.05751208031903672*epsb*etaone + -0.05751208031903672*epsc*etaone +
            0.06648677522367237*epsb*etapone + -0.06648677522367237*epsc*etapone + 0.06962650264965581*xipone +
            0.025700558133819457*as*xipone + -0.033243387611836185*epsb*xipone + 0.033243387611836185*epsc*xipone +
            0.06648677522367237*epsb*etaone*xipone + -0.06648677522367237*epsc*etaone*xipone + 0.060227992698537264*l1one*pow(epsc,2) +
            0.06962650264965581*l1pone*pow(epsc,2) + -0.02875604015951836*l4one*pow(epsc,2) + -0.033243387611836185*l4pone*pow(epsc,2) +
            0.06962650264965581*l1one*xipone*pow(epsc,2) + -0.033243387611836185*l4one*xipone*pow(epsc,2);
        }

        double V1_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();


            return 0.1306576264897831 + 0.015976186987051347*as + -0.062382885202461706*epsb + -2.4843077875504394*chi2one*epsb +
            -2.2280480847889863*chi2pone*epsb + 7.452923362651318*chi3pone*epsb + 3.3420721271834792*chi3ppone*epsb +
            0.062382885202461706*epsc + -2.4843077875504394*chi2one*epsc + -2.2280480847889863*chi2pone*epsc +
            7.452923362651318*chi3pone*epsc + 3.3420721271834792*chi3ppone*epsc + 0.12476577040492341*epsb*etaone +
            -0.12476577040492341*epsc*etaone + 0.5930701929996385*epsb*etapone + -0.5930701929996385*epsc*etapone +
            0.2659471008946895*epsb*etappone + -0.2659471008946895*epsc*etappone + 0.6210769468876098*xipone + 0.2853105709065347*as*xipone
            + -0.29653509649981924*epsb*xipone + -2.2280480847889863*chi2one*epsb*xipone + 6.6841442543669585*chi3pone*epsb*xipone +
            0.29653509649981924*epsc*xipone + -2.2280480847889863*chi2one*epsc*xipone + 6.6841442543669585*chi3pone*epsc*xipone +
            0.5930701929996385*epsb*etaone*xipone + -0.5930701929996385*epsc*etaone*xipone + 0.531894201789379*epsb*etapone*xipone +
            -0.531894201789379*epsc*etapone*xipone + 0.2785060105986233*xippone + 0.10280223253527784*as*xippone +
            -0.13297355044734474*epsb*xippone + 0.13297355044734474*epsc*xippone + 0.2659471008946895*epsb*etaone*xippone +
            -0.2659471008946895*epsc*etaone*xippone + 0.1306576264897831*l1one*pow(epsc,2) + 0.48182394158829817*l1pone*pow(epsc,2) +
            -0.062382885202461706*l4one*pow(epsc,2) + -0.23004832127614688*l4pone*pow(epsc,2) + 0.6210769468876098*l1one*xipone*pow(epsc,2)
            + 0.5570120211972466*l1pone*xipone*pow(epsc,2) + -0.29653509649981924*l4one*xipone*pow(epsc,2) +
            -0.2659471008946895*l4pone*xipone*pow(epsc,2) + 0.2785060105986233*l1one*xippone*pow(epsc,2) +
            -0.13297355044734474*l4one*xippone*pow(epsc,2);
        }

        double S1_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();

            return 0.04733475740152116 + 0.011600723452028892*as + 0.04733475740152116*l1one*pow(epsc,2);
        }

        double S1_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();


            return 0.2843488842022793 + 0.21589749983616183*as + -0.39656050100795237*epsb + -1.5147122368486774*chi2one*epsb +
            4.544136710546032*chi3pone*epsb + 0.39656050100795237*epsc + -1.5147122368486774*chi2one*epsc + 4.544136710546032*chi3pone*epsc
            + 0.7931210020159047*epsb*etaone + -0.7931210020159047*epsc*etaone + 0.37867805921216935*xipone + 0.09280578761623114*as*xipone
            + 0.2843488842022793*l1one*pow(epsc,2) + 0.37867805921216935*l1pone*pow(epsc,2) + -0.39656050100795237*l4one*pow(epsc,2) +
            0.37867805921216935*l1one*xipone*pow(epsc,2);
        }

        double S1_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();


            return 0.6559703043190044 + 0.3219759165015938*as + -1.5890932982243737*epsb + -12.128588768170292*chi2one*epsb +
            -12.117697894789417*chi2pone*epsb + 36.385766304510874*chi3pone*epsb + 18.176546842184127*chi3ppone*epsb +
            1.5890932982243737*epsc + -12.128588768170292*chi2one*epsc + -12.117697894789417*chi2pone*epsc +
            36.385766304510874*chi3pone*epsc + 18.176546842184127*chi3ppone*epsc + 3.1781865964487475*epsb*etaone +
            -3.1781865964487475*epsc*etaone + 6.344968016127237*epsb*etapone + -6.344968016127237*epsc*etapone + 3.032147192042573*xipone +
            1.9127915739217567*as*xipone + -3.1724840080636185*epsb*xipone + -12.117697894789417*chi2one*epsb*xipone +
            36.353093684368254*chi3pone*epsb*xipone + 3.1724840080636185*epsc*xipone + -12.117697894789417*chi2one*epsc*xipone +
            36.353093684368254*chi3pone*epsc*xipone + 6.344968016127237*epsb*etaone*xipone + -6.344968016127237*epsc*etaone*xipone +
            1.5147122368486772*xippone + 0.37122315046492455*as*xippone + 0.6559703043190044*l1one*pow(epsc,2) +
            2.274791073618234*l1pone*pow(epsc,2) + -1.5890932982243737*l4one*pow(epsc,2) + -3.1724840080636185*l4pone*pow(epsc,2) +
            3.032147192042573*l1one*xipone*pow(epsc,2) + 3.0294244736973543*l1pone*xipone*pow(epsc,2) +
            -3.1724840080636185*l4one*xipone*pow(epsc,2) + 1.5147122368486772*l1one*xippone*pow(epsc,2);
        }
        // }}}

        // B -> D^* form factors
        // {{{
        double A1_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();

            return 0.008363045281441904 + -0.003525762321913306*as + 0.008363045281441904*l2one*pow(epsc,2);
        }

        double A1_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.0812300174749896 + 0.033452181125767616*epsb + -0.26761744900614093*chi2one*epsb + 0.8028523470184228*chi3pone*epsb +
            0.033452181125767616*epsc + -0.26761744900614093*chi3pone*epsc + -0.06690436225153523*epsb*etaone + as*(0.00850018218230231 -
            0.028206098575306444*xipone) + 0.06690436225153523*xipone + (0.06690436225153523*l2pone - 0.033452181125767616*l5one +
            l2one*(0.0812300174749896 + 0.06690436225153523*xipone))*pow(epsc,2);
        }

        double A1_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.31721556277653973 + 0.25801570764842313*epsb + -3.134595457211949*chi2one*epsb + -2.1409395920491274*chi2pone*epsb +
            9.403786371635848*chi3pone*epsb + 3.211409388073691*chi3ppone*epsb + 0.25801570764842313*epsc + -3.134595457211949*chi3pone*epsc
            + -1.0704697960245637*chi3ppone*epsc + -0.5160314152968463*epsb*etaone + -0.5352348980122819*epsb*etapone +
            0.7836488643029873*xipone + 0.26761744900614093*epsb*xipone + -2.1409395920491274*chi2one*epsb*xipone +
            6.422818776147382*chi3pone*epsb*xipone + 0.26761744900614093*epsc*xipone + -2.1409395920491274*chi3pone*epsc*xipone +
            -0.5352348980122819*epsb*etaone*xipone + as*(0.17093446696648967 + 0.011589260307805627*xipone - 0.11282439430122577*xippone) +
            0.26761744900614093*xippone + (0.6498401397999168*l2pone - 0.25801570764842313*l5one - 0.26761744900614093*l5pone +
            0.5352348980122819*l2pone*xipone - 0.26761744900614093*l5one*xipone + l2one*(0.31721556277653973 + 0.7836488643029873*xipone +
            0.26761744900614093*xippone))*pow(epsc,2);
        }

        double A5_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();

            return 0.005602202787090022 + -0.0023618233360843633*as + 0.005602202787090022*l2one*pow(epsc,2);
        }

        double A5_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.043523865623533844 + -0.0499621450168773*epsb + -0.17927048918688068*chi2one*epsb + 0.537811467560642*chi3pone*epsb +
            0.04996214501687729*epsc + 0.17927048918688068*chi2one*epsc + -0.17927048918688068*chi3pone*epsc +
            0.0999242900337546*epsb*etaone + 0.0999242900337546*epsc*etaone + as*(0.012992295758625651 - 0.018894586688674903*xipone) +
            0.04481762229672017*xipone + (0.04481762229672017*l2pone + 0.04481762229672017*l3one + 0.0499621450168773*l5one -
            0.0999242900337546*l6one + l2one*(0.043523865623533844 + 0.04481762229672017*xipone))*pow(epsc,2);
        }

        double A5_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.11731273819106425 + -0.28823475536695886*epsb + -1.7513046783268447*chi2one*epsb + -1.4341639134950457*chi2pone*epsb +
            5.253914034980533*chi3pone*epsb + 2.1512458702425685*chi3ppone*epsb + 0.28823475536695875*epsc + 1.7513046783268447*chi2one*epsc
            + 1.4341639134950457*chi2pone*epsc + -1.7513046783268447*chi3pone*epsc + -0.7170819567475228*chi3ppone*epsc +
            0.5764695107339177*epsb*etaone + 0.5764695107339177*epsc*etaone + 0.7993943202700368*epsb*etapone +
            0.7993943202700368*epsc*etapone + 0.43782616958171117*xipone + -0.3996971601350184*epsb*xipone +
            -1.4341639134950457*chi2one*epsb*xipone + 4.302491740485137*chi3pone*epsb*xipone + 0.3996971601350183*epsc*xipone +
            1.4341639134950457*chi2one*epsc*xipone + -1.4341639134950457*chi3pone*epsc*xipone + 0.7993943202700368*epsb*etaone*xipone +
            0.7993943202700368*epsc*etaone*xipone + as*(0.11394596757274658 + 0.06614919269165541*xipone - 0.07557834675469961*xippone) +
            0.1792704891868807*xippone + (0.34819092498827076*l2pone + 0.43782616958171117*l3one + 0.3585409783737614*l3pone +
            0.28823475536695886*l5one + 0.3996971601350184*l5pone - 0.976166670868936*l6one - 0.7993943202700368*l6pone +
            0.3585409783737614*l2pone*xipone + 0.3585409783737614*l3one*xipone + 0.3996971601350184*l5one*xipone -
            0.7993943202700368*l6one*xipone + l2one*(0.11731273819106425 + 0.43782616958171117*xipone +
            0.1792704891868807*xippone))*pow(epsc,2);
        }

        double V4_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.007502586789045188*1 + 0.9117449999857681*as + epsb + epsc - 2*epsb*etaone + (l2one - l5one)*pow(epsc,2);
        }

        double V4_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.011671907229776847*as + -0.24008277724944604*chi2one*epsb + 0.720248331748338*chi3pone*epsb +
            -0.24008277724944604*chi3pone*epsc + -0.12004138862472302*epsb*etapone + (0.06002069431236151*l2pone -
            0.06002069431236151*l5pone)*pow(epsc,2) + 0.10292114090453137*(0.5793766430333188 + 0.5831716767309848*xipone)*(1 +
            0.9117449999857681*as + epsb + epsc - 2*epsb*etaone + (l2one - l5one)*pow(epsc,2));
        }

        double V4_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return -0.06713773492653477*as + -0.48016555449889203*chi2one*epsb + -1.9206622179955681*chi2pone*epsb +
            1.440496663496676*chi3pone*epsb + 2.880993326993352*chi3ppone*epsb + -0.48016555449889203*chi3pone*epsc +
            -0.9603311089977841*chi3ppone*epsc + -0.24008277724944602*epsb*etapone + -0.48016555449889203*epsb*etappone +
            0.10292114090453137*(1.4999940008253256 + 5.80135649772852*xipone + 2.332686706923939*xippone)*(1 + 0.9117449999857681*as + epsb
            + epsc - 2*epsb*etaone + (l2one - l5one)*pow(epsc,2)) + 0.09337525783821478*(0.9934924245310758 + 1.*xipone)*(1.*as -
            20.569284224343182*chi2one*epsb + 61.70785267302954*chi3pone*epsb - 20.569284224343182*chi3pone*epsc -
            10.284642112171591*epsb*etapone + (5.1423210560857955*l2pone - 5.1423210560857955*l5pone)*pow(epsc,2));
        }

        double P1_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.031488567512709106 + -0.0022676800849227453*as + -0.014123119863775365*epsb + 0.014123119863775365*epsc +
            0.02824623972755073*epsb*etaone + 0.02824623972755073*epsc*etaone + (0.031488567512709106*l2one + 0.014123119863775365*l5one -
            0.02824623972755073*l6one)*pow(epsc,2);
        }

        double P1_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.26142777033909864 + -0.11725448401958674*epsb + -1.0076341604066916*chi2one*epsb + 3.022902481220075*chi3pone*epsb +
            0.11725448401958674*epsc + 1.0076341604066916*chi2one*epsc + -1.0076341604066916*chi3pone*epsc + 0.23450896803917348*epsb*etaone
            + 0.23450896803917348*epsc*etaone + 0.22596991782040587*epsb*etapone + 0.22596991782040587*epsc*etapone +
            as*(-0.045688099445631854 - 0.018141440679381966*xipone) + 0.2519085401016729*xipone + -0.11298495891020294*epsb*xipone +
            0.11298495891020294*epsc*xipone + 0.22596991782040587*epsb*etaone*xipone + 0.22596991782040587*epsc*etaone*xipone +
            (0.2519085401016729*l2pone + 0.2519085401016729*l3one + 0.11725448401958674*l5one + 0.11298495891020294*l5pone -
            0.3474939269493764*l6one - 0.22596991782040587*l6pone + l2one*(0.26142777033909864 + 0.2519085401016729*xipone) +
            0.11298495891020294*l5one*xipone - 0.22596991782040587*l6one*xipone)*pow(epsc,2);
        }

        double P1_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.7662176967652994 + -0.34366073873619646*epsb + -10.38095697166454*chi2one*epsb + -8.061073283253533*chi2pone*epsb +
            31.142870914993626*chi3pone*epsb + 12.091609924880299*chi3ppone*epsb + 0.34366073873619657*epsc + 10.38095697166454*chi2one*epsc
            + 8.061073283253533*chi2pone*epsc + -10.38095697166454*chi3pone*epsc + -4.0305366416267665*chi3ppone*epsc +
            0.6873214774723929*epsb*etaone + 0.687321477472393*epsc*etaone + 2.3280115799542*epsb*etapone + 2.3280115799542*epsc*etapone +
            0.9038796712816234*epsb*etappone + 0.9038796712816234*epsc*etappone + 2.595239242916135*xipone + -1.1640057899771*epsb*xipone +
            -8.061073283253533*chi2one*epsb*xipone + 24.183219849760597*chi3pone*epsb*xipone + 1.1640057899771*epsc*xipone +
            8.061073283253533*chi2one*epsc*xipone + -8.061073283253533*chi3pone*epsc*xipone + 2.3280115799542*epsb*etaone*xipone +
            2.3280115799542*epsc*etaone*xipone + 1.8077593425632468*epsb*etapone*xipone + 1.8077593425632468*epsc*etapone*xipone +
            as*(-0.5977425693276348 - 0.401787676923819*xipone - 0.07256576271752786*xippone) + 1.0076341604066916*xippone +
            -0.4519398356408117*epsb*xippone + 0.4519398356408117*epsc*xippone + 0.9038796712816234*epsb*etaone*xippone +
            0.9038796712816234*epsc*etaone*xippone + (2.091422162712789*l2pone + 2.595239242916135*l3one + 2.0152683208133833*l3pone +
            0.34366073873619646*l5one + 0.9380358721566939*l5pone - 1.8513272674494927*l6one - 2.779951415595011*l6pone +
            2.0152683208133833*l2pone*xipone + 2.0152683208133833*l3one*xipone + 1.1640057899771*l5one*xipone +
            0.9038796712816234*l5pone*xipone - 3.2318912512358233*l6one*xipone - 1.8077593425632468*l6pone*xipone +
            0.4519398356408117*l5one*xippone - 0.9038796712816234*l6one*xippone + l2one*(0.7662176967652994 + 2.595239242916135*xipone +
            1.0076341604066916*xippone))*pow(epsc,2);
        }
        // }}}
    };

    BGLCoefficients::BGLCoefficients(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<BGLCoefficients>(new Implementation<BGLCoefficients>(p, o, *this))
    {
    }

    BGLCoefficients::~BGLCoefficients() = default;

    // B -> D form factors
    // {{{
    double BGLCoefficients::V1_a0() const { return _imp->V1_a0(); }
    double BGLCoefficients::V1_a1() const { return _imp->V1_a1(); }
    double BGLCoefficients::V1_a2() const { return _imp->V1_a2(); }
    double BGLCoefficients::S1_a0() const { return _imp->S1_a0(); }
    double BGLCoefficients::S1_a1() const { return _imp->S1_a1(); }
    double BGLCoefficients::S1_a2() const { return _imp->S1_a2(); }
    // }}}

    // B -> D^* form factors
    // {{{
    double BGLCoefficients::A1_a0() const { return _imp->A1_a0(); }
    double BGLCoefficients::A1_a1() const { return _imp->A1_a1(); }
    double BGLCoefficients::A1_a2() const { return _imp->A1_a2(); }
    double BGLCoefficients::A5_a0() const { return _imp->A5_a0(); }
    double BGLCoefficients::A5_a1() const { return _imp->A5_a1(); }
    double BGLCoefficients::A5_a2() const { return _imp->A5_a2(); }
    double BGLCoefficients::V4_a0() const { return _imp->V4_a0(); }
    double BGLCoefficients::V4_a1() const { return _imp->V4_a1(); }
    double BGLCoefficients::V4_a2() const { return _imp->V4_a2(); }
    double BGLCoefficients::P1_a0() const { return _imp->P1_a0(); }
    double BGLCoefficients::P1_a1() const { return _imp->P1_a1(); }
    double BGLCoefficients::P1_a2() const { return _imp->P1_a2(); }
    // }}}

    template <> struct Implementation<HQETUnitarityBounds>
    {
        // option to determine if we use z^3 terms in the leading-power IW function
        SwitchOption opt_zorder_bound;
        std::function<double ()> bound_0p;
        std::function<double ()> bound_0m;
        std::function<double ()> bound_1p;
        std::function<double ()> bound_1m;

        // parameters for the leading Isgur-Wise function xi
        UsedParameter xipone, xippone, xipppone;

        // parameters for the subleading Isgur-Wise function chi2
        UsedParameter chi2one, chi2pone, chi2ppone;

        // parameters for the subleading Isgur-Wise function chi3
        UsedParameter chi3pone, chi3ppone;

        // parameters for the subleading Isgur-Wise function eta
        UsedParameter etaone, etapone, etappone;

        // parameters for subsubleading 1/mc corrections in h+ (B->D), equal to delta{h+}
        UsedParameter l1one, l1pone, l1ppone;

        // parameters for subsubleading 1/mc corrections in hA1 (B->D^*), equal to delta{A1}
        UsedParameter l2one, l2pone, l2ppone;

        // parameters for subsubleading 1/m_c corrections
        UsedParameter l3one, l3pone, l3ppone;
        UsedParameter l4one, l4pone, l4ppone;
        UsedParameter l5one, l5pone, l5ppone;
        UsedParameter l6one, l6pone, l6ppone;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_zorder_bound(o, "z-order-bound", { "1", "2" }, "2"),
            xipone(p["B(*)->D(*)::xi'(1)@HQET"], u),
            xippone(p["B(*)->D(*)::xi''(1)@HQET"], u),
            xipppone(p["B(*)->D(*)::xi'''(1)@HQET"], u),
            chi2one(p["B(*)->D(*)::chi_2(1)@HQET"], u),
            chi2pone(p["B(*)->D(*)::chi_2'(1)@HQET"], u),
            chi2ppone(p["B(*)->D(*)::chi_2''(1)@HQET"], u),
            chi3pone(p["B(*)->D(*)::chi_3'(1)@HQET"], u),
            chi3ppone(p["B(*)->D(*)::chi_3''(1)@HQET"], u),
            etaone(p["B(*)->D(*)::eta(1)@HQET"], u),
            etapone(p["B(*)->D(*)::eta'(1)@HQET"], u),
            etappone(p["B(*)->D(*)::eta''(1)@HQET"], u),
            l1one(p["B(*)->D(*)::l_1(1)@HQET"], u),
            l1pone(p["B(*)->D(*)::l_1'(1)@HQET"], u),
            l1ppone(p["B(*)->D(*)::l_1''(1)@HQET"], u),
            l2one(p["B(*)->D(*)::l_2(1)@HQET"], u),
            l2pone(p["B(*)->D(*)::l_2'(1)@HQET"], u),
            l2ppone(p["B(*)->D(*)::l_2''(1)@HQET"], u),
            l3one(p["B(*)->D(*)::l_3(1)@HQET"], u),
            l3pone(p["B(*)->D(*)::l_3'(1)@HQET"], u),
            l3ppone(p["B(*)->D(*)::l_3''(1)@HQET"], u),
            l4one(p["B(*)->D(*)::l_4(1)@HQET"], u),
            l4pone(p["B(*)->D(*)::l_4'(1)@HQET"], u),
            l4ppone(p["B(*)->D(*)::l_4''(1)@HQET"], u),
            l5one(p["B(*)->D(*)::l_5(1)@HQET"], u),
            l5pone(p["B(*)->D(*)::l_5'(1)@HQET"], u),
            l5ppone(p["B(*)->D(*)::l_5''(1)@HQET"], u),
            l6one(p["B(*)->D(*)::l_6(1)@HQET"], u),
            l6pone(p["B(*)->D(*)::l_6'(1)@HQET"], u),
            l6ppone(p["B(*)->D(*)::l_6''(1)@HQET"], u)
        {
            if ("1" == opt_zorder_bound.value())
            {
                bound_0p = std::bind(&Implementation::bound_0p_z1, this);
                bound_0m = std::bind(&Implementation::bound_0m_z1, this);
                bound_1p = std::bind(&Implementation::bound_1p_z1, this);
                bound_1m = std::bind(&Implementation::bound_1m_z1, this);
            }
            else if ("2" == opt_zorder_bound.value())
            {
                bound_0p = std::bind(&Implementation::bound_0p_z2, this);
                bound_0m = std::bind(&Implementation::bound_0m_z2, this);
                bound_1p = std::bind(&Implementation::bound_1p_z2, this);
                bound_1m = std::bind(&Implementation::bound_1m_z2, this);
            }
            else
            {
                throw InternalError("Only z-order-bound=2 is presently supported");
            }
        }

        ~Implementation() = default;

        // bounds for z^1
        // {{{
        double bound_0p_z1() const
        {
            const std::array<double, 3> weak_bounds
            {
                pow(0.07787308942370388 + 0.0763250029248492*l1one*pow(0.18120460196719873,2),2) + 0.02239357695917938*pow(3.8071054412766747 - 3.811666973787887*chi2one + 11.43500092136366*chi3pone - 1.1013264558178975*etaone + 4.1630886290205975*xipone + 3.0639132699149516*l1one*pow(0.18120460196719873,2) + 4.080328058612812*l1pone*pow(0.18120460196719873,2) - 4.273014767654388*l4one*pow(0.18120460196719873,2) + 4.080328058612812*l1one*xipone*pow(0.18120460196719873,2),2),
                pow(0.06378188590605295 + 0.06251392701070964*l2one*pow(0.18120460196719873,2),2) + 0.015022540138963305*pow(4.4528226929609165 - 3.811666973787886*chi3pone + 4.163088629020597*xipone + 3.661246265806994*l2one*pow(0.18120460196719873,2) + 4.080328058612811*l2pone*pow(0.18120460196719873,2) - 4.514209096519603*l5one*pow(0.18120460196719873,2) + 4.080328058612811*l2one*xipone*pow(0.18120460196719873,2),2),
                pow(0.04510060404103671 + 0.04420402170787365*l2one*pow(0.18120460196719873,2),2) + 0.007511270069481652*pow(4.4528226929609165 + 3.811666973787886*chi2one - 3.811666973787886*chi3pone + 1.163491861232192*etaone + 4.163088629020597*xipone + 3.661246265806994*l2one*pow(0.18120460196719873,2) + 4.080328058612811*l2pone*pow(0.18120460196719873,2) + 4.080328058612811*l3one*pow(0.18120460196719873,2) + 4.514209096519603*l5one*pow(0.18120460196719873,2) - 9.028418193039206*l6one*pow(0.18120460196719873,2) + 4.080328058612811*l2one*xipone*pow(0.18120460196719873,2),2),
            };

            double result = 0.0;
            for (const auto & b : weak_bounds)
            {
                result += b;
            }

            if (result < 0.0)
            {
                throw InternalError("Contribution to 0^+ unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= result) && (result < 1.0))
            {
                return 0.0;
            }
            else
            {
                // add an r-fit like penalty
                static const double sigma = 0.0130561; // cf. [BG2016], eq. (2.8), p.5
                return -pow((result - 1.0) / sigma, 2) / 2.0;
            }
        }

        double bound_0m_z1() const
        {
            const std::array<double, 3> weak_bounds
            {
                pow(0.053405909191163586 + 0.010636706864722678*etaone + 0.05077378948243221*l2one*pow(0.18120460196719873,2) + 0.022772846510374423*l5one*pow(0.18120460196719873,2) - 0.045545693020748845*l6one*pow(0.18120460196719873,2),2) + 0.0008556908618011906*pow(-15.035034001006752 - 7.157853274953061*chi2one + 1.3441830910079258*chi3pone - 3.018891875414771*etaone - 2.9089665725799128*etapone - 14.60564877749723*xipone - 2.9089665725799128*etaone*xipone - 14.410529734621102*l2one*pow(0.18120460196719873,2) - 13.88580678644618*l2pone*pow(0.18120460196719873,2) - 13.88580678644618*l3one*pow(0.18120460196719873,2) - 6.463350187664441*l5one*pow(0.18120460196719873,2) - 6.228003657868921*l5pone*pow(0.18120460196719873,2) + 19.154704033197802*l6one*pow(0.18120460196719873,2) + 12.456007315737843*l6pone*pow(0.18120460196719873,2) - 13.88580678644618*l2one*xipone*pow(0.18120460196719873,2) - 6.228003657868921*l5one*xipone*pow(0.18120460196719873,2) + 12.456007315737843*l6one*xipone*pow(0.18120460196719873,2),2),
                pow(0.07480998702772482 - 0.015876875548846626*etaone + 0.07070378019068203*l1one*pow(0.18120460196719873,2) - 0.03399187873996447*l4one*pow(0.18120460196719873,2),2) + 0.0007031761385459602*pow(-21.02575811558329 + 10.99543036092901*chi2one - 41.916874484470014*chi3pone + 4.500537596414145*etaone + 4.789861601196494*etapone - 22.56926957370638*xipone + 4.789861601196494*etaone*xipone - 20.04204290559436*l1one*pow(0.18120460196719873,2) - 21.330476563405902*l1pone*pow(0.18120460196719873,2) + 9.63550591369812*l4one*pow(0.18120460196719873,2) + 10.254939281231*l4pone*pow(0.18120460196719873,2) - 21.330476563405902*l1one*xipone*pow(0.18120460196719873,2) + 10.254939281231*l4one*xipone*pow(0.18120460196719873,2),2),
                pow(0.06252207819688696 + 0.059403241258809435*l2one*pow(0.18120460196719873,2) - 0.02684686363639445*l5one*pow(0.18120460196719873,2),2) + 0.0010860026569941262*pow(-16.542981584936765 + 13.471146304272825*chi3pone - 15.177769336722374*xipone - 15.838121246690301*l2one*pow(0.18120460196719873,2) - 14.420644989449043*l2pone*pow(0.18120460196719873,2) + 7.1579239172158635*l5one*pow(0.18120460196719873,2) + 6.517305813227482*l5pone*pow(0.18120460196719873,2) - 14.420644989449043*l2one*xipone*pow(0.18120460196719873,2) + 6.517305813227482*l5one*xipone*pow(0.18120460196719873,2),2),
            };

            double result = 0.0;
            for (const auto & b : weak_bounds)
            {
                result += b;
            }

            if (result < 0.0)
            {
                throw InternalError("Contribution to 0^- unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= result) && (result < 1.0))
            {
                return 0.0;
            }
            else
            {
                static const double sigma = 0.0130561; // using the same relative uncertainty as for 0^+, cf. [BG2016], eq. (2.8), p.5
                return -pow((result - 1.0) / sigma, 2) / 2.0;
            }
        }

        double bound_1p_z1() const
        {
            const std::array<double, 3> weak_bounds
            {
                pow(0.00895448576377894 + 0.009033280565551251*l2one*pow(0.18120460196719873,2),2) + pow(0.0133673793613974 + 0.01348500532393532*l2one*pow(0.18120460196719873,2),2) + 0.000928214055953212*pow(-2.6585020549037592 - 1.2227089458123461*chi2one + 0.22961419116206505*chi3pone - 1.2350756258163746*etaone - 2.351292004199026*xipone - 2.3035098010173103*l2one*pow(0.18120460196719873,2) - 2.3719821468009723*l2pone*pow(0.18120460196719873,2) - 2.3719821468009723*l3one*pow(0.18120460196719873,2) - 2.6442571007294786*l5one*pow(0.18120460196719873,2) + 5.288514201458957*l6one*pow(0.18120460196719873,2) - 2.3719821468009723*l2one*xipone*pow(0.18120460196719873,2),2) + pow(0.1438601613520823 - 0.02258345505963029*chi2one - 0.010443075533085308*chi3pone - 0.005645863764907572*etaone + 0.1069390348911792*xipone + 0.13097946755643228*l2one*pow(0.18120460196719873,2) + 0.10788004259148257*l2pone*pow(0.18120460196719873,2) - 0.05394002129574128*l5one*pow(0.18120460196719873,2) + 0.10788004259148257*l2one*xipone*pow(0.18120460196719873,2),2),
                pow(0.01563097879202461 + 0.01576852324827991*l1one*pow(0.18120460196719873,2),2) + pow(0.021574370557675353 + 0.021764213759875837*l1one*pow(0.18120460196719873,2),2) + 0.00282838931341691*pow(-2.294683910554761 + 1.2227089458123461*chi2one - 4.661221592087319*chi3pone + 1.1522288822404556*etaone - 2.351292004199026*xipone - 1.962434781465052*l1one*pow(0.18120460196719873,2) - 2.3719821468009723*l1pone*pow(0.18120460196719873,2) + 2.46688489339753*l4one*pow(0.18120460196719873,2) - 2.3719821468009723*l1one*xipone*pow(0.18120460196719873,2),2) + pow(0.20716060223145838 - 0.1262008221275946*chi2one + 0.34215374934135023*chi3pone - 0.03155020553189865*etaone + 0.17259496446140282*xipone + 0.1861515276883168*l1one*pow(0.18120460196719873,2) + 0.1741137100790067*l1pone*pow(0.18120460196719873,2) - 0.08705685503950335*l4one*pow(0.18120460196719873,2) + 0.1741137100790067*l1one*xipone*pow(0.18120460196719873,2),2),
                pow(0.010249296744165842 + 0.01033948520686114*l2one*pow(0.18120460196719873,2),2) + 2*pow(0.010727105951260424 + 0.010821498885630722*l2one*pow(0.18120460196719873,2),2) + 0.0013320851299382112*pow(-3.322207246598387 - 0.49654737732514054*chi2one + 2.215803700462627*chi3pone - 0.12413684433128513*etaone - 2.351292004199026*xipone - 3.0403968018625096*l2one*pow(0.18120460196719873,2) - 2.3719821468009723*l2pone*pow(0.18120460196719873,2) + 1.1859910734004862*l5one*pow(0.18120460196719873,2) - 2.3719821468009723*l2one*xipone*pow(0.18120460196719873,2),2) + 0.0012160599192045114*pow(-2.8150006193039023 + 2.215803700462627*chi3pone - 2.351292004199026*xipone - 2.464318307667692*l2one*pow(0.18120460196719873,2) - 2.3719821468009723*l2pone*pow(0.18120460196719873,2) + 2.624206492727772*l5one*pow(0.18120460196719873,2) - 2.3719821468009723*l2one*xipone*pow(0.18120460196719873,2),2) + 0.0013320851299382112*pow(-3.322207246598387 - 1.7192563231374867*chi2one + 2.215803700462627*chi3pone - 0.42981408078437167*etaone - 2.351292004199026*xipone - 3.0403968018625096*l2one*pow(0.18120460196719873,2) - 2.3719821468009723*l2pone*pow(0.18120460196719873,2) - 2.3719821468009723*l3one*pow(0.18120460196719873,2) - 1.1859910734004862*l5one*pow(0.18120460196719873,2) + 2.3719821468009723*l6one*pow(0.18120460196719873,2) - 2.3719821468009723*l2one*xipone*pow(0.18120460196719873,2),2),
            };

            double result = 0.0;
            for (const auto & b : weak_bounds)
            {
                result += b;
            }

            if (result < 0.0)
            {
                throw InternalError("Contribution to 1^+ unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= result) && (result < 1.0))
            {
                return 0.0;
            }
            else
            {
                static const double sigma = 0.0093549; // same relative uncertainty as for 1^-, cf. [BG2016], eq. (2.8), p.5
                return -pow((result - 1.0) / sigma, 2) / 2.0;
            }
        }

        double bound_1m_z1() const
        {
            const std::array<double, 4> weak_bounds
            {
                pow(0.015004330244600596 - 0.0017269654080806617*etaone + 0.014033670261852926*l1one*pow(0.18120460196719873,2) - 0.006700418983830457*l4one*pow(0.18120460196719873,2),2) + 0.001717348589813567*pow(2.5111802828589873 - 2.530765430787031*chi2one + 7.592296292361093*chi3pone - 0.2883820366649611*etaone - 0.3333837264088947*etapone + 2.8965256083351667*xipone - 0.3333837264088947*etaone*xipone + 2.3434507680704035*l1one*pow(0.18120460196719873,2) + 2.709143602528773*l1pone*pow(0.18120460196719873,2) - 1.118887769276815*l4one*pow(0.18120460196719873,2) - 1.2934889366504103*l4pone*pow(0.18120460196719873,2) + 2.709143602528773*l1one*xipone*pow(0.18120460196719873,2) - 1.2934889366504103*l4one*xipone*pow(0.18120460196719873,2),2),
                0.00014635090217060402*pow(1.2524033813421684 - 0.10466929061730514*etaone + l2one*pow(0.18120460196719873,2) - l5one*pow(0.18120460196719873,2),2) + 0.02754117923723504*pow((0.5793766430333188 + 0.5831716767309848*xipone)*(1.2524033813421684 - 0.10466929061730514*etaone + l2one*pow(0.18120460196719873,2) - l5one*pow(0.18120460196719873,2)) + 0.0728964595913731*(0.03218802151990908 - 1.6747086498768822*chi2one - 0.7744213133197135*chi3pone - 0.8373543249384411*etapone + 8*l2pone*pow(0.18120460196719873,2) - 8*l5pone*pow(0.18120460196719873,2)),2),
                0.00021525910443458962*pow(1.2524033813421684 - 0.36240920393439746*etaone + l1one*pow(0.18120460196719873,2) - l4one*pow(0.18120460196719873,2),2) + 0.024885196861219643*pow((0.6843887459111011 + 0.744046728013827*xipone)*(1.2524033813421684 - 0.36240920393439746*etaone + l1one*pow(0.18120460196719873,2) - l4one*pow(0.18120460196719873,2)) + 0.09300584100172837*(0.03218802151990908 - 5.798547262950359*chi2one + 15.720933138974194*chi3pone - 2.8992736314751797*etapone + 8*l1pone*pow(0.18120460196719873,2) - 8*l4pone*pow(0.18120460196719873,2)),2),
                0.00011131839343464602*pow(1.2524033813421684 + 0.10466929061730514*etaone + l2one*pow(0.18120460196719873,2) - l5one*pow(0.18120460196719873,2),2) + pow(0.012605143100321134 + 0.011827576145421828*l2one*pow(0.18120460196719873,2) - 0.0053453871741067985*l5one*pow(0.18120460196719873,2),2) + 0.00011131839343464602*pow(1.2524033813421684 + 0.36240920393439746*etaone + l2one*pow(0.18120460196719873,2) + l5one*pow(0.18120460196719873,2) - 2*l6one*pow(0.18120460196719873,2),2) + pow(0.008913182164063897 + 0.0009741948907552025*etaone + 0.008363359297428022*l2one*pow(0.18120460196719873,2) + 0.0037797595188785132*l5one*pow(0.18120460196719873,2) - 0.0075595190377570265*l6one*pow(0.18120460196719873,2),2) + 0.0011764410635668663*pow(3.0667433429973205 - 2.5770362021716835*chi3pone + 2.9400362135836158*xipone + 2.872138398363049*l2one*pow(0.18120460196719873,2) + 2.758675717499152*l2pone*pow(0.18120460196719873,2) - 1.298042098237697*l5one*pow(0.18120460196719873,2) - 1.2467634633278375*l5pone*pow(0.18120460196719873,2) + 2.758675717499152*l2one*xipone*pow(0.18120460196719873,2) - 1.2467634633278375*l5one*xipone*pow(0.18120460196719873,2),2) + 0.0005882205317834331*pow(3.0667433429973214 + 2.577036202171683*chi2one - 2.577036202171683*chi3pone + 0.33455725788172064*etaone + 0.3213407069650347*etapone + 2.9400362135836153*xipone + 0.3213407069650347*etaone*xipone + 2.8721383983630497*l2one*pow(0.18120460196719873,2) + 2.7586757174991514*l2pone*pow(0.18120460196719873,2) + 2.7586757174991514*l3one*pow(0.18120460196719873,2) + 1.298042098237697*l5one*pow(0.18120460196719873,2) + 1.2467634633278375*l5pone*pow(0.18120460196719873,2) - 3.842847659803231*l6one*pow(0.18120460196719873,2) - 2.493526926655675*l6pone*pow(0.18120460196719873,2) + 2.7586757174991514*l2one*xipone*pow(0.18120460196719873,2) + 1.2467634633278375*l5one*xipone*pow(0.18120460196719873,2) - 2.493526926655675*l6one*xipone*pow(0.18120460196719873,2),2) + 0.027261049494251476*pow((0.5322392241805594 + 0.5112133261002769*xipone)*(1.2524033813421684 + 0.10466929061730514*etaone + l2one*pow(0.18120460196719873,2) - l5one*pow(0.18120460196719873,2)) + 0.06390166576253462*(0.03218802151990908 + 1.6747086498768822*chi2one - 7.473255912827241*chi3pone + 0.8373543249384411*etapone + 8*l2pone*pow(0.18120460196719873,2) - 8*l5pone*pow(0.18120460196719873,2)),2) + 0.027261049494251476*pow((0.5322392241805594 + 0.5112133261002769*xipone)*(1.2524033813421684 + 0.36240920393439746*etaone + l2one*pow(0.18120460196719873,2) + l5one*pow(0.18120460196719873,2) - 2*l6one*pow(0.18120460196719873,2)) + 0.06390166576253462*(0.03218802151990908 + 5.798547262950359*chi2one - 7.473255912827241*chi3pone + 2.8992736314751797*etapone + 8*l2pone*pow(0.18120460196719873,2) + 8*l3one*pow(0.18120460196719873,2) + 8*l5pone*pow(0.18120460196719873,2) - 8*l6one*pow(0.18120460196719873,2) - 16*l6pone*pow(0.18120460196719873,2)),2),
            };

            double result = 0.0;
            for (const auto & b : weak_bounds)
            {
                result += b;
            }

            if (result < 0.0)
            {
                throw InternalError("Contribution to 1^- unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= result) && (result < 1.0))
            {
                return 0.0;
            }
            else
            {
                static const double sigma = 0.0093549; // cf. [BG2016], eq. (2.8), p.5
                return -pow((result - 1.0) / sigma, 2) / 2.0;
            }
        }
        // }}}

        // bounds for z^2
        // {{{
        double bound_0p_z2() const
        {
            const std::array<double, 3> weak_bounds
            {
                pow(0.07787308942370388 + 0.0763250029248492*l1one*pow(0.18120460196719873,2),2) + 0.02239357695917938*pow(3.8071054412766747 - 3.811666973787887*chi2one + 11.43500092136366*chi3pone - 1.1013264558178975*etaone + 4.1630886290205975*xipone + 3.0639132699149516*l1one*pow(0.18120460196719873,2) + 4.080328058612812*l1pone*pow(0.18120460196719873,2) - 4.273014767654388*l4one*pow(0.18120460196719873,2) + 4.080328058612812*l1one*xipone*pow(0.18120460196719873,2),2) + 9.512806173691777*pow(0.46393115245451544 - 1.4808207700483622*chi2one - 1.479491066171462*chi2pone + 4.4424623101450855*chi3pone + 2.2192365992571927*chi3ppone - 0.2141230516283506*etaone - 0.4274777055618846*etapone + 1.881694189509151*xipone - 1.479491066171462*chi2one*xipone + 4.4384731985143855*chi3pone*xipone - 0.4274777055618846*etaone*xipone + 0.8079473464854077*xippone + 0.34293872540126374*l1one*pow(0.18120460196719873,2) + 1.1892519313822978*l1pone*pow(0.18120460196719873,2) - 0.8307717996510665*l4one*pow(0.18120460196719873,2) - 1.6585623082598282*l4pone*pow(0.18120460196719873,2) + 1.585194766320423*l1one*xipone*pow(0.18120460196719873,2) + 1.5837713397525006*l1pone*xipone*pow(0.18120460196719873,2) - 1.6585623082598282*l4one*xipone*pow(0.18120460196719873,2) + 0.7918856698762503*l1one*xippone*pow(0.18120460196719873,2),2),
                pow(0.06378188590605295 + 0.06251392701070964*l2one*pow(0.18120460196719873,2),2) + 0.015022540138963305*pow(4.4528226929609165 - 3.811666973787886*chi3pone + 4.163088629020597*xipone + 3.661246265806994*l2one*pow(0.18120460196719873,2) + 4.080328058612811*l2pone*pow(0.18120460196719873,2) - 4.514209096519603*l5one*pow(0.18120460196719873,2) + 4.080328058612811*l2one*xipone*pow(0.18120460196719873,2),2) + 7.459264833117849*pow(0.6094293123621641 - 1.570010746037894*chi3pone - 0.6842244888188334*chi3ppone + 1.972287796397628*xipone - 1.3684489776376667*chi3pone*xipone + 0.7473074664412248*xippone + 0.450786848478118*l2one*pow(0.18120460196719873,2) + 1.3144455546031966*l2pone*pow(0.18120460196719873,2) - 1.0490487252449183*l5one*pow(0.18120460196719873,2) - 1.6206727569475938*l5pone*pow(0.18120460196719873,2) + 1.6806711980443194*l2one*xipone*pow(0.18120460196719873,2) + 1.464902573764492*l2pone*xipone*pow(0.18120460196719873,2) - 1.6206727569475938*l5one*xipone*pow(0.18120460196719873,2) + 0.732451286882246*l2one*xippone*pow(0.18120460196719873,2),2),
                pow(0.04510060404103671 + 0.04420402170787365*l2one*pow(0.18120460196719873,2),2) + 0.007511270069481652*pow(4.4528226929609165 + 3.811666973787886*chi2one - 3.811666973787886*chi3pone + 1.163491861232192*etaone + 4.163088629020597*xipone + 3.661246265806994*l2one*pow(0.18120460196719873,2) + 4.080328058612811*l2pone*pow(0.18120460196719873,2) + 4.080328058612811*l3one*pow(0.18120460196719873,2) + 4.514209096519603*l5one*pow(0.18120460196719873,2) - 9.028418193039206*l6one*pow(0.18120460196719873,2) + 4.080328058612811*l2one*xipone*pow(0.18120460196719873,2),2) + 3.7296324165589216*pow(0.6094293123621642 + 1.5700107460378945*chi2one + 1.3684489776376672*chi2pone - 1.5700107460378945*chi3pone - 0.6842244888188336*chi3ppone + 0.2703817275100315*etaone + 0.41771205589104604*etapone + 1.9722877963976284*xipone + 1.3684489776376672*chi2one*xipone - 1.3684489776376672*chi3pone*xipone + 0.41771205589104604*etaone*xipone + 0.747307466441225*xippone + 0.450786848478118*l2one*pow(0.18120460196719873,2) + 1.3144455546031966*l2pone*pow(0.18120460196719873,2) + 1.6806711980443196*l3one*pow(0.18120460196719873,2) + 1.4649025737644925*l3pone*pow(0.18120460196719873,2) + 1.0490487252449183*l5one*pow(0.18120460196719873,2) + 1.6206727569475943*l5pone*pow(0.18120460196719873,2) - 3.718770207437431*l6one*pow(0.18120460196719873,2) - 3.2413455138951885*l6pone*pow(0.18120460196719873,2) + 1.6806711980443196*l2one*xipone*pow(0.18120460196719873,2) + 1.4649025737644925*l2pone*xipone*pow(0.18120460196719873,2) + 1.4649025737644925*l3one*xipone*pow(0.18120460196719873,2) + 1.6206727569475943*l5one*xipone*pow(0.18120460196719873,2) - 3.2413455138951885*l6one*xipone*pow(0.18120460196719873,2) + 0.7324512868822463*l2one*xippone*pow(0.18120460196719873,2),2),
            };

            double result = 0.0;
            for (const auto & b : weak_bounds)
            {
                result += b;
            }

            if (result < 0.0)
            {
                throw InternalError("Contribution to 0^+ unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= result) && (result < 1.0))
            {
                return 0.0;
            }
            else
            {
                // add an r-fit like penalty
                static const double sigma = 0.0130561; // cf. [BG2016], eq. (2.8), p.5
                return -pow((result - 1.0) / sigma, 2) / 2.0;
            }
        }

        double bound_0m_z2() const
        {
            const std::array<double, 3> weak_bounds
            {
                pow(0.053405909191163586 + 0.010636706864722678*etaone + 0.05077378948243221*l2one*pow(0.18120460196719873,2) + 0.022772846510374423*l5one*pow(0.18120460196719873,2) - 0.045545693020748845*l6one*pow(0.18120460196719873,2),2) + 0.0008556908618011906*pow(-15.035034001006752 - 7.157853274953061*chi2one + 1.3441830910079258*chi3pone - 3.018891875414771*etaone - 2.9089665725799128*etapone - 14.60564877749723*xipone - 2.9089665725799128*etaone*xipone - 14.410529734621102*l2one*pow(0.18120460196719873,2) - 13.88580678644618*l2pone*pow(0.18120460196719873,2) - 13.88580678644618*l3one*pow(0.18120460196719873,2) - 6.463350187664441*l5one*pow(0.18120460196719873,2) - 6.228003657868921*l5pone*pow(0.18120460196719873,2) + 19.154704033197802*l6one*pow(0.18120460196719873,2) + 12.456007315737843*l6pone*pow(0.18120460196719873,2) - 13.88580678644618*l2one*xipone*pow(0.18120460196719873,2) - 6.228003657868921*l5one*xipone*pow(0.18120460196719873,2) + 12.456007315737843*l6one*xipone*pow(0.18120460196719873,2),2) + 0.41972667705238187*pow(-1.8941241082941498 - 3.3296030511814187*chi2one - 2.585520224482183*chi2pone + 0.6252707270247599*chi3pone + 0.24276919585437917*chi3ppone - 0.399505882595457*etaone - 1.3531576582799327*etapone - 0.5253804190193674*etappone - 6.749814882475735*xipone - 2.585520224482183*chi2one*xipone + 0.48553839170875834*chi3pone*xipone - 1.3531576582799327*etaone*xipone - 1.0507608380387348*etapone*xipone - 2.637885889477819*xippone - 0.5253804190193674*etaone*xippone - 1.9070213965536285*l2one*pow(0.18120460196719873,2) - 5.205291956003253*l2pone*pow(0.18120460196719873,2) - 6.4592305637293235*l3one*pow(0.18120460196719873,2) - 5.015754430904282*l3pone*pow(0.18120460196719873,2) - 0.8553292160858292*l5one*pow(0.18120460196719873,2) - 2.3346556552916637*l5pone*pow(0.18120460196719873,2) + 4.607725357889744*l6one*pow(0.18120460196719873,2) + 6.918956392289015*l6pone*pow(0.18120460196719873,2) - 6.4592305637293235*l2one*xipone*pow(0.18120460196719873,2) - 5.015754430904282*l2pone*xipone*pow(0.18120460196719873,2) - 5.015754430904282*l3one*xipone*pow(0.18120460196719873,2) - 2.897066925718086*l5one*xipone*pow(0.18120460196719873,2) - 2.249645081705688*l5pone*xipone*pow(0.18120460196719873,2) + 8.04377893314186*l6one*xipone*pow(0.18120460196719873,2) + 4.499290163411376*l6pone*xipone*pow(0.18120460196719873,2) - 2.507877215452141*l2one*xippone*pow(0.18120460196719873,2) - 1.124822540852844*l5one*xippone*pow(0.18120460196719873,2) + 2.249645081705688*l6one*xippone*pow(0.18120460196719873,2),2),
                pow(0.07480998702772482 - 0.015876875548846626*etaone + 0.07070378019068203*l1one*pow(0.18120460196719873,2) - 0.03399187873996447*l4one*pow(0.18120460196719873,2),2) + 0.0007031761385459602*pow(-21.02575811558329 + 10.99543036092901*chi2one - 41.916874484470014*chi3pone + 4.500537596414145*etaone + 4.789861601196494*etapone - 22.56926957370638*xipone + 4.789861601196494*etaone*xipone - 20.04204290559436*l1one*pow(0.18120460196719873,2) - 21.330476563405902*l1pone*pow(0.18120460196719873,2) + 9.63550591369812*l4one*pow(0.18120460196719873,2) + 10.254939281231*l4pone*pow(0.18120460196719873,2) - 21.330476563405902*l1one*xipone*pow(0.18120460196719873,2) + 10.254939281231*l4one*xipone*pow(0.18120460196719873,2),2) + 0.5103614048190293*pow(-1.9037175574697134 + 3.8841418133342525*chi2one + 3.265091600692897*chi2pone - 14.807158931035932*chi3pone - 6.223605184789598*chi3ppone + 0.43230733217841577*etaone + 1.6920212410603224*etapone + 0.711174386503355*etappone - 7.919081517765644*xipone + 3.265091600692897*chi2one*xipone - 12.447210369579196*chi3pone*xipone + 1.6920212410603224*etaone*xipone + 1.42234877300671*etapone*xipone - 3.350970816964743*xippone + 0.711174386503355*etaone*xippone - 1.9251749183978024*l1one*pow(0.18120460196719873,2) - 5.951482006953827*l1pone*pow(0.18120460196719873,2) + 0.925556062248895*l4one*pow(0.18120460196719873,2) + 2.86126221480469*l4pone*pow(0.18120460196719873,2) - 7.535002560033591*l1one*xipone*pow(0.18120460196719873,2) - 6.3340822123190526*l1pone*xipone*pow(0.18120460196719873,2) + 3.622562932777087*l4one*xipone*pow(0.18120460196719873,2) + 3.0452028718895883*l4pone*xipone*pow(0.18120460196719873,2) - 3.1670411061595263*l1one*xippone*pow(0.18120460196719873,2) + 1.5226014359447941*l4one*xippone*pow(0.18120460196719873,2),2),
                pow(0.06252207819688696 + 0.059403241258809435*l2one*pow(0.18120460196719873,2) - 0.02684686363639445*l5one*pow(0.18120460196719873,2),2) + 0.0010860026569941262*pow(-16.542981584936765 + 13.471146304272825*chi3pone - 15.177769336722374*xipone - 15.838121246690301*l2one*pow(0.18120460196719873,2) - 14.420644989449043*l2pone*pow(0.18120460196719873,2) + 7.1579239172158635*l5one*pow(0.18120460196719873,2) + 6.517305813227482*l5pone*pow(0.18120460196719873,2) - 14.420644989449043*l2one*xipone*pow(0.18120460196719873,2) + 6.517305813227482*l5one*xipone*pow(0.18120460196719873,2),2) + 0.5761118796300714*pow(-2.152982304075243 + 6.308723146428653*chi3pone + 2.339518979349188*chi3ppone - 7.063955535061434*xipone + 4.679037958698376*chi3pone*xipone - 2.635906301172259*xippone - 2.153940696644513*l2one*pow(0.18120460196719873,2) - 5.501177764228245*l2pone*pow(0.18120460196719873,2) + 0.9734578608557002*l5one*pow(0.18120460196719873,2) + 2.486217353567367*l5pone*pow(0.18120460196719873,2) - 6.753386443625182*l2one*xipone*pow(0.18120460196719873,2) - 5.008834717587749*l2pone*xipone*pow(0.18120460196719873,2) + 3.0521439755443054*l5one*xipone*pow(0.18120460196719873,2) + 2.2637064879077555*l5pone*xipone*pow(0.18120460196719873,2) - 2.5044173587938743*l2one*xippone*pow(0.18120460196719873,2) + 1.1318532439538778*l5one*xippone*pow(0.18120460196719873,2),2),
            };

            double result = 0.0;
            for (const auto & b : weak_bounds)
            {
                result += b;
            }

            if (result < 0.0)
            {
                throw InternalError("Contribution to 0^- unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= result) && (result < 1.0))
            {
                return 0.0;
            }
            else
            {
                static const double sigma = 0.0130561; // using the same relative uncertainty as for 0^+, cf. [BG2016], eq. (2.8), p.5
                return -pow((result - 1.0) / sigma, 2) / 2.0;
            }
        }

        double bound_1p_z2() const
        {
            const std::array<double, 5> weak_bounds
            {
                pow(0.00895448576377894 + 0.009033280565551251*l2one*pow(0.18120460196719873,2),2) + pow(0.0133673793613974 + 0.01348500532393532*l2one*pow(0.18120460196719873,2),2) + 0.000928214055953212*pow(-2.6585020549037592 - 1.2227089458123461*chi2one + 0.22961419116206505*chi3pone - 1.2350756258163746*etaone - 2.351292004199026*xipone - 2.3035098010173103*l2one*pow(0.18120460196719873,2) - 2.3719821468009723*l2pone*pow(0.18120460196719873,2) - 2.3719821468009723*l3one*pow(0.18120460196719873,2) - 2.6442571007294786*l5one*pow(0.18120460196719873,2) + 5.288514201458957*l6one*pow(0.18120460196719873,2) - 2.3719821468009723*l2one*xipone*pow(0.18120460196719873,2),2) + pow(0.1438601613520823 - 0.02258345505963029*chi2one - 0.010443075533085308*chi3pone - 0.005645863764907572*etaone + 0.1069390348911792*xipone + 0.13097946755643228*l2one*pow(0.18120460196719873,2) + 0.10788004259148257*l2pone*pow(0.18120460196719873,2) - 0.05394002129574128*l5one*pow(0.18120460196719873,2) + 0.10788004259148257*l2one*xipone*pow(0.18120460196719873,2),2) + 0.011376860803770166*pow(2.370627524838556 + 3.4118414572935523*chi2one + 2.7939969310717263*chi2pone - 0.6407143901848045*chi3pone - 0.2623442592918539*chi3ppone + 2.0352214461909903*etaone + 2.822255876994335*etapone + 7.4181368138836*xipone + 2.7939969310717263*chi2one*xipone - 0.5246885185837078*chi3pone*xipone + 2.822255876994335*etaone*xipone + 2.68645398657033*xippone + 1.7734568392678665*l2one*pow(0.18120460196719873,2) + 5.263721457815983*l2pone*pow(0.18120460196719873,2) + 6.61876814767147*l3one*pow(0.18120460196719873,2) + 5.420186759421943*l3pone*pow(0.18120460196719873,2) + 4.357343508491814*l5one*pow(0.18120460196719873,2) + 6.042358853843418*l5pone*pow(0.18120460196719873,2) - 14.757045870827046*l6one*pow(0.18120460196719873,2) - 12.084717707686837*l6pone*pow(0.18120460196719873,2) + 6.61876814767147*l2one*xipone*pow(0.18120460196719873,2) + 5.420186759421943*l2pone*xipone*pow(0.18120460196719873,2) + 5.420186759421943*l3one*xipone*pow(0.18120460196719873,2) + 6.042358853843418*l5one*xipone*pow(0.18120460196719873,2) - 12.084717707686837*l6one*xipone*pow(0.18120460196719873,2) + 2.7100933797109716*l2one*xippone*pow(0.18120460196719873,2),2) + 0.04073136914413835*pow((2.5344087280025343 + 5.191928503443212*xipone + 0.06681694545011244*(16*xipone + 32*xippone))*(0.991277277263722 + l2one*pow(0.18120460196719873,2)) + (0.6489910629304015 + 0.5345355636008995*xipone)*(1.0399099620969015 - 1.6747086498768822*chi2one - 0.7744213133197135*chi3pone - 0.41867716246922054*etaone + 8*l2pone*pow(0.18120460196719873,2) - 4*l5one*pow(0.18120460196719873,2)) + 0.06681694545011244*(-2.141739969395875 - 3.3494172997537643*chi2one - 13.397669199015057*chi2pone - 1.548842626639427*chi3pone - 3.097685253278854*chi3ppone + 0.8373543249384411*etaone - 3.3494172997537643*etapone + 8*l5one*pow(0.18120460196719873,2) - 32*l5pone*pow(0.18120460196719873,2)),2),
                pow(0.01563097879202461 + 0.01576852324827991*l1one*pow(0.18120460196719873,2),2) + pow(0.021574370557675353 + 0.021764213759875837*l1one*pow(0.18120460196719873,2),2) + 0.00282838931341691*pow(-2.294683910554761 + 1.2227089458123461*chi2one - 4.661221592087319*chi3pone + 1.1522288822404556*etaone - 2.351292004199026*xipone - 1.962434781465052*l1one*pow(0.18120460196719873,2) - 2.3719821468009723*l1pone*pow(0.18120460196719873,2) + 2.46688489339753*l4one*pow(0.18120460196719873,2) - 2.3719821468009723*l1one*xipone*pow(0.18120460196719873,2),2) + pow(0.20716060223145838 - 0.1262008221275946*chi2one + 0.34215374934135023*chi3pone - 0.03155020553189865*etaone + 0.17259496446140282*xipone + 0.1861515276883168*l1one*pow(0.18120460196719873,2) + 0.1741137100790067*l1pone*pow(0.18120460196719873,2) - 0.08705685503950335*l4one*pow(0.18120460196719873,2) + 0.1741137100790067*l1one*xipone*pow(0.18120460196719873,2),2) + 0.03555160950637785*pow(1.6457602104526785 - 2.972389127976415*chi2one - 2.7590084905806913*chi2pone + 11.331367477813188*chi3pone + 5.258957985501109*chi3ppone - 1.501066933550325*etaone - 2.599972201137074*etapone + 6.5042961127520424*xipone - 2.7590084905806913*chi2one*xipone + 10.517915971002218*chi3pone*xipone - 2.599972201137074*etaone*xipone + 2.6528122762321*xippone + 1.2076947902829935*l1one*pow(0.18120460196719873,2) + 4.428179120482126*l1pone*pow(0.18120460196719873,2) - 3.2137359160391976*l4one*pow(0.18120460196719873,2) - 5.566456669413783*l4pone*pow(0.18120460196719873,2) + 5.766256940421067*l1one*xipone*pow(0.18120460196719873,2) + 5.352311279755763*l1pone*xipone*pow(0.18120460196719873,2) - 5.566456669413783*l4one*xipone*pow(0.18120460196719873,2) + 2.6761556398778814*l1one*xippone*pow(0.18120460196719873,2),2) + 0.037557428400381614*pow((3.2694898363708664 + 7.684377075048818*xipone + 0.11230395698986552*(16*xipone + 32*xippone))*(0.991277277263722 + l1one*pow(0.18120460196719873,2)) + (0.9605471343811023 + 0.8984316559189242*xipone)*(1.0399099620969015 - 5.798547262950359*chi2one + 15.720933138974194*chi3pone - 1.4496368157375898*etaone + 8*l1pone*pow(0.18120460196719873,2) - 4*l4one*pow(0.18120460196719873,2)) + 0.11230395698986552*(-2.141739969395875 - 11.597094525900719*chi2one - 46.388378103602875*chi2pone + 31.44186627794839*chi3pone + 62.88373255589678*chi3ppone + 2.8992736314751797*etaone - 11.597094525900719*etapone + 8*l4one*pow(0.18120460196719873,2) - 32*l4pone*pow(0.18120460196719873,2)),2),
                pow(0.010727105951260424 + 0.010821498885630722*l2one*pow(0.18120460196719873,2),2) + 0.0013320851299382112*pow(-3.322207246598387 - 0.49654737732514054*chi2one + 2.215803700462627*chi3pone - 0.12413684433128513*etaone - 2.351292004199026*xipone - 3.0403968018625096*l2one*pow(0.18120460196719873,2) - 2.3719821468009723*l2pone*pow(0.18120460196719873,2) + 1.1859910734004862*l5one*pow(0.18120460196719873,2) - 2.3719821468009723*l2one*xipone*pow(0.18120460196719873,2),2) + 0.00002269186265686789*pow(-114.85325433150433 - 46.62109573713025*chi2one - 30.435580980091824*chi2pone + 208.04298073316056*chi3pone + 67.90819571407988*chi3ppone - 7.850826311771085*etaone - 7.608895245022956*etapone - 239.66301640332577*xipone - 30.435580980091824*chi2one*xipone + 135.81639142815976*chi3pone*xipone - 7.608895245022956*etaone*xipone - 72.06053386803235*xippone - 96.3343409076623*l2one*pow(0.18120460196719873,2) - 186.35934313697047*l2pone*pow(0.18120460196719873,2) + 75.00601432826656*l5one*pow(0.18120460196719873,2) + 72.69462896087467*l5pone*pow(0.18120460196719873,2) - 222.70665761740776*l2one*xipone*pow(0.18120460196719873,2) - 145.38925792174933*l2pone*xipone*pow(0.18120460196719873,2) + 72.69462896087467*l5one*xipone*pow(0.18120460196719873,2) - 72.69462896087467*l2one*xippone*pow(0.18120460196719873,2),2),
                pow(0.010727105951260424 + 0.010821498885630722*l2one*pow(0.18120460196719873,2),2) + 0.0013320851299382112*pow(-3.322207246598387 - 1.7192563231374867*chi2one + 2.215803700462627*chi3pone - 0.42981408078437167*etaone - 2.351292004199026*xipone - 3.0403968018625096*l2one*pow(0.18120460196719873,2) - 2.3719821468009723*l2pone*pow(0.18120460196719873,2) - 2.3719821468009723*l3one*pow(0.18120460196719873,2) - 1.1859910734004862*l5one*pow(0.18120460196719873,2) + 2.3719821468009723*l6one*pow(0.18120460196719873,2) - 2.3719821468009723*l2one*xipone*pow(0.18120460196719873,2),2) + 0.00002269186265686789*pow(-114.85325433150433 - 161.42188499603031*chi2one - 105.38081044806793*chi2pone + 208.04298073316056*chi3pone + 67.90819571407988*chi3ppone - 27.18286994299909*etaone - 26.345202612016983*etapone - 239.66301640332577*xipone - 105.38081044806793*chi2one*xipone + 135.81639142815976*chi3pone*xipone - 26.345202612016983*etaone*xipone - 72.06053386803235*xippone - 96.3343409076623*l2one*pow(0.18120460196719873,2) - 186.35934313697047*l2pone*pow(0.18120460196719873,2) - 222.70665761740776*l3one*pow(0.18120460196719873,2) - 145.38925792174933*l3pone*pow(0.18120460196719873,2) - 75.00601432826656*l5one*pow(0.18120460196719873,2) - 72.69462896087467*l5pone*pow(0.18120460196719873,2) + 222.70665761740776*l6one*pow(0.18120460196719873,2) + 145.38925792174933*l6pone*pow(0.18120460196719873,2) - 222.70665761740776*l2one*xipone*pow(0.18120460196719873,2) - 145.38925792174933*l2pone*xipone*pow(0.18120460196719873,2) - 145.38925792174933*l3one*xipone*pow(0.18120460196719873,2) - 72.69462896087467*l5one*xipone*pow(0.18120460196719873,2) + 145.38925792174933*l6one*xipone*pow(0.18120460196719873,2) - 72.69462896087467*l2one*xippone*pow(0.18120460196719873,2),2),
                pow(0.010249296744165842 + 0.01033948520686114*l2one*pow(0.18120460196719873,2),2) + 0.0012160599192045114*pow(-2.8150006193039023 + 2.215803700462627*chi3pone - 2.351292004199026*xipone - 2.464318307667692*l2one*pow(0.18120460196719873,2) - 2.3719821468009723*l2pone*pow(0.18120460196719873,2) + 2.624206492727772*l5one*pow(0.18120460196719873,2) - 2.3719821468009723*l2one*xipone*pow(0.18120460196719873,2),2) + 0.009763964568486045*pow(3.3500180903719508 - 8.063326901525953*chi3pone - 3.127920198554829*chi3ppone + 9.607133278098072*xipone - 6.255840397109658*chi3pone*xipone + 3.3191810949223775*xippone + 2.5501773438200552*l2one*pow(0.18120460196719873,2) + 6.957467404366977*l2pone*pow(0.18120460196719873,2) - 5.845069358574001*l5one*pow(0.18120460196719873,2) - 7.408876961499929*l5pone*pow(0.18120460196719873,2) + 8.631661482579128*l2one*xipone*pow(0.18120460196719873,2) + 6.696776312848607*l2pone*xipone*pow(0.18120460196719873,2) - 7.408876961499929*l5one*xipone*pow(0.18120460196719873,2) + 3.3483881564243037*l2one*xippone*pow(0.18120460196719873,2),2)
            };

            double result = 0.0;
            for (const auto & b : weak_bounds)
            {
                result += b;
            }

            if (result < 0.0)
            {
                throw InternalError("Contribution to 1^+ unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= result) && (result < 1.0))
            {
                return 0.0;
            }
            else
            {
                static const double sigma = 0.0093549; // same relative uncertainty as for 1^-, cf. [BG2016], eq. (2.8), p.5
                return -pow((result - 1.0) / sigma, 2) / 2.0;
            }
        }

        double bound_1m_z2() const
        {
            const std::array<double, 4> weak_bounds
            {
                pow(0.015004330244600596 - 0.0017269654080806617*etaone + 0.014033670261852926*l1one*pow(0.18120460196719873,2) - 0.006700418983830457*l4one*pow(0.18120460196719873,2),2) + 0.001717348589813567*pow(2.5111802828589873 - 2.530765430787031*chi2one + 7.592296292361093*chi3pone - 0.2883820366649611*etaone - 0.3333837264088947*etapone + 2.8965256083351667*xipone - 0.3333837264088947*etaone*xipone + 2.3434507680704035*l1one*pow(0.18120460196719873,2) + 2.709143602528773*l1pone*pow(0.18120460196719873,2) - 1.118887769276815*l4one*pow(0.18120460196719873,2) - 1.2934889366504103*l4pone*pow(0.18120460196719873,2) + 2.709143602528773*l1one*xipone*pow(0.18120460196719873,2) - 1.2934889366504103*l4one*xipone*pow(0.18120460196719873,2),2) + 0.21164551686463315*pow(0.48728462933327854 - 2.033515158105124*chi2one - 1.823755323760786*chi2pone + 6.100545474315371*chi3pone + 2.735632985641179*chi3ppone - 0.056354587120502446*etaone - 0.26787976983992107*etapone - 0.12012380493603245*etappone + 2.331475227874015*xipone - 1.823755323760786*chi2one*xipone + 5.471265971282358*chi3pone*xipone - 0.26787976983992107*etaone*xipone - 0.2402476098720649*etapone*xipone + 1.0436672506957534*xippone - 0.12012380493603245*etaone*xippone + 0.4579487751702877*l1one*pow(0.18120460196719873,2) + 1.6887700306980196*l1pone*pow(0.18120460196719873,2) - 0.2186490497153637*l4one*pow(0.18120460196719873,2) - 0.8063084397651313*l4pone*pow(0.18120460196719873,2) + 2.176845200352108*l1one*xipone*pow(0.18120460196719873,2) + 1.952300678616352*l1pone*xipone*pow(0.18120460196719873,2) - 1.0393414290876786*l4one*xipone*pow(0.18120460196719873,2) - 0.932131957290189*l4pone*xipone*pow(0.18120460196719873,2) + 0.976150339308176*l1one*xippone*pow(0.18120460196719873,2) - 0.4660659786450945*l4one*xippone*pow(0.18120460196719873,2),2),
                0.00014635090217060402*pow(1.2524033813421684 - 0.10466929061730514*etaone + l2one*pow(0.18120460196719873,2) - l5one*pow(0.18120460196719873,2),2) + 0.02754117923723504*pow((0.5793766430333188 + 0.5831716767309848*xipone)*(1.2524033813421684 - 0.10466929061730514*etaone + l2one*pow(0.18120460196719873,2) - l5one*pow(0.18120460196719873,2)) + 0.0728964595913731*(0.03218802151990908 - 1.6747086498768822*chi2one - 0.7744213133197135*chi3pone - 0.8373543249384411*etapone + 8*l2pone*pow(0.18120460196719873,2) - 8*l5pone*pow(0.18120460196719873,2)),2) + 0.02754117923723504*pow(0.0728964595913731*(-0.1851480494207603 - 3.3494172997537643*chi2one - 13.397669199015057*chi2pone - 1.548842626639427*chi3pone - 3.097685253278854*chi3ppone - 1.6747086498768822*etapone - 3.3494172997537643*etappone) + (1.4999940008253256 + 4.6350131442665505*xipone + 0.0728964595913731*(16*xipone + 32*xippone))*(1.2524033813421684 - 0.10466929061730514*etaone + l2one*pow(0.18120460196719873,2) - l5one*pow(0.18120460196719873,2)) + (0.5793766430333188 + 0.5831716767309848*xipone)*(0.03218802151990908 - 1.6747086498768822*chi2one - 0.7744213133197135*chi3pone - 0.8373543249384411*etapone + 8*l2pone*pow(0.18120460196719873,2) - 8*l5pone*pow(0.18120460196719873,2)),2),
                0.00021525910443458962*pow(1.2524033813421684 - 0.36240920393439746*etaone + l1one*pow(0.18120460196719873,2) - l4one*pow(0.18120460196719873,2),2) + 0.024885196861219643*pow((0.6843887459111011 + 0.744046728013827*xipone)*(1.2524033813421684 - 0.36240920393439746*etaone + l1one*pow(0.18120460196719873,2) - l4one*pow(0.18120460196719873,2)) + 0.09300584100172837*(0.03218802151990908 - 5.798547262950359*chi2one + 15.720933138974194*chi3pone - 2.8992736314751797*etapone + 8*l1pone*pow(0.18120460196719873,2) - 8*l4pone*pow(0.18120460196719873,2)),2) + 0.024885196861219643*pow(0.09300584100172837*(-0.1851480494207603 - 11.597094525900719*chi2one - 46.388378103602875*chi2pone + 31.44186627794839*chi3pone + 62.88373255589678*chi3ppone - 5.798547262950359*etapone - 11.597094525900719*etappone) + (1.6094830164164176 + 5.475109967288809*xipone + 0.09300584100172837*(16*xipone + 32*xippone))*(1.2524033813421684 - 0.36240920393439746*etaone + l1one*pow(0.18120460196719873,2) - l4one*pow(0.18120460196719873,2)) + (0.6843887459111011 + 0.744046728013827*xipone)*(0.03218802151990908 - 5.798547262950359*chi2one + 15.720933138974194*chi3pone - 2.8992736314751797*etapone + 8*l1pone*pow(0.18120460196719873,2) - 8*l4pone*pow(0.18120460196719873,2)),2),
                0.00011131839343464602*pow(1.2524033813421684 + 0.10466929061730514*etaone + l2one*pow(0.18120460196719873,2) - l5one*pow(0.18120460196719873,2),2) + pow(0.012605143100321134 + 0.011827576145421828*l2one*pow(0.18120460196719873,2) - 0.0053453871741067985*l5one*pow(0.18120460196719873,2),2) + 0.00011131839343464602*pow(1.2524033813421684 + 0.36240920393439746*etaone + l2one*pow(0.18120460196719873,2) + l5one*pow(0.18120460196719873,2) - 2*l6one*pow(0.18120460196719873,2),2) + pow(0.008913182164063897 + 0.0009741948907552025*etaone + 0.008363359297428022*l2one*pow(0.18120460196719873,2) + 0.0037797595188785132*l5one*pow(0.18120460196719873,2) - 0.0075595190377570265*l6one*pow(0.18120460196719873,2),2) + 0.0011764410635668663*pow(3.0667433429973205 - 2.5770362021716835*chi3pone + 2.9400362135836158*xipone + 2.872138398363049*l2one*pow(0.18120460196719873,2) + 2.758675717499152*l2pone*pow(0.18120460196719873,2) - 1.298042098237697*l5one*pow(0.18120460196719873,2) - 1.2467634633278375*l5pone*pow(0.18120460196719873,2) + 2.758675717499152*l2one*xipone*pow(0.18120460196719873,2) - 1.2467634633278375*l5one*xipone*pow(0.18120460196719873,2),2) + 0.0005882205317834331*pow(3.0667433429973214 + 2.577036202171683*chi2one - 2.577036202171683*chi3pone + 0.33455725788172064*etaone + 0.3213407069650347*etapone + 2.9400362135836153*xipone + 0.3213407069650347*etaone*xipone + 2.8721383983630497*l2one*pow(0.18120460196719873,2) + 2.7586757174991514*l2pone*pow(0.18120460196719873,2) + 2.7586757174991514*l3one*pow(0.18120460196719873,2) + 1.298042098237697*l5one*pow(0.18120460196719873,2) + 1.2467634633278375*l5pone*pow(0.18120460196719873,2) - 3.842847659803231*l6one*pow(0.18120460196719873,2) - 2.493526926655675*l6pone*pow(0.18120460196719873,2) + 2.7586757174991514*l2one*xipone*pow(0.18120460196719873,2) + 1.2467634633278375*l5one*xipone*pow(0.18120460196719873,2) - 2.493526926655675*l6one*xipone*pow(0.18120460196719873,2),2) + 0.18576299474007757*pow(0.6660610642698668 - 2.1182916929684246*chi3pone - 0.8203250931716706*chi3ppone + 2.4203562884826924*xipone - 1.6406501863433411*chi3pone*xipone + 0.9358756694235182*xippone + 0.6263083175108694*l2one*pow(0.18120460196719873,2) + 1.828524719407211*l2pone*pow(0.18120460196719873,2) - 0.2830554972799635*l5one*pow(0.18120460196719873,2) - 0.8263884723701301*l5pone*pow(0.18120460196719873,2) + 2.2675971144866565*l2one*xipone*pow(0.18120460196719873,2) + 1.756289580317781*l2pone*xipone*pow(0.18120460196719873,2) - 1.024824053786403*l5one*xipone*pow(0.18120460196719873,2) - 0.7937423256650912*l5pone*xipone*pow(0.18120460196719873,2) + 0.8781447901588905*l2one*xippone*pow(0.18120460196719873,2) - 0.3968711628325456*l5one*xippone*pow(0.18120460196719873,2),2) + 0.09288149737003878*pow(0.6660610642698671 + 2.1182916929684237*chi2one + 1.6406501863433411*chi2pone - 2.1182916929684246*chi3pone - 0.8203250931716706*chi3ppone + 0.07295469933286425*etaone + 0.26413806278817864*etapone + 0.10228953910651392*etappone + 2.4203562884826924*xipone + 1.6406501863433411*chi2one*xipone - 1.6406501863433411*chi3pone*xipone + 0.26413806278817864*etaone*xipone + 0.20457907821302784*etapone*xipone + 0.9358756694235182*xippone + 0.10228953910651392*etaone*xippone + 0.6263083175108697*l2one*pow(0.18120460196719873,2) + 1.8285247194072114*l2pone*pow(0.18120460196719873,2) + 2.2675971144866556*l3one*pow(0.18120460196719873,2) + 1.756289580317781*l3pone*pow(0.18120460196719873,2) + 0.2830554972799635*l5one*pow(0.18120460196719873,2) + 0.82638847237013*l5pone*pow(0.18120460196719873,2) - 1.5909350483463298*l6one*pow(0.18120460196719873,2) - 2.4465192704053513*l6pone*pow(0.18120460196719873,2) + 2.2675971144866565*l2one*xipone*pow(0.18120460196719873,2) + 1.756289580317781*l2pone*xipone*pow(0.18120460196719873,2) + 1.756289580317781*l3one*xipone*pow(0.18120460196719873,2) + 1.0248240537864028*l5one*xipone*pow(0.18120460196719873,2) + 0.7937423256650912*l5pone*xipone*pow(0.18120460196719873,2) - 2.8433904332378965*l6one*xipone*pow(0.18120460196719873,2) - 1.5874846513301824*l6pone*xipone*pow(0.18120460196719873,2) + 0.8781447901588905*l2one*xippone*pow(0.18120460196719873,2) + 0.3968711628325456*l5one*xippone*pow(0.18120460196719873,2) - 0.7937423256650912*l6one*xippone*pow(0.18120460196719873,2),2) + 0.027261049494251476*pow((0.5322392241805594 + 0.5112133261002769*xipone)*(1.2524033813421684 + 0.10466929061730514*etaone + l2one*pow(0.18120460196719873,2) - l5one*pow(0.18120460196719873,2)) + 0.06390166576253462*(0.03218802151990908 + 1.6747086498768822*chi2one - 7.473255912827241*chi3pone + 0.8373543249384411*etapone + 8*l2pone*pow(0.18120460196719873,2) - 8*l5pone*pow(0.18120460196719873,2)),2) + 0.027261049494251476*pow(0.06390166576253462*(-0.1851480494207603 + 3.3494172997537643*chi2one + 13.397669199015057*chi2pone - 14.946511825654483*chi3pone - 29.893023651308965*chi3ppone + 1.6747086498768822*etapone + 3.3494172997537643*etappone) + (1.458425361043558 + 4.257913793444475*xipone + 0.06390166576253462*(16*xipone + 32*xippone))*(1.2524033813421684 + 0.10466929061730514*etaone + l2one*pow(0.18120460196719873,2) - l5one*pow(0.18120460196719873,2)) + (0.5322392241805594 + 0.5112133261002769*xipone)*(0.03218802151990908 + 1.6747086498768822*chi2one - 7.473255912827241*chi3pone + 0.8373543249384411*etapone + 8*l2pone*pow(0.18120460196719873,2) - 8*l5pone*pow(0.18120460196719873,2)),2) + 0.027261049494251476*pow((0.5322392241805594 + 0.5112133261002769*xipone)*(1.2524033813421684 + 0.36240920393439746*etaone + l2one*pow(0.18120460196719873,2) + l5one*pow(0.18120460196719873,2) - 2*l6one*pow(0.18120460196719873,2)) + 0.06390166576253462*(0.03218802151990908 + 5.798547262950359*chi2one - 7.473255912827241*chi3pone + 2.8992736314751797*etapone + 8*l2pone*pow(0.18120460196719873,2) + 8*l3one*pow(0.18120460196719873,2) + 8*l5pone*pow(0.18120460196719873,2) - 8*l6one*pow(0.18120460196719873,2) - 16*l6pone*pow(0.18120460196719873,2)),2) + 0.027261049494251476*pow((1.458425361043558 + 4.257913793444475*xipone + 0.06390166576253462*(16*xipone + 32*xippone))*(1.2524033813421684 + 0.36240920393439746*etaone + l2one*pow(0.18120460196719873,2) + l5one*pow(0.18120460196719873,2) - 2*l6one*pow(0.18120460196719873,2)) + 0.06390166576253462*(-0.1851480494207603 + 11.597094525900719*chi2one + 46.388378103602875*chi2pone - 14.946511825654483*chi3pone - 29.893023651308965*chi3ppone + 5.798547262950359*etapone + 11.597094525900719*etappone + 16*l3one*pow(0.18120460196719873,2) + 64*l3pone*pow(0.18120460196719873,2) - 16*l6one*pow(0.18120460196719873,2) - 64*l6pone*pow(0.18120460196719873,2)) + (0.5322392241805594 + 0.5112133261002769*xipone)*(0.03218802151990908 + 5.798547262950359*chi2one - 7.473255912827241*chi3pone + 2.8992736314751797*etapone + 8*l2pone*pow(0.18120460196719873,2) + 8*l3one*pow(0.18120460196719873,2) + 8*l5pone*pow(0.18120460196719873,2) - 8*l6one*pow(0.18120460196719873,2) - 16*l6pone*pow(0.18120460196719873,2)),2),
            };

            double result = 0.0;
            for (const auto & b : weak_bounds)
            {
                result += b;
            }

            if (result < 0.0)
            {
                throw InternalError("Contribution to 1^- unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= result) && (result < 1.0))
            {
                return 0.0;
            }
            else
            {
                static const double sigma = 0.0093549; // cf. [BG2016], eq. (2.8), p.5
                return -pow((result - 1.0) / sigma, 2) / 2.0;
            }
        }
        // }}}
    };

    HQETUnitarityBounds::HQETUnitarityBounds(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<HQETUnitarityBounds>(new Implementation<HQETUnitarityBounds>(p, o, *this))
    {
    }

    HQETUnitarityBounds::~HQETUnitarityBounds() = default;

    double
    HQETUnitarityBounds::bound_0p() const
    {
        return _imp->bound_0p();
    }

    double
    HQETUnitarityBounds::bound_0m() const
    {
        return _imp->bound_0m();
    }

    double
    HQETUnitarityBounds::bound_1p() const
    {
        return _imp->bound_1p();
    }

    double
    HQETUnitarityBounds::bound_1m() const
    {
        return _imp->bound_1m();
    }
}
