/* vim: set sw=4 sts=4 tw=120 et foldmethod=syntax : */

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
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <cmath>

#include <iostream>

namespace eos
{
    template <> struct Implementation<HQETUnitarityBounds>
    {
        // option to determine if we use z^3 terms in the leading-power IW function
        SwitchOption opt_zorder_bound;
        std::function<double ()> bound_0p;
        std::function<double ()> bound_0m;

        // parameters for the leading Isgur-Wise function xi
        UsedParameter xipone, xippone, xipppone, xippppone, xipppppone;

        // parameters for the subleading Isgur-Wise function chi2
        UsedParameter chi2one, chi2pone, chi2ppone, chi2pppone;

        // parameters for the subleading Isgur-Wise function chi3
        UsedParameter chi3pone, chi3ppone, chi3pppone;

        // parameters for the subleading Isgur-Wise function eta
        UsedParameter etaone, etapone, etappone, etapppone;

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
            opt_zorder_bound(o, "z-order-bound", { "2" }, "2"),
            xipone(p["B(*)->D(*)::xi'(1)@HQET"], u),
            xippone(p["B(*)->D(*)::xi''(1)@HQET"], u),
            xipppone(p["B(*)->D(*)::xi'''(1)@HQET"], u),
            xippppone(p["B(*)->D(*)::xi''''(1)@HQET"], u),
            xipppppone(p["B(*)->D(*)::xi'''''(1)@HQET"], u),
            chi2one(p["B(*)->D(*)::chi_2(1)@HQET"], u),
            chi2pone(p["B(*)->D(*)::chi_2'(1)@HQET"], u),
            chi2ppone(p["B(*)->D(*)::chi_2''(1)@HQET"], u),
            chi2pppone(p["B(*)->D(*)::chi_2'''(1)@HQET"], u),
            chi3pone(p["B(*)->D(*)::chi_3'(1)@HQET"], u),
            chi3ppone(p["B(*)->D(*)::chi_3''(1)@HQET"], u),
            chi3pppone(p["B(*)->D(*)::chi_3'''(1)@HQET"], u),
            etaone(p["B(*)->D(*)::eta(1)@HQET"], u),
            etapone(p["B(*)->D(*)::eta'(1)@HQET"], u),
            etappone(p["B(*)->D(*)::eta''(1)@HQET"], u),
            etapppone(p["B(*)->D(*)::eta'''(1)@HQET"], u),
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
            if ("2" == opt_zorder_bound.value())
            {
                bound_0p = std::bind(&Implementation::bound_0p_z2, this);
                bound_0m = std::bind(&Implementation::bound_0m_z2, this);
            }
            else
            {
                throw InternalError("Only z-order-bound=2 is presently supported");
            }
        }

        ~Implementation() = default;

        double bound_0p_z2() const
        {
            const std::array<double, 3> weak_bounds
            {
                pow(0.09389700723930866 + 0.0933887349513945*l1one*pow(0.17862041985779245,2),2) + 0.03862298802566855*pow(3.4148045551307606 - 2.7161430681278897*chi2one + 8.148429204383671*chi3pone - 1.4222042044725078*etaone + 3.8222474761508347*xipone - 6.530221765184699*pow(0.17862041985779245,2) - 24.942995562895145*chi2one*pow(0.17862041985779245,2) + 74.82898668868546*chi3pone*pow(0.17862041985779245,2) + 13.060443530369398*etaone*pow(0.17862041985779245,2) + 2.657502314742668*l1one*pow(0.17862041985779245,2) + 3.8015573335488884*l1pone*pow(0.17862041985779245,2) - 3.9810795585543546*l4one*pow(0.17862041985779245,2) + 3.8015573335488884*l1one*xipone*pow(0.17862041985779245,2),2) + 15.23923673813489*pow(0.41633604007414476 - 1.0381877434248652*chi2one - 1.0939159527318707*chi2pone + 3.114563230274595*chi3pone + 1.6408739290978063*chi3ppone - 0.25721363807054265*etaone - 0.5727871574847256*etapone + 1.760147912229388*xipone - 1.0939159527318707*chi2one*xipone + 3.2817478581956125*chi3pone*xipone - 0.5727871574847256*etaone*xipone + 0.7696975793569752*xippone - 1.1810273745137605*pow(0.17862041985779245,2) - 9.533927936847265*chi2one*pow(0.17862041985779245,2) - 10.045693496542476*chi2pone*pow(0.17862041985779245,2) + 28.601783810541797*chi3pone*pow(0.17862041985779245,2) + 15.068540244813715*chi3ppone*pow(0.17862041985779245,2) + 2.362054749027521*etaone*pow(0.17862041985779245,2) + 5.260042335498996*etapone*pow(0.17862041985779245,2) + 0.28486263672004786*l1one*pow(0.17862041985779245,2) + 1.0702986196241107*l1pone*pow(0.17862041985779245,2) - 0.7200006535515976*l4one*pow(0.17862041985779245,2) - 1.603364156071144*l4pone*pow(0.17862041985779245,2) - 2.630021167749498*xipone*pow(0.17862041985779245,2) - 10.045693496542476*chi2one*xipone*pow(0.17862041985779245,2) + 30.13708048962743*chi3pone*xipone*pow(0.17862041985779245,2) + 5.260042335498996*etaone*xipone*pow(0.17862041985779245,2) + 1.4530641908850788*l1one*xipone*pow(0.17862041985779245,2) + 1.5310622850438729*l1pone*xipone*pow(0.17862041985779245,2) - 1.603364156071144*l4one*xipone*pow(0.17862041985779245,2) + 0.7655311425219364*l1one*xippone*pow(0.17862041985779245,2),2),
                pow(0.12784586698325978 + 0.12715382670186573*l2one*pow(0.17862041985779245,2),2) + 0.07160045039259076*pow(3.4417932031169736 - 2.7161430681278897*chi3pone + 3.822247476150835*xipone - 6.8988262614568585*pow(0.17862041985779245,2) - 24.94299556289515*chi3pone*pow(0.17862041985779245,2) + 2.6431819892228363*l2one*pow(0.17862041985779245,2) + 3.801557333548889*l2pone*pow(0.17862041985779245,2) - 4.205795330555218*l5one*pow(0.17862041985779245,2) + 3.801557333548889*l2one*xipone*pow(0.17862041985779245,2),2) + 29.552589768397787*pow(0.40916995066113737 - 1.0110379501803748*chi3pone - 0.5347770210962225*chi3ppone + 1.7315762735650717*xipone - 1.069554042192445*chi3pone*xipone + 0.75255613121932*xippone - 1.2096717008122218*pow(0.17862041985779245,2) - 9.284604850601443*chi3pone*pow(0.17862041985779245,2) - 4.910986104106538*chi3ppone*pow(0.17862041985779245,2) + 0.2745444030911792*l2one*pow(0.17862041985779245,2) + 1.0408236642601045*l2pone*pow(0.17862041985779245,2) - 0.7374633594130332*l5one*pow(0.17862041985779245,2) - 1.656144497399369*l5pone*pow(0.17862041985779245,2) - 2.7165975168643626*xipone*pow(0.17862041985779245,2) - 9.821972208213076*chi3pone*xipone*pow(0.17862041985779245,2) + 1.4150649054924773*l2one*xipone*pow(0.17862041985779245,2) + 1.4969649649294912*l2pone*xipone*pow(0.17862041985779245,2) - 1.656144497399369*l5one*xipone*pow(0.17862041985779245,2) + 0.7484824824647456*l2one*xippone*pow(0.17862041985779245,2),2),
                pow(0.09040067949053633 + 0.08991133311470835*l2one*pow(0.17862041985779245,2),2) + 0.03580022519629537*pow(3.4417932031169736 + 2.7161430681278897*chi2one - 2.7161430681278897*chi3pone + 1.5024818555594324*etaone + 3.822247476150835*xipone - 6.8988262614568585*pow(0.17862041985779245,2) + 24.94299556289515*chi2one*pow(0.17862041985779245,2) - 24.94299556289515*chi3pone*pow(0.17862041985779245,2) - 13.797652522913717*etaone*pow(0.17862041985779245,2) + 2.6431819892228363*l2one*pow(0.17862041985779245,2) + 3.801557333548889*l2pone*pow(0.17862041985779245,2) + 3.801557333548889*l3one*pow(0.17862041985779245,2) + 4.205795330555218*l5one*pow(0.17862041985779245,2) - 8.411590661110436*l6one*pow(0.17862041985779245,2) + 3.801557333548889*l2one*xipone*pow(0.17862041985779245,2),2) + 14.776294884198888*pow(0.4091699506611374 + 1.0110379501803748*chi2one + 1.0695540421924452*chi2pone - 1.0110379501803748*chi3pone - 0.5347770210962226*chi3ppone + 0.2634520297761883*etaone + 0.5916424509412961*etapone + 1.7315762735650717*xipone + 1.0695540421924452*chi2one*xipone - 1.0695540421924452*chi3pone*xipone + 0.5916424509412961*etaone*xipone + 0.7525561312193201*xippone - 1.2096717008122222*pow(0.17862041985779245,2) + 9.284604850601443*chi2one*pow(0.17862041985779245,2) + 9.821972208213076*chi2pone*pow(0.17862041985779245,2) - 9.284604850601443*chi3pone*pow(0.17862041985779245,2) - 4.910986104106538*chi3ppone*pow(0.17862041985779245,2) - 2.4193434016244444*etaone*pow(0.17862041985779245,2) - 5.433195033728726*etapone*pow(0.17862041985779245,2) + 0.2745444030911792*l2one*pow(0.17862041985779245,2) + 1.0408236642601045*l2pone*pow(0.17862041985779245,2) + 1.4150649054924773*l3one*pow(0.17862041985779245,2) + 1.4969649649294914*l3pone*pow(0.17862041985779245,2) + 0.7374633594130333*l5one*pow(0.17862041985779245,2) + 1.6561444973993693*l5pone*pow(0.17862041985779245,2) - 3.131071216225436*l6one*pow(0.17862041985779245,2) - 3.3122889947987386*l6pone*pow(0.17862041985779245,2) - 2.716597516864363*xipone*pow(0.17862041985779245,2) + 9.821972208213076*chi2one*xipone*pow(0.17862041985779245,2) - 9.821972208213076*chi3pone*xipone*pow(0.17862041985779245,2) - 5.433195033728726*etaone*xipone*pow(0.17862041985779245,2) + 1.4150649054924773*l2one*xipone*pow(0.17862041985779245,2) + 1.4969649649294914*l2pone*xipone*pow(0.17862041985779245,2) + 1.4969649649294914*l3one*xipone*pow(0.17862041985779245,2) + 1.6561444973993693*l5one*xipone*pow(0.17862041985779245,2) - 3.3122889947987386*l6one*xipone*pow(0.17862041985779245,2) + 0.7484824824647457*l2one*xippone*pow(0.17862041985779245,2),2),
            };

            double result = 0.0;
            for (const auto & b : weak_bounds)
            {
                result += b;
            }

            //static_assert(std::numeric_limits<double>::is_iec559, "Unexpected absence of IEC559");
            //return (0.0 < result) && (result < 1.0) ? std::log(1.0) : std::log(1.0e-100);

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
                static const double sigma2 = 0.71 / 5.69; // cf. [CLN1997], eq. (10), p.6
                return -pow((result - 1.0), 2) / sigma2 / 2.0;
            }
        }

        double bound_0m_z2() const
        {
            const std::array<double, 3> weak_bounds
            {
                pow(0.08235327828970843 + 0.012226527430372875*etaone - 0.05613957214256857*pow(0.17862041985779245,2) + 0.11227914428513713*etaone*pow(0.17862041985779245,2) + 0.07630698319762817*l2one*pow(0.17862041985779245,2) + 0.03422488716605567*l5one*pow(0.17862041985779245,2) - 0.06844977433211134*l6one*pow(0.17862041985779245,2),2) + 0.0006704953120457946*pow(-23.157364603462458 - 16.844087755078398*chi2one + 16.844087755078398*chi3pone - 3.4452898308770963*etaone - 3.777419697877775*etapone - 25.443274663873034*xipone - 3.777419697877775*etaone*xipone + 15.81946289443559*pow(0.17862041985779245,2) + 154.68331218116623*chi2one*pow(0.17862041985779245,2) - 464.0499365434987*chi3pone*pow(0.17862041985779245,2) - 31.63892578887118*etaone*pow(0.17862041985779245,2) - 34.688954300297446*etapone*pow(0.17862041985779245,2) - 21.50239916712641*l2one*pow(0.17862041985779245,2) - 23.575254957536092*l2pone*pow(0.17862041985779245,2) - 23.575254957536092*l3one*pow(0.17862041985779245,2) - 9.644165638005001*l5one*pow(0.17862041985779245,2) - 10.573874198944175*l5pone*pow(0.17862041985779245,2) + 29.862205474954177*l6one*pow(0.17862041985779245,2) + 21.14774839788835*l6pone*pow(0.17862041985779245,2) + 17.344477150148723*xipone*pow(0.17862041985779245,2) - 34.688954300297446*etaone*xipone*pow(0.17862041985779245,2) - 23.575254957536092*l2one*xipone*pow(0.17862041985779245,2) - 10.573874198944175*l5one*xipone*pow(0.17862041985779245,2) + 21.14774839788835*l6one*xipone*pow(0.17862041985779245,2),2) + 0.5320248102904437*pow(-2.0351618624215435 - 5.559087973749947*chi2one - 4.783760323361458*chi2pone + 5.559087973749947*chi3pone + 2.391880161680729*chi3ppone - 0.30731004820015934*etaone - 1.2466693785745335*etapone - 0.5363980150822842*etappone - 8.383233096101703*xipone - 4.783760323361458*chi2one*xipone + 4.783760323361458*chi3pone*xipone - 1.2466693785745335*etaone*xipone - 1.0727960301645685*etapone*xipone - 3.6129747601418005*xippone - 0.5363980150822842*etaone*xippone + 1.4110510706589834*pow(0.17862041985779245,2) + 51.05044292035792*chi2one*pow(0.17862041985779245,2) + 43.93042248757584*chi2pone*pow(0.17862041985779245,2) - 153.15132876107376*chi3pone*pow(0.17862041985779245,2) - 65.89563373136377*chi3ppone*pow(0.17862041985779245,2) - 2.8221021413179668*etaone*pow(0.17862041985779245,2) - 11.44846497339134*etapone*pow(0.17862041985779245,2) - 4.925872082049402*etappone*pow(0.17862041985779245,2) - 1.917952813504359*l2one*pow(0.17862041985779245,2) - 6.106731660891956*l2pone*pow(0.17862041985779245,2) - 7.7805885494163824*l3one*pow(0.17862041985779245,2) - 6.695427554097704*l3pone*pow(0.17862041985779245,2) - 0.8602321292403814*l5one*pow(0.17862041985779245,2) - 2.738965600384307*l5pone*pow(0.17862041985779245,2) + 5.210181110470757*l6one*pow(0.17862041985779245,2) + 8.480936207191363*l6pone*pow(0.17862041985779245,2) + 5.72423248669567*xipone*pow(0.17862041985779245,2) + 43.93042248757584*chi2one*xipone*pow(0.17862041985779245,2) - 131.79126746272755*chi3pone*xipone*pow(0.17862041985779245,2) - 11.44846497339134*etaone*xipone*pow(0.17862041985779245,2) - 9.851744164098804*etapone*xipone*pow(0.17862041985779245,2) - 7.7805885494163824*l2one*xipone*pow(0.17862041985779245,2) - 6.695427554097704*l2pone*xipone*pow(0.17862041985779245,2) - 6.695427554097704*l3one*xipone*pow(0.17862041985779245,2) - 3.489716851989995*l5one*xipone*pow(0.17862041985779245,2) - 3.0030050064227494*l5pone*xipone*pow(0.17862041985779245,2) + 9.98243871040274*l6one*xipone*pow(0.17862041985779245,2) + 6.006010012845499*l6pone*xipone*pow(0.17862041985779245,2) + 2.462936041024701*xippone*pow(0.17862041985779245,2) - 4.925872082049402*etaone*xippone*pow(0.17862041985779245,2) - 3.347713777048852*l2one*xippone*pow(0.17862041985779245,2) - 1.5015025032113747*l5one*xippone*pow(0.17862041985779245,2) + 3.0030050064227494*l6one*xippone*pow(0.17862041985779245,2),2),
                pow(0.08107635004175542 - 0.012827568848485876*etaone - 0.058899326148354775*pow(0.17862041985779245,2) - 0.11779865229670955*etaone*pow(0.17862041985779245,2) + 0.07468798073339414*l1one*pow(0.17862041985779245,2) - 0.035907341553385846*l4one*pow(0.17862041985779245,2),2) + 0.00009656888000882487*pow(-60.39427414853287 + 43.44233946393831*chi2one - 130.32701839181496*chi3pone + 9.574316510664799*etaone + 10.442770749004271*etapone - 66.00328921656983*xipone + 10.442770749004271*etaone*xipone + 43.96162651473802*pow(0.17862041985779245,2) - 398.94145975074144*chi2one*pow(0.17862041985779245,2) + 398.94145975074144*chi3pone*pow(0.17862041985779245,2) + 87.92325302947604*etaone*pow(0.17862041985779245,2) + 95.89847733472939*etapone*pow(0.17862041985779245,2) - 55.74605566575123*l1one*pow(0.17862041985779245,2) - 60.80259398466964*l1pone*pow(0.17862041985779245,2) + 26.800733416390276*l4one*pow(0.17862041985779245,2) + 29.231738334615436*l4pone*pow(0.17862041985779245,2) + 47.949238667364696*xipone*pow(0.17862041985779245,2) + 95.89847733472939*etaone*xipone*pow(0.17862041985779245,2) - 60.80259398466964*l1one*xipone*pow(0.17862041985779245,2) + 29.231738334615436*l4one*xipone*pow(0.17862041985779245,2),2) + 0.2998419192092904*pow(-2.711516247045916 + 7.277558949884631*chi2one + 6.236998141865022*chi2pone - 21.832676849653886*chi3pone - 9.355497212797532*chi3ppone + 0.4360577422754174*etaone + 1.7493965717268902*etapone + 0.7496320704773283*etappone - 11.039796189889389*xipone + 6.236998141865022*chi2one*xipone - 18.710994425595064*chi3pone*xipone + 1.7493965717268902*etaone*xipone + 1.4992641409546565*etapone*xipone - 4.738032036033057*xippone + 0.7496320704773283*etaone*xippone + 2.0022116026160823*pow(0.17862041985779245,2) - 66.83157552551025*chi2one*pow(0.17862041985779245,2) - 57.275855165299234*chi2pone*pow(0.17862041985779245,2) + 66.83157552551025*chi3pone*pow(0.17862041985779245,2) + 28.637927582649617*chi3ppone*pow(0.17862041985779245,2) + 4.004423205232165*etaone*pow(0.17862041985779245,2) + 16.065129793182617*etapone*pow(0.17862041985779245,2) + 6.884051737601438*etappone*pow(0.17862041985779245,2) - 2.5389278855874524*l1one*pow(0.17862041985779245,2) - 8.003437427494251*l1pone*pow(0.17862041985779245,2) + 1.2206267979399612*l4one*pow(0.17862041985779245,2) + 3.8477698618741063*l4pone*pow(0.17862041985779245,2) + 8.032564896591309*xipone*pow(0.17862041985779245,2) - 57.275855165299234*chi2one*xipone*pow(0.17862041985779245,2) + 57.275855165299234*chi3pone*xipone*pow(0.17862041985779245,2) + 16.065129793182617*etaone*xipone*pow(0.17862041985779245,2) + 13.768103475202876*etapone*xipone*pow(0.17862041985779245,2) - 10.18578804662789*l1one*xipone*pow(0.17862041985779245,2) - 8.729402476534554*l1pone*xipone*pow(0.17862041985779245,2) + 4.896966912068792*l4one*xipone*pow(0.17862041985779245,2) + 4.196788200778742*l4pone*xipone*pow(0.17862041985779245,2) + 3.442025868800719*xippone*pow(0.17862041985779245,2) + 6.884051737601438*etaone*xippone*pow(0.17862041985779245,2) - 4.364701238267277*l1one*xippone*pow(0.17862041985779245,2) + 2.098394100389371*l4one*xippone*pow(0.17862041985779245,2),2),
                pow(0.1146648004474738 - 0.07982781284440091*pow(0.17862041985779245,2) + 0.10768215603064406*l2one*pow(0.17862041985779245,2) - 0.04866616867642605*l5one*pow(0.17862041985779245,2),2) + 0.43722981259908034*pow(-1.259741831367107 + 0.930829081491412*chi3pone - 1.387283507292596*xipone + 0.8813520597575941*pow(0.17862041985779245,2) + 8.54802750337334*chi3pone*pow(0.17862041985779245,2) - 1.1888825039179742*l2one*pow(0.17862041985779245,2) - 1.3028032884376906*l2pone*pow(0.17862041985779245,2) + 0.5373068166991254*l5one*pow(0.17862041985779245,2) + 0.588792488230531*l5pone*pow(0.17862041985779245,2) + 0.9658047434792989*xipone*pow(0.17862041985779245,2) - 1.3028032884376906*l2one*xipone*pow(0.17862041985779245,2) + 0.588792488230531*l5one*xipone*pow(0.17862041985779245,2),2) + 59.96258889104367*pow(-0.26491321571785836 + 0.7392459455764908*chi3pone + 0.3179396050091804*chi3ppone - 1.097494530458977*xipone + 0.6358792100183608*chi3pone*xipone - 0.47384893651759824*xippone + 0.18928896486929*pow(0.17862041985779245,2) + 6.788673452724926*chi3pone*pow(0.17862041985779245,2) + 2.9197159200007277*chi3ppone*pow(0.17862041985779245,2) - 0.2553376213584142*l2one*pow(0.17862041985779245,2) - 0.8121637821894662*l2pone*pow(0.17862041985779245,2) + 0.11539798429490701*l5one*pow(0.17862041985779245,2) + 0.36705152528399193*l5pone*pow(0.17862041985779245,2) + 0.7670229207833379*xipone*pow(0.17862041985779245,2) + 5.839431840001455*chi3pone*xipone*pow(0.17862041985779245,2) - 1.0346604634635794*l2one*xipone*pow(0.17862041985779245,2) - 0.8899867250964532*l2pone*xipone*pow(0.17862041985779245,2) + 0.4676072851236219*l5one*xipone*pow(0.17862041985779245,2) + 0.4022230393585197*l5pone*xipone*pow(0.17862041985779245,2) + 0.32988610343566505*xippone*pow(0.17862041985779245,2) - 0.4449933625482266*l2one*xippone*pow(0.17862041985779245,2) + 0.20111151967925986*l5one*xippone*pow(0.17862041985779245,2),2),
            };

            double result = 0.0;
            for (const auto & b : weak_bounds)
            {
                result += b;
            }

            //static_assert(std::numeric_limits<double>::is_iec559, "Unexpected absence of IEC559");
            //return (0.0 < result) && (result < 1.0) ? std::log(1.0) : std::log(1.0e-100);
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
                static const double sigma2 = 0.06 / 1.46; // cf. [CLN1997], eq. (10), p.6
                return -pow((result - 1.0), 2) / sigma2 / 2.0;
            }
        }

        double bound_1m_z2() const
        {
            const std::array<double, 4> weak_bounds
            {
                pow(0.01519421454368591 - 0.0017149979611038946*etaone + 0.014221600928947439*l1one*pow(0.17862041985779245,2) - 0.006790147058236399*l4one*pow(0.17862041985779245,2),2) + 0.0019811350500335083*pow(2.3573183205788975 - 2.3613965081232338*chi2one + 7.084189524369702*chi3pone - 0.2654272510182673*etaone - 0.3082453645276031*etapone + 2.730933975988174*xipone - 0.3082453645276031*etaone*xipone + 2.201052435782271*l1one*pow(0.17862041985779245,2) + 2.556121151122426*l1pone*pow(0.17862041985779245,2) - 1.050899248018569*l4one*pow(0.17862041985779245,2) - 1.2204278970774323*l4pone*pow(0.17862041985779245,2) + 2.556121151122426*l1one*xipone*pow(0.17862041985779245,2) - 1.2204278970774323*l4one*xipone*pow(0.17862041985779245,2),2) + 0.24438861020842015*pow(0.45421911755301403 - 1.889839924681759*chi2one - 1.7008870071413291*chi2pone + 5.669519774045277*chi3pone + 2.5513305107119937*chi3ppone - 0.05149303022644889*etaone - 0.24669063178437883*etapone - 0.11101281248891848*etappone + 2.189714694106252*xipone - 1.7008870071413291*chi2one*xipone + 5.102661021423987*chi3pone*xipone - 0.24669063178437883*etaone*xipone - 0.22202562497783696*etapone*xipone + 0.9835303180004296*xippone - 0.11101281248891848*etaone*xippone + 0.4270053627535599*l1one*pow(0.17862041985779245,2) + 1.5853929982450305*l1pone*pow(0.17862041985779245,2) - 0.20387502238588281*l4one*pow(0.17862041985779245,2) - 0.7569507580029399*l4pone*pow(0.17862041985779245,2) + 2.0456792356121163*l1one*xipone*pow(0.17862041985779245,2) + 1.841144949468343*l1pone*xipone*pow(0.17862041985779245,2) - 0.9767158362258269*l4one*xipone*pow(0.17862041985779245,2) - 0.8790603128915485*l4pone*xipone*pow(0.17862041985779245,2) + 0.9205724747341715*l1one*xippone*pow(0.17862041985779245,2) - 0.43953015644577426*l4one*xippone*pow(0.17862041985779245,2),2),
                pow(0.023115252451343042 - 0.0019352698939998842*etaone + 0.018489376230471206*l2one*pow(0.17862041985779245,2) - 0.018489376230471206*l5one*pow(0.17862041985779245,2),2) + 0.02754117923723504*pow(0.7627509471178091*(1.2501910374482086 - 0.10466929061730514*etaone + l2one*pow(0.17862041985779245,2) - l5one*pow(0.17862041985779245,2)) + 0.8912939958929453*xipone*(1.2501910374482086 - 0.10466929061730514*etaone + l2one*pow(0.17862041985779245,2) - l5one*pow(0.17862041985779245,2)) + 0.11141174948661817*(0.03322062512449487 - 1.6747086498768822*chi2one - 0.6917274858187126*chi3pone - 0.8373543249384411*etapone + 8*l2pone*pow(0.17862041985779245,2) - 8*l5pone*pow(0.17862041985779245,2)),2) + 0.02754117923723504*pow(0.11141174948661817*(-0.18629478987192996 - 3.3494172997537643*chi2one - 13.397669199015057*chi2pone - 1.3834549716374251*chi3pone - 2.7669099432748503*chi3ppone - 1.6747086498768822*etapone - 3.3494172997537643*etappone) + 1.620901280884261*(1.2501910374482086 - 0.10466929061730514*etaone + l2one*pow(0.17862041985779245,2) - l5one*pow(0.17862041985779245,2)) + 0.11141174948661817*(16*xipone + 32*xippone)*(1.2501910374482086 - 0.10466929061730514*etaone + l2one*pow(0.17862041985779245,2) - l5one*pow(0.17862041985779245,2)) + 0.7627509471178091*(0.03322062512449487 - 1.6747086498768822*chi2one - 0.6917274858187126*chi3pone - 0.8373543249384411*etapone + 8*l2pone*pow(0.17862041985779245,2) - 8*l5pone*pow(0.17862041985779245,2)) + 8*xipone*(0.7627509471178091*(1.2501910374482086 - 0.10466929061730514*etaone + l2one*pow(0.17862041985779245,2) - l5one*pow(0.17862041985779245,2)) + 0.11141174948661817*(0.03322062512449487 - 1.6747086498768822*chi2one - 0.6917274858187126*chi3pone - 0.8373543249384411*etapone + 8*l2pone*pow(0.17862041985779245,2) - 8*l5pone*pow(0.17862041985779245,2))),2),
                pow(0.02197242119121575 - 0.006278597399767078*etaone + 0.017575250928101457*l1one*pow(0.17862041985779245,2) - 0.017575250928101457*l4one*pow(0.17862041985779245,2),2) + 0.024885196861219643*pow(0.768056221971285*(1.2501910374482086 - 0.3572408397155849*etaone + l1one*pow(0.17862041985779245,2) - l4one*pow(0.17862041985779245,2)) + 0.8912939958929453*xipone*(1.2501910374482086 - 0.3572408397155849*etaone + l1one*pow(0.17862041985779245,2) - l4one*pow(0.17862041985779245,2)) + 0.11141174948661817*(0.03322062512449487 - 5.715853435449358*chi2one + 15.472851656471194*chi3pone - 2.857926717724679*etapone + 8*l1pone*pow(0.17862041985779245,2) - 8*l4pone*pow(0.17862041985779245,2)),2) + 0.024885196861219643*pow(0.11141174948661817*(-0.18629478987192996 - 11.431706870898717*chi2one - 45.72682748359487*chi2pone + 30.94570331294239*chi3pone + 61.89140662588478*chi3ppone - 5.715853435449358*etapone - 11.431706870898717*etappone) + 1.657671595809546*(1.2501910374482086 - 0.3572408397155849*etaone + l1one*pow(0.17862041985779245,2) - l4one*pow(0.17862041985779245,2)) + 0.11141174948661817*(16*xipone + 32*xippone)*(1.2501910374482086 - 0.3572408397155849*etaone + l1one*pow(0.17862041985779245,2) - l4one*pow(0.17862041985779245,2)) + 0.768056221971285*(0.03322062512449487 - 5.715853435449358*chi2one + 15.472851656471194*chi3pone - 2.857926717724679*etapone + 8*l1pone*pow(0.17862041985779245,2) - 8*l4pone*pow(0.17862041985779245,2)) + 8*xipone*(0.768056221971285*(1.2501910374482086 - 0.3572408397155849*etaone + l1one*pow(0.17862041985779245,2) - l4one*pow(0.17862041985779245,2)) + 0.11141174948661817*(0.03322062512449487 - 5.715853435449358*chi2one + 15.472851656471194*chi3pone - 2.857926717724679*etapone + 8*l1pone*pow(0.17862041985779245,2) - 8*l4pone*pow(0.17862041985779245,2))),2),
                pow(0.02299739585136725 + 0.0019254026286424128*etaone + 0.018395105357904118*l2one*pow(0.17862041985779245,2) - 0.018395105357904118*l5one*pow(0.17862041985779245,2),2) + pow(0.021962158618685525 + 0.020621230054384904*l2one*pow(0.17862041985779245,2) - 0.009319615218852892*l5one*pow(0.17862041985779245,2),2) + pow(0.02299739585136725 + 0.006571482884714323*etaone + 0.018395105357904118*l2one*pow(0.17862041985779245,2) + 0.018395105357904118*l5one*pow(0.17862041985779245,2) - 0.036790210715808236*l6one*pow(0.17862041985779245,2),2) + pow(0.015529591288767113 + 0.0016644371935421911*etaone + 0.014581411607863404*l2one*pow(0.17862041985779245,2) + 0.00658996311930023*l5one*pow(0.17862041985779245,2) - 0.01317992623860046*l6one*pow(0.17862041985779245,2),2) + 0.004020754305663588*pow(2.378775055680936 - 2.4034675845265463*chi3pone + 2.7708386888405045*xipone + 2.228016339671317*l2one*pow(0.17862041985779245,2) + 2.6016614777363856*l2pone*pow(0.17862041985779245,2) - 1.0069358099537038*l5one*pow(0.17862041985779245,2) - 1.1758020175454833*l5pone*pow(0.17862041985779245,2) + 2.6016614777363856*l2one*xipone*pow(0.17862041985779245,2) - 1.1758020175454833*l5one*xipone*pow(0.17862041985779245,2),2) + 0.0020103771528317935*pow(2.3787750556809373 + 2.403467584526546*chi2one - 2.4034675845265463*chi3pone + 0.2543233373625381*etaone + 0.2969741370043454*etapone + 2.7708386888405045*xipone + 0.2969741370043454*etaone*xipone + 2.228016339671318*l2one*pow(0.17862041985779245,2) + 2.6016614777363856*l2pone*pow(0.17862041985779245,2) + 2.601661477736385*l3one*pow(0.17862041985779245,2) + 1.0069358099537038*l5one*pow(0.17862041985779245,2) + 1.1758020175454829*l5pone*pow(0.17862041985779245,2) - 3.189673637452891*l6one*pow(0.17862041985779245,2) - 2.3516040350909657*l6pone*pow(0.17862041985779245,2) + 2.6016614777363856*l2one*xipone*pow(0.17862041985779245,2) + 1.1758020175454829*l5one*xipone*pow(0.17862041985779245,2) - 2.3516040350909657*l6one*xipone*pow(0.17862041985779245,2),2) + 0.49753962588006134*pow(0.45202058718255755 - 1.9123762000367492*chi3pone - 0.864247616998666*chi3ppone + 2.208911302833781*xipone - 1.728495233997332*chi3pone*xipone + 0.9963482550524346*xippone + 0.4263038125528885*l2one*pow(0.17862041985779245,2) + 1.6023164402896142*l2pone*pow(0.17862041985779245,2) - 0.19266491323964904*l5one*pow(0.17862041985779245,2) - 0.7241552828302755*l5pone*pow(0.17862041985779245,2) + 2.0700738893283583*l2one*xipone*pow(0.17862041985779245,2) + 1.8710297961549769*l2pone*xipone*pow(0.17862041985779245,2) - 0.9355548661381743*l5one*xipone*pow(0.17862041985779245,2) - 0.8455983332315951*l5pone*xipone*pow(0.17862041985779245,2) + 0.9355148980774884*l2one*xippone*pow(0.17862041985779245,2) - 0.42279916661579753*l5one*xippone*pow(0.17862041985779245,2),2) + 0.24876981294003087*pow(0.45202058718255755 + 1.9123762000367492*chi2one + 1.7284952339973318*chi2pone - 1.9123762000367492*chi3pone - 0.8642476169986659*chi3ppone + 0.048661675593823825*etaone + 0.23629454180695242*etapone + 0.10678704046961368*etappone + 2.208911302833781*xipone + 1.7284952339973318*chi2one*xipone - 1.7284952339973318*chi3pone*xipone + 0.23629454180695242*etaone*xipone + 0.21357408093922736*etapone*xipone + 0.9963482550524344*xippone + 0.10678704046961368*etaone*xippone + 0.4263038125528885*l2one*pow(0.17862041985779245,2) + 1.602316440289614*l2pone*pow(0.17862041985779245,2) + 2.0700738893283583*l3one*pow(0.17862041985779245,2) + 1.8710297961549764*l3pone*pow(0.17862041985779245,2) + 0.19266491323964902*l5one*pow(0.17862041985779245,2) + 0.7241552828302756*l5pone*pow(0.17862041985779245,2) - 1.3208846926174722*l6one*pow(0.17862041985779245,2) - 2.2939088988921457*l6pone*pow(0.17862041985779245,2) + 2.0700738893283583*l2one*xipone*pow(0.17862041985779245,2) + 1.8710297961549764*l2pone*xipone*pow(0.17862041985779245,2) + 1.8710297961549764*l3one*xipone*pow(0.17862041985779245,2) + 0.9355548661381743*l5one*xipone*pow(0.17862041985779245,2) + 0.845598333231595*l5pone*xipone*pow(0.17862041985779245,2) - 2.7167080655079436*l6one*xipone*pow(0.17862041985779245,2) - 1.69119666646319*l6pone*xipone*pow(0.17862041985779245,2) + 0.9355148980774882*l2one*xippone*pow(0.17862041985779245,2) + 0.4227991666157975*l5one*xippone*pow(0.17862041985779245,2) - 0.845598333231595*l6one*xippone*pow(0.17862041985779245,2),2) + 0.027261049494251476*pow(0.763288230730238*(1.2501910374482086 + 0.10466929061730514*etaone + l2one*pow(0.17862041985779245,2) - l5one*pow(0.17862041985779245,2)) + 0.8912939958929453*xipone*(1.2501910374482086 + 0.10466929061730514*etaone + l2one*pow(0.17862041985779245,2) - l5one*pow(0.17862041985779245,2)) + 0.11141174948661817*(0.03322062512449487 + 1.6747086498768822*chi2one - 7.39056208532624*chi3pone + 0.8373543249384411*etapone + 8*l2pone*pow(0.17862041985779245,2) - 8*l5pone*pow(0.17862041985779245,2)),2) + 0.027261049494251476*pow(0.763288230730238*(1.2501910374482086 + 0.3572408397155849*etaone + l2one*pow(0.17862041985779245,2) +
                l5one*pow(0.17862041985779245,2) - 2*l6one*pow(0.17862041985779245,2)) + 0.8912939958929453*xipone*(1.2501910374482086 + 0.3572408397155849*etaone + l2one*pow(0.17862041985779245,2) + l5one*pow(0.17862041985779245,2) - 2*l6one*pow(0.17862041985779245,2)) + 0.11141174948661817*(0.03322062512449487 + 5.715853435449358*chi2one - 7.39056208532624*chi3pone + 2.857926717724679*etapone + 8*l2pone*pow(0.17862041985779245,2) + 8*l3one*pow(0.17862041985779245,2) + 8*l5pone*pow(0.17862041985779245,2) - 8*l6one*pow(0.17862041985779245,2) - 16*l6pone*pow(0.17862041985779245,2)),2) + 0.027261049494251476*pow(0.11141174948661817*(-0.18629478987192996 + 3.3494172997537643*chi2one + 13.397669199015057*chi2pone - 14.78112417065248*chi3pone - 29.56224834130496*chi3ppone + 1.6747086498768822*etapone + 3.3494172997537643*etappone) + 1.6246113422052164*(1.2501910374482086 + 0.10466929061730514*etaone + l2one*pow(0.17862041985779245,2) - l5one*pow(0.17862041985779245,2)) + 0.11141174948661817*(16*xipone + 32*xippone)*(1.2501910374482086 + 0.10466929061730514*etaone + l2one*pow(0.17862041985779245,2) - l5one*pow(0.17862041985779245,2)) + 0.763288230730238*(0.03322062512449487 + 1.6747086498768822*chi2one - 7.39056208532624*chi3pone + 0.8373543249384411*etapone + 8*l2pone*pow(0.17862041985779245,2) - 8*l5pone*pow(0.17862041985779245,2)) + 8*xipone*(0.763288230730238*(1.2501910374482086 + 0.10466929061730514*etaone + l2one*pow(0.17862041985779245,2) - l5one*pow(0.17862041985779245,2)) + 0.11141174948661817*(0.03322062512449487 + 1.6747086498768822*chi2one - 7.39056208532624*chi3pone + 0.8373543249384411*etapone + 8*l2pone*pow(0.17862041985779245,2) - 8*l5pone*pow(0.17862041985779245,2))),2) + 0.027261049494251476*pow(1.6246113422052164*(1.2501910374482086 + 0.3572408397155849*etaone + l2one*pow(0.17862041985779245,2) + l5one*pow(0.17862041985779245,2) - 2*l6one*pow(0.17862041985779245,2)) + 0.11141174948661817*(16*xipone + 32*xippone)*(1.2501910374482086 + 0.3572408397155849*etaone + l2one*pow(0.17862041985779245,2) + l5one*pow(0.17862041985779245,2) - 2*l6one*pow(0.17862041985779245,2)) + 0.11141174948661817*(-0.18629478987192996 + 11.431706870898717*chi2one + 45.72682748359487*chi2pone - 14.78112417065248*chi3pone - 29.56224834130496*chi3ppone + 5.715853435449358*etapone + 11.431706870898717*etappone + 16*l3one*pow(0.17862041985779245,2) + 64*l3pone*pow(0.17862041985779245,2) - 16*l6one*pow(0.17862041985779245,2) - 64*l6pone*pow(0.17862041985779245,2)) + 0.763288230730238*(0.03322062512449487 + 5.715853435449358*chi2one - 7.39056208532624*chi3pone + 2.857926717724679*etapone + 8*l2pone*pow(0.17862041985779245,2) + 8*l3one*pow(0.17862041985779245,2) + 8*l5pone*pow(0.17862041985779245,2) - 8*l6one*pow(0.17862041985779245,2) - 16*l6pone*pow(0.17862041985779245,2)) + 8*xipone*(0.763288230730238*(1.2501910374482086 + 0.3572408397155849*etaone + l2one*pow(0.17862041985779245,2) + l5one*pow(0.17862041985779245,2) - 2*l6one*pow(0.17862041985779245,2)) + 0.11141174948661817*(0.03322062512449487 + 5.715853435449358*chi2one - 7.39056208532624*chi3pone + 2.857926717724679*etapone + 8*l2pone*pow(0.17862041985779245,2) + 8*l3one*pow(0.17862041985779245,2) + 8*l5pone*pow(0.17862041985779245,2) - 8*l6one*pow(0.17862041985779245,2) - 16*l6pone*pow(0.17862041985779245,2))),2),
            };

            double result = 0.0;
            for (const auto & b : weak_bounds)
            {
                result += b;
            }

            //static_assert(std::numeric_limits<double>::is_iec559, "Unexpected absence of IEC559");
            //return (0.0 < result) && (result < 1.0) ? std::log(1.0) : std::log(1.0e-100);
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
                static const double sigma2 = 0.06 / 1.46; // cf. [CLN1997], eq. (10), p.6
                return -pow((result - 1.0), 2) / sigma2 / 2.0;
            }
        }
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
    HQETUnitarityBounds::bound_1m() const
    {
        return _imp->bound_0m();
    }
}
