$OFFUPPER
$OFFSYMXREF
$OFFSYMLIST

$Title    Chemical Equilibrium Calculation Using GAMS
$Stitle   Steam Reactors for Chemical Looping

Sets    j  Components    / C,CH4,CO,H2,CO2,H2O,Fe2O3,Fe3O4,FeO,Fe,O2,N2 /
        i  Reactions     / 4,5,6,7,8,9,10,11,7a,8a,ox1,ox2,ox3,Fe1,Fe2 /
        k  Streams       / GasIn,SolidIn,GasOut,SolidOut /
        e  Elements      / C,O,H,Fe,N /
        iter Iterations  / 1*20 /
        ignore(j) Not included in EQ calcs - may be zero   / O2,N2,CH4,CO,CO2 /
        gas(j)              / CH4,CO,H2,CO2,H2O,O2,N2 /
        solid(j)            / C,Fe2O3,Fe3O4,FeO,Fe /
        gas_streams(k)   / GasIn,GasOut /
;

* 4-11 regressed from Fig 2 of Huang, et al, J Sust Bioenergy Syst, 2013,3,33-39
* 7a-8a regressed from Fig 3 of Xin, et al, J Nat Gas Chem, 2015, 14(4), 248-253
* ox1-3 regressed from Fig 5 of Fan, et al, AiChE J, 2015, 61, 1, 2-22
* Fe1 and Fe2 are regressed from Fig 4 and 6 of Xin, et al.

* Rxn 4  : CO + 3 Fe2O3 <-> CO2 + 2 Fe3O4
* Rxn 5  : CO + Fe2O3 <-> CO2 + 2 FeO
* Rxn 6  : H2 + 3 Fe2O3 <-> H2O + 2 Fe3O4
* Rxn 7  : H2 + Fe2O3 <-> H2O + 2 FeO
* Rxn 8  : CH4 + 3 Fe2O3 <-> 2 H2 + CO + 2 Fe3O4
* Rxn 9  : CH4 + 4 Fe2O3 <-> 2 H2O + CO2 + 8 FeO
* Rxn 10 : C + 3 Fe2O3 <-> CO + 2 Fe3O4
* Rxn 11 : C + 2 Fe2O3 <-> CO2 + 4 FeO
* Rxn 7a : CH4 + Fe3O4 <-> 3 Fe + CO2 + 2 H2O
* Rxn 8a : CH4 + 4 FeO <-> 4 Fe + CO2 + 2 H2O
* Rxn ox1: 2 Fe + O2 <-> 2 FeO
* Rxn ox2: 6 FeO + O2 <-> 2 Fe3O4
* Rxn ox3: 4 Fe3O4 + O2 <-> 6 Fe2O3
* Rxn Fe1: Fe3O4 + Fe <-> FeO
* Rxn Fe2: Fe2O3 + FeO <-> Fe3O4

Table  EqConst(i,*)    Constants in the equilibrium equations
         A        B
4        7.515672 4830.116457
5        3.295034 818.8999531
6        11.14578 870.2915707
7        7.651818 -3891.91264
8        41.68227 -26442.10669
9        55.84036 -37059.91603
10       28.49731 -15572.83047
11       28.11654 -19250.97118
7a       36.52848 -36362.43546
8a       30.85680 -30715.77616
ox1      -24.0225 70418.51336
ox2      -36.0141 59803.18789
ox3      -15.9728 63988.6020
Fe1      8.272653 -4573.15
Fe2      -0.63963 2697.507
;

Table ElemComp(j,e)    Number of elements in each component
         C       H       O       Fe      N
CH4      1       4       0       0       0
CO2      1       0       2       0       0
CO       1       0       1       0       0
H2       0       2       0       0       0
H2O      0       2       1       0       0
Fe2O3    0       0       3       2       0
FeO      0       0       1       1       0
Fe3O4    0       0       4       3       0
O2       0       0       2       0       0
C        1       0       0       0       0
Fe       0       0       0       1       0
N2       0       0       0       0       2;

Variables
         T               Temperature
         objective       Objective
         y(k,j)          Component mole fractions
         logYeq(j)       Log of either partial pressure (gas) or mole fraction (solid) used for equilibrium conditions
         logYeqSolidout(j) Yeq at solid outlet for countercurrent reactors
         F(k)            Molar flow rates
         p(k,j)          Partial pressures
         Pt              Total pressure;

Parameters
         Fc(k,j,iter)      Outlet component flow rates at each iteration
         ystore(k,j,iter)  Outlet component mole fractions at each iteration
         Ti(iter)          Temperature at each iteration
         iter1
         statstore(iter)   Model status at each iteration;

Equations
         EQ_rxn4
         EQ_rxn5
         EQ_rxn6
         EQ_rxn7
         EQ_rxn8
         EQ_rxn9
         EQ_rxn10
         EQ_rxn11
         EQ_rxn7a
         EQ_rxn8a

         EQ_ox1
         EQ_ox2
         EQ_ox3

         EQ_fe1
         EQ_fe2
         EQ_rxn10_counter

         EqYGasIn(j)
         EqYGasOut(j)
         EqYSolidIn(j)
         EqYSolidOut1(j)
         EqYSolidOut2(j)

         ElemBal(e)

         SumFrac_GasOut
         SumFrac_SolidOut

         CH4conversion
         COtoH2ratio
         Fullyreduced
         No_C_Rxn

         Partialpressures(k,j)

         obj
         obj1
;


SumFrac_GasOut .. 1 =e= sum(j,y('GasOut',j));
SumFrac_SolidOut .. 1 =e= sum(j,y('SolidOut',j));

ElemBal(e) .. F('SolidIn')*sum(j$ElemComp(j,e),y('SolidIn',j)*ElemComp(j,e)) +
              F('GasIn')*sum(j$ElemComp(j,e),y('GasIn',j)*ElemComp(j,e))   =e=
              F('SolidOut')*sum(j$ElemComp(j,e),y('SolidOut',j)*ElemComp(j,e)) +
              F('GasOut')*sum(j$ElemComp(j,e),y('GasOut',j)*ElemComp(j,e)) ;


EQ_rxn4 .. EqConst('4','A') + EqConst('4','B')/T =e=
                 logYeq('CO2') + 2*logYeq('Fe3O4')
                 -logYeq('CO') - 3*logYeq('Fe2O3');
EQ_rxn5 .. EqConst('5','A') + EqConst('5','B')/T =e=
                 logYeq('CO2') + 2*logYeq('FeO')
                 -logYeq('CO') - logYeq('Fe2O3');
EQ_rxn6 .. EqConst('6','A') + EqConst('6','B')/T =e=
                 logYeq('H2O') + 2*logYeq('Fe3O4')
                 -logYeq('H2') - 3*logYeq('Fe2O3');
EQ_rxn7 .. EqConst('7','A') + EqConst('7','B')/T =e=
                 logYeq('H2O') + 2*logYeq('FeO')
                 -logYeq('H2') - logYeq('Fe2O3');
EQ_rxn8 .. EqConst('8','A') + EqConst('8','B')/T =e=
                 2*logYeq('H2') + logYeq('CO') + 2*logYeq('Fe3O4')
                 -logYeq('CH4') - 3*logYeq('Fe2O3');
EQ_rxn9 .. EqConst('9','A') + EqConst('9','B')/T =e=
                 2*logYeq('H2O') + logYeq('CO2') + 8*logYeq('FeO')
                 -logYeq('CH4') - 4*logYeq('Fe2O3');
EQ_rxn10 .. EqConst('10','A') + EqConst('10','B')/T =e=
                 logYeq('CO') + 2*logYeq('Fe3O4')
                 -logYeq('C') - 3*logYeq('Fe2O3');
EQ_rxn11 .. EqConst('11','A') + EqConst('11','B')/T =e=
                 logYeq('CO2') + 4*logYeq('FeO')
                 -logYeq('C') - 2*logYeq('Fe2O3');
EQ_rxn7a .. EqConst('7a','A') + EqConst('7a','B')/T =e=
                 3*logYeq('Fe') + logYeq('CO2') + 2*logYeq('H2O')
                 -logYeq('Fe3O4') - logYeq('CH4');
EQ_rxn8a .. EqConst('8a','A') + EqConst('8a','B')/T =e=
                 4*logYeq('Fe') + logYeq('CO2') + 2*logYeq('H2O')
                 -4*logYeq('FeO') - logYeq('CH4');
EQ_rxn10_counter .. EqConst('10','A') + EqConst('10','B')/T =e=
                 logYeq('CO') + 2*logYeq('Fe3O4')
                 -logYeqSolidout('C') - 3*logYeq('Fe2O3');

EQ_fe1 .. EqConst('Fe1','A') + EqConst('Fe1','B')/T =e=
                 4*logYeqSolidout('FeO') - logYeqSolidout('Fe3O4')
                 - logYeqSolidout('Fe');
EQ_fe2 .. EqConst('Fe2','A') + EqConst('Fe2','B')/T =e=
                 logYeqSolidout('Fe3O4') - logYeqSolidout('Fe2O3')
                 - logYeqSolidout('FeO');

EQ_ox1 .. EqConst('ox1','A') + EqConst('ox1','B')/T =e=
                 2*logYeq('Fe3O4') - 6*logYeq('FeO') - logYeq('O2');
EQ_ox2 .. EqConst('ox2','A') + EqConst('ox2','B')/T =e=
                 6*logYeq('Fe2O3') - 4*logYeq('Fe3O4') - logYeq('O2');
EQ_ox3 .. EqConst('ox3','A') + EqConst('ox3','B')/T =e=
                 2*logYeq('FeO') - 2*logYeq('Fe') - logYeq('O2');

EqYGasIn(j)$(not ignore(j) and gas(j)) .. logYeq(j) + log(Pt)
         =e= log(p('GasIn',j));
EqYGasOut(j)$(not ignore(j) and gas(j)) .. logYeq(j) + log(Pt)
         =e= log(p('GasOut',j));
EqYSolidIn(j)$(not ignore(j) and solid(j)) .. logYeq(j) + log(F('SolidIn'))
         =e= log(y('SolidIn',j)*F('SolidIn'));
EqYSolidOut1(j)$(not ignore(j) and solid(j)) .. logYeq(j) + log(F('SolidOut'))
         =e= log(y('SolidOut',j)*F('SolidOut'));
EqYSolidOut2(j)$(not ignore(j) and solid(j)) .. logYeqSolidout(j) +
          log(F('SolidOut')) =e= log(y('SolidOut',j)*F('SolidOut'));

CH4conversion .. F('GasOut')*y('GasOut','CH4') =e=
         0.05*F('GasIn')*y('GasIn','CH4');
COtoH2ratio .. y('GasOut','CO') =e= y('GasOut','H2')*1;
Fullyreduced .. y('SolidOut','Fe3O4') =e= 0;
No_C_Rxn .. y('SolidIn','C')*F('SolidIn') =e= y('SolidOut','C')*F('SolidOut');

Partialpressures(k,j)$(gas_streams(k) and gas(j)) .. y(k,j)*Pt =e= p(k,j);

obj .. objective =g= 1;
obj1 .. objective =e= -y('GasOut','H2');
* obj .. objective =e= y('SolidOut','Fe2O3');
* obj .. objective =e= -y('GasOut','CH4');
* obj .. objective =e= y('GasOut','H2');



****************************************************************************
****************************************************************************
* Concurrent Steam Reactor
Model  EQsteam_concur  /
* EQ_rxn4
* EQ_rxn5
         EQ_rxn6
         EQ_rxn7
* EQ_rxn8
* EQ_rxn9
* EQ_rxn10
* EQ_rxn11
* EQ_rxn7a
* EQ_rxn8a
         EQ_fe1
* EQ_fe2
         EqYGasOut
         EqYSolidOut1
         EqYSolidOut2
         ElemBal
         SumFrac_GasOut
         SumFrac_SolidOut
         Partialpressures
         obj
         No_C_rxn
* COtoH2ratio
* Fullyreduced
* CH4conversion
/;

* option nlp = conopt;
option nlp = ipopth;
* option nlp = knitro;
* option nlp = snopt;
* option nlp = minos;
* option nlp = baron;
option iterlim = 1000;

iter1 = 0;
loop(iter,
         iter1 = iter1 + 1;
         T.fx = 773.15 + 50*iter1;
         Pt.fx = 100000;
         y.l(k,j) = 0;
         logYeq.l(j) = 0;
         F.l(k) = 0;

         y.lo(k,j) = 0;
         F.lo(k) = 0;
         F.fx('GasIn') = 1;
         F.fx('SolidIn') = 1;

* Example inlet conditions for steam reactor
         y.fx('GasIn','CH4') = 0;
         y.fx('GasIn','CO2') = 0;
         y.fx('GasIn','CO') = 0;
         y.fx('GasIn','H2') = 1e-5;
         y.fx('GasIn','H2O') = 1-1e-5;
         y.fx('GasIn','O2') = 0;
         y.fx('GasIn','N2') = 0;

         y.fx('SolidIn','Fe2O3') = 1e-5;
         y.fx('SolidIn','Fe3O4') = 1e-5;
         y.fx('SolidIn','FeO') = 0.2;
         y.fx('SolidIn','Fe') = 0.8-3e-5;
         y.fx('SolidIn','C') = 1e-5;

         y.fx('SolidIn',j)$gas(j) = 0;
         y.fx('SolidOut',j)$gas(j) = 0;
         y.fx('GasIn',j)$solid(j) = 0;
         y.fx('GasOut',j)$solid(j) = 0;
         y.fx('GasOut','O2') = 0;

         y.l('SolidOut',j) = y.l('SolidIn',j);
         y.l('GasOut',j) = y.l('GasIn',j);
         F.l('SolidOut') = F.l('SolidIn');
         F.l('GasOut') = F.l('GasIn');
         p.l(k,j) = y.l(k,j)*Pt.l;
         logYeq.l(j)$(not ignore(j)) = log((p.l('GasIn',j)+p.l('SolidIn',j)));
         p.fx(k,j)$((ignore(j) and gas(j)) or solid(j)) = 0;

         EQsteam_concur.optfile = 1;
         solve EQsteam_concur using nlp  minimizing  objective;
         statstore(iter) = EQsteam_concur.modelstat;
         loop(j,
                 loop(k,
                         Fc(k,j,iter) = y.l(k,j)*F.l(k);
                         ystore(k,j,iter) = y.l(k,j);
                 );
         );
         Ti(iter) = T.l;
);

display Ti,Fc,ystore,statstore;


****************************************************************************
****************************************************************************
* Countercurrent Steam Reactor
* This is pretty sensitive to whether reaction 6 or 7 is chosen
* EQ in rxn 6 or rxn 7 essentially sets the H2/H2O ratio in the outlet gas
Model  EQsteam_counter  /
* EQ_rxn4
* EQ_rxn5
* EQ_rxn6
         EQ_rxn7
* EQ_rxn8
* EQ_rxn9
* EQ_rxn10
* EQ_rxn10_counter
* EQ_rxn11
* EQ_rxn7a
* EQ_rxn8a
         EQ_fe1
         EQ_fe2
         EqYSolidIn
         EqYGasOut
         EqYSolidOut2
         ElemBal
         SumFrac_GasOut
         SumFrac_SolidOut
         obj
         Partialpressures
         No_C_rxn
* COtoH2ratio
* Fullyreduced
* CH4conversion
/;

* option nlp = conopt;
option nlp = ipopth;
* option nlp = knitro;
* option nlp = snopt;
* option nlp = minos;
* option nlp = baron;
option iterlim = 1000;

iter1 = 0;
loop(iter,
         iter1 = iter1 + 1;
         T.fx = 773.15 + 50*iter1;
         Pt.fx = 100000;
         y.l(k,j) = 0;
         logYeq.l(j) = 0;
         F.l(k) = 0;

         y.lo(k,j) = 0;
         F.lo(k) = 0;
         F.fx('GasIn') = 1;
         F.fx('SolidIn') = 1;

* Example inlet conditions for steam reactor
         y.fx('GasIn','CH4') = 0;
         y.fx('GasIn','CO2') = 0;
         y.fx('GasIn','CO') = 0;
         y.fx('GasIn','H2') = 1e-5;
         y.fx('GasIn','H2O') = 1-1e-5;
         y.fx('GasIn','O2') = 0;
         y.fx('GasIn','N2') = 0;

         y.fx('SolidIn','Fe2O3') = 1e-5;
         y.fx('SolidIn','Fe3O4') = 1e-5;
         y.fx('SolidIn','FeO') = 0.2;
         y.fx('SolidIn','Fe') = 0.8-3e-5;
         y.fx('SolidIn','C') = 1e-5;

         y.fx('SolidIn',j)$gas(j) = 0;
         y.fx('SolidOut',j)$gas(j) = 0;
         y.fx('GasIn',j)$solid(j) = 0;
         y.fx('GasOut',j)$solid(j) = 0;
         y.fx('GasOut','O2') = 0;

         y.l('SolidOut',j) = y.l('SolidIn',j);
         y.l('GasOut',j) = y.l('GasIn',j);
         F.l('SolidOut') = F.l('SolidIn');
         F.l('GasOut') = F.l('GasIn');
         p.l(k,j) = y.l(k,j)*Pt.l;
         logYeq.l(j)$(not ignore(j)) = log((p.l('GasIn',j)+p.l('SolidIn',j)));
         p.fx(k,j)$((ignore(j) and gas(j)) or solid(j)) = 0;

         EQsteam_counter.optfile = 1;
         solve EQsteam_counter using nlp  minimizing  objective;
         statstore(iter) = EQsteam_counter.modelstat;
         loop(j,
                 loop(k,
                         Fc(k,j,iter) = y.l(k,j)*F.l(k);
                         ystore(k,j,iter) = y.l(k,j);
                 );
         );
         Ti(iter) = T.l;
);

display Ti,Fc,ystore,statstore;


























