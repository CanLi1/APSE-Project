$OFFUPPER
$OFFSYMXREF
$OFFSYMLIST

$Title    Chemical Equilibrium Calculation Using GAMS
$Stitle   Water Gas Shift Equilibrium Reactor

Sets    j  Components    / CH4,CO,H2,CO2,H2O,O2,N2 /
        i  Reactions     / WGS /
        k  Streams       / GasIn,GasOut /
        e  Elements      / C,O,H,N /
        iter Iterations  / 1*20 /
        ignore(j) Not included in EQ calcs - may be zero / O2,CH4 /
        gas(j)           / CH4,CO,H2,CO2,H2O,O2,N2 /
        gas_streams(k)   / GasIn,GasOut /
;

* Callaghan, Caitlin (2006). "Kinetics and Catalysis of the
*                                                Water-Gas-Shift Reaction".
* WGS  : CO + H2O <-> CO2 + H2

Table  EqConst(i,*)    Constants in the equilibrium equations
         A        B      C
WGS      -2.4198  2180.6 0.0003855
;

Table ElemComp(j,e)    Number of elements in each component
         C       H       O       N
CH4      1       4       0       0
CO2      1       0       2       0
CO       1       0       1       0
H2       0       2       0       0
H2O      0       2       1       0
O2       0       0       2       0
N2       0       0       0       2;

Variables
         T               Temperature
         objective       Objective
         y(k,j)          Component mole fractions
         logPeq(j)       Log of partial pressure (gas) used for equilibrium
         F(k)            Molar flow rates
         p(k,j)          Partial pressures
         Pt              Total pressure;

Parameters
         Fc(k,j,iter)      Outlet component flow rates at each iteration
         ystore(k,j,iter)  Outlet component mole fractions at each iteration
         Ti(iter)          Temperature at each iteration
         statstore(iter)   Model status at each iteration
         iter1;

Equations
         EQ_WGS

         EqYGasOut(j)

         ElemBal(e)

         SumFrac_GasOut

         COtoH2ratio
         Fullyreduced
         No_C_Rxn
         No_CH4_Rxn

         PartialPressures(k,j)

         obj
         obj1
;


SumFrac_GasOut .. 1 =e= sum(j,y('GasOut',j));

ElemBal(e) .. F('GasIn')*sum(j$ElemComp(j,e),y('GasIn',j)*ElemComp(j,e))   =e=
              F('GasOut')*sum(j$ElemComp(j,e),y('GasOut',j)*ElemComp(j,e)) ;


EQ_WGS .. EqConst('WGS','A') + EqConst('WGS','B')/T + EqConst('WGS','C')*T =e=
                 logPeq('CO2') + logPeq('H2')
                 -logPeq('CO') - 3*logPeq('H2O');


EqYGasOut(j)$(not ignore(j) and gas(j)) .. logPeq(j) + log10(Pt)
         =e= log10(p('GasOut',j));

COtoH2ratio .. y('GasOut','CO') =e= y('GasOut','H2')*1;
No_CH4_Rxn .. y('GasIn','CH4')*F('GasIn') =e= y('GasOut','CH4')*F('GasOut');


obj .. objective =g= 1;
obj1 .. objective =e= y('GasOut','H2');

Partialpressures(k,j)$(gas_streams(k) and gas(j)) .. y(k,j)*Pt =e= p(k,j);

****************************************************************************
****************************************************************************
* WGS Reactor
Model  WGS_reactor  /
         SumFrac_GasOut
         ElemBal
         EQ_WGS
         EqYGasOut
*         COtoH2ratio
         No_CH4_Rxn
         PartialPressures
         obj
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
         T.fx = 250 + 50*iter1;
         Pt.fx = 100000;
         y.l(k,j) = 0;
         logPeq.l(j) = 0;
         F.l(k) = 0;

         y.lo(k,j) = 0;
         F.lo(k) = 0;

* Example inlet conditions for air reactor
         F.fx('GasIn') = 10;

         y.fx('GasIn','CH4') = 0;
         y.fx('GasIn','CO2') = 0.20;
         y.fx('GasIn','CO') = 0.20;
         y.fx('GasIn','H2') = 0.20;
         y.fx('GasIn','H2O') = 0.20;
         y.fx('GasIn','N2') = 0.20;
         y.fx('GasIn','O2') = 0;

         y.l('GasOut',j) = y.l('GasIn',j);
         F.l('GasOut') = F.l('GasIn');
         y.fx('GasOut',j)$ignore(j) = 0;
         p.l(k,j) = y.l(k,j)*Pt.l;
         logPeq.l(j)$(not ignore(j)) = log10(p.l('GasIn',j));
         p.fx(k,j)$ignore(j) = 0;

         WGS_reactor.optfile = 1;
         solve WGS_reactor using nlp  minimizing  objective;
         statstore(iter) = WGS_reactor.modelstat;
         loop(j,
                 loop(k,
                         Fc(k,j,iter) = y.l(k,j)*F.l(k);
                         ystore(k,j,iter) = y.l(k,j);
                 );
         );
         Ti(iter) = T.l;
);

display Ti,Fc,ystore,statstore;