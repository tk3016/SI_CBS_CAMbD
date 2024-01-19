SETS
s                   'candidate solvents'
/water, n-pentane, n-heptane, ethanol, 1-propanol, 1-butanol, 1-pentanol, acetone/
i                   'components in the mixture'
/ibuprofen, s1, s2/

ii(i)               'selected solvents in the mixture'
/s1,s2/

jj(i)               'solid solute in the mixture'
/ibuprofen/

k                   'groups'
/CH3, CH2, CH, aCH, aCCH2, aCCH, CH3CO, COOH, OH, H2O/

st                  'heated(initial) state or cooled(final) state'
/heated, cooled/

prp                 'melting/boiling T in K (p1, p2); molar mass in g/mol (p3); density in g/mL (p4)'
/p1*p4/

pr                  'UNIFAC group contributions'
/Q,R/

ic                  'Integer cuts'
/1*5/

;

set dyn(ic);

alias(k,m);
alias(s,ss);
alias(ic,ics);

Parameter
yv(ii,s,ic)
yieldk(ic,jj)
xi_sk(ic)
v_sk(ic)
Tk(ic,st)
yk(ic,ii,s)
xk(ic,i,st)
nmk(ic,i,st)
nmtk(ic,st)
;

Table GS(k, pr)     'Group specifications (R = volume, Q = surface area)'
$onDelim
,Q,R
CH3,0.848,0.9011
CH2,0.54,0.6744
CH,0.228,0.4469
aCH,0.400,0.5313
aCCH2,0.660,1.0396
aCCH,0.348,0.8121
CH3CO,1.488,1.6724
COOH,1.224,1.3013
OH,1.2,1
H2O,1.4,0.92
$offDelim
;

Table v(s, k)       'number of groups in molecule i'
$ondelim
,CH3,CH2,CH,aCH,aCCH2,aCCH,CH3CO,COOH,OH,H2O
water,0,0,0,0,0,0,0,0,0,1
n-pentane,2,3,0,0,0,0,0,0,0,0
n-heptane,2,5,0,0,0,0,0,0,0,0
ethanol,1,1,0,0,0,0,0,0,1,0
1-propanol,1,2,0,0,0,0,0,0,1,0
1-butanol,1,3,0,0,0,0,0,0,1,0
1-pentanol,1,4,0,0,0,0,0,0,1,0
acetone,1,0,0,0,0,0,1,0,0,0
$offdelim
;

Table vib(jj,k)     'identity of ibuprofen'
$onDelim
,CH3,CH2,CH,aCH,aCCH2,aCCH,CH3CO,COOH,OH,H2O
ibuprofen,3,0,1,4,1,1,0,1,0,0
$offDelim
;

Table a(k,m)        'group interaction parameters'
$onDelim
,CH3,CH2,CH,aCH,aCCH2,aCCH,CH3CO,COOH,OH,H2O
CH3,0,0,0,61.13,76.5,76.5,476.4,663.5,986.5,1318
CH2,0,0,0,61.13,76.5,76.5,476.4,663.5,986.5,1318
CH,0,0,0,61.13,76.5,76.5,476.4,663.5,986.5,1318
aCH,-11.12,-11.12,-11.12,0,167,167,25.77,537.4,636.1,903.8
aCCH2,-69.7,-69.7,-69.7,-146.8,0,0,-52.1,872.3,803.2,5695
aCCH,-69.7,-69.7,-69.7,-146.8,0,0,-52.1,872.3,803.2,5695
CH3CO,26.76,26.76,26.76,140.1,365.8,365.8,0,669.4,164.5,472.5
COOH,315.3,315.3,315.3,62.32,89.86,89.86,-297.8,0,-151,-66.17
OH,156.4,156.4,156.4,89.6,25.82,25.82,84,199,0,353.5
H2O,300,300,300,362.3,377.6,377.6,-195.4,-14.09,-229.1,0
$offDelim
;

Table prop(s,prp)   'All relevant properties'
$onDelim
,p1,p2,p3,p4
water,373.15,273.15,18.015,0.998
n-pentane,309,143,72.151,0.626
n-heptane,371.53,182.6,100.205,0.684
ethanol,351.39,159.01,46.069,0.789
1-propanol,370,147,60.096,0.804
1-butanol,390.8,183.3,74.123,0.81
1-pentanol,411,194.65,88.1482,0.811
acetone,329.3,178.7,58.08,0.791
$offDelim
;

Parameters
Hm(jj)              'Enthalpy of fusion in J/mol'
/
ibuprofen  25500
/

Tm(jj)              'Melting point in K'
/
ibuprofen  347.15
/

Mwib(jj)              'Molar mass of solids in g/mol'
/
ibuprofen  206.28
/

Rg                  'gas constant'
/8.3144/ ;

Parameter nib(jj, k);
nib(jj,k) = vib(jj,k);

Parameter qib(jj);
qib(jj) = sum(k, nib(jj,k)*GS(k,'Q'));

Parameter rib(jj);
rib(jj) = sum(k, nib(jj,k)*GS(k,'R'));

Parameter qs(s);
qs(s) = sum(k, v(s,k)*GS(k,'Q'));

Parameter rs(s);
rs(s) = sum(k, v(s,k)*GS(k,'R'));

option decimals=4;
display qs,rs;

Parameter eib(jj,k);
eib(jj,k) = vib(jj,k)*GS(k,'Q')/qib(jj);

Scalar  Tos         'Minimum temperature offset'
/10/
;

Positive Variables

yield(jj)           'Crystal yield of ibuprofen'
xi_s                'Solvent use in terms of mass'
v_s                 'Solvent use in terms of volume'
x(i,st)             'Liquid phase mole fraction of component i'
rc(i)               'van der Waals volume of component i'
qc(i)               'van der Waals area of component i'
J(i,st)             'UNIFAC intermediate variable to estimate solubility'
L(i,st)             'UNIFAC intermediate variable to estimate solubility'
T(st)               'Process temperature'
Mw(i)               'Molecular weight of component i'
den(ii)             'Density of selected solvents ii in g/mL'
vol(ii)             'Volume of selected solvents ii in mL at cooling state'
mass(i,st)          'Mass in g of each component'

th(k,st)            'UNIFAC intermediate to estimate solubility'
w(k,st)             'UNIFAC intermediate to estimate solubility'

x12(ii,st)          'Miscibility mole fraction of selected solvents'
th12(k,st)          'UNIFAC intermediate to estimate miscibility'
w12(k,st)           'UNIFAC intermediate to estimate miscibility'
sumr12(st)          'UNIFAC intermediate to estimate miscibility'
sumq12(st)          'UNIFAC intermediate to estimate miscibility'
b12(ii,k,st)        'UNIFAC intermediate to estimate miscibility'
nmt(st)             'Total mass at state st'
nm(i,st)            'Mass of each component at each state'
;

Free Variables
z                   'Objective function'
ps(k,m,st)          'UNIFAC intermediate'
bib(i,k,st)         'UNIFAC intermediate'

lng(i,st)           'Natural log of activity coefficient of ibuprofen'
lngc(i,st)          'Natural log of combinatorial activity coefficient of ibuprofen'
lngr(i,st)          'Natural log of residual activity coefficient of ibuprofen'

d_lngc12(ii,st)     'Natural log of combinatorial activity coefficient for miscibility'
d_lngr12(ii,st)     'Natural log of residual activity coefficient for miscibility'
d_lng12(ii,st)      'Natural log of activity coefficient for miscibility'

d_th12(k,st)        'UNIFAC intermediate to estimate miscibility'
d_w12(k,st)         'UNIFAC intermediate to estimate miscibility'
a12(k,st)           'UNIFAC intermediate to estimate miscibility'
c12(k,st)           'UNIFAC intermediate to estimate miscibility'
;

Binary Variables
y(ii,s)             'Selected solvents'
;

Integer Variables
n(i,k)              'Number of groups k in component i'
;

EQUATIONS
eq_z                'Objective function to maximise crystal yield'
eq_yield(jj)        'Equation to estimate yield of ibuprofen'
eq_soluse           'Equation to estimate solvent use'
eq_volsoluse        'Equation to estimate solvent use'
eq_x(st)            'Mass balance enforcing sum of all mole fractions to be one'
eq_s1               'Mass balance on first solvent'
eq_s2               'Mass balance on second solvent'
eq_API              'Mass balance on API'
eq_T1(s)            'Operational constraint on temperature'
eq_T2(s)            'Operational constraint on temperature'
eq_T(s)             'Temperature of heated state cannot be less than the temperature of cooled state'
eq_melt(ii,s)
eq_boil(ii,s)

eq_n1(ii,k)         'Equation to estimate number of gorups k in s1'
eq_n2(jj,k)         'Equation to estimate number of gorups k in s1'
eq_nmt(st)          'Equation to estimate total number of moles at each state'
eq_nm(i,st)         'Equation to estimate number of moles of i at each state'
eq_Mw(ii)           'Equation to assign molecular weight'
eq_Mwib(jj)         'Equation to assign molecular weight of ibuprofen'
eq_den(ii)          'Equation to assign molar density'
eq_vol(ii)          'Equation to assign volume'
eq_mass(i,st)       'Equation to estimate mass'

*Equations for solubility calculations
eq_qc1(ii)
eq_qc2(jj)
eq_rc1(ii)
eq_rc2(jj)
eq_J(jj,st)
eq_L(jj,st)
eq_lngc(jj,st)
eq_ps(k,m,st)
eq_bib(jj,k,st)
eq_th(k,st)
eq_w(k,st)
eq_lngr(jj,st)
eq_lng(jj,st)
eq_solub1(st)

*Equations for logical conditions regarding solvent selection
eq_logic1_s1
eq_logic1_s2
eq_logic2(s)
*eq_logic3(s,ss)

*Equations for miscibility constraints
eq_x12a(st)
eq_x12b(st)
eq_sumr12(st)
eq_sumq12(st)
eq_dlngc12(st)
eq_th12(k,st)
eq_dth12(k,st)
eq_w12(k,st)
eq_dw12(k,st)
eq_a12(k,st)
eq_b12(k,st)
eq_c12(k,st)
eq_dlngr12(st)
eq_dlng12(st)
eq_misc12(st)

*Equations for interger cuts
eq_cut1(ic)

*eq_thp
;

*******************Choose only 2 solvents*************************************
eq_logic1_s1..                         sum(s,y('s1',s)) =e= 1;
eq_logic1_s2..                         sum(s,y('s2',s)) =e= 1;

******************Each solvent at most once***********************************
eq_logic2(s)..                          y('s1',s) + y('s2',s) =l= 1;


******************UNIFAC******************************************************
*Process equations
eq_z..                                  z =e= yield('ibuprofen');
eq_yield(jj)..                          yield(jj) =e= 1 - (nm(jj,'cooled') / nm(jj,'heated'));
eq_soluse..                             xi_s * (mass('ibuprofen','heated') - mass('ibuprofen','cooled')) =e= mass('s1','cooled') + mass('s2','cooled');
eq_volsoluse..                          v_s * (mass('ibuprofen','heated') - mass('ibuprofen','cooled')) =e= sum(ii,vol(ii));
eq_x(st)..                              sum(i, x(i,st)) =e= 1;
eq_melt(ii,s)..                         y(ii,s)*prop(s,'p2') =l= T('cooled') - Tos;
eq_boil(ii,s)..                         prop(s,'p1') =g= (y(ii,s) * (T('heated') + Tos));
eq_T1(s)..                              T('heated') =l= prop(s,'p1')-Tos + 600*(1 - sum(ii, y(ii,s)));
eq_T2(s)..                              T('cooled') =g= prop(s,'p2')+Tos - 600*(1 - sum(ii, y(ii,s)));
eq_T(s)..                               T('heated') =g= T('cooled');

eq_s1..                                 nm('s1','heated') =e= nm('s1','cooled');
eq_s2..                                 nm('s2','heated') =l= nm('s2','cooled');
eq_API(jj)..                            nm(jj,'heated') =g= nm(jj,'cooled');
eq_nmt(st)..                            nmt(st) =e= sum(i, nm(i,st));
eq_nm(i,st)..                           x(i,st) =e= nm(i,st)/nmt(st);

eq_Mw(ii)..                             Mw(ii) =e= sum(s, prop(s, 'p3') * y(ii,s));
eq_Mwib(jj)..                           Mw(jj) =e= Mwib(jj);
eq_mass(i,st)..                         mass(i,st) =e= nm(i,st) * Mw(i);
eq_den(ii)..                            den(ii) =e= sum(s, prop(s,'p4')*y(ii,s));
eq_vol(ii)..                            vol(ii) * den(ii) =e= mass(ii,'cooled');
eq_n1(ii,k)..                           n(ii,k) =e= sum(s, v(s,k)*y(ii,s));
eq_n2(jj,k)..                           n(jj,k) =e= vib(jj,k);

*****Solubility calculations using UNIFAC
eq_qc1(ii)..                            qc(ii) =e= sum(s,qs(s) * y(ii,s));
eq_rc1(ii)..                            rc(ii) =e= sum(s,rs(s) * y(ii,s));
eq_qc2(jj)..                            qc(jj) =e= qib(jj);
eq_rc2(jj)..                            rc(jj) =e= rib(jj);

*Combinatorial part of the activity coefficient
eq_J(jj,st)..                           J(jj,st) * sum(i,x(i,st)*rc(i)) =e= rib(jj);
eq_L(jj,st)..                           L(jj,st) * sum(i,x(i,st)*qc(i)) =e= qib(jj);
eq_lngc(jj,st)..                        lngc(jj,st) =e= 1 - J(jj,st) + log(J(jj,st)) - 5*qib(jj)*(1 - J(jj,st)/L(jj,st) + log(J(jj,st)) - log(L(jj,st)));

*Residual part of the activity coefficient
eq_ps(k,m,st)..                         ps(k,m,st) =e= exp(-a(k,m)/T(st));
eq_bib(jj,k,st)..                       bib(jj,k,st) =e= sum(m,eib(jj,m)*ps(m,k,st));
eq_th(k,st)..                           th(k,st)*sum(i,x(i,st)*qc(i)) =e= sum(i,x(i,st)*GS(k,'Q')*n(i,k));
eq_w(k,st)..                            w(k,st) =e= sum(m,th(m,st)*ps(m,k,st));
eq_lngr(jj,st)..                        lngr(jj,st) =e= qib(jj)*(1-sum(k,th(k,st)*bib(jj,k,st)/w(k,st)-eib(jj,k)*(log(bib(jj,k,st)) - log(w(k,st)))));

*Overall activity coefficient
eq_lng(jj,st)..                         lng(jj,st) =e= lngc(jj,st) + lngr(jj,st);

*Solubility oonstraint
eq_solub1(st)..                         log(x('ibuprofen',st)) + lng('ibuprofen',st) =e= Hm('ibuprofen')/Rg * (1/Tm('ibuprofen') - 1/T(st));

***********************************Miscibility between s1 & s2**********
eq_x12a(st)..x12('s1',st)*(x('s1',st)+x('s2',st))=e=x('s1',st);
eq_x12b(st)..x12('s1',st)+x12('s2',st)=e=1;

eq_sumr12(st)..sumr12(st) =e= x12('s1',st)*rc('s1')+(1-x12('s1',st))*rc('s2');

eq_sumq12(st)..sumq12(st)=e=x12('s1',st)*qc('s1')+(1-x12('s1',st))*qc('s2');

*differ of combinatorial activity coefficient, i.e. dlngc/dx(s1)
eq_dlngc12(st)..d_lngc12('s1',st)*power(sumr12(st),2)*sumq12(st)=e=5*qc('s1')*(qc('s2')-qc('s1'))*power(sumr12(st),2)
                                        +(rc('s1')-rc('s2'))*(5*qc('s1')-1)*sumr12(st)*sumq12(st)
                                        +5*rc('s1')*(qc('s1')-qc('s2'))*sumr12(st)*sumq12(st)
                                        +((rc('s1')-rc('s2'))*(rc('s1')-5*rc('s1')*(sumq12(st))))*sumq12(st);

eq_th12(k,st)..th12(k,st)*sumq12(st)=e=x12('s1',st)*GS(k,'Q')*n('s1',k)+(1-x12('s1',st))*GS(k,'Q')*n('s2',k);
eq_w12(k,st)..w12(k,st)=e=sum(m,th12(m,st)*ps(m,k,st));
eq_dth12(k,st)..d_th12(k,st)*power(sumq12(st),2)=e=((GS(k,'Q')*n('s1',k)- GS(k,'Q')*n('s2',k))*sumq12(st))-(x12('s1',st)*GS(k,'Q')*n('s1',k)+ x12('s2',st)*GS(k,'Q')*n('s2',k))*(qc('s1')-qc('s2'));
eq_dw12(k,st)..d_w12(k,st)=e=sum(m,ps(m,k,st)*d_th12(m,st));
eq_a12(k,st)..a12(k,st)*power(w12(k,st),2)=e=d_th12(k,st)*w12(k,st)-th12(k,st)*d_w12(k,st);
eq_b12(k,st)..b12('s1',k,st)=e=sum(m,GS(m,'Q')*n('s1',m)*ps(m,k,st));
eq_c12(k,st)..c12(k,st)*w12(k,st)=e=d_w12(k,st);

*differ of residual activity coefficient, i.e. dlngr/dx(s1)
*simplier form, not the general form of unifac
eq_dlngr12(st)..d_lngr12('s1',st)=e=-sum(k,a12(k,st)*b12('s1',k,st)+GS(k,'Q')*n('s1',k)*c12(k,st));

*differ of activity coefficient,i.e. dlng/dx(s1)
eq_dlng12(st)..d_lng12('s1',st)=e=d_lngc12('s1',st)+d_lngr12('s1',st);

eq_misc12(st)..d_lng12('s1',st)*x12('s1',st)+1=g=0;

************************************INTEGER CUTS********************************
eq_cut1(ic)$dyn(ic)..                                           sum((ii,s),yv(ii,s,ic)*y(ii,s))-sum((ii,s),(1-yv(ii,s,ic))*y(ii,s))=l=sum((ii,s),yv(ii,s,ic))-1;


************************************BOUNDS************************************
*Bounds on optimisation variables
x.up(i,st)=1;
x.lo(jj,st)=1e-10;
x.lo('s1',st)=1e-8;
x.lo('s2',st)=0;

x12.up(ii,st)=0.90;
x12.lo('s1',st)=0.1;
x12.lo('s2',st)=0.1;

xi_s.lo=3.5;
xi_s.up = 200;

v_s.lo = 4;
v_s.up=500;

yield.lo(jj)=0.5;
yield.up(jj)=1;

T.up(st) =318.15;
T.lo(st) = 293.15;

mass.lo(jj,st)=1e-5;
mass.lo('s1',st)=1e-5;
mass.lo('s2',st)=0;
mass.up(i,st)=1000;

den.lo('s1')=0.5;
den.lo('s2')=0;
den.up(ii)=2;
vol.lo('s1')=1e-7;
vol.lo('s2')=0;
vol.up(ii)=2000;

*Bounds on other variables
*variable bounds
rc.lo(ii)= 0.000001;
rc.up(i)= 10;
qc.lo(ii)= 0.000001;
qc.up(i)= 10;

J.lo(jj,st)=0.000001;
J.up(jj,st)=10;
L.lo(jj,st)=0.000001;
L.up(jj,st)=10;

bib.up(i,k,st) = 5;
bib.lo(i,k,st) = -5;
th.up(k,st)=1;
th.lo(k,st)=0;
w.up(k,st)=3;
w.lo(k,st)=0;

lngc.lo(jj,st)=-30;
lngc.up(jj,st)=30;
lngr.lo(jj,st)=-30;
lngr.up(jj,st)=30;
lng.lo(jj,st) = -30;
lng.up(jj,st) = 30;

n.up(ii,k) = 10;

ps.lo(k,m,st) = 0;
ps.up(k,m,st) = 5;

nm.up(i,st) = 34;
nmt.up(st) = 45;

Mw.up(ii) = 110;
Mw.up(jj) = Mwib(jj);

*bounds for imiscibility calculations
*s1 & s2 bounds
*================================================================================
th12.up(k,st)=1;
th12.lo(k,st)=0;
w12.up(k,st)=5;
w12.lo(k,st)=0;

sumr12.lo(st)=0.000001;
sumr12.up(st)=10;
sumq12.lo(st)=0.000001;
sumq12.up(st)=10;

d_lngc12.lo('s1',st)=-30;
d_lngc12.up('s1',st)=30;

d_th12.lo(k,st)=-5;
d_th12.up(k,st)=5;

d_w12.lo(k,st)=-5;
d_w12.up(k,st)= 5;

a12.lo(k,st)=-5;
a12.up(k,st)= 5;

b12.up(ii,k,st)=5 ;

c12.lo(k,st)=-5;
c12.up(k,st)= 5;

d_lngr12.lo('s1',st)=-30;
d_lngr12.up('s1',st)=30;

d_lng12.lo('s1',st)=-30;
d_lng12.up('s1',st)=30;

*********************INITIAL POINTS*********************************************
qc.fx(jj)=qib(jj);
rc.fx(jj)=rib(jj);
n.fx(jj,k)=vib(jj,k);

*fixing initial amount of ibuprofen to be crystallised
nm.fx('ibuprofen','heated')=1;

*initilaization  of mole fraction
x.l('ibuprofen','heated')=0.13892;
x.l('ibuprofen','cooled')=3.67181e-5;


x.l('s2','heated')=0.342767;
x.l('s2','cooled')=0.899967;

x.l('s1',st)=1-x.l('ibuprofen',st)-x.l('s2',st);

nmt.l('heated')=7.19837;
nmt.l('cooled')=37.3115;
nm.l(i,st)=x.l(i,st)*nmt.l(st);

Yield.l(jj)=1-nm.l(jj,'cooled')/nm.l(jj,'heated');

y.l('s1','ethanol')=1;
y.l('s2','water')=1;

T.l('heated')=308.498;
T.l('cooled')=293.15;
n.l(ii,k)=sum(s,v(s,k)*y.l(ii,s));
Mw.l(ii)=sum(s,prop(s,'p3')*y.l(ii,s));
mass.l(i,st)=nm.l(i,st)*Mw.l(i);
den.l(ii)=sum(s,prop(s,'p4')*y.l(ii,s));
vol.l(ii)= mass.l(ii,'cooled')/den.l(ii);

rc.l(ii)= sum(s,rs(s)*y.l(ii,s));
qc.l(ii)= sum(s,qs(s)*y.l(ii,s));

J.l(jj,st)=rib(jj)/sum(i,x.l(i,st)*rc.l(i));
L.l(jj,st)=qib(jj)/sum(i,x.l(i,st)*qc.l(i));

lngc.l(jj,st)=1-J.l(jj,st)+log(J.l(jj,st))-5*qib(jj)*(1-J.l(jj,st)/L.l(jj,st)+log(J.l(jj,st))-log(L.l(jj,st)));

ps.l(k,m,st)=exp(-a(k,m)/T.l(st));
bib.l(jj,k,st) = sum(m,eib(jj,m)*ps.l(m,k,st));
th.l(k,st)=sum(i,x.l(i,st)*GS(k,'Q')*n.l(i,k))/sum(i,x.l(i,st)*qc.l(i));
w.l(k,st)=sum(m,th.l(m,st)*ps.l(m,k,st));

lngr.l(jj,st)=qib(jj)*(1-sum(k,th.l(k,st)*bib.l(jj,k,st)/w.l(k,st)-eib(jj,k)*(log(bib.l(jj,k,st))-log(w.l(k,st)))));

lng.l(jj,st) = lngc.l(jj,st) + lngr.l(jj,st);

*initialisation of imiscibility calculations(s1,s2)
x12.l('s1',st)=x.l('s1',st)/(x.l('s1',st)+x.l('s2',st));
x12.l('s2',st)=1-x12.l('s1',st);

sumr12.l(st)=x12.l('s1',st)*rc.l('s1')+(1-x12.l('s1',st))*rc.l('s2');
sumq12.l(st)=x12.l('s1',st)*qc.l('s1')+(1-x12.l('s1',st))*qc.l('s2');

d_lngc12.l('s1',st)=5*qc.l('s1')*(qc.l('s2')-qc.l('s1'))/sumq12.l(st)
                  +(rc.l('s1')-rc.l('s2'))*(5*qc.l('s1')-1)/sumr12.l(st)
                  +5*rc.l('s1')*(qc.l('s1')-qc.l('s2'))/sumr12.l(st)
                  +((rc.l('s1')-rc.l('s2'))*(rc.l('s1')-5*rc.l('s1')*(sumq12.l(st))))/power(sumr12.l(st),2);

th12.l(k,st)=(x12.l('s1',st)*GS(k,'Q')*n.l('s1',k)+(1-x12.l('s1',st))*GS(k,'Q')*n.l('s2',k))/sumq12.l(st);

w12.l(k,st)=sum(m,th12.l(m,st)*ps.l(m,k,st));

d_th12.l(k,st)=((GS(k,'Q')*n.l('s1',k)- GS(k,'Q')*n.l('s2',k))*sumq12.l(st)-(x12.l('s1',st)*GS(k,'Q')*n.l('s1',k)+ x12.l('s2',st)*GS(k,'Q')*n.l('s2',k))*(qc.l('s1')-qc.l('s2')))/power(sumq12.l(st),2);

d_w12.l(k,st)=sum(m,ps.l(m,k,st)*d_th12.l(m,st));

a12.l(k,st)=(d_th12.l(k,st)*w12.l(k,st)-th12.l(k,st)*d_w12.l(k,st))/power(w12.l(k,st),2);
b12.l(ii,k,st)=sum(m,GS(m,'Q')*n.l(ii,m)*ps.l(m,k,st));
c12.l(k,st)=d_w12.l(k,st)/w12.l(k,st);

d_lngr12.l('s1',st)=-sum(k,a12.l(k,st)*b12.l('s1',k,st)+GS(k,'Q')*n.l('s1',k)*c12.l(k,st));

d_lng12.l('s1',st)=d_lngc12.l('s1',st)+d_lngr12.l('s1',st);

model N2    /ALL/;
N2.optfile=1;
option nlp=snopt;
option mip=cplex;
option rminlp=snopt;
option minlp=sbb;

option decimals=7;
option reslim=3600;
*option reslim=300;
option iterlim=10000000;
option domlim=100;
option optca =1e-9;
option optcr =1e-6;
option  limrow=0;
option limcol=0;
option sysout=off;
N2.nodlim=10000000;

dyn(ic) = no;
Parameter stab12(st,ic);
yv(ii,s,ic)=0;
loop(ics,
        Solve N2 using MINLP maximising z;
        yv(ii,s,ics) = y.l(ii,s);
        yieldk(ics,jj) = yield.l(jj);
        xi_sk(ics) = xi_s.l;
        v_sk(ics) = v_s.l;
        Tk(ics,st) = T.l(st);
        yk(ics,ii,s) = y.l(ii,s);
        xk(ics,i,st) = x.l(i,st);
        nmk(ics,i,st) = nm.l(i,st);
        nmtk(ics,st) = nmt.l(st);
        stab12(st,ics) = d_lng12.l('s1',st)*x12.l('s1',st)+1;
        dyn(ics) = yes;
);

*=== export results to GDX file (occurs during execution phase)===*
execute_unload "Crystallisation_ibuprofen.gdx" yieldk,xi_sk,Tk,yk,xk,nmk,nmtk,stab12,v_sk
