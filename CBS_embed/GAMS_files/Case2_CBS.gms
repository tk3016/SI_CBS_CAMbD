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

prp                 'melting/boiling T (pi, p2); molar mass (p3); density (p4)'
/p1*p4/

pr                  'UNIFAC contributions'
/Q,R/

ic                  'Integer cuts'
/1*1/

f                   'Inputs to NN'
*/1*7/

hl1                 'Hidden layer'
hl2                 'Hidden layer'
hl3                 'Hidden layer'
*/1*5/
;

*Load the appropriate gdx file here
$gdxin Case2_CBS3.gdx
$load f=f
$load hl1=hl1
$load hl2=hl2
$load hl3=hl3
$gdxin


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

Parameters
input_offset(f)  'vector of input offsets'
$gdxin Case2_slle.gdx
$load input_offset=input_offset
$gdxin

input_gain(f)    'vector of input gains'
$gdxin Case2_slle.gdx
$load input_gain=input_gain
$gdxIn

bias1(hl1)         'bias for hidden neurons in layer 1'
$gdxin Case2_slle.gdx
$load bias1=bias1
$gdxIn

bias2(hl2)
$gdxin Case2_slle.gdx
$load bias2=bias2
$gdxIn

bias3(hl3)
$gdxin Case2_slle.gdx
$load bias3=bias3
$gdxIn

bias4
$gdxin Case2_slle.gdx
$load bias4=bias4
$gdxIn

wt4(hl3)              'weight vector for neurons in layer 4'
$gdxin Case2_slle.gdx
$load wt4=wt4
$gdxIn

wt3(hl2,hl3)              'weight vector for neurons in layer 3'
$gdxin Case2_slle.gdx
$load wt3=wt3
$gdxIn

wt2(hl1,hl2)              'weight vector for neurons in layer 2'
$gdxin Case2_slle.gdx
$load wt2=wt2
$gdxIn

table wt1(f,hl1)    'weight matrix for layer 1'
$gdxin Case2_slle.gdx
$load wt1=wt1
$gdxin

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

Parameter eib(jj,k);
eib(jj,k) = vib(jj,k)*GS(k,'Q')/qib(jj);

Scalar  Tos         'Minimum temperature offset'
/10/
;
Positive Variables

yield(jj)           'Crystal yield of lovastatin'
xi_s                'Solvent use'
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

misc12

inp(f,st)
inter1v(hl1,st)
inter2v(hl2,st)
inter3v(hl3,st)
inter4v(st)
inter1z(hl1,st)
inter2z(hl2,st)
inter3z(hl3,st)
inter4z(st)
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

*Equations for miscibility constraints
eq_x12a(st)
eq_x12b(st)
eq_dlng12(st)
eq_misc12(st)

*Equations for interger cuts
eq_cut1(ic)

*Equations for neural network
eq_inp1(st)
eq_inp2(st)
eq_inp3(st)
eq_inp4(st)
eq_inp5(st)
eq_inp6(st)
eq_inp7(st)
eq_inter1v(hl1,st)
eq_inter1z(hl1,st)
eq_inter2v(hl2,st)
eq_inter2z(hl2,st)
eq_inter3v(hl3,st)
eq_inter3z(hl3,st)
eq_inter4v(st)
eq_inter4z(st)
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
eq_x12a(st)..                           x12('s1',st)*(x('s1',st)+x('s2',st))=e=x('s1',st);
eq_x12b(st)..                           x12('s1',st)+x12('s2',st)=e=1;

*Thresholding constraints
eq_misc12(st)..                         misc12(st)=g=0.95;

************************************INTEGER CUTS********************************
eq_cut1(ic)$dyn(ic)..                                           sum((ii,s),yv(ii,s,ic)*y(ii,s))-sum((ii,s),(1-yv(ii,s,ic))*y(ii,s))=l=sum((ii,s),yv(ii,s,ic))-1;



*******Neural Network Implementation*******************************************
eq_inp1(st)..                                                   inp('1',st) =e= input_gain('1')*(qc('s1') - input_offset('1'));
eq_inp2(st)..                                                   inp('2',st) =e= input_gain('2')*(qc('s2') - input_offset('2'));
eq_inp3(st)..                                                   inp('3',st) =e= input_gain('3')*(rc('s1') - input_offset('3'));
eq_inp4(st)..                                                   inp('4',st) =e= input_gain('4')*(rc('s2') - input_offset('4'));
eq_inp5(st)..                                                   inp('5',st) =e= input_gain('5')*(x('s1',st) - input_offset('5'));
eq_inp6(st)..                                                   inp('6',st) =e= input_gain('6')*(x('s2',st) - input_offset('6'));
eq_inp7(st)..                                                   inp('7',st) =e= input_gain('7')*(T(st) - input_offset('7'));

eq_inter1v(hl1,st)..                                            inter1v(hl1,st) =e= sum(f, wt1(f,hl1)*inp(f,st)) + bias1(hl1);
eq_inter1z(hl1,st)..                                            inter1z(hl1,st) =e= 1 - (2/(exp(2*inter1v(hl1,st)) + 1));
eq_inter2v(hl2,st)..                                            inter2v(hl2,st) =e= sum(hl1, wt2(hl1,hl2)*inter1z(hl1,st)) + bias2(hl2);
eq_inter2z(hl2,st)..                                            inter2z(hl2,st) =e= 1 - (2/(exp(2*inter2v(hl2,st)) + 1));
eq_inter3v(hl3,st)..                                            inter3v(hl3,st) =e= sum(hl2, wt3(hl2,hl3)*inter2z(hl2,st)) + bias3(hl3);
eq_inter3z(hl3,st)..                                            inter3z(hl3,st) =e= 1 - (2/(exp(2*inter3v(hl3,st)) + 1));
eq_inter4v(st)..                                                inter4v(st) =e= sum(hl3, wt4(hl3)*inter3z(hl3,st)) + bias4;
eq_inter4z(st)..                                                inter4z(st) =e= 1/(1 + exp(-inter4v(st)));
eq_dlng12(st)..                                                 misc12(st) =e= inter4z(st);


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
inp.lo(f,st) = -15;
inp.up(f,st) = 15;

inter4v.lo(st) = -70;
inter4v.up(st) = 70;

inter3v.lo(hl3,st) = -70;
inter3v.up(hl3,st) = 70;

inter2v.lo(hl2,st) = -70;
inter2v.up(hl2,st) = 70;

inter1v.lo(hl1,st) = -50;
inter1v.up(hl1,st) = 50;

inter4z.lo(st) = 0;
inter4z.up(st) = 1;

inter3z.lo(hl3,st) = -1;
inter3z.up(hl3,st) = 1;

inter2z.lo(hl2,st) = -1;
inter2z.up(hl2,st) = 1;

inter1z.lo(hl1,st) = -1;
inter1z.up(hl1,st) = 1;

misc12.lo(st) = 0;
misc12.up(st) = 1;

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
Mw.l(jj) = Mwib(jj);
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

x12.l('s1',st) = x.l('s1',st)/(x.l('s1',st) + x.l('s2',st));
x12.l('s2',st) = 1 - x12.l('s1',st);

*initialisation of imiscibility calculations(s1,s2)


inp.l('1',st) = input_gain('1')*(qc.l('s1') - input_offset('1'));
inp.l('2',st) = input_gain('2')*(qc.l('s2') - input_offset('2'));
inp.l('3',st) = input_gain('3')*(rc.l('s1') - input_offset('3'));
inp.l('4',st) = input_gain('4')*(rc.l('s2') - input_offset('4'));
inp.l('5',st) = input_gain('5')*(x.l('s1',st) - input_offset('5'));
inp.l('6',st) = input_gain('6')*(x.l('s2',st) - input_offset('6'));
inp.l('7',st) = input_gain('7')*(T.l(st) - input_offset('7'));

inter1v.l(hl1,st) = sum(f, wt1(f,hl1)*inp.l(f,st)) + bias1(hl1);
inter1z.l(hl1,st) = 1 - (2/(exp(2*inter1v.l(hl1,st))+1));
inter2v.l(hl2,st) = sum(hl1, wt2(hl1,hl2)*inter1z.l(hl1,st)) + bias2(hl2);
inter2z.l(hl2,st) = 1 - (2/(exp(2*inter2v.l(hl2,st))+1));
inter3v.l(hl3,st) = sum(hl2, wt3(hl2,hl3)*inter2z.l(hl2,st)) + bias3(hl3);
inter3z.l(hl3,st) = 1 - (2/(exp(2*inter3v.l(hl3,st))+1));
inter4v.l(st) = sum(hl3, wt4(hl3)*inter3z.l(hl3,st)) + bias4;
inter4z.l(st) = 1/(1 + exp(-inter4v.l(st)));
misc12.l(st) = inter4z.l(st);

model N2    /ALL/;
N2.optfile=0;
option nlp=snopt;
option mip=cplex;
option rminlp=snopt;
option minlp=antigone;

option decimals=7;
option reslim=5759;
option iterlim=10000000;
option domlim=10000;
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
        stab12(st,ics) = misc12.l(st);
        dyn(ics) = yes;
);

*=== export results to GDX file (occurs during execution phase)===*
execute_unload "Crystallisation_ibuprofen_classify.gdx" yieldk,xi_sk,Tk,yk,xk,nmk,nmtk,stab12,v_sk
