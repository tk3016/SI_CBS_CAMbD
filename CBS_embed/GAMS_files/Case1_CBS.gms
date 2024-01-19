*Title: 'Select 2 solvents to maximise the solubility of Ibuprofen
*N=2
*Solved with DICOPT
*Miscibility constraints are included
*==============================================================================
SETS

s 'solvents'
/methanol, ethanol, 2-propanol, acetone, mibk, ethylacetate,
chloroform, toluene, water/

i 'components in mixture'
/ibuprofen, s1, s2/

ii(i) 'selected solvent'
/s1, s2/


k 'groups'
/ch3, ch2, ch, ach, acch3, acch2, acch, oh, ch3oh, h2o, ch3co, ch3coo, cooh, chcl3/

f 'inputs to NN'

hl1 'hidden layer 1'

c 'integer cuts'
/1*5/

dyn(c) 'dynamic set of c'

*Load the appropriate gdx file for the neural network model
$gdxin Case1_CBS3.gdx
$load f=f
$load hl1=hl1
$gdxin

;

alias(k,m);
alias(k,p);
alias(s,ss);


PARAMETERS
Hm 'enthalpy of fusion of ibuprofen' /25500/
Tm 'melting point of ibuprofen' /347.15/
T 'system temperature' /300/
Rg 'gas constant' /8.3144/
temp 'NN calibration parameter' /1/


R(k) 'Van der Waals volume of group k'
/ch3             0.9011
ch2              0.6744
ch               0.4469
ach              0.5313
acch3            1.2663
acch2            1.0396
acch             0.8121
oh               1.0000
ch3oh            1.4311
h2o              0.9200
ch3co            1.6724
ch3coo           1.9031
cooh             1.3013
chcl3            2.8700 /

Q(k) 'Van der Waals surface area of group k'
/ch3             0.848
ch2              0.540
ch               0.228
ach              0.400
acch3            0.968
acch2            0.660
acch             0.348
oh               1.200
ch3oh            1.432
h2o              1.400
ch3co            1.488
ch3coo           1.728
cooh             1.224
chcl3            2.410 /;

*Change the gdx file to change the neural network model
Parameters
input_offset(f)  'vector of input offsets'
$gdxin Case1_N2_slle.gdx
$load input_offset=input_offset
$gdxin

input_gain(f)    'vector of input gains'
$gdxin Case1_N2_slle.gdx
$load input_gain=input_gain
$gdxIn

bias1(hl1)         'bias for hidden neurons in layer 1'
$gdxin Case1_N2_slle.gdx
$load bias1=bias1
$gdxIn

bias2
$gdxin Case1_N2_slle.gdx
$load bias2=bias2
$gdxIn

wt2(hl1)              'weight vector for neurons in layer 2'
$gdxin Case1_N2_slle.gdx
$load wt2=wt2
$gdxIn

table wt1(f,hl1)    'weight matrix for layer 1'
$gdxin Case1_N2_slle.gdx
$load wt1=wt1
$gdxin
;

table v(s,k) 'number of groups in molecule i'
                 ch3     ch2     ch      ach     acch3   acch2   acch    oh      ch3oh   h2o     ch3co   ch3coo  cooh    chcl3

methanol                                                                         1
ethanol          1       1                                               1
2-propanol       2               1                                       1
acetone          1                                                                               1
mibk             2       1       1                                                               1
ethylacetate     1       1                                                                               1
chloroform                                                                                                               1
toluene                                  5       1
water                                                                                    1
;

table vib(i,k)'identity of ibuprofen'
             ch3     ch ach      acch2 acch   cooh
ibuprofen     3       1  4         1    1       1      ;


                                                                                                                          ;
*===============================================================================
*group interaction parameters
*===============================================================================
table a(k,m) 'group interaction parameters'
         ch3     ch2     ch      ach     acch3   acch2   acch    oh      ch3oh   h2o     ch3co   ch3coo  cooh    chcl3
ch3      0       0       0       61.13   76.5    76.5    76.5    986.5   697.2   1318    476.4   232.1   663.5   24.9
ch2      0       0       0       61.13   76.5    76.5    76.5    986.5   697.2   1318    476.4   232.1   663.5   24.9
ch       0       0       0       61.13   76.5    76.5    76.5    986.5   697.2   1318    476.4   232.1   663.5   24.9
ach      -11.12  -11.12  -11.12  0       167     167     167     636.1   637.4   903.8   25.77   5.994   537.4   -231.9
acch3    -69.7   -69.7   -69.7   -146.8  0       0       0       803.2   603.3   5695    -52.1   5688    872.3   -80.25
acch2    -69.7   -69.7   -69.7   -146.8  0       0       0       803.2   603.3   5695    -52.1   5688    872.3   -80.25
acch     -69.7   -69.7   -69.7   -146.8  0       0       0       803.2   603.3   5695    -52.1   5688    872.3   -80.25
oh       156.4   156.4   156.4   89.6    25.82   25.82   25.82   0       -137.1  353.5   84      101.1   199     -98.12
ch3oh    16.51   16.51   16.51   -50     -44.5   -44.5   -44.5   249.1   0       -181    23.39   -10.72  -202    -139.4
h2o      300     300     300     362.3   377.6   377.6   377.6   -229.1  289.6   0       -195.4  72.87   -14.09  353.7
ch3co    26.76   26.76   26.76   140.1   365.8   365.8   365.8   164.5   108.7   472.5   0       -213.7  669.4   -354.6
ch3coo   114.8   114.8   114.8   85.84   -170    -170    -170    245.4   249.6   200.8   372.2   0       660.2   -209.7
cooh     315.3   315.3   315.3   62.32   89.86   89.86   89.86   -151    339.8   -66.17  -297.8  -256.3  0       39.63
chcl3    36.7    36.7    36.7    288.5   69.9    69.9    69.9    742.1   649.1   826.8   552.1   176.5   504.2   0
                                                                                                                    ;



parameter ps(k,m);
ps(k,m)= exp(-a(k,m)/T);
parameter nib(i,k);
nib('ibuprofen',k)= vib('ibuprofen',k);
parameter qib;
qib = sum(k,nib('ibuprofen',k)*Q(k));
parameter rib;
rib = sum(k,nib('ibuprofen',k)*R(k));
parameter qs(s);
qs(s) = sum(k,v(s,k)*Q(k));
parameter rs(s);
rs(s) = sum(k,v(s,k)*R(k));
parameter eib(i,k);
eib('ibuprofen',k)=vib('ibuprofen',k)*Q(k)/qib;
parameter bib(i,k);
bib('ibuprofen',k)=sum(m,eib('ibuprofen',m)*ps(m,k));
parameter yv(ii,s,c) 'store y values from previous iterations';
parameter zv(c) 'store objective values from previous iterations';
parameter xv(i,c) 'store mole fractions from previous iterations';


POSITIVE VARIABLES

x(i)             'liquid phase mole fraction of component i'
rc(i)            'van der waals volume of component i'
qc(i)            'van der waals area of component i'
J(i),L(i)

b(i,k)          'intermediates'
th(k)
w(k)

*xm(ii)
th12(k)
*w12(k)
*sumr12
*sumq12
*b12(ii,k)

FREE VARIABLES
z            'objective function'

lng(i)       'natural log of activity coefficient of ibuprofen'
lngc(i)      'natural log of combinatorial activity coefficient of ibuprofen'
lngr(i)      'natural log of residual activity coefficient of ibuprofen'

*d_lngc12(ii)  'differention of combinatorial act_coeff'
*d_lngr12(ii)  'differention of residual act_coeff'
*d_lng12(ii)   'differentiation of act_coef wrt to x(ii)'
*d_th12(k)     'diff of intermediate'
*d_w12(k)      'diff of intermediate'
*a12(k)     'intermediate'
*c12(k)     'intermediate'
misc12
inp(f)           'inputs to the neural network'
inter1v(hl1)     'pre-activated values'
inter1z(hl1)	 'activated values'
inter2v			 'pre-activated values for 2nd layer'	
inter2z			 'activated values for 2nd layer'


BINARY VARIABLES
y(ii,s)     'selected solvents'

INTEGER VARIABLE
n(i,k)    'number of groups k in component i' ;


EQUATIONS
eq_z                   'objective function maximise solubility of ibuprofen'
eq_x                   'mole fraction constrain'
eq_n(ii,k)
eq_qc(ii)
eq_rc(ii)
eq_J,eq_L
eq_lngc
eq_th(k)
eq_w(k)
eq_lngr
eq_lng
eq_solub
eq_inp1
eq_inp2
eq_inp3
eq_inp4
eq_inp5
eq_inp6
*eq_inp7

eq_inter1v(hl1)
eq_inter1z(hl1)
eq_inter2v
eq_inter2z
*equations for logical conditions
logic1(ii),logic2(s)
logic3(s,ss)
*equations for miscibility constraints
*eq_xm1
eq_dlng12,eq_misc12
IntCut(c)
;



*********************choose only 2 solvents*************************************
logic1(ii)..
sum(s,y(ii,s))=e=1;

**********************each solvent at most once*********************************
logic2(s)..y('s1',s)+y('s2',s)=l=1;

*****************************solvent ordering***********************************
logic3(s,ss)$(ord(ss) <ord(s))..y('s1',s)+y('s2',s-ord(ss))=l=1;



*******************************UNIFAC*******************************************
eq_z..
z=e=x('ibuprofen');

eq_x..
sum(i,x(i))=e=1;

eq_n(ii,k)..
n(ii,k)=e=sum(s,v(s,k)*y(ii,s));


eq_qc(ii)..
qc(ii)=e= sum(s,qs(s)*y(ii,s));

eq_rc(ii)..
rc(ii)=e= sum(s,rs(s)*y(ii,s));

*combinatorial part of activity coefficient
eq_J..J('ibuprofen')*sum(i,x(i)*rc(i))=e=rib;
eq_L..L('ibuprofen')*sum(i,x(i)*qc(i))=e=qib;

eq_lngc..
lngc('ibuprofen')=e=1-J('ibuprofen')+log(J('ibuprofen'))-5*qib*(1-J('ibuprofen')/L('ibuprofen')+log(J('ibuprofen'))-log(L('ibuprofen')));

*residual part of activity coefficient
eq_th(k)..th(k)*sum(i,x(i)*qc(i))=e=sum(i,x(i)*Q(k)*n(i,k));
eq_w(k)..w(k)=e=sum(m,th(m)*ps(m,k));

eq_lngr..
lngr('ibuprofen')=e=qib*(1-sum(k,th(k)*bib('ibuprofen',k)/w(k)-eib('ibuprofen',k)*(log(bib('ibuprofen',k))-log(w(k)))));

*activity coefficient
eq_lng..
lng('ibuprofen') =e=  lngc('ibuprofen') + lngr('ibuprofen');

*ibuprofen solubility constraint
eq_solub..
log(x('ibuprofen'))+ lng('ibuprofen') =e= Hm/Rg *(1/Tm-1/T);

*******Neural Network Implementation*******************************************
eq_inp1..                                                   inp('1') =e= input_gain('1')*(qc('s1') - input_offset('1'));
eq_inp2..                                                   inp('2') =e= input_gain('2')*(qc('s2') - input_offset('2'));
eq_inp3..                                                   inp('3') =e= input_gain('3')*(rc('s1') - input_offset('3'));
eq_inp4..                                                   inp('4') =e= input_gain('4')*(rc('s2') - input_offset('4'));
eq_inp5..                                                   inp('5') =e= input_gain('5')*(x('s1') - input_offset('5'));
eq_inp6..                                                   inp('6') =e= input_gain('6')*(x('s2') - input_offset('6'));
*eq_inp7..                                                   inp('7') =e= input_gain('7')*(xm('s1') - input_offset('7'));

eq_inter1v(hl1)..  inter1v(hl1) =e= sum(f, wt1(f,hl1)*inp(f)) + bias1(hl1);
eq_inter1z(hl1)..  inter1z(hl1) =e= 1 - (2/(exp(2*inter1v(hl1))+1));
eq_inter2v..     inter2v =e= sum(hl1, wt2(hl1)*inter1z(hl1)) + bias2;
eq_inter2z..     inter2z =e= 1/(1 + exp(-inter2v/temp));
eq_dlng12..     misc12 =e= inter2z;

*Thresholding constraint
eq_misc12..misc12=g=0.5;


*Integer cuts equation
IntCut(c)$(dyn(c)).. sum((ii,s),yv(ii,s,c)*y(ii,s))-sum((ii,s),(1-yv(ii,s,c))*y(ii,s))=l=sum((ii,s),yv(ii,s,c))-1;

*******************************BOUNDS*******************************************
*bounds on optmization variables
x.up(i)=1;
x.lo(i)=0.01;

n.up(ii,k)=6;

*variable bounds
rc.lo(i)= 0.1;
rc.up(i)= 10;
qc.lo(i)= 0.1;
qc.up(i)= 10;

J.lo('ibuprofen')=0.1;
J.up('ibuprofen')=10;
L.lo('ibuprofen')=0.1;
L.up('ibuprofen')=10;

th.up(k)=1;
th.lo(k)=0;
w.up(k)=3;
w.lo(k)=0;

lngc.lo('ibuprofen')=-30;
lngc.up('ibuprofen')=30;
lngr.lo('ibuprofen')=-30;
lngr.up('ibuprofen')=30;
lng.lo('ibuprofen') = -30;
lng.up('ibuprofen') = 30;

*bounds for imiscibility calculations
*s1 & s2 bounds
*================================================================================
inp.lo(f) = -10;
inp.up(f) = 10;

inter2v.lo = -20;
inter2v.up = 20;

inter2z.lo = 0;
inter2z.up = 1;

inter1v.lo(hl1) = -20;
inter1v.up(hl1) = 20;

inter1z.lo(hl1) = -1;
inter1z.up(hl1) = 1;

misc12.lo = 0;
misc12.up = 1;

*********************INITIAL POINTS*********************************************

*fixing of ibuprofen identity
qc.fx('ibuprofen')=qib;
rc.fx('ibuprofen')=rib;
n.fx('ibuprofen',k)=vib('ibuprofen',k);

*initilaization  of optimal solvent
x.l('ibuprofen')=0.31623;
x.l('s2')=0.01;
x.l('s1')=1-x.l('ibuprofen')-x.l('s2');

y.l('s1','methanol')=1;
y.l('s2','chloroform')=1;

z.l=x.l('ibuprofen');

n.l(ii,k)=sum(s,v(s,k)*y.l(ii,s));

rc.l(ii)= sum(s,rs(s)*y.l(ii,s));
qc.l(ii)= sum(s,qs(s)*y.l(ii,s));

J.l('ibuprofen')=rib/sum(i,x.l(i)*rc.l(i));
L.l('ibuprofen')=qib/sum(i,x.l(i)*rc.l(i));

lngc.l('ibuprofen')=1-J.l('ibuprofen')+log(J.l('ibuprofen'))-5*qib*(1-J.l('ibuprofen')/L.l('ibuprofen')+log(J.l('ibuprofen'))-log(L.l('ibuprofen')));

th.l(k)=sum(i,x.l(i)*Q(k)*n.l(i,k))/sum(i,x.l(i)*qc.l(i));
w.l(k)=sum(m,th.l(m)*ps(m,k));

lngr.l('ibuprofen')=qib*(1-sum(k,th.l(k)*bib('ibuprofen',k)/w.l(k)-eib('ibuprofen',k)*(log(bib('ibuprofen',k))-log(w.l(k)))));

lng.l('ibuprofen') = lngc.l('ibuprofen') + lngr.l('ibuprofen');

*initialisation of imiscibility calculations(s1,s2)


inp.l('1') = input_gain('1')*(qc.l('s1') - input_offset('1'));
inp.l('2') = input_gain('2')*(qc.l('s2') - input_offset('2'));
inp.l('3') = input_gain('3')*(rc.l('s1') - input_offset('3'));
inp.l('4') = input_gain('4')*(rc.l('s2') - input_offset('4'));
inp.l('5') = input_gain('5')*(x.l('s1') - input_offset('5'));
inp.l('6') = input_gain('6')*(x.l('s2') - input_offset('6'));

inter1v.l(hl1) = sum(f, wt1(f,hl1)*inp.l(f)) + bias1(hl1);
inter1z.l(hl1) = 1 - (2/(exp(2*inter1v.l(hl1))+1));
inter2v.l = sum(hl1, wt2(hl1)*inter1z.l(hl1)) + bias2;
inter2z.l = 1/(1 + exp(-inter2v.l));
misc12.l = inter2z.l;

MODEL N2 / ALL/;
N2.optfile=0;
option nlp=conopt3;
option mip=cplex;
option rminlp=conopt3;
option minlp=baron;
*option threads=4;

option decimals=5;
OPTION OPTCA = 1e-10;
option iterlim = 1000000;
option optcr  = 1e-5;
OPTION reslim = 3600;

option sysout=on;
dyn(c) = no;
yv(ii,s,c) = 0;
alias (c,cc);
parameter
stab12(c)     'stabbility criterion'
inv1(c);


loop(cc,
    solve N2 using minlp maximising z;
    yv(ii,s,cc) = y.l(ii,s);
    zv(cc) = z.l;
    xv(i,cc) = x.l(i);
    stab12(cc)= misc12.l;
    dyn(cc) = yes;


);


display
zv
yv
xv
stab12
;
