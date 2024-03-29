grammar tucan;

/*
 * TUCAN
 */
tucan : sum_formula '/' tuples ('/' node_attributes)? EOF ;

/*
 * Hill system formula
 */
sum_formula_start : sum_formula EOF ; // for testing only

sum_formula : with_carbon | without_carbon ;

with_carbon : c h? ac? ag? al? am? ar? as? at? au? b? ba? be? bh? bi? bk? br? ca? cd? ce? cf? cl? cm? cn? co? cr? cs? cu? db? ds? dy? er? es? eu? f? fe? fl? fm? fr? ga? gd? ge? he? hf? hg? ho? hs? i? in? ir? k? kr? la? li? lr? lu? lv? mc? md? mg? mn? mo? mt? n? na? nb? nd? ne? nh? ni? no? np? o? og? os? p? pa? pb? pd? pm? po? pr? pt? pu? ra? rb? re? rf? rg? rh? rn? ru? s? sb? sc? se? sg? si? sm? sn? sr? ta? tb? tc? te? th? ti? tl? tm? ts? u? v? w? xe? y? yb? zn? zr? ;

without_carbon : ac? ag? al? am? ar? as? at? au? b? ba? be? bh? bi? bk? br? ca? cd? ce? cf? cl? cm? cn? co? cr? cs? cu? db? ds? dy? er? es? eu? f? fe? fl? fm? fr? ga? gd? ge? h? he? hf? hg? ho? hs? i? in? ir? k? kr? la? li? lr? lu? lv? mc? md? mg? mn? mo? mt? n? na? nb? nd? ne? nh? ni? no? np? o? og? os? p? pa? pb? pd? pm? po? pr? pt? pu? ra? rb? re? rf? rg? rh? rn? ru? s? sb? sc? se? sg? si? sm? sn? sr? ta? tb? tc? te? th? ti? tl? tm? ts? u? v? w? xe? y? yb? zn? zr? ;

// chemical elements ordered by atomic number
h  : 'H'  count? ;
he : 'He' count? ;
li : 'Li' count? ;
be : 'Be' count? ;
b  : 'B'  count? ;
c  : 'C'  count? ;
n  : 'N'  count? ;
o  : 'O'  count? ;
f  : 'F'  count? ;
ne : 'Ne' count? ;
na : 'Na' count? ;
mg : 'Mg' count? ;
al : 'Al' count? ;
si : 'Si' count? ;
p  : 'P'  count? ;
s  : 'S'  count? ;
cl : 'Cl' count? ;
ar : 'Ar' count? ;
k  : 'K'  count? ;
ca : 'Ca' count? ;
sc : 'Sc' count? ;
ti : 'Ti' count? ;
v  : 'V'  count? ;
cr : 'Cr' count? ;
mn : 'Mn' count? ;
fe : 'Fe' count? ;
co : 'Co' count? ;
ni : 'Ni' count? ;
cu : 'Cu' count? ;
zn : 'Zn' count? ;
ga : 'Ga' count? ;
ge : 'Ge' count? ;
as : 'As' count? ;
se : 'Se' count? ;
br : 'Br' count? ;
kr : 'Kr' count? ;
rb : 'Rb' count? ;
sr : 'Sr' count? ;
y  : 'Y'  count? ;
zr : 'Zr' count? ;
nb : 'Nb' count? ;
mo : 'Mo' count? ;
tc : 'Tc' count? ;
ru : 'Ru' count? ;
rh : 'Rh' count? ;
pd : 'Pd' count? ;
ag : 'Ag' count? ;
cd : 'Cd' count? ;
in : 'In' count? ;
sn : 'Sn' count? ;
sb : 'Sb' count? ;
te : 'Te' count? ;
i  : 'I'  count? ;
xe : 'Xe' count? ;
cs : 'Cs' count? ;
ba : 'Ba' count? ;
la : 'La' count? ;
ce : 'Ce' count? ;
pr : 'Pr' count? ;
nd : 'Nd' count? ;
pm : 'Pm' count? ;
sm : 'Sm' count? ;
eu : 'Eu' count? ;
gd : 'Gd' count? ;
tb : 'Tb' count? ;
dy : 'Dy' count? ;
ho : 'Ho' count? ;
er : 'Er' count? ;
tm : 'Tm' count? ;
yb : 'Yb' count? ;
lu : 'Lu' count? ;
hf : 'Hf' count? ;
ta : 'Ta' count? ;
w  : 'W'  count? ;
re : 'Re' count? ;
os : 'Os' count? ;
ir : 'Ir' count? ;
pt : 'Pt' count? ;
au : 'Au' count? ;
hg : 'Hg' count? ;
tl : 'Tl' count? ;
pb : 'Pb' count? ;
bi : 'Bi' count? ;
po : 'Po' count? ;
at : 'At' count? ;
rn : 'Rn' count? ;
fr : 'Fr' count? ;
ra : 'Ra' count? ;
ac : 'Ac' count? ;
th : 'Th' count? ;
pa : 'Pa' count? ;
u  : 'U'  count? ;
np : 'Np' count? ;
pu : 'Pu' count? ;
am : 'Am' count? ;
cm : 'Cm' count? ;
bk : 'Bk' count? ;
cf : 'Cf' count? ;
es : 'Es' count? ;
fm : 'Fm' count? ;
md : 'Md' count? ;
no : 'No' count? ;
lr : 'Lr' count? ;
rf : 'Rf' count? ;
db : 'Db' count? ;
sg : 'Sg' count? ;
bh : 'Bh' count? ;
hs : 'Hs' count? ;
mt : 'Mt' count? ;
ds : 'Ds' count? ;
rg : 'Rg' count? ;
cn : 'Cn' count? ;
nh : 'Nh' count? ;
fl : 'Fl' count? ;
mc : 'Mc' count? ;
lv : 'Lv' count? ;
ts : 'Ts' count? ;
og : 'Og' count? ;

count : greater_than_one ;

/*
 * Tuples
 */
tuples_start : tuples EOF ; // for testing only

tuples : tuple* ;

tuple : '(' node_index '-' node_index ')' ;

node_index : greater_than_zero ;

/*
 * Node attributes
 */
node_attributes_start : node_attributes EOF ; // for testing only

node_attributes : node_attribute* ;

node_attribute : '(' node_index ':' node_property (',' node_property)* ')' ;

node_property : node_property_key '=' node_property_value ;

node_property_key : 'mass' | 'rad' ;

node_property_value : greater_than_zero ;

/*
 * Numbers
 */
greater_than_zero : '1' | greater_than_one ;
greater_than_one : '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9' | GREATER_THAN_NINE ;
GREATER_THAN_NINE : [1-9] [0-9]+ ;
