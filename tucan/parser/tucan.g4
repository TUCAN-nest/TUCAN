grammar tucan;

/*
 * Tuples
 */
tuples : tuple* ;

// for testing only
tuples_start : tuples EOF ;

tuple : '(' node_index '-' node_index ')' ;

node_index : '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9' | GREATER_THAN_NINE ;
//node_index : ONE_TO_NINE | GREATER_THAN_NINE ;

/*
 * Hill system formula
 */
sum_formula : (with_carbon | without_carbon) ;

// for testing only
sum_formula_start : sum_formula EOF ;

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

//count : TWO_TO_NINE | GREATER_THAN_NINE ;
count : '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9' | GREATER_THAN_NINE ;

//TWO_TO_NINE : [2-9] ;
//ONE_TO_NINE : [1-9] ;
GREATER_THAN_NINE : [1-9] [0-9]+ ;
