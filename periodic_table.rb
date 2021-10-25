# (c) CC BY-SA Ulrich Schatzschneider, Universität Würzburg and NFDI4Chem 2021

# Periodic table of the elements containing all 118 elements from hydrogen to oganesson.
# List of element symbols organized by increasing atomic number (number of protons).
# Colors assigned according to Jmol as taken from http://jmol.sourceforge.net/jscolors with the exception of H, which was set to lightgrey for better visibility.
# Colors beyond meitnerium are not offically assigned so far and were thus set to the same value as for Mt.

module PeriodicTable
  ELEMENTS = %w[H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar
                K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
                Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe
                Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
                Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og].freeze
  ELEMENT_COLORS = { 'H' => 'lightgrey', 'He' => '"#D9FFFF"', 'Li' => '"#CC80FF"', 'Be' => '"#C2FF00"', 'B' => '"#FFB5B5"', 'C' => '"#909090"',
                     'N' => '"#3050F8"', 'O' => '"#FF0D0D"', 'F' => '"#90E050"', 'Ne' => '"#B3E3F5"',
                     'Na' => '"#AB5CF2"', 'Mg' => '"#8AFF00"', 'Al' => '"#BFA6A6"', 'Si' => '"#F0C8A0"', 'P' => '"#FF8000"', 'S' => '"#FFFF30"', 'Cl' => '"#1FF01F"', 'Ar' => '"#80D1E3"',
                     'K' => '"#8F40D4"', 'Ca' => '"#3DFF00"', 'Sc' => '"#E6E6E6"', 'Ti' => '"#BFC2C7"', 'V' => '"#A6A6AB"', 'Cr' => '"#8A99C7"', 'Mn' => '"#9C7AC7"', 'Fe' => '"#E06633"', 'Co' => '"#F090A0"', 'Ni' => '"#50D050"', 'Cu' => '"#C88033"', 'Zn' => '"#7D80B0"',
                     'Ga' => '"#C28F8F"', 'Ge' => '"#668F8F"', 'As' => '"#BD80E3"', 'Se' => '"#FFA100"', 'Br' => '"#A62929"', 'Kr' => '"#5CB8D1"',
                     'Rb' => '"#702EB0"', 'Sr' => '"#00FF00"', 'Y' => '"#94FFFF"', 'Zr' => '"#94E0E0"', 'Nb' => '"#73C2C9"', 'Mo' => '"#54B5B5"', 'Tc' => '"#3B9E9E"', 'Ru' => '"#248F8F"', 'Rh' => '"#0A7D8C"', 'Pd' => '"#006985"', 'Ag' => '"#C0C0C0"', 'Cd' => '"#FFD98F"',
                     'In' => '"#A67573"', 'Sn' => '"#668080"', 'Sb' => '"#9E63B5"', 'Te' => '"#D47A00"', 'I' => '"#940094"', 'Xe' => '"#429EB0"',
                     'Cs' => '"#57178F"', 'Ba' => '"#00C900"',
                     'La' => '"#70D4FF"', 'Ce' => '"#FFFFC7"', 'Pr' => '"#D9FFC7"', 'Nd' => '"#C7FFC7"', 'Pm' => '"#A3FFC7"', 'Sm' => '"#8FFFC7"', 'Eu' => '"#61FFC7"', 'Gd' => '"#45FFC7"', 'Tb' => '"#30FFC7"', 'Dy' => '"#1FFFC7"', 'Ho' => '"#00FF9C"', 'Er' => '"#00E675"', 'Tm' => '"#00D452"', 'Yb' => '"#00BF38"',
                     'Lu' => '"#00AB24"', 'Hf' => '"#4DC2FF"', 'Ta' => '"#4DA6FF"',  'W' => '"#2194D6"', 'Re' => '"#267DAB"', 'Os' => '"#266696"', 'Ir' => '"#175487"', 'Pt' => '"#D0D0E0"', 'Au' => '"#FFD123"', 'Hg' => '"#B8B8D0"', 'Tl' => '"#A6544D"', 'Pb' => '"#575961"', 'Bi' => '"#9E4FB5"', 'Po' => '"#AB5C00"', 'At' => '"#754F45"', 'Rn' => '"#428296"',
                     'Fr' => '"#420066"', 'Ra' => '"#007D00"', 'Ac' => '"#70ABFA"', 'Th' => '"#00BAFF"', 'Pa' => '"#00A1FF"',  'U' => '"#008FFF"', 'Np' => '"#0080FF"', 'Pu' => '"#006BFF"', 'Am' => '"#545CF2"', 'Cm' => '"#785CE3"', 'Bk' => '"#8A4FE3"', 'Cf' => '"#A136D4"', 'Es' => '"#B31FD4"', 'Fm' => '"#B31FBA"', 'Md' => '"#B30DA6"', 'No' => '"#BD0D87"', 'Lr' => '"#C70066"',
                     'Rf' => '"#CC0059"', 'Db' => '"#D1004F"', 'Sg' => '"#D90045"', 'Bh' => '"#E00038"', 'Hs' => '"#E6002E"', 'Mt' => '"#EB0026"', 'Ds' => '"#EB0026"', 'Rg' => '"#EB0026"', 'Cn' => '"#EB0026"', 'Nh' => '"#EB0026"', 'Fl' => '"#EB0026"', 'Mc' => '"#EB0026"', 'Lv' => '"#EB0026"', 'Ts' => '"#EB0026"', 'Og' => '"#EB0026"' }.freeze
end
