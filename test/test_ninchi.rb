require './inchi'
require 'minitest/autorun'
require './periodic_table'

class TestnInChI < Minitest::Test
  @@periodic_table_elements = PeriodicTable::Elements
  @@periodic_table_colors = PeriodicTable::ElementColor

  hydrogen_molfile = read_molfile('test/testfiles/hydrogen.mol')
  hydrogen_array = create_molecule_array(hydrogen_molfile, @@periodic_table_elements, )
  @@canonicalized_hydrogen = canonicalize_molecule(hydrogen_array)
  @@ninchi_string_hydrogen = 'nInChI=1S/H2/c(0-1)'
  @@dot_file_string_hydrogen = "graph test\n{\n  bgcolor=grey\n"\
  "  0 [label=\"H 0\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
  "  1 [label=\"H 1\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
  "  0 -- 1 [color=black,style=bold];\n"\
  "}\n"

  cisplatin_molfile = read_molfile('test/testfiles/cisplatin.mol')
  cisplatin_array = create_molecule_array(cisplatin_molfile, @@periodic_table_elements, )
  @@canonicalized_cisplatin = canonicalize_molecule(cisplatin_array)
  @@ninchi_string_cisplatin = 'nInChI=1S/H6Cl2N2Pt/c(0-6)(1-6)(2-6)(3-7)(4-7)(5-7)(6-10)(7-10)(8-10)(9-10)'
  @@dot_file_string_cisplatin = "graph test\n{\n  bgcolor=grey\n"\
  "  0 [label=\"H 0\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
  "  1 [label=\"H 1\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
  "  2 [label=\"H 2\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
  "  3 [label=\"H 3\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
  "  4 [label=\"H 4\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
  "  5 [label=\"H 5\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
  "  6 [label=\"N 6\" color=\"#3050F8\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  7 [label=\"N 7\" color=\"#3050F8\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  8 [label=\"Cl 8\" color=\"#1FF01F\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  9 [label=\"Cl 9\" color=\"#1FF01F\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  10 [label=\"Pt 10\" color=\"#D0D0E0\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  0 -- 6 [color=black,style=bold];\n"\
  "  1 -- 6 [color=black,style=bold];\n"\
  "  2 -- 6 [color=black,style=bold];\n"\
  "  3 -- 7 [color=black,style=bold];\n"\
  "  4 -- 7 [color=black,style=bold];\n"\
  "  5 -- 7 [color=black,style=bold];\n"\
  "  6 -- 10 [color=black,style=bold];\n"\
  "  7 -- 10 [color=black,style=bold];\n"\
  "  8 -- 10 [color=black,style=bold];\n"\
  "  9 -- 10 [color=black,style=bold];\n"\
  "}\n"

  tpp_molfile = read_molfile('test/testfiles/tpp.mol')
  tpp_array = create_molecule_array(tpp_molfile, @@periodic_table_elements, )
  @@canonicalized_tpp = canonicalize_molecule(tpp_array)
  @@ninchi_string_tpp = 'nInChI=1S/C20N4/c(0-1)(0-12)(1-13)(2-3)(2-14)(3-15)(4-6)(4-16)(5-12)(5-16)(6-17)(7-14)(7-17)(8-10)(8-18)(9-13)(9-18)(10-19)(11-15)(11-19)(12-20)(13-20)(14-21)(15-21)(16-22)(17-22)(18-23)(19-23)'
  @@dot_file_string_tpp =  "graph test\n{\n  bgcolor=grey\n"\
  "  0 [label=\"C 0\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  1 [label=\"C 1\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  2 [label=\"C 2\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  3 [label=\"C 3\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  4 [label=\"C 4\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  5 [label=\"C 5\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  6 [label=\"C 6\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  7 [label=\"C 7\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  8 [label=\"C 8\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  9 [label=\"C 9\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  10 [label=\"C 10\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  11 [label=\"C 11\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  12 [label=\"C 12\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  13 [label=\"C 13\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  14 [label=\"C 14\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  15 [label=\"C 15\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  16 [label=\"C 16\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  17 [label=\"C 17\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  18 [label=\"C 18\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  19 [label=\"C 19\" color=\"#909090\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  20 [label=\"N 20\" color=\"#3050F8\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  21 [label=\"N 21\" color=\"#3050F8\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  22 [label=\"N 22\" color=\"#3050F8\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  23 [label=\"N 23\" color=\"#3050F8\",style=filled,shape=circle,fontname=Calibri];\n"\
  "  0 -- 1 [color=black,style=bold];\n"\
  "  0 -- 12 [color=black,style=bold];\n"\
  "  1 -- 13 [color=black,style=bold];\n"\
  "  2 -- 3 [color=black,style=bold];\n"\
  "  2 -- 14 [color=black,style=bold];\n"\
  "  3 -- 15 [color=black,style=bold];\n"\
  "  4 -- 6 [color=black,style=bold];\n"\
  "  4 -- 16 [color=black,style=bold];\n"\
  "  5 -- 12 [color=black,style=bold];\n"\
  "  5 -- 16 [color=black,style=bold];\n"\
  "  6 -- 17 [color=black,style=bold];\n"\
  "  7 -- 14 [color=black,style=bold];\n"\
  "  7 -- 17 [color=black,style=bold];\n"\
  "  8 -- 10 [color=black,style=bold];\n"\
  "  8 -- 18 [color=black,style=bold];\n"\
  "  9 -- 13 [color=black,style=bold];\n"\
  "  9 -- 18 [color=black,style=bold];\n"\
  "  10 -- 19 [color=black,style=bold];\n"\
  "  11 -- 15 [color=black,style=bold];\n"\
  "  11 -- 19 [color=black,style=bold];\n"\
  "  12 -- 20 [color=black,style=bold];\n"\
  "  13 -- 20 [color=black,style=bold];\n"\
  "  14 -- 21 [color=black,style=bold];\n"\
  "  15 -- 21 [color=black,style=bold];\n"\
  "  16 -- 22 [color=black,style=bold];\n"\
  "  17 -- 22 [color=black,style=bold];\n"\
  "  18 -- 23 [color=black,style=bold];\n"\
  "  19 -- 23 [color=black,style=bold];\n"\
  "}\n"

  def test_hydrogen_ninchi_string
    assert_equal @@ninchi_string_hydrogen, write_ninchi_string(@@canonicalized_hydrogen, @@periodic_table_elements)
  end

  def test_cisplatin_ninchi_string
    assert_equal @@ninchi_string_cisplatin, write_ninchi_string(@@canonicalized_cisplatin, @@periodic_table_elements)
  end

  def test_tpp_ninchi_string
    assert_equal @@ninchi_string_tpp, write_ninchi_string(@@canonicalized_tpp, @@periodic_table_elements)
  end

  def test_hydrogen_dot_file_string
    assert_equal @@dot_file_string_hydrogen, write_dot_file(@@canonicalized_hydrogen, @@periodic_table_elements, @@periodic_table_colors)
  end

  def test_cisplatin_dot_file_string
    assert_equal @@dot_file_string_cisplatin, write_dot_file(@@canonicalized_cisplatin, @@periodic_table_elements, @@periodic_table_colors)
  end

  def test_tpp_dot_file_string
    assert_equal @@dot_file_string_tpp, write_dot_file(@@canonicalized_tpp, @@periodic_table_elements, @@periodic_table_colors)
  end
end
