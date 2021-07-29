require "./inchi_v1.4"
require 'minitest/autorun'


class TestnInChI < Minitest::Test

    hydrogen_molfile = read_molfile("test/testfiles/hydrogen.mol")
    hydrogen_array = create_molecule_array(hydrogen_molfile, "hydrogen")
    @@canonicalized_hydrogen = canonicalize_molecule(hydrogen_array)
    @@ninchi_string_hydrogen = "nInChI=1S/H2/c(0-1)"
    @@dot_file_string_hydrogen = "graph test\n{\n  bgcolor=grey\n"\
    "  0 [label=\"H 0\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
    "  1 [label=\"H 1\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
    "  0 -- 1 [color=black,style=bold];\n"\
    "}\n"

    cisplatin_molfile = read_molfile("test/testfiles/cisplatin.mol")
    cisplatin_array = create_molecule_array(cisplatin_molfile, "cisplation")
    @@canonicalized_cisplatin = canonicalize_molecule(cisplatin_array)
    @@ninchi_string_cisplatin = "nInChI=1S/H6Cl2N2Pt/c(0-6)(1-6)(2-6)(3-7)(4-7)(5-7)(6-10)(7-10)(8-10)(9-10)"
    @@dot_file_string_cisplatin = "graph test\n{\n  bgcolor=grey\n"\
    "  0 [label=\"H 0\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
    "  1 [label=\"H 1\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
    "  2 [label=\"H 2\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
    "  3 [label=\"H 3\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
    "  4 [label=\"H 4\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
    "  5 [label=\"H 5\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
    "  6 [label=\"N 6\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
    "  7 [label=\"N 7\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
    "  8 [label=\"Cl 8\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
    "  9 [label=\"Cl 9\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
    "  10 [label=\"Pt 10\" color=lightgrey,style=filled,shape=circle,fontname=Calibri];\n"\
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


    def test_hydrogen_ninchi_string
        assert_equal @@ninchi_string_hydrogen, create_ninchi_string(@@canonicalized_hydrogen)
    end


    def test_cisplatin_ninchi_string
        assert_equal @@ninchi_string_cisplatin, create_ninchi_string(@@canonicalized_cisplatin)
    end


    def test_hydrogen_dot_file_string
        assert_equal @@dot_file_string_hydrogen, create_dot_file(@@canonicalized_hydrogen)
    end


    def test_cisplatin_dot_file_string
        assert_equal @@dot_file_string_cisplatin, create_dot_file(@@canonicalized_cisplatin)
    end

end