require 'minitest/autorun'
require './test/utils'

class MoleculeCollection
  def initialize
    @hydrogen = Molecule.new(
      'test/testfiles/hydrogen/hydrogen.mol',
      "nInChI=1S/H2/c(0-1)"
    )

    @cisplatin = Molecule.new(
      'test/testfiles/cisplatin/cisplatin.mol',
      "nInChI=1S/H6Cl2N2Pt/c(0-6)(1-6)(2-6)(3-7)(4-7)(5-7)(6-10)(7-10)(8-10)(9-10)"    )

    @tpp = Molecule.new(
      'test/testfiles/tpp/tpp.mol',
      "nInChI=1S/C20H14N4/c(0-14)(1-15)(2-16)(3-17)(4-18)(5-19)(6-20)(7-21)(8-22)(9-23)(10-24)(11-25)(12-36)(13-37)(14-15)(14-27)(15-29)(16-17)(16-26)(17-28)(18-19)(18-30)(19-32)(20-26)(20-27)(21-22)(21-31)(22-33)(23-28)(23-30)(24-29)(24-31)(25-32)(25-33)(26-34)(27-36)(28-34)(29-36)(30-37)(31-35)(32-37)(33-35)"    )
  end

  def molecules
    [@hydrogen, @cisplatin, @tpp]
  end
end

class RegressionTests < Minitest::Test
  # Metaprogramming shenanigans for test parameterization (as in pytest)
  # inspired by https://stackoverflow.com/questions/18770988/.
  MoleculeCollection.new.molecules.each do |molecule|
    define_method("test_regression_ninchi_string_with_#{molecule.name}") do
      assert_equal molecule.reference_ninchi_string, molecule.ninchi_string(permute_molfile: false)
    end
  end
end
