require 'minitest/autorun'
require './inchi'
require './periodic_table'


class Molecule
  include Inchi
  include PeriodicTable
  attr_reader :reference_ninchi_string

  def initialize(name, ninchi_string)
    @reference_ninchi_string = ninchi_string
    molfile = read_molfile("test/testfiles/#{name}/#{name}.mol")
    @molecule_array = create_molecule_array(molfile, PeriodicTable::ELEMENTS)
    @canonicalized_molecule = canonicalize_molecule(@molecule_array)
  end

  def ninchi_string
    write_ninchi_string(@canonicalized_molecule, PeriodicTable::ELEMENTS)
  end
end

class MoleculeCollection
  def initialize
    @hydrogen = Molecule.new(
      'hydrogen',
      'nInChI=1S/H2/c(0-1)'
    )

    @cisplatin = Molecule.new(
      'cisplatin',
      'nInChI=1S/H6Cl2N2Pt/c(0-6)(1-6)(2-6)(3-7)(4-7)(5-7)(6-10)(7-10)(8-10)(9-10)'
    )

    @tpp = Molecule.new(
      'tpp',
      'nInChI=1S/C20N4/c(0-1)(0-12)(1-13)(2-3)(2-14)(3-15)(4-6)(4-16)(5-12)(5-16)(6-17)(7-14)(7-17)(8-10)(8-18)(9-13)(9-18)(10-19)(11-15)(11-19)(12-20)(13-20)(14-21)(15-21)(16-22)(17-22)(18-23)(19-23)',
      )
  end

  def molecules
    [@hydrogen, @cisplatin, @tpp]
  end
end

class TestSuite < Minitest::Test
  # Metaprogramming shenanigans for test parameterization (as in pytest)
  # inspired by https://stackoverflow.com/questions/18770988/.
  MoleculeCollection.new.molecules.each do |molecule|
    define_method("test_ninchi_string_#{molecule}") do
      assert_equal molecule.reference_ninchi_string, molecule.ninchi_string
    end
  end
end
