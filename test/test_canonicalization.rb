require 'minitest/autorun'
require './test/utils'

class MoleculeCollection
  def initialize
    @molecule_paths = Dir['test/testfiles/*/*.mol']
  end

  def molecules
    @molecule_paths.map { |path| Molecule.new(path) }
  end
end

class CanonicalizationValidationTests < Minitest::Test
  # Metaprogramming shenanigans for test parameterization (as in pytest)
  # inspired by https://stackoverflow.com/questions/18770988/.
  MoleculeCollection.new.molecules.each do |molecule|
    define_method("test_canonicalization_#{molecule.name}") do
      ninchi_strings = []
      10.times do
        ninchi_strings.push(molecule.ninchi_string(permute_molfile: false))
      end
      assert ninchi_strings.uniq.size == 1
    end
  end
end
