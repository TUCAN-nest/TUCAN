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
      puts "\n#{'-' * 75}\n"
      puts "\nNow working on: ",molecule.name,"\n"
      puts "\nProcessing original molfile"\n"
      ninchi_string_original_molfile = molecule.ninchi_string(permute_molfile: false)
      puts "\nProcessing permuted molfile"\n"
      ninchi_string_permuted_molfile = molecule.ninchi_string(permute_molfile: true, random_seed: 181)
      assert_equal(ninchi_string_original_molfile, ninchi_string_permuted_molfile)
    end
  end
end
