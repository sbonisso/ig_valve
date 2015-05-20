require_relative 'test_helper'
require_relative 'load_helper'
require 'ig_valve/label_prediction'

class TestLabelPrediction < MiniTest::Test
  #
  #
  #
  def test_get_genes()
    v_genes = ["IGHV4-28", "IGHV6-1", "IGHV4-59", "IGHV3-11", "IGHV4-59", "IGHV2-70", "IGHV2-5", "IGHV2-5", "IGHV4-4", "IGHV3-74"]
    d_genes = ["IGHD2-2", "IGHD2-2", "IGHD2-15", "IGHD4-4", "IGHD2-2", "IGHD3-3", "IGHD3-9", "IGHD6-19", "IGHD3-3", "IGHD3-16"]
    j_genes = ["IGHJ4", "IGHJ6", "IGHJ4", "IGHJ3", "IGHJ3", "IGHJ5", "IGHJ6", "IGHJ6", "IGHJ1", "IGHJ3"]
    i = 0
    IO.foreach("#{get_data_dir}/truth.tab") do |l|
      lp = IgValve::LabelPrediction.new(l.chomp)
      assert_equal(v_genes[i], lp.get_v_genes[0], "err V-gene #{i}")
      assert_equal(d_genes[i], lp.get_d_genes[0], "err D-gene #{i}")
      assert_equal(j_genes[i], lp.get_j_genes[0], "err J-gene #{i}")
      i+=1
    end
  end
  #
  #
  #
  def test_get_alleles()
    v_alleles = ["IGHV4-28*05", "IGHV6-1*02", "IGHV4-59*10", "IGHV3-11*04", "IGHV4-59*03", "IGHV2-70*13", "IGHV2-5*04", "IGHV2-5*07", "IGHV4-4*01", "IGHV3-74*03"]    
    d_alleles = ["IGHD2-2*02", "IGHD2-2*01", "IGHD2-15*01", "IGHD4-4*01", "IGHD2-2*03", "IGHD3-3*01", "IGHD3-9*01", "IGHD6-19*01", "IGHD3-3*01", "IGHD3-16*02"]
    j_alleles = ["IGHJ4*02", "IGHJ6*01", "IGHJ4*02", "IGHJ3*01", "IGHJ3*01", "IGHJ5*01", "IGHJ6*04", "IGHJ6*02", "IGHJ1*01", "IGHJ3*01"]
    
    i = 0
    IO.foreach("#{get_data_dir}/truth.tab") do |l|
      lp = IgValve::LabelPrediction.new(l.chomp)
      assert_equal(v_alleles[i], lp.get_v_alleles[0], "err V-allele #{i}")
      assert_equal(d_alleles[i], lp.get_d_alleles[0], "err D-allele #{i}")
      assert_equal(j_alleles[i], lp.get_j_alleles[0], "err J-allele #{i}")
      i+=1
    end
  end  

end
