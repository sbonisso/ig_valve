require_relative 'test_helper'
require_relative 'load_helper'
require 'ig_valve/validate_labels'

class TestValidateLabels < MiniTest::Test
  #
  #
  #
  def test_validate_genes_perfect()
    vl = IgValve::ValidateLabels.new("#{get_data_dir}/perf_preds.tab",
                            "#{get_data_dir}/truth.tab")
    h = vl.get_gene_accuracy
    [:V, :D, :J].each{|seg| assert_equal(1.0, h[seg], "identity gene error")}
  end
  #
  def test_validate_alleles_perfect()
    vl = IgValve::ValidateLabels.new("#{get_data_dir}/perf_preds.tab",
                            "#{get_data_dir}/truth.tab")
    h = vl.get_allele_accuracy
    [:V, :D, :J].each{|seg| assert_equal(1.0, h[seg], "identity gene error")}
  end
  #
  def test_validate_partitions_perfect()
    vl = IgValve::ValidateLabels.new("#{get_data_dir}/perf_preds.tab",
                            "#{get_data_dir}/truth.tab")
    h = vl.get_partition_accuracy
    [:V, :D, :J].each{|seg| assert_equal(1.0, h[seg], "identity gene error")}
  end
  #
  #
  #
  def test_validate_genes_one_err()
    vl = IgValve::ValidateLabels.new("#{get_data_dir}/one_err_preds.tab",
                            "#{get_data_dir}/truth.tab")
    h = vl.get_gene_accuracy
    [:V, :D, :J].each{|seg| assert_equal(0.9, h[seg], "one err gene error")}
  end
  #
  def test_validate_allele_one_err()
    vl = IgValve::ValidateLabels.new("#{get_data_dir}/one_err_preds.tab",
                            "#{get_data_dir}/truth.tab")
    h = vl.get_allele_accuracy
    [:V, :D, :J].each{|seg| assert_equal(0.9, h[seg], "one err gene error")}
  end
  #
  def test_validate_partitions_one_err()
    vl = IgValve::ValidateLabels.new("#{get_data_dir}/one_err_preds.tab",
                            "#{get_data_dir}/truth.tab")
    h = vl.get_partition_accuracy
    assert_in_delta(0.999652777777, h[:V], 0.0001, "one err partitions")
    assert_in_delta(0.995454545454, h[:D], 0.0001, "one err partitions")
    assert_in_delta(0.998305084745, h[:J], 0.0001, "one err partitions")
  end
  #
  def test_validate_clone_partition()
    vl = IgValve::ValidateLabels.new("#{get_data_dir}/ten_abs_preds1.tab",
                            "#{get_data_dir}/ten_abs_preds2.tab")
    sim = vl.get_clone_accuracy_2[:clone]
    assert_in_delta( 1.0, sim, 0.001, 'clone similarity too different')
  end
  
end
