require_relative 'test_helper'
require_relative 'load_helper'
require 'ig_valve/compare_abs'

class TestCompareAbs < MiniTest::Test

  include IgValve::CompareAbs
  #
  # test for identity for comparing gene predictions
  #
  def test_cmp_genes_identity()
    truth_p = IO.readlines("#{get_data_dir}/truth.tab").map do |l|
      IgValve::LabelPrediction.new(l)
    end
    pred_p = IO.readlines("#{get_data_dir}/perf_preds.tab").map do |l|
      IgValve::LabelPrediction.new(l)
    end
    (0..truth_p.size-1).each do |i|
      retval = cmp_genes(truth_p[i], pred_p[i])
      [:V, :D, :J].each{|seg| assert_equal(retval[seg], true)}
    end
  end
  #
  #
  #
  def test_cmp_alleles_identity()
    truth_p = IO.readlines("#{get_data_dir}/truth.tab").map do |l|
      IgValve::LabelPrediction.new(l)
    end
    pred_p = IO.readlines("#{get_data_dir}/perf_preds.tab").map do |l|
      IgValve::LabelPrediction.new(l)
    end
    (0..truth_p.size-1).each do |i|
      retval = cmp_alleles(truth_p[i], pred_p[i])
      [:V, :D, :J].each{|seg| assert_equal(retval[seg], true)}
    end
  end
  #
  #
  #
  def test_cmp_partitions_identity()
    truth_p = IO.readlines("#{get_data_dir}/truth.tab").map do |l|
      IgValve::LabelPrediction.new(l)
    end
    pred_p = IO.readlines("#{get_data_dir}/perf_preds.tab").map do |l|
      IgValve::LabelPrediction.new(l)
    end
    (0..truth_p.size-1).each do |i|
      retval = cmp_partitions(truth_p[i], pred_p[i])
      [:V, :D, :J].each{|seg| assert_equal(retval[seg], 1.0)}
    end
  end
  #
  #
  #
  def test_jaccard_num()
    a = (1..10).to_a
    b = (7..15).to_a
    assert_in_delta((4.0/15), jaccard(a,b), 0.0001, "jaccard for nums err 1")
    a.push(11)
    assert_in_delta((5.0/15), jaccard(a,b), 0.0001, "jaccard for nums err 2")
  end
  #
  #
  #
  def test_jaccard_str()
    a = ["a", "b", "c", "d", "e"]
    b = ["c", "d", "e", "f", "g", "h"]
    assert_in_delta((3.0/8), jaccard(a,b), 0.0001, "jaccard for str err 1")
    a.push("f")
    assert_in_delta((4.0/8), jaccard(a,b), 0.0001, "jaccard for str err 2")
  end
  
end
