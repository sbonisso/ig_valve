require_relative 'test_helper'
require_relative 'load_helper'
require 'ig_valve/evaluate_labels'

class TestEvaluateLabels < MiniTest::Test
  #
  #
  #
  def test_validate_genes_one_err()
    vl = IgValve::EvaluateLabels.new(["#{get_data_dir}/one_err_preds.tab",
                                    "#{get_data_dir}/perf_preds.tab",
                                    "#{get_data_dir}/perf_preds.tab"],
                                   ["one_err", "perf1", "perf2"],
                                   "\t")
    h = vl.get_gene_similarity
    #
    v1 = {:V=>0.9, :D=>0.9, :J=>0.9, :total => 0.7, :CDR3 => 0}
    v2 = {:V=>1.0, :D=>1.0, :J=>1.0, :total => 1.0, :CDR3 => 0}
    assert_equal(v1, h['one_err']['perf1'], "one_err - perf1 incorrect")
    assert_equal(v1, h['one_err']['perf2'], "one_err - perf2 incorrect")
    assert_equal(v1, h['perf1']['one_err'], "perf1 - one_err incorrect")
    assert_equal(v2, h['perf1']['perf2'], "perf1 - perf2 incorrect")
  end
  #
  #
  #
  def test_validate_clones()
    vl = IgValve::EvaluateLabels.new(["#{get_data_dir}/ten_abs_preds1.tab",
                                    "#{get_data_dir}/ten_abs_preds2.tab"],
                                   ['pred1', 'pred2'])
    vh = vl.get_clone_similarity('rand')
    sim = vh['pred1']['pred2'][:clone]
    #assert_in_delta( (9.0/11), sim, 0.001, 'clone similarity too different')
    assert_in_delta( 1.0, sim, 0.001, 'clone similarity too different for rand')
  end
  #
  #
  #
  def test_validate_cdr3()
    vl = IgValve::EvaluateLabels.new(["#{get_data_dir}/ten_abs_preds1.tab",
                                    "#{get_data_dir}/ten_abs_preds2.tab"],
                                   ['pred1', 'pred2'])
    vh = vl.get_cdr3_similarity
    sim = vh['pred1']['pred2'][:CDR3]
    assert_in_delta( (9.0/10), sim, 0.001, 'cdr3 similarity too different')
  end

end
