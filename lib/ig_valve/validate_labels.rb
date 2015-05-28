require 'ig_valve/label_prediction'
require 'ig_valve/compare_abs'
require 'multiset'

module IgValve

  class ValidateLabels
    #
    include CompareAbs
    #
    #
    def initialize(pred_f, truth_f, delim=",")
      @pred_file = pred_f
      @truth_file = truth_f
      @delim = delim
      @truth_h = load_truth()
    end

    #
    #
    #
    def load_truth()
      h = {}
      IO.readlines(@truth_file).each do |l|
        t_lab = LabelPrediction.new(l, @delim)
        id = t_lab.get_id
        h[id] = t_lab
      end
      h
    end

    #
    private :load_truth
    #
    # yields, for each prediction, a predicted and truth LabelPrediction object
    #
    def each_pred()
      num = 0
      IO.readlines(@pred_file).each do |l|
        num += 1
        pred_lab = LabelPrediction.new(l, @delim)
        id = pred_lab.get_id
        truth_lab = @truth_h[id]
        yield [pred_lab, truth_lab]
      end
      num
    end

    #
    # yields, for each pair of predictions, the LabelPrediction objects
    #
    def each_pair_pred()
      i = -1
      each_pred do |pred_lab_i, truth_lab_i|
        i += 1
        j = -1
        each_pred do |pred_lab_j, truth_lab_j|
          j+=1
          next if i >= j # skip lower triangle of pairwise matrix
          yield pred_lab_i, truth_lab_i, pred_lab_j, truth_lab_j
        end

      end
    end

    #
    #
    #
    def validate_all(cmp_method)
      #acc_h = {V: 0, D: 0, J: 0, total: 0}
      acc_h = {V: 0, D: 0, J: 0, total: 0, CDR3: 0}
      num = each_pred do |pred_lab, truth_lab|
        next if pred_lab.nil? || truth_lab.nil?
        h = cmp_method.call(pred_lab, truth_lab)
        # puts [pred_lab.get_id,
        #        [:V, :D, :J,:total].map{|seg| h[seg]},
        #        pred_lab.get_d_alleles.to_s
        #       ].flatten.join("\t")
        #[:V, :D, :J,:total, :CDR3].each do |seg|
        acc_h.each_key do |seg|
          #acc_h[seg] += (h[seg] ? 1 : 0)
          acc_h[seg] += h[seg] if h.has_key?(seg)
        end
      end
      #[:V, :D, :J, :total, :CDR3]
      acc_h.keys.each { |seg| acc_h[seg] /= num.to_f }
      acc_h
    end

    #
    #
    #
    def get_gene_accuracy()
      self.validate_all(method(:cmp_genes_i))
    end

    #
    #
    #
    def get_allele_accuracy()
      self.validate_all(method(:cmp_alleles_i))
    end

    #
    #
    #
    def get_partition_accuracy()
      self.validate_all(method(:cmp_partitions))
    end

    #
    # compare CDR3 sequence at ID level
    #
    def get_cdr3_accuracy()
      self.validate_all(method(:cmp_cdr3))
    end

    #
    # compare clones as comparing multisets across the dataset using Tanimoto similarity,
    # T(A,B) = \frac{A \dot B}{|A|^2+|B|^2-A \dot B}
    #
    # def get_clone_accuracy()
    #   truth_ms = Multiset.new
    #   @truth_h.each_key do |k|
    #     truth_ms.add(@truth_h[k].get_cdr3)
    #   end
    #   pred_ms = Multiset.new
    #   IO.readlines(@pred_file).each do |l|
    #     num += 1
    #     pred_lab = LabelPrediction.new(l, @delim)
    #     pred_ms.add(pred_lab.get_cdr3)
    #   end
    #   self.compare_clones(pred_ms, truth_ms)
    # end
    #
    # compute FM index on CDR3 partitions
    #
    def get_clone_accuracy()
      fowlkes_mallows_index(@truth_file, @pred_file)
    end

    #
    # compute Rand index over cdr3 partitions
    #
    # def get_clone_accuracy_2()
    #   a = 0
    #   b = 0
    #   c = 0
    #   d = 0
    #   num = 0
    #   each_pair_pred do |p1, t1, p2, t2|
    #     p_cdr3 = (p1.get_cdr3 == p2.get_cdr3)
    #     t_cdr3 = (t1.get_cdr3 == t2.get_cdr3)
    #     num += 1
    #     if p_cdr3 && t_cdr3 then
    #       a += 1
    #     elsif !p_cdr3 && !t_cdr3 then
    #       b += 1
    #     elsif p_cdr3 && !t_cdr3 then
    #       c += 0
    #     elsif !p_cdr3 && t_cdr3 then
    #       d += 0
    #     else
    #     end
    #   end
    #   denom = (a + b + c + d).to_f
    #   #puts [a, b, c, d, denom, num].join("\t")
    #   {:clone => (a+b)/denom}
    # end

  end

end
