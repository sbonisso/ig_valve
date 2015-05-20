require 'ig_valve/label_prediction'
require 'ig_valve/compare_abs'
require 'ig_valve/validate_labels'

module IgValve

  class EvaluateLabels
    #
    include CompareAbs
    #
    #
    def initialize(pred_lst, name_lst, delim=",")
      @pred_lst = pred_lst
      @name_lst = name_lst
      @h = {}
      @delim = delim
      #each_pair(@name_lst.size){|i,j| @h[[i,j]] = {}}
    end

    #
    # yields a block for each pair of index in list_len
    #
    def each_pair(lst_len)
      (0..lst_len-1).each do |i|
        (0..lst_len-1).each do |j|
          next if i >= j
          yield [i, j]
        end
      end
    end

    #
    # return a hash of the names of pairs of tools
    #
    def get_name_hash()
      p_h = {}
      @name_lst.each do |ns|
        p_h[ns] = {}
        @name_lst.each { |ns_nest| p_h[ns][ns_nest] = nil }
      end
      p_h
    end

    #
    private :get_name_hash
    #
    # yields a ValidateLabels object for each pair of tools
    #
    def generic_pair_similarity()
      # init
      p_h = get_name_hash
      #
      each_pair(@name_lst.size) do |i, j|
        n1 = @name_lst[i]
        n2 = @name_lst[j]
        #
        #puts [@pred_lst[i], @pred_lst[j]].join("\t")
        vl = ValidateLabels.new(@pred_lst[i], @pred_lst[j], @delim)
        val = yield vl
        p_h[n1][n2] = val
        p_h[n2][n1] = val
      end
      p_h
    end

    #
    #
    #
    def get_gene_similarity()
      generic_pair_similarity { |vl| vl.get_gene_accuracy() }
    end

    #
    #
    #
    def get_allele_similarity()
      generic_pair_similarity { |vl| vl.get_allele_accuracy() }
    end

    #
    #
    #
    def get_partition_similarity()
      generic_pair_similarity { |vl| vl.get_partition_accuracy() }
    end

    #
    #
    #
    def get_cdr3_similarity()
      generic_pair_similarity { |vl| vl.get_cdr3_accuracy() }
    end

    #
    # return Multiset of CDR3 from a prediction file
    #
    def get_cdr3_multiset(pred_f)
      pred_ms = Multiset.new
      IO.readlines(pred_f).each do |l|
        pred_lab = LabelPrediction.new(l, @delim)
        pred_ms.add(pred_lab.get_cdr3)
      end
      pred_ms
    end

    #
    private :get_cdr3_multiset
    #
    #
    #
    def get_clone_similarity_old()
      p_h = get_name_hash
      each_pair(@name_lst.size) do |i, j|
        n1 = @name_lst[i]
        n2 = @name_lst[j]
        #
        pred1_ms = get_cdr3_multiset(@pred_lst[i])
        pred2_ms = get_cdr3_multiset(@pred_lst[j])
        #
        val = self.jaccard_multi(pred1_ms, pred2_ms)
        p_h[n1][n2] = val
        p_h[n2][n1] = val
      end
      p_h
    end

    #
    #
    #
    def get_clone_similarity(index='fm')
      p_h = get_name_hash
      #
      each_pair(@name_lst.size) do |i, j|
        n1 = @name_lst[i]
        n2 = @name_lst[j]
        val = if index == 'fm'
                fowlkes_mallows_index(@pred_lst[i], @pred_lst[j])
              elsif index == 'rand'
                rand_index(@pred_lst[i], @pred_lst[j])
              elsif index == 'jaccard'
                jaccard_index(@pred_lst[i], @pred_lst[j])
              else
                fowlkes_mallows_index(@pred_lst[i], @pred_lst[j])
              end
        p_h[n1][n2] = val
      end
      p_h
    end

  end
end
