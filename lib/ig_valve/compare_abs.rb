require 'ig_valve/label_prediction'
require 'ig_valve/utils'
require 'multiset'
require 'tempfile'
require 'fileutils'
require 'open3'


module IgValve
#
# module for mix-ins for comparing two VDJ predictions
#
  module CompareAbs
    #
    #
    def is_label_pred?(obj)
      obj.class == LabelPrediction
    end

    #
    # check if gene predictions are the same for
    #
    def cmp_genes(pred1, pred2)
      raise 'not a LabelPrediction object' if !is_label_pred?(pred1)
      {:V => (pred1.get_v_genes & pred2.get_v_genes).size > 0,
       :D => (pred1.get_d_genes & pred2.get_d_genes).size > 0,
       :J => (pred1.get_j_genes & pred2.get_j_genes).size > 0,
       :total => (((pred1.get_v_genes & pred2.get_v_genes).size > 0) &
           ((pred1.get_d_genes & pred2.get_d_genes).size > 0) &
           ((pred1.get_j_genes & pred2.get_j_genes).size > 0))}
    end

    #
    def cmp_genes_i(pred1, pred2)
      raise 'not a LabelPrediction object' if !is_label_pred?(pred1)
      # {:V => ((pred1.get_v_genes & pred2.get_v_genes).size > 0 ? 1 : 0),
      #   :D => ((pred1.get_d_genes & pred2.get_d_genes).size > 0 ? 1 : 0),
      #   :J => ((pred1.get_j_genes & pred2.get_j_genes).size > 0 ? 1 : 0),
      #   :total => (((pred1.get_v_genes & pred2.get_v_genes).size > 0 ? 1 : 0) &
      #              ((pred1.get_d_genes & pred2.get_d_genes).size > 0 ? 1 : 0) &
      #              ((pred1.get_j_genes & pred2.get_j_genes).size > 0 ? 1 : 0))
      # }
      th = cmp_genes(pred1, pred2)
      {:V => (th[:V] ? 1 : 0),
       :D => (th[:D] ? 1 : 0),
       :J => (th[:J] ? 1 : 0),
       :total => (th[:total] ? 1 : 0)
      }
    end

    #
    #
    #
    def cmp_alleles(pred1, pred2)
      raise 'not a LabelPrediction object' if !is_label_pred?(pred1)
      {:V => (pred1.get_v_alleles & pred2.get_v_alleles).size > 0,
       :D => (pred1.get_d_alleles & pred2.get_d_alleles).size > 0,
       :J => (pred1.get_j_alleles & pred2.get_j_alleles).size > 0,
       :total => (((pred1.get_v_alleles & pred2.get_v_alleles).size > 0) &
           ((pred1.get_d_alleles & pred2.get_d_alleles).size > 0) &
           ((pred1.get_j_alleles & pred2.get_j_alleles).size > 0))
      }
    end

    #
    def cmp_alleles_i(pred1, pred2)
      raise 'not a LabelPrediction object' if !is_label_pred?(pred1)
      # {:V => ((pred1.get_v_alleles & pred2.get_v_alleles).size > 0 ? 1 : 0),
      #   :D => ((pred1.get_d_alleles & pred2.get_d_alleles).size > 0 ? 1 : 0),
      #   :J => ((pred1.get_j_alleles & pred2.get_j_alleles).size > 0 ? 1 : 0),
      #   :total => (((pred1.get_v_alleles & pred2.get_v_alleles).size > 0 ? 1 : 0) &
      #              ((pred1.get_d_alleles & pred2.get_d_alleles).size > 0 ? 1 : 0) &
      #              (((pred1.get_j_alleles & pred2.get_j_alleles).size > 0 ? 1 : 0)))
      # }
      th = cmp_alleles(pred1, pred2)
      #puts th
      {:V => (th[:V] ? 1 : 0),
       :D => (th[:D] ? 1 : 0),
       :J => (th[:J] ? 1 : 0),
       :total => (th[:total] ? 1 : 0)
      }
    end

    #
    #
    #
    def cmp_partitions(pred1, pred2)
      {:V => jaccard(pred1.get_v_partition,
                     pred2.get_v_partition),
       :D => jaccard(pred1.get_d_partition,
                     pred2.get_d_partition),
       :J => jaccard(pred1.get_j_partition,
                     pred2.get_j_partition),
       :total => jaccard(pred1.get_v_partition |
                             pred1.get_d_partition |
                             pred1.get_j_partition,
                         pred2.get_v_partition |
                             pred2.get_d_partition |
                             pred2.get_j_partition)
      }
    end

    #
    #
    #
    def cmp_cdr3(pred1, pred2)
      raise 'not a LabelPrediction object' if !is_label_pred?(pred1)
      {:CDR3 => ((pred1.get_cdr3 == pred2.get_cdr3) ? 1 : 0)}
    end

    #
    # compare clones as comparing multisets across the dataset using Tanimoto similarity,
    #
    def compare_clones(pred_lbls, truth_lbls)
      truth_ms = Multiset.new
      truth_lbls.each_key do |k|
        truth_ms.add(truth_ms[k].get_cdr3)
      end
      pred_ms = Multiset.new
      pred_lbls.each_key do |k|
        pred_ms.add(pred_lbls[k].get_cdr3)
      end
      # compute T(A,B)
      {:clone => self.jaccard_multi(pred_ms, truth_ms)}
    end

    #
    # compute Jaccard coefficient given two sets, e.g., two arrays
    # of partitions
    # J(A,B) = \frac{ |A \wedge B| }{ |A \cup B| }
    #
    def jaccard(set_a, set_b)
      (set_a & set_b).size.to_f / (set_a | set_b).size.to_f
    end

    #
    # compute Jaccard similarity given two multiset objects
    #
    def jaccard_multi(mset_a, mset_b)
      raise 'not a multiset object' if mset_a.class != Multiset
      raise 'not a multiset object' if mset_b.class != Multiset
      numerator = (mset_a & mset_b).size
      denom = (mset_a | mset_b).size
      (numerator / denom.to_f)
    end

    #
    # given a prediction file, returh a hash of:
    # CDR3_seq.to_sym => [id_x, id_y, ...]
    #
    def get_cdr3_hash(pred_file)
      h1 = {}
      IO.foreach(pred_file) do |l|
        lp = LabelPrediction.new(l)
        id = lp.get_id
        cdr3_seq = lp.get_cdr3.to_sym
        if h1.has_key?(cdr3_seq) then
          h1[cdr3_seq].push(id)
        else
          h1[cdr3_seq] = [id]
        end
      end
      h1
    end

    #

    # a = two ids in same partition in both A and B
    # b = two ids in different partitions in both A and B
    # c = two ids in same partition in A, but different in B
    # d = two ids in different partition in A, but same in B
    # (a+b+c+d) = (n choose 2)
    #
    def compute_pair_confusion(lbl_f1, lbl_f2)
      # create temp files for each with: [CDR3_seq, [id_x, id_y, ...]]
      h1 = get_cdr3_hash(lbl_f1)
      h2 = get_cdr3_hash(lbl_f2)
      num_f1 = `wc -l #{lbl_f1}`.match(/^(\d+)/)[1].to_i
      num_f2 = `wc -l #{lbl_f2}`.match(/^(\d+)/)[1].to_i
      #puts [h1.size, h2.size].join("\t")
      #puts [num_f1, num_f2].join("\t")
      #
      # create adjacency file, a node for each id,
      # an edge if two ids are in the same CDR3_seq partition
      adj_f1 = Tempfile.new('ad1')
      h1.each_key do |k|
        id_lst = h1[k]
        id_lst.each_with_index do |id1, i|
          id_lst.each_with_index do |id2, j|
            next if i >= j
            adj_f1.puts [id1, id2].join("\t")
          end
        end
      end
      adj_f1.close
      #
      adj_f2 = Tempfile.new('ad1')
      h2.each_key do |k|
        id_lst = h2[k]
        id_lst.each_with_index do |id1, i|
          id_lst.each_with_index do |id2, j|
            next if i >= j
            adj_f2.puts [id1, id2].join("\t")
          end
        end
      end
      adj_f2.close
      #
      # find intersection of two adjacency files, this is computing a in Rand index
      #FileUtils.cp(adj_f1.path, "/tmp/adj_f1.tab")
      #FileUtils.cp(adj_f2.path, "/tmp/adj_f2.tab")
      adj1_srt_f = Tempfile.new('adj1_sort')
      `sort #{adj_f1.path} > #{adj1_srt_f.path}`
      adj1_srt_f.close
      adj2_srt_f = Tempfile.new('adj2_sort')
      `sort #{adj_f2.path} > #{adj2_srt_f.path}`
      adj2_srt_f.close
      cmd = "bash -c 'comm -12 #{adj1_srt_f.path} #{adj2_srt_f.path} | wc -l'"
      out, err, pip = Open3.capture3(cmd)
      num_int = out.to_i
      a = num_int
      cmd_c = "bash -c 'comm -23 #{adj1_srt_f.path} #{adj2_srt_f.path} | wc -l'"
      c_val = `#{cmd_c}`.to_i
      cmd_d = "bash -c 'comm -13 #{adj1_srt_f.path} #{adj2_srt_f.path} | wc -l'"
      d_val = `#{cmd_d}`.to_i
      #
      # find number of remaining entries in files, this correspond to c+d in Rand index
      num_adj1 = `wc -l #{adj_f1.path}`.match(/^(\d+)/)[1].to_i
      num_adj2 = `wc -l #{adj_f2.path}`.match(/^(\d+)/)[1].to_i
      cd = num_adj1+num_adj2-(2*num_int)
      # explicitly rm tempfiles
      adj1_srt_f.unlink
      adj2_srt_f.unlink
      adj_f1.unlink
      adj_f2.unlink
      #
      # b = total - a + (c+d)
      total = (num_f1*(num_f1-1))/2.0
      b = total - (a + cd)
      #puts [a, b, cd, total].join("\t")
      {:a => a, :b => b, :c => c_val, :d => d_val, :total => total, :n_a => num_f1, :n_b => num_f2}
    end

    #
    # compute Rand index given two output files, RI(A,B) = (a+b)/(a+b+c+d)
    #
    def rand_index(lbl_f1, lbl_f2)
      h = compute_pair_confusion(lbl_f1, lbl_f2)
      {:clone => (h[:a]+h[:b])/(h[:c]+h[:d]+h[:a]+h[:b]).to_f}
    end
    #
    #
    #
    def jaccard_index(lbl_f1, lbl_f2)
      h = compute_pair_confusion(lbl_f1, lbl_f2)
      {:clone => h[:a]/(h[:c]+h[:d]+h[:a]).to_f}
    end
    #
    # compute Fowlkes-Mallows index
    #
    def fowlkes_mallows_index(lbl_f1, lbl_f2)
      h = compute_pair_confusion(lbl_f1, lbl_f2)
      a = h[:a]
      c = h[:c]
      d = h[:d]
      {:clone => (a/Math.sqrt((a+c)*(a+d)))}
    end

  end

end
