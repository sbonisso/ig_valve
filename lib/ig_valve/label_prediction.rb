
module IgValve
  #
  #
  class LabelPrediction
    #
    #
    def initialize(str, delim=";")
      #ary = str.chomp.split("\t")
      ary = str.split("\t")
      ary[-1] = ary[-1].chomp
      @id = ary[0]
      #
      # labels, within each sep by ';'
      @v_preds = ary[1].split(delim)
      @d_preds = ary[2].split(delim)
      ary[3] = "?,?" if ary[3].nil?
      @j_preds = ary[3].split(delim)
      if ary.size == 5 then
        #
        # cdr3 sequence
        @cdr3_seq = ary[4]
      end
      if ary.size > 5 then
        #
        # partitions, shown as range x..y
        @v_part = Range.new(*ary[4].split("..").map{|v| v.to_i})
        @d_part = Range.new(*ary[5].split("..").map{|v| v.to_i})
        @j_part = Range.new(*ary[6].split("..").map{|v| v.to_i})
      end

    end
    
    def get_id() @id end

    def get_v_alleles() @v_preds end
    def get_d_alleles() @d_preds end
    def get_j_alleles() @j_preds end

    def get_v_genes() @v_preds.map{|v| v.split("*")[0]}.uniq end
    def get_d_genes() @d_preds.map{|v| v.split("*")[0]}.uniq end 
    def get_j_genes() @j_preds.map{|v| v.split("*")[0]}.uniq end

    def get_v_partition() @v_part.to_a end
    def get_d_partition() @d_part.to_a end
    def get_j_partition() @j_part.to_a end

    def get_cdr3() @cdr3_seq end

  end

end
