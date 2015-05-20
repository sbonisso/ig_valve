class Array
  # sums elements of an array together
  def sum() self.inject(0){|s,val| !val.is_a?(Numeric) ? raise(TypeError,"not a number") : s += val} end
end
#
class Range
  # checks for overlap between two ranges
  def overlap?(range)
    self.include?(range.first) || range.include?(self.first)
  end
  #
  # intersection of two ranges
  def intersection(range)
    int_a = (self.to_a & range.to_a)    
    int_a.empty? ? nil : Range.new(int_a.min, int_a.max)
  end  
end
