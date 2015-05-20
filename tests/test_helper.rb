# minitest stuff
require 'rubygems'
gem 'minitest'
require 'minitest/autorun'

# minitest-reporters
require "minitest/reporters"
#Minitest::Reporters.use!
Minitest::Reporters.use! Minitest::Reporters::SpecReporter.new

def get_data_dir()
  "#{File.dirname(__FILE__)}/data/"
end