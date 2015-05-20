require "bundler/gem_tasks"

require 'rake'
require 'rake/testtask'

Rake::TestTask.new(:all) do |t|
  t.pattern = "tests/test_*.rb"
end


task :default => :all
