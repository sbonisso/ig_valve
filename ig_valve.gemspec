# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'ig_valve/version'

Gem::Specification.new do |spec|
  spec.name          = "ig_valve"
  spec.version       = IgValve::VERSION
  spec.authors       = ["sbonisso"]
  spec.email         = ["sbonisso@ucsd.edu"]

  # if spec.respond_to?(:metadata)
  #   spec.metadata['allowed_push_host'] = "TODO: Set to 'http://mygemserver.com' to prevent pushes to rubygems.org, or delete to allow pushes to any server."
  # end

  spec.summary       = %q{Ig Vdj Antibody Labeling Validation and Evaluation}
  spec.description   = %q{Validates antibody VDJ labeling of reads in either supervised (to simulated ground truth) or unsupervised (to real data) modes.}
  spec.homepage      = "https://bitbucket.org/sbonisso/ig_valve"
  spec.license       = "MIT"

  spec.files         = `git ls-files -z`.split("\x0").reject { |f| f.match(%r{^(test|spec|features)/}) }
  spec.bindir        = "bin"
  #spec.executables   = spec.files.grep(%r{^exe/}) { |f| File.basename(f) }
  spec.executables   = "ig_valve"
  spec.require_paths = ["lib"]

  spec.add_dependency "multiset", "~>0.1"
  spec.add_dependency "cluster_eval", "~>0.1"

  spec.add_development_dependency "minitest", "~> 0.0"
  spec.add_development_dependency "bundler", "~> 1.8"
  spec.add_development_dependency "rake", "~> 10.0"
end
