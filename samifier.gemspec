# -*- encoding: utf-8 -*-
require File.expand_path('../lib/samifier/version', __FILE__)

Gem::Specification.new do |gem|
  gem.authors       = ["Sean McCarthy"]
  gem.email         = ["sean@intersect.org.au"]
  gem.description   = %q{Turn Mascot search for peptides into a SAM file}
  gem.summary       = %q{Turn Mascot search for peptides into a SAM file}
  gem.homepage      = "https://github.com/IntersectAustralia/samifier"

  gem.files         = `git ls-files`.split($\)
  gem.executables   = gem.files.grep(%r{^bin/}).map{ |f| File.basename(f) }
  gem.test_files    = gem.files.grep(%r{^(test|spec|features)/})
  gem.name          = "samifier"
  gem.require_paths = ["lib"]
  gem.version       = Samifier::VERSION
end
