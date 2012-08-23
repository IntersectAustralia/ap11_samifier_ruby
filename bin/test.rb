#!/usr/bin/env ruby
require 'samifier'

def usage
  $stderr.puts "read_n.rb mascot_results_file protein_table_file chromosome_directory genome_file"
end

unless ARGV.size == 4
  usage
  exit
end

mascot_results_file = ARGV[0]
protein_table_file = ARGV[1]
chromosome_dir = ARGV[2]
genome_file = ARGV[3]

@samifier = Samifier::Base.new(genome_file, protein_table_file)

results = @samifier.parse_mascot_search_results(mascot_results_file)
results.each do |result|
  next unless result[:protein] == 'RL36A_YEAST'
  protein_location = @samifier.get_sequence_locations(result[:protein])
  if protein_location.nil?
    $stderr.puts "#{result[:protein]} not found in #{genome_file}"
    next
  end
  unless protein_location[:direction] == '+'
    #$stderr.puts "#{result[:protein]}: #{protein_location[:direction]} non-Watson encoding not supported"
    next
  end
  seq_parts = @samifier.extract_sequence_parts(chromosome_dir, protein_location)
  puts seq_parts.inspect

  peptide_sequence = @samifier.get_peptide_sequence(result, seq_parts)
  puts result[:protein]
  puts result[:peptide]
  puts peptide_sequence.inspect
  puts @samifier.nucleotide_to_amino_acid_sequence(peptide_sequence[:sequence])
  puts '----'
  #puts @samifier.to_sam_entry(protein_name, seq_parts)
end
