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

@samifier = Samifier::Base.new(protein_table_file)

results = @samifier.parse_mascot_search_results(mascot_results_file)
results.each do |result|
  next unless result[:protein] == 'RS4A_YEAST'
  puts result.inspect
  peptide = result[:peptide]
  protein_location = @samifier.get_sequence_locations(genome_file, result[:protein])
  puts protein_location.inspect
  seq_parts = @samifier.extract_sequence_parts(chromosome_dir, protein_location)
  puts seq_parts.inspect
  full_seq = @samifier.nucleotide_sequence(seq_parts)
  amino_seq = @samifier.nucleotide_to_amino_acid_sequence(full_seq).split(/\s+/).join
  puts amino_seq
  puts peptide


  peptide_start = (result[:start] - 1) * 3
  peptide_stop = result[:stop] * 3 - 1
  part_seq = ''
  cigar = ''
  start_index = 0
  stop_index = 0
  seq_parts.each do |part|
    if part[:type] == 'intron'
      # add to cigar
      next
    end

    seq_size = part[:sequence].size
    if peptide_start >= seq_size 
      peptide_start -= seq_size
      peptide_stop -= seq_size
      next
    else
      #if peptide_start < 
      #start_index = 
    end

    if peptide_stop > seq_size
      stop_index = peptide_stop
    else
      # possibly spanning
      stop_index = seq_size - 1
    end
    puts part[:sequence][start_index..stop_index]
  end
  #puts @samifier.to_sam_entry(protein_name, seq_parts)
  exit
end
