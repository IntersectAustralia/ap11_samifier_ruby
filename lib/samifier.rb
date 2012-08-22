require "samifier/version"
require 'stringio'
require 'set'

module Samifier

  class Base
    LINE_LENGTH = 70

    @@codons = {
      # Isoleucine
      ATT: 'I', ATC: 'I', ATA: 'I',

      # Leucine
      CTT: 'L', CTC: 'L', CTA: 'L', CTG: 'L', TTA: 'L', TTG: 'L',

      # Valine
      GTT: 'V', GTC: 'V', GTA: 'V', GTG: 'V',

      TTT: 'F', TTC: 'F',

      ATG: 'M', # Also the start codon

      TGT: 'C', TGC: 'C',

      GCT: 'A', GCC: 'A', GCA: 'A', GCG: 'A',

      GGT: 'G', GGC: 'G', GGA: 'G', GGG: 'G',

      CCT: 'P', CCC: 'P', CCA: 'P', CCG: 'P',

      ACT: 'T', ACC: 'T', ACA: 'T', ACG: 'T',

      TCT: 'S', TCC: 'S', TCA: 'S', TCG: 'S', AGT: 'S', AGC: 'S',

      TAT: 'Y', TAC: 'Y',

      TGG: 'W',

      CAA: 'Q', CAG: 'Q',

      AAT: 'N', AAC: 'N',

      CAT: 'H', CAC: 'H',

      GAA: 'E', GAG: 'E',

      GAT: 'D', GAC: 'D',

      AAA: 'K', AAG: 'K',

      CGT: 'R', CGC: 'R', CGA: 'R', CGG: 'R', AGA: 'R', AGG: 'R'
    }

    @@start_codon = 'ATG'
    @@stop_codons = Set.new [ 'TGA', 'TAA', 'TAG' ]

    def initialize(protein_table_file)
      @protein_to_oln = protein_to_ordered_locus_name_map(protein_table_file)
    end

    def parse_mascot_search_results(search_results_file)
      peptides_start = false
      results = []
      IO.foreach(search_results_file) do |line|
        if peptides_start
          break if line.match(/^--/)
          results << get_proteins_from_query(line)
        elsif line.match(/Content-Type: application\/x-Mascot; name="peptides"/)
          peptides_start = true
        end
      end
      results.flatten!
    end

    def get_proteins_from_query(line)
      results_list = []
      line.match(/^(q\d+_p\d+)=([^;]+)\;(.+)$/) do |m|
        id = m[1]
        peptide_part= m[2]
        proteins_part = m[3]
        parts = peptide_part.split(/,/)
        peptide = parts[4]
        proteins_part.split(/,/).each do |protein_part|
          # "14311_ARATH":0:192:197:1
          protein_part.match(/^"([^"]+)":\d:(\d+):(\d+):\d/) do |m|
            protein = m[1]
            next unless @protein_to_oln.has_key?(protein)
            results_list << {
              protein: protein,
              id: "#{protein}.#{id}",
              peptide: peptide,
              start: m[2],
              stop: m[3]
            }
          end
        end
      end
      results_list
    end

    def protein_to_ordered_locus_name_map(protein_table_file)
      protein_map = {}
      IO.foreach(protein_table_file) do |line|
        (oln, acc, protein_id, db_id) = line.split(/\s+/)
        protein_map[protein_id] = oln
      end
      protein_map
    end

    def get_sequence_locations(genome_file, protein)
      sequence_locations = {
        chromosome: nil,
        locations: []
      }
      ordered_locus_name = @protein_to_oln[protein]
      IO.foreach(genome_file) do |line|
        (chromosome, db, type, start, stop, a, b, c, desc) = line.split(/\s+/)
        next if type.nil?
        type.match(/(CDS|intron)/) do |m|
          seq_type = m[1]
          if desc =~ /Name=#{ordered_locus_name}_#{seq_type}/
            sequence_locations[:chromosome] = chromosome
            sequence_locations[:locations] << {
              type: seq_type,
              start: start.to_i,
              stop: stop.to_i
            }
          end
        end
      end
      sequence_locations[:locations].sort!{|a,b| a[:start] <=> b[:start]}
      sequence_locations
    end

    def roll_cigar(sequences)
      cigar = ''
      sequences.each do |seq|
        cigar << case seq[:type]
          when 'CDS' then
            "#{seq[:stop]-seq[:start]+1}M"
          when 'intron' then
            "#{seq[:stop]-seq[:start]+1}N"
        end
      end 
      cigar
    end

    def to_sam_entry(name, sequence_parts)
      qname = name
      flag = 0
      rname = sequence_parts[:TODO]
      pos = sequence_parts.first[:start]
      mapq = 255
      cigar = roll_cigar(sequence_parts)
      rnext = '='
      pnext = 0
      tlen = 0
      seq = nucleotide_sequence(sequence_parts)
      qual = '*'
      [qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual].join("\t")
    end

    def nucleotide_sequence(sequence_parts)
      sequence_parts.select{|s| s[:type] == 'CDS'}.map{|s|s[:sequence]}.join
    end

    def extract_sequence_parts(chromosome_dir, protein_location)
      # {:chromosome=>"chrV", :locations=>[{:type=>"CDS", :start=>545611, :stop=>546414}]}
      nucleotide_sequence_parts = []
      read_stop = 0

      chromosome = protein_location[:chromosome]
      locations = protein_location[:locations]
      chromosome_file = File.join(chromosome_dir, chromosome + ".fa")
      file = File.open(chromosome_file)
      file.readline # Skip header line
      read_cursor = 0
      line = nil

      locations.each do |location|
        start_index = location[:start]
        stop_index = location[:stop]

        if location[:type] == 'intron'
          nucleotide_sequence_parts << {
            type: location[:type],
            start: start_index,
            stop: stop_index
          }
          next
        end

        nucleotide_sequence = ''

        while read_cursor < start_index
          line = file.gets
          line.chomp!
          read_cursor += line.size
        end

        read_start = start_index % line.size - 1
        read_stop  = line.size

        while read_cursor < stop_index
          nucleotide_sequence << line[read_start..read_stop]
          read_start = 0
          line = file.gets
          line.chomp!
          read_cursor += line.size
        end

        read_stop = stop_index % line.size - 1
        nucleotide_sequence << line[read_start..read_stop]

        nucleotide_sequence_parts << {
          sequence: nucleotide_sequence,
          type: location[:type],
          start: start_index,
          stop: stop_index
        }
      end
      nucleotide_sequence_parts
    end

    def nucleotide_to_amino_acid_sequence(nucleotide_seq, line_length=60)
      aminos = ''
      count = 0
      seq_io = StringIO.new(nucleotide_seq)
      while codon = seq_io.read(3)
        if @@codons.has_key?(codon.to_sym)
          count += 1
          aminos << @@codons[codon.to_sym]
          aminos << "\n" if count % line_length == 0
        else
          $stderr.puts "codon #{codon} not found" unless @@stop_codons.include?(codon)
        end
      end
      aminos
    end

  end
end
