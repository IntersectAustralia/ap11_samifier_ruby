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

    def initialize(genome_file, protein_table_file)
      @protein_to_oln = protein_to_ordered_locus_name_map(protein_table_file)
      @genome = load_genome(genome_file)
    end

    def load_genome(genome_file)
      genome = {}
      IO.foreach(genome_file) do |line|
        next if line =~ /^\#/
        (chromosome, source, type, start, stop, score, strand, phase, attributes) = line.split(/\s+/)
        next if type.nil?
        type.match(/(gene|CDS|intron)/) do |type_match|
          seq_type = type_match[1]
          if seq_type == 'gene'
            attributes.match(/ID=(.+?);/) do |attr_match|
              ordered_locus_name = attr_match[1]
              if !genome.has_key?(ordered_locus_name)
                genome[ordered_locus_name] = {
                  chromosome: chromosome,
                  start: start.to_i,
                  direction: strand,
                  locations: []
                }
              end
            end
          else
            attributes.match(/Name=([^_]+)_#{seq_type}/) do |attr_match|
              ordered_locus_name = attr_match[1]
              next unless @protein_to_oln.has_value?(ordered_locus_name)
              genome[ordered_locus_name][:locations] << {
                type: seq_type,
                start: start.to_i,
                stop: stop.to_i,
                direction: strand
              }
            end
          end
        end
      end
      genome.each { |k,v| v[:locations].sort!{|a,b| a[:start] <=> b[:start]} }
      genome
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
              start: m[2].to_i,
              stop: m[3].to_i
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

    def get_sequence_locations(protein)
      ordered_locus_name = @protein_to_oln[protein]
      @genome[ordered_locus_name]
    end

    def get_peptide_sequence(peptide, seq_parts)
      part_seq = ''
      cigar = ''
      peptide_start = (peptide[:start] - 1) * 3
      peptide_stop = peptide[:stop] * 3 - 1
      remaining = peptide_stop - peptide_start + 1

      seq_parts.each do |part|
        if part[:type] == 'intron'
          cigar << "#{part[:stop]-part[:start]+1}N" unless cigar.empty?
          next
        end

        seq_size = part[:sequence].size

        if peptide_start >= seq_size 
          peptide_start -= seq_size
          next
        end

        if (peptide_start + remaining) < seq_size
          part_seq << part[:sequence][peptide_start, remaining]
          cigar << "#{remaining}M"
          break
        else
          part_seq << part[:sequence][peptide_start..seq_size-1]
          cigar << "#{seq_size-peptide_start}M"
          remaining -= (seq_size - peptide_start)
          peptide_start = 0
        end
      end
      {
        sequence: part_seq,
        cigar: cigar
      }
    end

    def to_sam_entry(qname, rname, pos, cigar, peptide_sequence)
      flag = 0
      mapq = 255
      rnext = '='
      pnext = 0
      tlen = 0
      qual = '*'
      [qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, peptide_sequence, qual]
    end

    def to_sam_file(mascot_results_file, chromosome_dir, outfile)
      results = parse_mascot_search_results(mascot_results_file)
      output_files = {}
      sam_entries = []
      results.each do |result|
        protein_location = get_sequence_locations(result[:protein])
        if protein_location.nil?
          $stderr.puts "#{result[:protein]} not found in genome file"
          next
        end
        # TODO: Directionality supported in a future release
        next unless protein_location[:direction] == '+'
        seq_parts = extract_sequence_parts(chromosome_dir, protein_location)

        peptide_sequence = get_peptide_sequence(result, seq_parts)
        protein_name = result[:protein]
        peptide_start = result[:start] + @genome[@protein_to_oln[protein_name]][:start]
        sam_entries << to_sam_entry(result[:id], protein_location[:chromosome], peptide_start, peptide_sequence[:cigar], peptide_sequence[:sequence])
      end
      sam_entries.sort_by!{|e| [e[2],e[3]]}.each do |sam_entry|
        chromosome = sam_entry[2]
        output_files[chromosome] ||= File.open("#{outfile}.#{chromosome}.sam", 'w')
        output_files[chromosome].puts sam_entry.join("\t")
      end
      output_files.values(&:close)
    end

    def nucleotide_sequence(sequence_parts)
      sequence_parts.select{|s| s[:type] == 'CDS'}.map{|s|s[:sequence]}.join
    end

    def extract_sequence_parts(chromosome_dir, protein_location)
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
      file.close
      nucleotide_sequence_parts
    end

    def self.nucleotide_to_amino_acid_sequence(nucleotide_seq)
      aminos = ''
      count = 0
      seq_io = StringIO.new(nucleotide_seq)
      while codon = seq_io.read(3)
        if @@codons.has_key?(codon.to_sym)
          count += 1
          aminos << @@codons[codon.to_sym]
        else
          $stderr.puts "codon #{codon} not found" unless @@stop_codons.include?(codon)
        end
      end
      aminos
    end

  end
end
