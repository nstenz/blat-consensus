#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# Not yet implemented
#my $overhang_length = 50;
#my $overhang_length = 0;
my $target_perc_cov_threshold = 0;
my $contig_perc_cov_threshold = 0;

my $blat = check_path_for_exec("blat");
my $mafft = check_path_for_exec("mafft");

my $targets;
GetOptions(
	't|targets=s' => \$targets,
);

my $input = shift;

die "You must specify an assembly as input.\n" if (!defined($input));
die "You must specify a fasta file containing target sequences (-t).\n" if (!defined($targets));
die "You must you have specified an incorrect argument.\n" if (@ARGV);

# Output file names
#my $output_fasta = $input.".hits";
my $blat_out_name = $input.".psl";
(my $coverage_output = $input) =~ s/\.fa(sta)?/.cov/i;
(my $consensus_output = $input) =~ s/\.fa(sta)?/.con.fasta/i;

# Run blat 
print "Running blat against targets and contigs...\n";
my $return = system("$blat $input $targets -t=dna -q=dna -noHead $blat_out_name");
die "Error running blat: '$return'.\n" if ($return);
#my $return = system("$blat $input $targets -t=dna -q=dna -noHead $blat_out_name") if (!-e $blat_out_name);

# Parse blat output
my %targets;
open(my $blat_out, "<", $blat_out_name);
while (my $line = <$blat_out>) {
	chomp($line);

	my @line = split(/\s+/, $line);

	# Blat headers
	my ($match, $mismatch, $rep_matches, $n_count, $q_num_inserts, $q_base_inserts, 
	    $t_num_inserts, $t_base_inserts, $strand, $q_name, $q_size, $q_start, $q_end,
		$t_name, $t_size, $t_start, $t_end, $block_count, $block_sizes, $q_starts, $t_starts) = @line; 

	# Check that contig meets target coverage threshold
	if (($match + $mismatch) / $q_size >= $target_perc_cov_threshold && 
		($match + $mismatch) / $q_size >= $contig_perc_cov_threshold) {
		#push(@{$targets{$q_name}}, {'NAME' => $t_name, 'STRAND' => $strand, 'T_STARTS' => $t_starts, 'Q_STARTS' => $q_starts, 'SIZES' => $block_sizes});
		if (!exists($targets{$q_name}->{$t_name})) {
			$targets{$q_name}->{$t_name} = {'NAME' => $t_name, 'STRAND' => $strand, 'T_STARTS' => $t_starts, 'Q_STARTS' => $q_starts, 'SIZES' => $block_sizes};
		}
		else {
			my $hit = $targets{$q_name}->{$t_name};

			my $current_size = $hit->{SIZES};
			my @current_size = split(",", $current_size);

			$current_size = 0;
			foreach my $size (@current_size) {
				$current_size += $size;
			}

			my $new_size = $block_sizes;
			my @new_size = split(",", $new_size);

			$new_size = 0;
			foreach my $size (@new_size) {
				$new_size += $size;
			}

			if ($current_size > $new_size) {
				$targets{$q_name}->{$t_name} = {'NAME' => $t_name, 'STRAND' => $strand, 'T_STARTS' => $t_starts, 'Q_STARTS' => $q_starts, 'SIZES' => $block_sizes};
			}
		}
	}
}
close($blat_out);

# Load contig and target sequences
my %contig_seqs = parse_fasta($input);
my %target_seqs = parse_fasta($targets);

# Create consensus sequence for each target
my %consensuses;
my %coverage_maps;
foreach my $target (sort { $a cmp $b } keys %targets) {
	#my @hits = @{$targets{$target}};
	my @hits = values %{$targets{$target}};

	foreach my $hit (@hits) {
		print "$hit->{'NAME'}\n";
	}

	# Holds bases/introns in contigs which map to specific sites in target
	my %introns;
	my %overhangs;
	my %alignment;

	# Add target bases to alignment
	my $target_seq = $target_seqs{$target};
	my $target_length = length($target_seq);
	foreach my $index (0 .. length($target_seq) - 1) {
		my $base = substr($target_seq, $index, 1);
		$alignment{$index}->{"BASES"} = {};
		$alignment{$index}->{"TARGET"} = $base;
		$alignment{$index}->{"INTRONS"} = {};
	}

	print "Building consensus sequence for '$target'...\n";

	# Extract homologous sequence from each contig
	foreach my $hit (@hits) {
		print "  ",$hit->{'NAME'},"\n";

		my $seq = $contig_seqs{$hit->{'NAME'}};
		my $contig_length = length($seq);

		print "    Contig length: $contig_length\n";
		print "    Target length: $target_length\n";

		my @sizes = split(",", $hit->{'SIZES'});
		my @t_starts = split(",", $hit->{'T_STARTS'});
		my @q_starts = split(",", $hit->{'Q_STARTS'});

		# Output strandedness
		if ($hit->{'STRAND'} eq "-") {
			print "    Contig will be reverse complemented (- hit).\n\n";

			$seq = rev_comp($seq);
			foreach my $index (0 .. $#t_starts) {
				$t_starts[$index] = $contig_length - $sizes[$index] - $t_starts[$index];
				$q_starts[$index] = $target_length - $sizes[$index] - $q_starts[$index];
			}
			@sizes = reverse(@sizes);
			@t_starts = reverse(@t_starts);
			@q_starts = reverse(@q_starts);
		}
		else {
			print "    Contig will not be reverse complemented (+ hit).\n\n";
		}

		my $prematch = substr($seq, 0, $t_starts[0] + 1);
		my $postmatch = substr($seq, $t_starts[$#t_starts], length($seq) - $t_starts[$#t_starts]);

		$overhangs{$q_starts[0]}->{$prematch}++;
		$overhangs{$q_starts[$#q_starts] + $sizes[$#q_starts] - 1}->{$postmatch}++;

		# Iterate through aligned blocks of target and contig
		foreach my $index (0 .. $#t_starts) {
			my $size = $sizes[$index];
			my $q_start = $q_starts[$index];
			my $t_start = $t_starts[$index];

			# Where the intron starts 
			my $intron_start = $t_start + $size - 1;

			# If this is the last (or only) block, we set the next block 
			# start to a value which makes the intron length 0
			my $intron_length = 0;
			if ($index + 1 < scalar(@t_starts)) {
				$intron_length = $t_starts[$index + 1] - $intron_start - 1;
			}

			print "    match_start : $t_start ($q_start)\n";
			print "    match_length : $size\n";
			print "    match_end : ",($t_start + $size - 1),"\n";

			if ($intron_length > 0) {
				print "    intron_start : $intron_start (",($q_start + $intron_start - $t_start),"), intron_length : $intron_length\n\n";
			}
			else {
				print "\n\n";
			}

			# Place aligned characters into %alignment
			foreach my $offset (0 .. $size - 1) {
				my $t_index = $t_start + $offset;
				my $q_index = $q_start + $offset;

				my $char = substr($seq, $t_index, 1);
				die "Could not get character at $t_index\n" if (!$char);
				$alignment{$q_index}->{"BASES"}{$char}++;
			}

			# Place introns into %alignment
			if ($intron_length > 0) {
				my $intron = substr($seq, $intron_start + 1, $intron_length);
				die "Could not get intron at ",($intron_start + 1),"\n" if (!$intron);
				$alignment{$q_start + $intron_start - $t_start}->{"INTRONS"}{$intron}++;
				$introns{$q_start + $intron_start - $t_start}->{$intron}++;
			}
		}
	}

	use Data::Dumper;
	$Data::Dumper::Sortkeys = sub { [sort { $a <=> $b } keys %{$_[0]}] };
	print "Overhangs:\n",Dumper(\%overhangs),"\n";

	# Handle overhangs
	foreach my $site (keys %overhangs) {
		if (exists($introns{$site})) {
			$alignment{$site}->{INTRONS} = align_overhangs_to_introns($overhangs{$site}, $introns{$site});
		}
		else {
			# Total up how many times an overhang occurs at this site
			my $count = 0;
			foreach my $overhang (values %{$overhangs{$site}}) {
				$count += $overhang;
			}

			# Ignore sites with only one overhanging sequence
			if ($count > 1) {
				$alignment{$site}->{INTRONS} = align_overhangs($overhangs{$site});
			}
		}
	}

	my ($consensus, $coverage) = consense_alignment(\%alignment);
	#print "\n\n",$consensus,"\n$coverage\n\n";

#	#die if ($alignment =~ /a|t|c|g/);
	$consensuses{$target} = $consensus;	
	$coverage_maps{$target} = $coverage;	
}

# Output consensus sequences
open(my $output, ">", $consensus_output);
foreach my $target (sort { $a cmp $b } keys %consensuses) {
	print {$output} ">$target\n";
	print {$output} "$consensuses{$target}\n";
}
close($output);

# Output coverage maps
open($output, ">", $coverage_output);
foreach my $target (sort { $a cmp $b } keys %coverage_maps) {
	print {$output} ">$target\n";
	print {$output} "$coverage_maps{$target}\n\n";
}
close($output);

sub align_overhangs {
	my ($overhangs) = @_;

	# Too few hits to include as an intron
	return {} if (scalar(values %{$overhangs}) == 1);

	my $count = 0;
	my $blat_out_name = "overhangs.psl";
	my $overhang_fasta = "overhangs.fasta";

	# Write overhang sequences to file
	my %overhangs;
	open(my $out, ">", $overhang_fasta);
	foreach my $sequence (keys %{$overhangs}) {
		$overhangs{"OVERHANG_$count"} = $sequence;
		print {$out} ">OVERHANG_$count\n";
		print {$out} "$sequence\n";
		$count++;
	}
	close($out);

	# Initialize hash used later to determine which overhangs align to eachother
	my %hits;
	foreach my $overhang (keys %overhangs) {
		foreach my $index (0 .. $count - 1) {
			next if ($overhang eq "OVERHANG_$index");

			$hits{$overhang}->{"OVERHANG_$index"}++;
		}
	}

	my $return = system("$blat $overhang_fasta $overhang_fasta -t=dna -q=dna -noHead $blat_out_name >/dev/null");
	die "Error running blat: '$return'.\n" if ($return);

	# Parse blat output
	my %sequence_indices;
	open(my $blat_out, "<", $blat_out_name);
	while (my $line = <$blat_out>) {
		chomp($line);

		my @line = split(/\s+/, $line);

		# Blat headers
		my ($match, $mismatch, $rep_matches, $n_count, $q_num_inserts, $q_base_inserts, 
			$t_num_inserts, $t_base_inserts, $strand, $q_name, $q_size, $q_start, $q_end,
			$t_name, $t_size, $t_start, $t_end, $block_count, $block_sizes, $q_starts, $t_starts) = @line; 

		# All hits should be in the + strand
		next if ($strand eq "-");

		# Exclude self-matches
		next if ($t_name eq $q_name);

		delete($hits{$q_name}->{$t_name});

		#print "$line\n";

		my @sizes = split(",", $block_sizes);
		my @starts = split(",", $q_starts);

		# Add aligned part of sequence to final sequences
		foreach my $index (0 .. $#sizes) {
			my $size = $sizes[$index];
			my $start = $starts[$index];
			my $end = $start + $size - 1;

			my $range = "$start-$end";
			if (exists($sequence_indices{$q_name})) {
				$sequence_indices{$q_name} = intersection($range, $sequence_indices{$q_name});
			}
			else {
				$sequence_indices{$q_name} = $range; 
			}
		}
	}
	close($blat_out);

	unlink($blat_out_name, $overhang_fasta);

	# Determine if any overhangs are missing homologous sequence with others
	my $fewest_hits;
	my @fewest_hits;
	foreach my $overhang (keys %hits) {
		my $num_unmatched = scalar(keys %{$hits{$overhang}});

		if ($num_unmatched) {
			if (!defined($fewest_hits) || $num_unmatched > $fewest_hits) {
				undef(@fewest_hits);
				$fewest_hits = $num_unmatched;
				push(@fewest_hits, $overhang);
			}
			elsif ($num_unmatched == $fewest_hits) {
				push(@fewest_hits, $overhang);
			}
		}
	}

	# Remove any overhangs which don't align to other overhangs
	if (@fewest_hits) {
		my $overhang_to_delete;
		if (scalar(@fewest_hits) == 1) {
			$overhang_to_delete = $fewest_hits[0];
		}
		else {
			my $shortest_overhang;
			my $shortest_overhang_length;
			foreach my $overhang (@fewest_hits) {
				my $sequence = $overhangs{$overhang};
				if (!defined($shortest_overhang) || length($sequence) < $shortest_overhang_length) {
					$shortest_overhang = $overhang;
					$shortest_overhang_length = length($sequence);
				}
			}
			$overhang_to_delete = $shortest_overhang;
		}
		delete($overhangs{$overhang_to_delete});

		my $overhangs;
		foreach my $overhang (values %overhangs) {
			$overhangs->{$overhang}++;
		}
		return align_overhangs($overhangs);
	}
	else {
		# Translate sequence indices for overhangs to actual sequence
		my %final_sequences;
		foreach my $overhang (keys %sequence_indices) {
			my $range = $sequence_indices{$overhang};
			my $full_sequence = $overhangs{$overhang};

			my $sequence = get_partial_seq($range, $full_sequence);
			$final_sequences{$sequence}++;
		}
		return \%final_sequences;
	}
}

sub align_overhangs_to_introns {
	my ($overhangs, $introns) = @_;

	my $count = 0;
	my $intron_fasta = "introns.fasta";
	my $overhang_fasta = "overhangs.fasta";
	my $blat_out_name = "overhangs-to-introns.psl";

	# For now just align to the largest intron if there are multiple
	open(my $out, ">", $intron_fasta);
#	foreach my $sequence (keys %{$introns}) {
#		#$introns{"INTRON_$count"} = $sequence;
#		print {$out} ">INTRON_$count\n";
#		print {$out} "$sequence\n";
#		$count++;
#	}
	print {$out} ">INTRON_$count\n";
	print {$out} (sort { length($a) <=> length($b) } keys(%{$introns}))[0],"\n";
	$count++;
	close($out);

	# Add all overhangs to fasta file
	my %overhangs;
	open($out, ">", $overhang_fasta);
	foreach my $sequence (keys %{$overhangs}) {
		$overhangs{"OVERHANG_$count"} = $sequence;
		print {$out} ">OVERHANG_$count\n";
		print {$out} "$sequence\n";
		$count++;
	}
	close($out);

	# Run blat against overhangs and the longest intron
	my $return = system("$blat $overhang_fasta $intron_fasta -t=dna -q=dna -noHead $blat_out_name >/dev/null");
	die "Error running blat: '$return'.\n" if ($return);

	# Parse blat output
	my %final_sequences = %{$introns};
	open(my $blat_out, "<", $blat_out_name);
	while (my $line = <$blat_out>) {
		chomp($line);

		my @line = split(/\s+/, $line);

		# Blat headers
		my ($match, $mismatch, $rep_matches, $n_count, $q_num_inserts, $q_base_inserts, 
			$t_num_inserts, $t_base_inserts, $strand, $q_name, $q_size, $q_start, $q_end,
			$t_name, $t_size, $t_start, $t_end, $block_count, $block_sizes, $q_starts, $t_starts) = @line; 

		# All hits should be in the + strand
		next if ($strand eq "-");

		my @sizes = split(",", $block_sizes);
		my @starts = split(",", $t_starts);

		# Add aligned part of sequence to final sequences
		my $sequence;
		foreach my $index (0 .. $#sizes) {
			my $size = $sizes[$index];
			my $start = $starts[$index];
			$sequence .= substr($overhangs{$t_name}, $start, $size);
		}
		$final_sequences{$sequence}++;
	}
	close($blat_out);

	# Clean up
	unlink($overhang_fasta, $intron_fasta, $blat_out_name);

	return \%final_sequences;
}

sub consense_alignment {
	my $align = shift;

	my $coverage;
	my $alignment;

	my %align = %{$align};
	SITE: foreach my $site (sort { $a <=> $b } keys %align) {

		my %site = %{$align{$site}};
		my %bases = %{$site{"BASES"}};

		if (scalar(keys %bases) == 0) {
			$coverage .= "0 ";
			$alignment .= "?";
		}
		else {
			# An intron at site 0 isn't really an intron so it needs to be handled differently
			if ($site == 0) {
				consense_introns(\%site, \$coverage, \$alignment);
				consense_bases(\%site, \$coverage, \$alignment);
			}
			else {
				consense_bases(\%site, \$coverage, \$alignment);
				consense_introns(\%site, \$coverage, \$alignment);
			}
		}

#		# Only the target sequence's nucleotide present, no contigs mapped to this site
#		if (scalar(keys %bases) == 0) {
#			$coverage .= "0 ";
#			$alignment .= "?";
#		}
#		else {
#			my $target_base = $site{"TARGET"};
#
#			my %introns = %{$site{"INTRONS"}};
#
#			my $num_single_bases = 0;
#			foreach my $base_present (keys %bases) {
#				$num_single_bases += $bases{$base_present};
#			}
#
#			foreach my $base_present (keys %bases) {
#				# Perhaps filter by proportion here? (remove base if it is present in < x% of contigs)
#
#				# Ambiguous, can't determine consensus
#				if ($base_present eq "N") {
#					$coverage .= "$num_single_bases ";
#					$alignment .= "N";
#					next SITE;
#				}
#			}
#
#			# If only one nucleotide is found at the site output it
#			if (scalar(keys %bases) == 1) {
#				$coverage .= "$num_single_bases ";
#				$alignment .= (keys %bases)[0];
#			}
#			# Handle ambiguous non-introns with IUPAC codes
#			else {
#				$coverage .= "$num_single_bases ";
#				#$alignment .= get_IUPAC(\%single_bases);
#				$alignment .= get_IUPAC(\%bases);
#			}
#
#			# We have introns
#			if (scalar(keys %introns) > 0) {
#
#				# We only have one intron, output it
#				if (scalar(keys %introns) == 1) {
#					my $intron = lc((keys %introns)[0]);
#					$alignment .= $intron;
#
#					$coverage .= "(";
#					foreach my $index (0 .. length($intron) - 1) {
#						$coverage .= $introns{(keys %introns)[0]}." ";
#					}
#					chop($coverage);
#					$coverage .= ") ";
#				}
#				# We have multiple introns, align them
#				else {
#
#					# Create a temporary fasta file to hold introns
#					open(my $intron_fasta, ">", "introns.fa");
#					my $count = 0;
#					foreach my $intron (keys %introns) {
#						print {$intron_fasta} ">$count\n";	
#						print {$intron_fasta} "$intron\n";	
#						$count++;
#					}
#					close($intron_fasta);
#
#					# Align introns with mafft
#					my $return = system("$mafft --maxiterate 1000 --genafpair --thread 8 introns.fa > introns-aligned.fa 2> /dev/null");
#					die "Error running mafft: '$return'.\n" if ($return);
#
#					# Load aligned introns
#					my %intron_sequences = parse_fasta("introns-aligned.fa");
#					unlink("introns.fa", "introns-aligned.fa");
#
#					# Determine the consensus base at each site of the intron
#					$coverage .= "(";
#					foreach my $index (0 .. length((values %intron_sequences)[0]) - 1) {
#
#						my %bases;
#						foreach my $intron (values %intron_sequences) {
#							my $base = substr(uc($intron), $index, 1); 
#							die "$index" if ($base eq '');
#
#							# Not quite sure what to do with gaps, ignore for now
#							$bases{$base}++ if ($base ne "-");
#						}
#
#						# Get the IUPAC code for the current assortment of nucleotides
#						$alignment .= lc(get_IUPAC(\%bases));
#						$coverage .= scalar(keys %intron_sequences)." ";
#					}
#					chop($coverage);
#					$coverage .= ") ";
#				}
#			}
#		}
	}

	return ($alignment, $coverage);
}

sub consense_bases {
	my ($site, $coverage, $alignment) = @_;

	my %bases = %{$site->{BASES}};
	my $target_base = $site->{TARGET};

	my $num_single_bases = 0;
	foreach my $base_present (keys %bases) {
		$num_single_bases += $bases{$base_present};
	}

	foreach my $base_present (keys %bases) {
		# Perhaps filter by proportion here? (remove base if it is present in < x% of contigs)

		# Ambiguous, can't determine consensus
		if ($base_present eq "N") {
			$$coverage .= "$num_single_bases ";
			$$alignment .= "N";
			next SITE;
		}
	}

	# If only one nucleotide is found at the site output it
	if (scalar(keys %bases) == 1) {
		$$coverage .= "$num_single_bases ";
		$$alignment .= (keys %bases)[0];
	}
	# Handle ambiguous non-introns with IUPAC codes
	else {
		$$coverage .= "$num_single_bases ";
		$$alignment .= get_IUPAC(\%bases);
	}

	return;
}

sub consense_introns {
	my ($site, $coverage, $alignment) = @_;

	my %bases = %{$site->{BASES}};
	my %introns = %{$site->{INTRONS}};

	if (scalar(keys %introns) > 0) {

		# We only have one intron, output it
		if (scalar(keys %introns) == 1) {
			my $intron = lc((keys %introns)[0]);
			$$alignment .= $intron;

			$$coverage .= "(";
			foreach my $index (0 .. length($intron) - 1) {
				$$coverage .= $introns{(keys %introns)[0]}." ";
			}
			chop($$coverage);
			$$coverage .= ") ";
		}
		# We have multiple introns, align them
		else {

			# Create a temporary fasta file to hold introns
			open(my $intron_fasta, ">", "introns.fa");
			my $count = 0;
			foreach my $intron (keys %introns) {
				print {$intron_fasta} ">$count\n";	
				print {$intron_fasta} "$intron\n";	
				$count++;
			}
			close($intron_fasta);

			# Align introns with mafft
			my $return = system("$mafft --maxiterate 1000 --genafpair --thread 8 introns.fa > introns-aligned.fa 2> /dev/null");
			die "Error running mafft: '$return'.\n" if ($return);

			# Load aligned introns
			my %intron_sequences = parse_fasta("introns-aligned.fa");
			unlink("introns.fa", "introns-aligned.fa");

			# Determine the consensus base at each site of the intron
			$$coverage .= "(";
			foreach my $index (0 .. length((values %intron_sequences)[0]) - 1) {

				my %bases;
				my $base_count = 0;
				foreach my $intron (values %intron_sequences) {
					my $base = substr(uc($intron), $index, 1); 
					die "$index" if ($base eq '');

					# Not quite sure what to do with gaps, ignore for now
					#$bases{$base}++ if ($base ne "-");
					if ($base ne "-") {
						$bases{$base}++;
						$base_count++;
					}
				}

				# Get the IUPAC code for the current assortment of nucleotides
				$$alignment .= lc(get_IUPAC(\%bases));
				$$coverage .= scalar(keys %intron_sequences)." ";
			}
			chop($$coverage);
			$$coverage .= ") ";
		}
	}

	return;
}

sub get_IUPAC {
	my $bases = shift;

	my %bases = %{$bases};

	my $code;
	if (scalar(keys %bases) == 1) {
		$code = "A" if (exists($bases{"A"}));
		$code = "T" if (exists($bases{"T"}));
		$code = "C" if (exists($bases{"C"}));
		$code = "G" if (exists($bases{"G"}));
	}
	elsif (scalar(keys %bases) == 2) {
		$code = "R" if (exists($bases{"A"}) && exists($bases{"G"}));
		$code = "Y" if (exists($bases{"C"}) && exists($bases{"T"}));
		$code = "S" if (exists($bases{"G"}) && exists($bases{"C"}));
		$code = "W" if (exists($bases{"A"}) && exists($bases{"T"}));
		$code = "K" if (exists($bases{"G"}) && exists($bases{"T"}));
		$code = "M" if (exists($bases{"A"}) && exists($bases{"C"}));
	}
	elsif (scalar(keys %bases) == 3) {
		$code = "B" if (exists($bases{"C"}) && exists($bases{"G"}) && exists($bases{"T"}));
		$code = "D" if (exists($bases{"A"}) && exists($bases{"G"}) && exists($bases{"T"}));
		$code = "H" if (exists($bases{"A"}) && exists($bases{"C"}) && exists($bases{"T"}));
		$code = "V" if (exists($bases{"A"}) && exists($bases{"C"}) && exists($bases{"G"}));
	}
#	else {
#		$code = "N";
#	}
#
#	if (!defined($code)) {
#		use Data::Dumper;
#		print Dumper(\%bases);
#		die "IUPAC code not defined\n";
#	}

	$code = "N" if (!defined($code));

	return $code;
}

sub get_partial_seq {
	my ($range, $seq) = (@_);

	my $partial_seq;

	my @range = split(",", $range);
	foreach my $segment (@range) {
		my ($start, $end) = split("-", $segment);
		$partial_seq .= substr($seq, $start, $end - $start + 1);	
	}
	
	return $partial_seq;
}

sub intersection {
	my ($new, $old) = (@_);

	my @new = split(",", $new);
	my @old = split(",", $old);

	# Populate a hash with all members of @new and @old

	my %old;
	foreach my $segment (@old) {
		my ($start, $end) = split("-", $segment);
		foreach my $site ($start .. $end) {
			$old{$site}++;
		}
	}
	foreach my $segment (@new) {
		my ($start, $end) = split("-", $segment);
		foreach my $site ($start .. $end) {
			$old{$site}++;
		}
	}

	# Remove site if both arrays don't have it
	foreach my $site (keys %old) {
		if ($old{$site} != 2) {
			delete($old{$site});
		}
	}

	my @sites = sort { $a <=> $b } keys %old;

	return "0-0" if (scalar(@sites) == 0);

	# Convert sites to a string (i.e. 1,2,3,4,6,7,8 => 1-4,6-8)

	my @ends;
	my @starts = ($sites[0]);
	my $previous = $sites[0] - 1;
	foreach my $index (0 .. $#sites) {
		my $site = $sites[$index];
		if ($previous != $site - 1) {
			push(@ends, $sites[$index - 1]);
			push(@starts, $site);
		}
		$previous = $site;
	}
	push(@ends, $sites[$#sites]);

	my $intersection;
	foreach my $index (0 .. $#starts) {
		$intersection .= $starts[$index]."-".$ends[$index].",";
	}
	chop($intersection);

	return $intersection;
}

sub rev_comp {
	my $seq = shift;
	
	my %comp = ('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C');

	my $rev_comp;
	my $rev = reverse($seq);
	foreach my $index (0 .. length($rev) - 1) {
		my $char = substr($rev, $index, 1);
		#$rev_comp .= $comp{$char};
		if (exists($comp{$char})) {
			$rev_comp .= $comp{$char};
		}
		else {
			$rev_comp .= 'N';
		}
	}

	return $rev_comp;
}

sub check_path_for_exec {
	my ($exec, $continue) = @_;
	
	my $path = $ENV{PATH}.":."; # include current directory as well
	my @path_dirs = split(":", $path);

	my $exec_path;
	foreach my $dir (@path_dirs) {
		$dir .= "/" if ($dir !~ /\/$/);
		$exec_path = $dir.$exec if (-e $dir.$exec);
	}

	die "Could not locate: '$exec'. This script requires this program in your path.\n" if (!defined($exec_path) && !defined($continue));
	return $exec_path;
}

sub parse_fasta {
	my $filename = shift;

	my $taxon;
	my %align;
	open(my $alignment_file, '<', $filename) 
		or die "Could not open '$filename': $!\n";

	while (my $line = <$alignment_file>) {
		$line =~ s/^\s+|\s+$//g;

		# Taxon name
		#if ($line =~ /^>(.*)/) {
		if ($line =~ /^>(\S+)/) {
			$taxon = $1;
		}
		else {
			# Taxon sequence
			$taxon =~ s/-/_/g;
			$align{$taxon} .= uc($line);
		}
	}
	close($alignment_file);
	
	return %align;
}
