#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# Not yet implemented
#my $overhang_length = 50;
my $overhang_length = 0;
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
#my $return = system("$blat $input $targets -t=dna -q=dna -noHead $blat_out_name") if (!-e $blat_out_name);

# Check that things worked correctly
if ($return) {
	die "Error running blat: '$return'.\n";
}
else {
	print "Successfully completed blat.\n";
}

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
		push(@{$targets{$q_name}}, {'NAME' => $t_name, 'STRAND' => $strand, 'T_STARTS' => $t_starts, 'Q_STARTS' => $q_starts, 'SIZES' => $block_sizes});
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
	my @hits = @{$targets{$target}};

	# Holds bases/introns in contigs which map to specific sites in target
	my %alignment;

	# Add target bases to alignment
	my $target_seq = $target_seqs{$target};
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

		my @sizes = split(",", $hit->{'SIZES'});
		my @t_starts = split(",", $hit->{'T_STARTS'});
		my @q_starts = split(",", $hit->{'Q_STARTS'});

		# Output strandedness
		if ($hit->{'STRAND'} eq "-") {
			print "    Contig will be reverse complemented (- hit).\n";
		}
		else {
			print "    Contig will not be reverse complemented (+ hit).\n";
		}

		# Iterate through aligned blocks of target and contig
		foreach my $index (0 .. $#t_starts) {
			my $size = $sizes[$index];
			my $q_start = $q_starts[$index];
			my $t_start = $t_starts[$index];

			# Where the intron starts 
			#my $intron_start = $t_start + $size - 1;
			my $intron_start = $t_start + $size;

			# If this is the last (or only) block, we set the next block 
			# "start" to a value which makes the intron length 0
			my $next_t_start;
			if ($index + 1 < scalar(@t_starts)) {
				$next_t_start = $t_starts[$index + 1];
			}
			else {
				# only/last block
				$next_t_start = $intron_start + 1;
			}

			# Determine how long the intron is
			#my $intron_length = $next_t_start - $intron_start;
			my $intron_length = $next_t_start - $intron_start + 1;

			print "    match_start : $t_start ($q_start)\n";
			print "    match_length : $size\n";
			print "    match_end : ",($t_start + $size - 1),"\n";
			#print "    intron_start : $intron_start, intron_length : $intron_length\n\n";
			print "    intron_start : $intron_start (",($q_start + $intron_start - $t_start),"), intron_length : $intron_length\n\n";

			# Place aligned characters into %alignment
			#foreach my $offset (0 .. $size - 2) {
			foreach my $offset (0 .. $size - 1) {
				my $t_index = $t_start + $offset;
				my $q_index = $q_start + $offset;

				my $char = substr($seq, $t_index, 1);

				# Reverse complement character if needed
				if ($hit->{'STRAND'} eq "-") {
					#$alignment{$q_index}->{rev_comp($char)}++;
					#$alignment{$q_index}->{"BASES"}->{rev_comp($char)}++;
					$alignment{$q_index}->{"BASES"}{rev_comp($char)}++;
				}
				else {
					#$alignment{$q_index}->{$char}++;
					$alignment{$q_index}->{"BASES"}{$char}++;
				}
			}

			# Place introns into %alignment
			if ($intron_length > 0) {

				# Reverse complement intron if needed
				my $intron;
				if ($hit->{'STRAND'} eq "-") {
					#$intron = rev_comp(substr($seq, $intron_start, $intron_length));
					$intron = rev_comp(substr($seq, $intron_start + 1, $intron_length - 1));
				}
				else {
					#$intron = substr($seq, $intron_start, $intron_length);
					$intron = substr($seq, $intron_start + 1, $intron_length - 1);
				}
				#$alignment{$q_start + $intron_start - $t_start}->{$intron}++;
				#$alignment{$q_start + $intron_start - $t_start}->{"INTRONS"}->{$intron}++;
				$alignment{$q_start + $intron_start - $t_start}->{"INTRONS"}{$intron}++;

				#my $intron_first_char = substr($intron, 0, 1);
				#$alignment{$q_start + $intron_start - $t_start}->{$intron_first_char}--;

				#delete($alignment{$q_start + $intron_start - $t_start}->{$intron_first_char}) if ($alignment{$q_start + $intron_start - $t_start}->{$intron_first_char} == 0);
			}
		}

#		my $start = $t_starts[0];
#		my $end = $t_starts[$#t_starts] + $sizes[$#sizes];
#
#		print "  length:",length($seq),"\n";
#		print "  num blocks:",scalar(@t_starts),"\n";
#		print "  init_start: $start\n";
#		print "  init_end: $end\n";
#
##		# Add overhangs
##		$start = ($start - $overhang_length < 0) ? (0) : ($start - $overhang_length);
##		$end = ($end + $overhang_length > $contig_length) ? ($contig_length) : ($end + $overhang_length);
#
#		print "  final_start: $start\n";
#		print "  final_end: $end\n\n";
	}

#	use Data::Dumper;
#	#$Data::Dumper::Sortkeys = sub { [sort { $a <=> $b } keys %{$_[0]}] };
#	$Data::Dumper::Sortkeys = sub { [sort { $b <=> $a } keys %{$_[0]}] };
#
#	#print Dumper(@alignment),"\n";
#	print Dumper(\%alignment),"\n";
	#my $alignment = consense_alignment(\%alignment);
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

sub consense_alignment {
	my $align = shift;

	my $coverage;
	my $alignment;

	my %align = %{$align};
	SITE: foreach my $site (sort { $a <=> $b } keys %align) {
		my %bases_present = %{$align{$site}};

		# Only the target sequence's nucleotide present, no contigs mapped to this site
		if (scalar(keys %bases_present) == 1) {
			$coverage .= "0 ";
			$alignment .= "?";
		}
		else {
			my $target_base = $bases_present{"TARGET"};

		#	my %introns;
		#	my %single_bases;
			my %introns = %{$bases_present{"INTRONS"}};
			my %single_bases = %{$bases_present{"BASES"}};

			my $num_single_bases = 0;
			#my $site_coverage = 0;
			#foreach my $base_present (keys %bases_present) {
			foreach my $base_present (keys %single_bases) {
				#next if ($base_present eq "TARGET");

			#	if (length($base_present) == 1) {
			#		$num_single_bases += $bases_present{$base_present};
			#	}
				$num_single_bases += $single_bases{$base_present};
			}

			#foreach my $base_present (keys %bases_present) {
			foreach my $base_present (keys %single_bases) {
				#next if ($base_present eq "TARGET");

				# Perhaps filter by proportion here? (remove base if it is present in < x% of contigs)

			#	# Intron located here
			#	if (length($base_present) > 1) {
			#		$introns{$base_present} = $bases_present{$base_present};
			#	}
			#	else {
			#		$single_bases{$base_present} = $bases_present{$base_present};
			#	}

				# Ambiguous, can't determine consensus
				if ($base_present eq "N") {
					$coverage .= "$num_single_bases ";
					$alignment .= "N";
					next SITE;
				}
			}

			# If only one nucleotide is found at the site output it
			#if (scalar(keys %single_bases) == 1 && scalar(keys %introns) == 0) {
			if (scalar(keys %single_bases) == 1) {
				$coverage .= "$num_single_bases ";
				$alignment .= (keys %single_bases)[0];
			}
			# Handle ambiguous non-introns with IUPAC codes
			#elsif (scalar(keys %introns) == 0) {
			else {
				$coverage .= "$num_single_bases ";
				$alignment .= get_IUPAC(\%single_bases);
			}

			# We have introns
			#else {
			if (scalar(keys %introns) > 0) {

				# We only have one intron, output it
				if (scalar(keys %introns) == 1) {
					#$alignment .= lc((keys %introns)[0]);
					my $intron = lc((keys %introns)[0]);
					#$alignment .= ucfirst($intron);
					$alignment .= $intron;

					$coverage .= "(";
					foreach my $index (0 .. length($intron) - 1) {
						#$coverage .= $bases_present{(keys %introns)[0]}." ";
						$coverage .= $introns{(keys %introns)[0]}." ";
					}
					chop($coverage);
					$coverage .= ") ";
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
					my $return = system("$mafft --maxiterate 1000 --genafpair --thread 40 introns.fa > introns-aligned.fa 2> /dev/null");
					die "Error running mafft: '$return'.\n" if ($return);

					# Load aligned introns
					my %intron_sequences = parse_fasta("introns-aligned.fa");
					#unlink("introns.fa", "introns-aligned.fa");

					# Determine the consensus base at each site of the intron
					$coverage .= "(";
					foreach my $index (0 .. length((values %introns)[0]) - 1) {

						my %bases;
						foreach my $intron (values %intron_sequences) {
							my $base = substr(uc($intron), $index, 1); 
							die "$index" if ($base eq '');

							# Not quite sure what to do with gaps, ignore for now
							$bases{$base}++ if ($base ne "-");
						}

						# Get the IUPAC code for the current assortment of nucleotides
					#	if ($index == 0) {
					#		$coverage .= scalar(keys %intron_sequences)." ";
					#		$alignment .= get_IUPAC(\%bases);
					#	}
					#	else {
					#		$coverage .= scalar(keys %intron_sequences)." ";
					#		$alignment .= lc(get_IUPAC(\%bases));
					#	}
						$alignment .= lc(get_IUPAC(\%bases));
						$coverage .= scalar(keys %intron_sequences)." ";
					}
					chop($coverage);
					$coverage .= ") ";
				}
			}
		}
	}

	#return $alignment;
	return ($alignment, $coverage);
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
	else {
		$code = "N";
	}

	if (!defined($code)) {
		use Data::Dumper;
		print Dumper(\%bases);
		die "IUPAC code not defined\n";
	}

	return $code;
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
		chomp($line);

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
