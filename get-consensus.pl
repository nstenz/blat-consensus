#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#my $overhang_length = 50;
my $overhang_length = 0;
#my $target_perc_cov_threshold = 1;
my $target_perc_cov_threshold = 0.40;

my $blat = check_path_for_exec("blat");
my $cap3 = check_path_for_exec("cap3");
my $muscle = check_path_for_exec("muscle");
my $mafft = check_path_for_exec("mafft");
my $blastn = check_path_for_exec("blastn");
my $makeblastdb = check_path_for_exec("makeblastdb");

my $targets;
GetOptions(
	't|targets=s' => \$targets,
);

my $input = shift;

die "You must specify an assembly as input.\n" if (!defined($input));
die "You must specify a fasta file containing target sequences (-t).\n" if (!defined($targets));
die "You must you have specified an incorrect argument.\n" if (@ARGV);

my $output_fasta = $input.".hits";
#my $output_psl = $input.".hits.psl";
my $blat_out_name = $input.".psl";
#my $blat_out_name = $input.".blat";
my $blast_out_name = $input.".blast";

system("$blat $input $targets -t=dna -q=dna -noHead $blat_out_name");
#system("$blat $input $targets -t=dna -q=dna -noHead -tileSize=5 $blat_out_name");
#system("$blat $input $targets -t=dna -q=dna -noHead $blat_out_name") if (!-e $blat_out_name);
#system("$blat $input $targets -t=dnax -q=dnax -noHead $blat_out_name");

my %targets;
open(my $blat_out, "<", $blat_out_name);
while (my $line = <$blat_out>) {
	chomp($line);

	my @line = split(/\s+/, $line);

	my ($match, $mismatch, $rep_matches, $n_count, $q_num_inserts, $q_base_inserts, 
	    $t_num_inserts, $t_base_inserts, $strand, $q_name, $q_size, $q_start, $q_end,
		$t_name, $t_size, $t_start, $t_end, $block_count, $block_sizes, $q_starts, $t_starts) = @line; 

	if (($match + $mismatch) / $q_size >= $target_perc_cov_threshold) {
#		print "threshold: $target_perc_cov_threshold, ".(($match + $mismatch) / $q_size)."\n";
#		$contigs{$t_name}++;
#		print {$output} "$line\n";
		#push(@{$targets{$q_name}}, {'NAME' => $t_name, 'STRAND' => $strand, 'STARTS' => $t_starts, 'SIZES' => $block_sizes});
		push(@{$targets{$q_name}}, {'NAME' => $t_name, 'STRAND' => $strand, 'T_STARTS' => $t_starts, 'Q_STARTS' => $q_starts, 'SIZES' => $block_sizes});
	}
}
close($blat_out);

my %contig_seqs = parse_fasta($input);
my %target_seqs = parse_fasta($targets);

#mkdir("sequences");
my $count = 0;
foreach my $target (sort { $a cmp $b } keys %targets) {
	my @hits = @{$targets{$target}};

	my @alignment;
	my %alignment;

	# Add target bases to alignment
	my $target_seq = $target_seqs{$target};
	foreach my $index (0 .. length($target_seq) - 1) {
		my $base = substr($target_seq, $index, 1);
		$alignment{$index}->{"TARGET"} = $base;
	}

	foreach my $hit (@hits) {
		print $hit->{'NAME'},"\n";

		my $seq = $contig_seqs{$hit->{'NAME'}};
		my $contig_length = length($seq);

		my @sizes = split(",", $hit->{'SIZES'});
		my @t_starts = split(",", $hit->{'T_STARTS'});
		my @q_starts = split(",", $hit->{'Q_STARTS'});

		if ($hit->{'STRAND'} eq "-") {
			print "  contig will be reverse complemented\n";
		}

		#foreach my $index (0 .. $#starts - 1) {
		#foreach my $index (0 .. $#t_starts - 1) {
		foreach my $index (0 .. $#t_starts) {
			my $size = $sizes[$index];
			#my $start = $starts[$index];
			my $q_start = $q_starts[$index];
			my $t_start = $t_starts[$index];
			#my $next_start = $starts[$index + 1];

		#	my $intron_start = $start + $size;
		#	my $intron_length = $next_start - $start - $size - 1;
			#my $intron_start = $t_start + $size;
			my $intron_start = $t_start + $size - 1;

			# If this is the last (or only) block, we set the next block 
			# "start" to a value which makes the intron length 0
			my $next_t_start;
			if ($index + 1 < scalar(@t_starts)) {
				$next_t_start = $t_starts[$index + 1];
			}
			else {
				# only one block
				$next_t_start = $intron_start + 1;
			}

			#my $intron_length = $next_t_start - $t_start - $size - 1;
			#my $intron_length = $next_t_start - $intron_start - 1;
			#my $intron_length = $next_t_start - $intron_start + 1;
			my $intron_length = $next_t_start - $intron_start;

		#	print "  match_start : $start\n";
		#	print "  match_length : $size\n";
		#	print "  match_end : ",($start + $size - 1),"\n";
			print "  match_start : $t_start ($q_start)\n";
			print "  match_length : $size\n";
			print "  match_end : ",($t_start + $size - 1),"\n";

			print "  intron_start : $intron_start, intron_length : $intron_length\n\n";

			#foreach my $index ($start .. $start + $size - 1) {
			#foreach my $offset (0 .. $size - 1) {
			foreach my $offset (0 .. $size - 2) {
				my $t_index = $t_start + $offset;
				my $q_index = $q_start + $offset;

				my $char = substr($seq, $t_index, 1);
				if ($hit->{'STRAND'} eq "-") {
					$alignment{$q_index}->{rev_comp($char)}++;
				}
				else {
					$alignment{$q_index}->{$char}++;
				}
			}

			#$alignment{$intron_start}->{substr($seq, $intron_start, $intron_length)}++;

			if ($intron_length > 0) {
				if ($hit->{'STRAND'} eq "-") {
					#$alignment{$intron_start}->{rev_comp(substr($seq, $intron_start, $intron_length))}++;
					#$alignment{$intron_start - $t_starts[0] + $q_starts[0]}->{rev_comp(substr($seq, $intron_start, $intron_length))}++;

					$alignment{$q_start + $intron_start - $t_start}->{rev_comp(substr($seq, $intron_start, $intron_length))}++;
				}
				else {
					#$alignment{$intron_start}->{substr($seq, $intron_start, $intron_length)}++;
					#$alignment{$intron_start - $t_starts[0] + $q_starts[0]}->{substr($seq, $intron_start, $intron_length)}++;
					$alignment{$q_start + $intron_start - $t_start}->{substr($seq, $intron_start, $intron_length)}++;
				}
			}
		}

	#	my $start = $starts[0];
	#	my $end = $starts[$#starts] + $sizes[$#sizes];
		my $start = $t_starts[0];
		my $end = $t_starts[$#t_starts] + $sizes[$#sizes];

		print "  length:",length($seq),"\n";
		#print "  num gaps:",scalar(@starts),"\n";
		print "  num blocks:",scalar(@t_starts),"\n";
		print "  init_start: $start\n";
		print "  init_end: $end\n";

#		# Add overhangs
#		$start = ($start - $overhang_length < 0) ? (0) : ($start - $overhang_length);
#		$end = ($end + $overhang_length > $contig_length) ? ($contig_length) : ($end + $overhang_length);

		print "  final_start: $start\n";
		print "  final_end: $end\n\n";

#		print {$fasta} ">",$hit->{'NAME'},"\n";
#		my $subseq = substr($seq, $start, $end - $start + 1);

#		if ($hit->{'STRAND'} eq "-") {
#			print "  reverse complementing strand\n";
##			$subseq = rev_comp($subseq);
#		}
#		print "\n";

#		for (my $i = 0; $i < scalar(@sizes); $i++) {
#			my $size = $sizes[$i];
#			my $start = $starts[$i];
#
#			print "  match from $start to ",($start + $size - 1),"\n";
#
#			my $subseq = substr($seq, $start, $size);
#			if ($hit->{'STRAND'} eq "-") {
#				$subseq = rev_comp($subseq);
#			}
#
#			print "  ",length($subseq)," == $size\n";
#			print {$fasta} $subseq;
#		}
#		print {$fasta} "\n";
		#print {$fasta} "$subseq\n";
	}
#	close($fasta);

	use Data::Dumper;
	#$Data::Dumper::Sortkeys = sub { [sort { $a <=> $b } keys %{$_[0]}] };
	$Data::Dumper::Sortkeys = sub { [sort { $b <=> $a } keys %{$_[0]}] };

	#print Dumper(@alignment),"\n";
	print Dumper(\%alignment),"\n";
	my $alignment = consense_alignment(\%alignment);
	print "\n\n",$alignment,"\n\n";

	$count ++;
	die if ($alignment =~ /a|t|c|g/);
	#die if $count > 1;
}

sub consense_alignment {
	my $align = shift;

	my $alignment;

	my %align = %{$align};
	SITE: foreach my $site (sort { $a <=> $b } keys %align) {
		my %bases_present = %{$align{$site}};
		#print "$site(",scalar(keys %bases_present),") ";	

		# Only the target sequence's nucleotide present, no contigs mapped to this site
		if (scalar(keys %bases_present) == 1) {
			$alignment .= "?";
		}
		else {
			my $target_base = $bases_present{"TARGET"};

			my %introns;
			my %single_bases;

			#print "target base = $target_base\n";

			# Check whether all mapped contigs have the same bases at the site
			#my $all_equal_target = 1;
#			my $total = 0;
#			my $support = 0;
#			my $against = 0;
			foreach my $base_present (keys %bases_present) {
				next if ($base_present eq "TARGET");
				#$all_equal_target = 0 if ($base_present ne $target_base);

				# Perhaps filter by proportion here? (remove base if it is present in < x% of contigs)
			#	$support++ if ($base_present eq $target_base);
			#	$against++ if ($base_present ne $target_base);
			#	$total++;

				# Intron located here
				if (length($base_present) > 1) {
					$introns{$base_present} = $bases_present{$base_present};
				}
				else {
					$single_bases{$base_present} = $bases_present{$base_present};
				}

				# Ambiguous, can't determine consensus
				if ($base_present eq "N") {
					$alignment .= "N";
					next SITE;
				}
			}

			# If only one nucleotide is found at the site output it
			#if ($all_equal_target) {
			#if (($support / $total == 1 || $against / $total == 1) && scalar(keys %bases_present) == 2) {
			#if (scalar(keys %bases_present) == 2) {
			if (scalar(keys %single_bases) == 1 && scalar(keys %introns) == 0) {
				$alignment .= (keys %single_bases)[0];
			}
			# Handle ambiguous non-introns with IUPAC codes
			elsif (scalar(keys %introns) == 0) {
				$alignment .= get_IUPAC(\%single_bases);
			}
			else {
				#print "introns: ".scalar(keys %introns),"\n";

			#	# Check if we have single bases too
			#	if (scalar(keys %single_bases)) {
			#		print "what to do, what to do?\n";
			#		die;
			#	}
			#	# Only introns (> 1 base)
			#	else {
					if (scalar(keys %introns) == 1) {
						$alignment .= lc((keys %introns)[0]);
					}
					else {
						open(my $intron_fasta, ">", "introns.fa");
						my $count = 0;
						foreach my $intron (keys %introns) {
							print {$intron_fasta} ">$count\n";	
							print {$intron_fasta} "$intron\n";	
							$count++;
						}
						close($intron_fasta);

						my $return = system("$mafft --maxiterate 1000 --genafpair --thread 40 introns.fa > introns-aligned.fa 2> /dev/null");
						die if ($return);

						my %intron_sequences = parse_fasta("introns-aligned.fa");
						unlink("introns.fa", "introns-aligned.fa");

						foreach my $index (0 .. length((values %introns)[0]) - 1) {

							my %bases;
							foreach my $intron (values %intron_sequences) {
								my $base = substr(uc($intron), $index, 1); 

								# Not quite sure what to do with gaps, ignore for now
								$bases{$base}++ if ($base ne "-");
							}
							$alignment .= lc(get_IUPAC(\%bases));
						}
					}
			#	}
			}
		}
	}

	return $alignment;
}

sub get_IUPAC {
	my $bases = shift;

	my %bases = %{$bases};

	my $code;
	if (scalar(keys %bases) == 1) {
		$code .= "A" if (exists($bases{"A"}));
		$code .= "T" if (exists($bases{"T"}));
		$code .= "C" if (exists($bases{"C"}));
		$code .= "G" if (exists($bases{"G"}));
	}
	elsif (scalar(keys %bases) == 2) {
		$code .= "R" if (exists($bases{"A"}) && exists($bases{"G"}));
		$code .= "Y" if (exists($bases{"C"}) && exists($bases{"T"}));
		$code .= "S" if (exists($bases{"G"}) && exists($bases{"C"}));
		$code .= "W" if (exists($bases{"A"}) && exists($bases{"T"}));
		$code .= "K" if (exists($bases{"G"}) && exists($bases{"T"}));
		$code .= "M" if (exists($bases{"A"}) && exists($bases{"C"}));
	}
	elsif (scalar(keys %bases) == 3) {
		$code .= "B" if (exists($bases{"C"}) && exists($bases{"G"}) && exists($bases{"T"}));
		$code .= "D" if (exists($bases{"A"}) && exists($bases{"G"}) && exists($bases{"T"}));
		$code .= "H" if (exists($bases{"A"}) && exists($bases{"C"}) && exists($bases{"T"}));
		$code .= "V" if (exists($bases{"A"}) && exists($bases{"C"}) && exists($bases{"G"}));
	}
	else {
		$code .= "N";
	}

	#die "IUPAC code not defined\n" if (!defined($code));
	use Data::Dumper;
	if (!defined($code)) {
		print Dumper(\%bases);
		die "IUPAC code not defined\n" if (!defined($code));
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
		$rev_comp .= $comp{$char};
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
