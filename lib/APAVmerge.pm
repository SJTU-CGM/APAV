#!/usr/bin/perl

package APAVmerge;

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw/sum max min/;


sub mergeElePAV{

	my $usage = "\n\tUsage: apav mergeElePAV --pav <elecov_file> [option] 

'apav mergeElePAV' is used to merge neighboring elements with the same PAV.

Necessary input description:
  -i, --pav  	<file>          Element PAV file.

Options:
  -o, --out     <string>        Output file name.

  -h, --help                    Print usage page.
  \n";

	my $cmdline = $0. " mergeElePAV";
	$cmdline .= " ".join(" ", @ARGV);

	my $stime = `date +"%Y-%m-%d %H:%M:%S"`;
        chomp($stime);

	my ($pav, $out, $help);
	GetOptions(
                'pav|i=s'    	=> \$pav,
                'out|o=s'       => \$out,
                'help|h!'       => \$help
        ) or die $!."\n";

	die $usage if !defined($pav);
        die $usage if $help;
	APAVutils::check_file('--pav/-i', $pav);
	$out = APAVutils::check_out($out, $pav, ".mpav");

	open(PAV, "<$pav") or die "Could not open file '$pav'\n";
        open(OUT, ">$out") or die "Could not open file '$out'\n";

	print STDOUT "[".$stime."] [mergeElePAV] Merge neighboring elements with the same PAV...\n";
        print OUT "##Date: $stime\n";
        print OUT "##CMDline: $cmdline\n";
	print OUT "##Args: --pav $pav --out $out\n";

	my $cur_gene;
	my @tmp_blocks;
	my @tmp_blocks_pavs;

	while(<PAV>){
		if($_ =~ /^#/){
                        next;
                }
                if($_ =~ /^Chr\tStart\tEnd\tLength\tAnnotation/){
                        print OUT $_;
                        next;
                }
                chomp $_;
		my @info = split(/\t/, $_);
                my $chr = $info[0];
                my $start = $info[1];
                my $end = $info[2];
                my $len = $info[3];
                my $anno = $info[4];
                my @pavs = @info[(5)..($#info)];

		die "It doesn't look like a PAV profile. It should only contain 0 and 1 in sample columns.\n" if(!isAll01(\@pavs));
		if($anno =~ /^([^:]*):(\[.*\])$/){
			$cur_gene = $1 if !$cur_gene;
			my $block = $2;
			if($1 ne $cur_gene){
				printBlock(\@tmp_blocks, \@tmp_blocks_pavs, $cur_gene) if(@tmp_blocks);
				@tmp_blocks = ();
				@tmp_blocks_pavs = ();
				$cur_gene = $1;
			}
			if(@tmp_blocks){
				if(!isSame(\@pavs, \@tmp_blocks_pavs)){
					printBlock(\@tmp_blocks, \@tmp_blocks_pavs, $cur_gene);
					@tmp_blocks = ();
					@tmp_blocks_pavs = ();
				}	
			}
			push(@tmp_blocks, [$chr, $start, $end, $len, $block]);
			@tmp_blocks_pavs = @pavs;

		}else{
			die "Unexpected format in column 'Annotation'\n";
		}
	}

	printBlock(\@tmp_blocks, \@tmp_blocks_pavs, $cur_gene);

	close(PAV);
	close(OUT);

	my $etime = `date +"%Y-%m-%d %H:%M:%S"`;
        chomp($etime);
        print STDOUT "[-- ".$etime." --] apav mergeElePAV is finished\n";

}


sub printBlock{
	my ($blocks, $pavs, $g_id) = @_;
        my @starts = map{@$_[1]} @$blocks;
        my @ends = map{@$_[2]} @$blocks;
        my @lens = map{@$_[3]} @$blocks;
        my @ids = map{@$_[4]} @$blocks;
	print OUT @$blocks[0]->[0]."\t".join(";", @starts)."\t".join(";", @ends)."\t".join(";", @lens)."\t".$g_id.":".join(";", @ids)."\t".join("\t", @$pavs)."\n";
}



sub mergeElecov{

	my $usage = "\n\tUsaege: apav mergeElecov --elecov <elecov_file> [options]

'apav merge' is used to merge neighboring elements with the same coverage.

Necessary input description:
  -i, --elecov  <file>	        Element coverage file.

Options:
  --asgene			Treating elements as gene elements.
  -o, --out     <string>        Output file name.
  -h, --help                    Print usage page.
  \n";

  	my $cmdline = $0. " mergeElecov";
        $cmdline .= " ".join(" ", @ARGV);

        my $stime = `date +"%Y-%m-%d %H:%M:%S"`;
        chomp($stime);

	my ($cov, $out, $help);
	my $asgene = 0;
	GetOptions(
		'elecov|i=s'	=> \$cov,
		'out|o=s'	=> \$out,
		'asgene!'	=> \$asgene,
		'help|h!'	=> \$help
	) or die $!."\n";

	die $usage if !defined($cov);
	die $usage if $help;
	APAVutils::check_file('--elecov/-i', $cov);
	$out = APAVutils::check_out($out, $cov, ".mcov");

	print STDOUT "[".$stime."] [mergeElecov] Merge neighboring elements with the same coverage...\n";

	mergecov($cov, $out, $asgene);

	my $etime = `date +"%Y-%m-%d %H:%M:%S"`;
        chomp($etime);
        print STDOUT "[".$etime."] [mergeElecov] Finished\n";
}

sub mergecov{
	my ($cov, $out, $asgene) = @_;

	open(COV, "<$cov") or die "Could not open file '$cov'\n";
	open(OUT, ">$out") or die "Could not open file '$out'\n";

	if($asgene){
		mergeGene();
	}else{
		mergeItem();
	}

	close(COV);
	close(OUT);

}


sub mergeGene{
        my $cur_gene;
        my @tmp_streams;
        my @tmp_streams_covs;
	my @tmp_eles;
	my @tmp_eles_covs;
	my $tag;
	
	while(<COV>){
		if($_ =~ /^#/){
			print OUT $_;
			next;
		}
		if($_ =~ /^Chr\tStart\tEnd\tLength\tAnnotation/){
			print OUT $_;
			next;
		}
		chomp $_;
                my @info = split(/\t/, $_);
                my $chr = $info[0];
                my $start = $info[1];
                my $end = $info[2];
                my $len = $info[3];
                my $anno = $info[4];
                my @covs = @info[(5)..($#info)];

		die "It seems to have been merged\n" if($start =~ /,/ | $end =~ /,/);

		if(!isAll01(\@covs)){
			printStream(\@tmp_streams, \@tmp_streams_covs) if(@tmp_streams);
			printEle(\@tmp_eles, \@tmp_eles_covs, $cur_gene) if(@tmp_eles);
			@tmp_streams = ();
			@tmp_streams_covs = ();
			@tmp_eles = ();
			@tmp_eles_covs = ();
			print OUT "$_\n";	
		}else{
                	if($anno =~ /^([^:]+):\[(DOWN|UP):(\d+)-(\d+)\]/){
				$cur_gene = $1 if !$cur_gene;
				$tag = $2 if !$tag;
				if($1 ne $cur_gene || $2 ne $tag){
					printStream(\@tmp_streams, \@tmp_streams_covs) if(@tmp_streams);
					printEle(\@tmp_eles, \@tmp_eles_covs, $cur_gene) if(@tmp_eles);
					@tmp_streams = ();
					@tmp_streams_covs = ();
                	                @tmp_eles = ();
					@tmp_eles_covs = ();
        	                        $cur_gene = $1;
					$tag = $2;
                        	}	
				if(@tmp_streams){
					if(!isSame(\@covs, \@tmp_streams_covs)){
						printStream(\@tmp_streams, \@tmp_streams_covs);
						@tmp_streams = (); 
						@tmp_streams_covs = ();
					}
				}
				push(@tmp_streams, [$chr, $start, $end, $len, $1, $2, $3, $4]);
				@tmp_streams_covs = @covs;
			}elsif($anno =~ /^([^:]+):/){
				$cur_gene = $1 if !$cur_gene;
				$tag = "ele" if !$tag;
				if($1 ne $cur_gene || $tag ne "ele"){
					printStream(\@tmp_streams, \@tmp_streams_covs) if(@tmp_streams);
        	                        printEle(\@tmp_eles, \@tmp_eles_covs, $cur_gene) if(@tmp_eles);
                                	@tmp_streams = ();
					@tmp_streams_covs = ();
					@tmp_eles = ();
					@tmp_eles_covs = ();
					$cur_gene = $1;
                        	        $tag = "ele";
				}
				if(@tmp_eles){
					if(!isSame(\@covs, \@tmp_eles_covs)){
						printEle(\@tmp_eles, \@tmp_eles_covs, $cur_gene);
						@tmp_eles = ();
						@tmp_eles_covs = ();
					}
				}
				my $item;
				if($anno =~ /^.*\[(.*)\]$/){
					$item = $1;
				}else{
					die "Unexpected format in column 'Annotation'\n";
				}
				push(@tmp_eles, [$chr, $start, $end, $len, $item]);
				@tmp_eles_covs = @covs;
			}
		}
	}
	printEle(\@tmp_eles, \@tmp_eles_covs, $cur_gene) if(@tmp_eles);
	printStream(\@tmp_streams, \@tmp_streams_covs) if(@tmp_streams);
}

sub printStream{
	my ($info, $covs) = @_;
        my $start = @$info[0]->[1];
        my $end = @$info[-1]->[2];
        my @lens = map{@$_[3]} @$info;
	my $len = sum @lens;
	my @range = (@$info[0]->[6], @$info[0]->[7], @$info[-1]->[6], @$info[-1]->[7]);
	my $min = min @range;
	my $max = max @range;
        print OUT @$info[0]->[0]."\t".$start."\t".$end."\t".$len."\t".@$info[0]->[4].":"."[".@$info[0]->[5].":".$min."-".$max."]"."\t".join("\t",@$covs)."\n";
}

sub printEle{
	my ($eles, $covs, $g_id) = @_;
	my @starts = map{@$_[1]} @$eles;
        my @ends = map{@$_[2]} @$eles;
        my @lens = map{@$_[3]} @$eles;
	my @ids = map{split(",", @$_[4])} @$eles;
	my ($fs, $fe, $lens) = mergeElePos(\@starts, \@ends);
	print OUT @$eles[0]->[0]."\t".$fs."\t".$fe."\t".$lens."\t".$g_id.":[".mergeEleID(\@ids)."]"."\t".join("\t", @$covs)."\n";
}


sub mergeElePos{
	my ($starts, $ends) = @_;
	my $ts = @$starts[0];
	my $te = @$ends[0];
	my (@fs, @fe, @lens);
	for my $i (1..(@$starts-1)){
		if(@$starts[$i] == $te + 1){
			$te = @$ends[$i];
		}elsif(@$starts[$i] > $te){
			push(@fs, $ts);
			push(@fe, $te);
			$ts = @$starts[$i];
			$te = @$ends[$i];
		}else{
			$te = @$ends[$i] if(@$ends[$i] > $te);
		}
	}
	push(@fs, $ts);
	push(@fe, $te);
	for my $i (0..$#fs){
		push(@lens, $fe[$i] - $fs[$i] + 1);
	}
	return(join(",", @fs), join(",", @fe), join(",", @lens));
}

sub mergeEleID{
        my ($d) = @_;
        my %h;
        foreach my $i (@$d){
                if($i =~ /^T(\d+):exon(\d+)(.*)/){
                        if(exists($h{$1})){  ## transcript
                                if(exists($h{$1}{$2})){
                                        push(@{$h{$1}{$2}}, $3);
                                }else{
                                        $h{$1}{$2} = [$3];
                                }
                        }else{
                                $h{$1}= {$2 => [$3]};
                        }
                }
        }
        my @trans;
        foreach my $tran (sort{$a <=> $b} keys %h){
                my @exons;
                foreach my $exon (sort{$a <=> $b} keys %{$h{$tran}}){
                        if($#{$h{$tran}{$exon}}){
                                map{$_ =~ s/://}@{$h{$tran}{$exon}};
                                $h{$tran}{$exon} = "(".join("&", @{$h{$tran}{$exon}}).")";
                        }else{
                                $h{$tran}{$exon} = $h{$tran}{$exon}[0];
                        }
                        push(@exons, "exon".$exon.$h{$tran}{$exon});
                }
		push(@trans, "T".$tran.":".join("+", @exons))
        }
	return(join(",", @trans))
}

sub mergeItem{
	my $cur_item;
	my @tmp_parts;
	my @tmp_parts_covs;

	while(<COV>){
		if($_ =~ /^#/){
                        next;
                }
                if($_ =~ /^Chr\tStart\tEnd\tLength\tAnnotation/){
                        print OUT $_;
                        next;
                }
			chomp $_;
			my @info = split(/\t/, $_);
			my $chr = $info[0];
			my $start = $info[1];
			my $end = $info[2];
			my $len = $info[3];
			my $part = $info[4];
			my @covs = @info[5..$#info];
			
			if($part =~ /^(.*)\:\[(.*)\]$/){
				$cur_item = $1 if !$cur_item;
				if($1 ne $cur_item){
					printParts(\@tmp_parts, \@tmp_parts_covs) if @tmp_parts;
					$cur_item = $1;
					@tmp_parts = (); 
					@tmp_parts_covs = ();
				}
				if(isAll01(\@covs)){
					if(@tmp_parts){
						if(!isSame(\@covs, \@tmp_parts_covs)){
							printParts(\@tmp_parts, \@tmp_parts_covs);
							@tmp_parts= ();
							@tmp_parts_covs = ();
						}
					}
					push(@tmp_parts, [$chr, $start, $end, $len, $1, $2]);
					@tmp_parts_covs = @covs;
				}else{
					if (@tmp_parts){
						printParts(\@tmp_parts, \@tmp_parts_covs);
						@tmp_parts = ();
						@tmp_parts_covs = ();
					}
					print OUT "$_\n";
				}
			}
	
	}
	printParts(\@tmp_parts, \@tmp_parts_covs) if @tmp_parts;
	
}


sub printParts {
	my ($parts, $covs) = @_;
	my @starts = map{@$_[1]} @$parts;
	my @ends = map{@$_[2]} @$parts;
	my @lens = map{@$_[3]} @$parts;
	if(scalar(@$parts) eq 1){
		print OUT @$parts[0]->[0]."\t".join(",", @starts)."\t".join(",", @ends)."\t".join(",", @lens)."\t".@$parts[0]->[4].":[".@$parts[0]->[5]."]"."\t".join("\t",@$covs)."\n";
	}else{
		print OUT @$parts[0]->[0]."\t".join(",", @starts)."\t".join(",", @ends)."\t".join(",", @lens)."\t".@$parts[0]->[4].":[".@$parts[0]->[5]."-".@$parts[-1]->[5]."]"."\t".join("\t",@$covs)."\n";
	}
}


sub isAll01 {
	my ($data) = @_;
	my $f = 1;
	foreach my $d (@$data){
		if ($d != 0 && $d != 1){
			$f = 0;
			last;
		}	
	}
	return $f;
}

sub isSame {
	my ($a, $b) = @_;
	my $f = 1;
	for my $i (0..(@$a-1)){
		if($a->[$i] ne $b->[$i]){
			$f = 0;
			last;
		}	
	}

	return $f;
}


1;



