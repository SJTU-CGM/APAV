#!/usr/bin/perl

package APAVstaCov;

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use APAVmerge;

sub staCov{
	my $usage = "\n\tUsaege: apav staCov --bed <BED_file> --bamdir <bam_dir> [options]

'apav staCov' is used to calculate coverage of target regions and element regions.
The script will call samtools program, so please make sure samtools is in your PATH.

Necessary input description:
  -i, --bed	<file>		Annotation of target regions in a BED file.
  -b, --bamdir	<dir>		The directory contains mapping results (sorted '.bam' files).

Options:
  -o, --out	<string>	Output file prefix.
  --asgene			Treating target regions as genes.
  --rep		<type>		Representative transcript selection methods:
  				\"cdslen\": using the transcript with the longest CDS region;
				\"exonlen\": using the transcript with the longest exon region;
				\"len\": using the longest transcript;
				\"highcov\": using the transcript with the highest coverage;
				\"none\": using gene body region;
				\"all\": considering all transcript reigon.
				(Default:\"cdslen\")
  --mincov	<int>		Depth threshold.
  				(Default:1)
  --rmele			Remove analysis on elements.
  --merge			Merge neighboring elements with the same coverage(all 0 or all 1).
  --thread	<int>		Thread number.
  -h, --help                    Print usage page.

Warning: '--rep XXX' only works with '--asgene' ;
	 '--merge' only works without '--rmele'
  \n";

  	my $cmdline = $0. " staCov";
  	$cmdline .= " ".join(" ", @ARGV);

	my $stime = `date +"%Y-%m-%d %H:%M:%S"`;
	chomp($stime);

	my ($bed, $out, $bam_dir, $help);
	my $asgene = 0;
	my $rep = "cdslen";
	my $mincov = 1;
	my $rmele = 0;
	my $merge = 0;
	my $thread = 1;

	GetOptions(
		'bed|i=s'	=> \$bed,
		'out|o=s'	=> \$out,
		'bamdir|b=s'	=> \$bam_dir,
		'asgene!'	=> \$asgene,
		'rep=s'		=> \$rep,
		'mincov=i'	=> \$mincov,
		'rmele!'	=> \$rmele,
		'merge!'	=> \$merge,
		'thread=i'	=> \$thread,
		'help|h!'	=> \$help	
	) or die $!."\n";

	die $usage if !defined($bed) & !defined($bam_dir);
	die $usage if $help;
	APAVutils::check_file('--bed/-i', $bed);
	APAVutils::check_file('--bamdir/-b', $bam_dir);
	$out = APAVutils::check_out($out, $bed, "");
	
	die "Please install 'samtools' first\n" if(system("command -v samtools > /dev/null 2>&1") != 0);

	print STDOUT "[".$stime."] [staCov] Calculate coverage of target regions and element regions...\n";

	my $header = "##Date: $stime\\n##CMDline: $cmdline\\n##Args: --bed $bed --bamdir $bam_dir --out $out";
	$header .= " --asgene --rep $rep" if $asgene;
	$header .= " --mincov $mincov --thread $thread";
	$header .= " --merge" if $merge;

	$bam_dir.="/" unless($bam_dir=~/\/$/);
	my @bams = <$bam_dir*.bam>;
	my @samples = @bams;
	map{ $_ =~ s/.*\///g }@samples;
	map{ $_ =~ s/.bam$//g }@samples;

	my $region = get_region($bed, $asgene, $rep);

        if(scalar(keys(%$region)) eq 0){
		die "No valid regions extracted\n";
	}

	my $tmpdir_init = "${out}_tmp";
	my $tmpdir = $tmpdir_init;
	my $counter = 0;
	while (-d $tmpdir){
		$tmpdir = $tmpdir_init.time().$counter;
		$counter ++;
	}
	system("mkdir $tmpdir");
	foreach my $i(0..$#bams){
		print "-- $bams[$i]\n";
		system("samtools depth -@ $thread -aa -b $bed $bams[$i] > $tmpdir/$samples[$i].depth");
		depth2cov($region, $samples[$i], $out, $tmpdir, $asgene, $rep, $mincov, $header, $rmele);
		system("rm $tmpdir/$samples[$i].depth");
	}

	paste_cov($asgene, $rep, \@samples, $out, $tmpdir, $header, $rmele);

	if(!$rmele && $merge){
		if(!$asgene || ($asgene && ($rep eq "cdslen" || $rep eq "exonlen" || $rep eq "all"))){
			APAVmerge::mergecov("${out}_ele.cov", "${out}_ele.mcov", $asgene);
		}
	}

	system("rm -r $tmpdir");	

	my $etime = `date +"%Y-%m-%d %H:%M:%S"`;
        chomp($etime);
        print STDOUT "[".$etime."] [staCov] Finished\n";
}

sub depth2cov {

	my ($region, $sample, $out, $tmpdir, $asgene, $rep, $mincov, $header, $rmele) = @_;
	my %region = %$region;

		open(DEPTH, "<$tmpdir/$sample.depth");
		open(OUT1, ">$tmpdir/$sample.cov");
		if(!$rmele){
			open(OUT2, ">$tmpdir/${sample}_ele.cov");
		}
	
		my @chr_region;
		my @tmp;
		my $cur_chr;
		my $min_pos;
		my %chr_res;
		my @chr_res_list;

		while(<DEPTH>){
			next if($_ =~ /^[#\s]+/);
			$_ =~ s/\r\n$//;
			chomp $_;
			my @arr = split(/\t/, $_);
			my $chr = $arr[0];
			my $pos = $arr[1];
			my $depth = $arr[2];

			if(!defined($cur_chr)){
				$cur_chr = $chr;
				@chr_region = @{$region{$chr}};
				$min_pos = $chr_region[0][0];
			}

			if($chr ne $cur_chr){
				foreach(@tmp){
					if($asgene){
						if($_->[2] =~ /^([^:]*)/){
							push(@chr_res_list, $1) if !(exists $chr_res{$1});
							push(@{$chr_res{$1}}, [$_->[0], $_->[1], $_->[3], $_->[2]]);
						}
					}else{
						push(@chr_res_list, $_->[2]) if !(exists $chr_res{$_->[2]});
						push(@{$chr_res{$_->[2]}}, [$_->[0], $_->[1], $_->[3]]);
					}
		                }

				print_chr_res($cur_chr, \@chr_res_list, \%chr_res, $asgene, $rep, $rmele);

				$cur_chr = $chr;
				@chr_region = @{$region{$chr}};
				$min_pos = $chr_region[0][0];
				@tmp = ();
				@chr_res_list =();
				%chr_res = ();
			}

			if($pos < $min_pos){
				next;
			}

			my @n;
			foreach(@chr_region){
				if($pos >= $_->[0]){
					push(@tmp, [$_->[0], $_->[1], $_->[2], 0]);
					push(@n, 1);
				}else{
					last;
				}
			}
			foreach(@n){
				shift(@chr_region);
			}

			@tmp = sort {$a->[1] <=> $b->[1]} @tmp;

			my @m;
			foreach(@tmp){
				if($pos > $_->[1]){
					if($asgene){
						if($_->[2] =~ /^([^:]*)/){
							push(@chr_res_list, $1) if !(exists $chr_res{$1});
                        		                push(@{$chr_res{$1}}, [$_->[0], $_->[1], $_->[3], $_->[2]]);
						}
					}else{
						push(@chr_res_list, $_->[2]) if !(exists $chr_res{$_->[2]});
						push(@{$chr_res{$_->[2]}}, [$_->[0], $_->[1], $_->[3]]);
					}

					push(@m, 1);
				}else{
					$_->[3]++ if($depth >= $mincov);
				}
			}
			foreach(@m){
				shift(@tmp);
			}

			if(@chr_region){
				$min_pos = $chr_region[0][0];
			}else{
				$min_pos = 9**9**9;
			}

			if(@tmp){
				if($tmp[0][0] < $min_pos){
					$min_pos = $tmp[0][0];
				}
			}
		}

		foreach(@tmp){
			if($asgene){
				if($_->[2] =~ /^([^:]*)/){
					push(@chr_res_list, $1) if !(exists $chr_res{$1});
					push(@{$chr_res{$1}}, [$_->[0], $_->[1],  $_->[3], $_->[2]]);
				}
			}else{
				push(@chr_res_list, $_->[2]) if !(exists $chr_res{$_->[2]});
				push(@{$chr_res{$_->[2]}}, [$_->[0], $_->[1], $_->[3]]);
			}
		}

		print_chr_res($cur_chr, \@chr_res_list, \%chr_res, $asgene, $rep, $rmele);
		@chr_res_list = ();
		%chr_res = ();


		close(DEPTH);
		close(OUT1);
		close(OUT2);
}


sub paste_cov{
	my ($asgene, $rep, $samples, $out, $tmpdir, $header, $rmele) = @_;
	my @samples = @$samples;

	if(!($asgene && $rep eq "all")){
		system("cut -f 1,2,3,4,5,6 $tmpdir/$samples[0].cov > ${out}.cov");
		foreach (1..$#samples){
			system("cut -f 6 $tmpdir/$samples[$_].cov | paste ${out}.cov - > ${out}.cov.tmp ");
			system("mv ${out}.cov.tmp ${out}.cov");
		}
		system("sed -i '1i\\Chr\tStart\tEnd\tLength\tAnnotation\t".join("\t", @samples)."' ${out}.cov");
		system("sed -i '1i\\".$header."' ${out}.cov");
	}
	if($asgene && ($rep eq "highcov")){
		system("cut -f 1,2,3,4,5 $tmpdir/$samples[0].cov > ${out}.cov.pos");
		foreach(1..$#samples){
			system("cut -f 7 $tmpdir/$samples[$_].cov | paste ${out}.cov.pos - > ${out}.cov.pos.tmp ");
			system("mv ${out}.cov.pos.tmp ${out}.cov.pos");
		}
		system("sed -i '1i\\Chr\tStart\tEnd\tLength\tAnnotation\t".join("\t", @samples)."' ${out}.cov.pos");
		system("sed -i '1i\\".$header."' ${out}.cov.pos");
	}
	if(!$rmele && (!$asgene || ($asgene && ($rep eq "cdslen" || $rep eq "exonlen" || $rep eq "all")))){
		system("cut -f 1,2,3,4,5,6 $tmpdir/$samples[0]_ele.cov > ${out}_ele.cov");
		foreach (1..$#samples){
			system("cut -f 6 $tmpdir/$samples[$_]_ele.cov | paste ${out}_ele.cov - > ${out}_ele.cov.tmp ");
			system("mv ${out}_ele.cov.tmp ${out}_ele.cov");
		}
		system("sed -i '1i\\Chr\tStart\tEnd\tLength\tAnnotation\t".join("\t", @samples)."' ${out}_ele.cov");
		system("sed -i '1i\\".$header."' ${out}_ele.cov");
	}


}





sub print_chr_res {
	my ($chr, $chr_res_list, $chr_res, $asgene, $rep, $rmele) = @_;
	my %chr_res = %$chr_res;
	foreach my $target (@$chr_res_list){
			my $flag = 0;
			my %trans;
			my @parts;
                        foreach(@{$chr_res{$target}}){
				my $start = $_->[0];
				my $end = $_->[1];
				my $covered = $_->[2];
				if($asgene){
					my $anno = $_->[3];
					if($rep eq "none" || $rep eq "len"){
						print_gene_cov($chr, $rep, $_, $rmele);
					}else{
						if( $anno =~ /^[^:]*:\[[UP|DOWN].*\]$/){
							if($flag){
								print_gene_cov($chr, $rep, \%trans, $rmele);
								$flag = 0;
							}
							my $cur_len = $end - $start + 1;
							my $cov = covround($covered/$cur_len);
							print OUT2 "$chr\t$start\t$end\t$cur_len\t$anno\t$cov\n" if(!$rmele);
						}elsif($anno =~ /^([^:]*):\[(T.*)\]$/){
							if($rep eq "cdslen" || $rep eq "exonlen" || $rep eq "highcov"){
								my @exons = split(/,/, $2);
								foreach(@exons){
									my @a = split(/:/, $_, 2);
									push(@{$trans{$a[0]}}, [$start, $end, $covered, $a[1], $1]);
								}
							}elsif($rep eq "all"){
								push(@{$trans{"all"}}, [$start, $end, $covered, $anno]);
							}
							$flag = 1;
						}
					}

				}else{
					push(@parts, [$start, $end, $covered]);
				}                             
                        }
			if($flag){
				print_gene_cov($chr, $rep, \%trans, $rmele);
			}
			if(@parts){
				print_general_cov($chr, $target, \@parts, $rmele);
			}              
                }
}

sub print_general_cov {
	my ($chr, $target, $parts, $rmele) = @_;
	my @parts = @$parts;
	my @starts;
	my @ends;
	my $sum_len = 0;
	my $sum_covered = 0;
	my $n = 1;
	foreach(@parts){
		my $l = $_->[1] - $_->[0] + 1;
		my $cov = covround($_->[2]/$l);
		print OUT2 "$chr\t$_->[0]\t$_->[1]\t$l\t$target:[$n]\t$cov\n" if(!$rmele);
		$sum_len += $l;
		$sum_covered += $_->[2];
		$n ++;
		push(@starts, $_->[0]);
		push(@ends, $_->[1]);
	}
	my $whole_cov = covround($sum_covered/$sum_len);
	print OUT1 "$chr\t".join(",", @starts)."\t".join(",", @ends)."\t$sum_len\t$target\t$whole_cov\n";
}

sub print_gene_cov {
	my ($chr, $rep, $t, $rmele) = @_;

	if($rep eq "none" || $rep eq "len"){
		my @info = @$t;
		my $l = $info[1] - $info[0] + 1;
		my $cov = covround($info[2]/$l);
		print OUT1 "$chr\t$info[0]\t$info[1]\t$l\t$info[3]\t$cov\n";
	}elsif($rep eq "all"){
		my %trans = %$t;
		foreach(@{$trans{"all"}}){
			my $l = $_->[1] - $_->[0] + 1;
			my $cov = covround($_->[2]/$l);
			print OUT2 "$chr\t$_->[0]\t$_->[1]\t$l\t$_->[3]\t$cov\n" if(!$rmele);	
		}

	}else{
		my %trans = %$t;
		my $len = 0;
		my $rep_trans;

		my $min_pos;
		my $max_pos;

		if($rep eq "cdslen" || $rep eq "exonlen" || $rep eq "len"){
			my @keys = keys %trans;
			$rep_trans = $keys[0];
		}elsif($rep eq "highcov"){
			foreach my $key (keys %trans){
				my $sum = 0;
                		foreach(@{$trans{$key}}){
					if(!defined($min_pos)){
						$min_pos = $_->[0];
					}else{
						$min_pos = $_->[0] if($_->[0] < $min_pos);
					}
					if(!defined($max_pos)){
						$max_pos = $_->[1];
					}else{
						$max_pos = $_->[1] if($_->[1] > $max_pos);
					}
                        		$sum += $_->[2];
                		}
                		if($sum >= $len){
                        		$len = $sum;
                        		$rep_trans = $key;
                		}
			}
		}

		my @starts;
		my @ends;
		my $sum_len = 0;
		my $sum_covered = 0;

		foreach(@{$trans{$rep_trans}}){
			my $l = $_->[1] - $_->[0] + 1;
			my $cov = covround($_->[2]/$l);
			if( $rep eq "cdslen" || $rep eq "exonlen" || $rep eq "len"){
               			print OUT2 "$chr\t$_->[0]\t$_->[1]\t$l\t$_->[4]:[${rep_trans}:$_->[3]]\t$cov\n" if(!$rmele);
			}
			push(@starts, $_->[0]);
			push(@ends, $_->[1]);
			$sum_len += $l;
			$sum_covered += $_->[2];
        	}
		my $whole_cov = covround($sum_covered/$sum_len);
		my $gene = $trans{$rep_trans}[0]->[4];
		if($rep eq "cdslen" || $rep eq "exonlen" || $rep eq "len"){
			print OUT1 "$chr\t".join(",", @starts)."\t".join(",", @ends)."\t$sum_len\t".$gene."\t$whole_cov\n";
		}elsif($rep eq "highcov"){
			print OUT1 "$chr\t$min_pos\t$max_pos\t-\t$gene\t$whole_cov\t".$rep_trans."($chr:(".join(",", @starts).")-(".join(",", @ends)."))\n";
		}
	}
}


sub get_region {

	my ($bed, $asgene, $rep) = @_;

	open(BED, "<$bed") or die "Could not open file '$bed'.\n";

	my %region;
	my $cur_gene;
	my %cur_trans;
	my $flag;
	my $cur_chr;

	while(<BED>){
        	next if($_ =~ /^[#\s]+/);
        	$_ =~ s/\r\n$//;
        	chomp $_;
        	my @arr = split(/\t/, $_);
       		my $chr = $arr[0];
        	my $start = $arr[1] + 1;
        	my $end = $arr[2];
        	my $anno = $arr[3];

		$cur_chr = $chr if(!defined($cur_chr));

		if($asgene){
			if($rep eq "cdslen" || $rep eq "exonlen" || $rep eq "len"){
				if($chr ne $cur_chr){
					if(%cur_trans){
						map { push(@{$region{$cur_chr}}, $_) } @{get_longest($rep, \%cur_trans)};
						%cur_trans = ();
					}
					$cur_chr = $chr;
				}
				if($anno =~ /.*:\[UP.*\]$/){
					$flag = "up" if !defined($flag);
					if($flag ne "up"){
						if(%cur_trans){
							map { push(@{$region{$cur_chr}}, $_) } @{get_longest($rep, \%cur_trans)};
							%cur_trans = ();
						}
					}
					push(@{$region{$chr}}, [$start, $end, $anno]);
					$flag = "up";
				}elsif($anno =~ /.*:\[DOWN.*\]$/){
					$flag = "down" if !defined($flag);
					if($flag ne "down"){
						if(%cur_trans){
							map { push(@{$region{$cur_chr}}, $_) } @{get_longest($rep, \%cur_trans)};
							%cur_trans = ();
						}
					}
					push(@{$region{$chr}}, [$start, $end, $anno]);
					$flag = "down";
				}elsif($anno =~ /(.*):\[(T.*)\]$/){
					$flag = "gene";
					$cur_gene = $1 if !defined($cur_gene);
					if($1 ne $cur_gene){
						if(%cur_trans){
							map { push(@{$region{$cur_chr}}, $_) } @{get_longest($rep, \%cur_trans)};
							%cur_trans = ();
                                        	}
						$cur_gene = $1;
					}
					my @exons = split(/,/, $2);
					foreach(@exons){
						my @a = split(/:/, $_, 2);
						if($rep eq "cdslen"){
							if($a[1] =~ /CDS/){
								push(@{$cur_trans{$a[0]}}, [$start, $end, ($end - $start + 1), $cur_gene, $a[1]]);
							}
						}else{
							push(@{$cur_trans{$a[0]}}, [$start, $end, ($end - $start + 1), $cur_gene, $a[1]]);
						}
					}				
							
				}
		
			}elsif($rep eq "highcov" || $rep eq "all"){
				if($anno =~ /.*:\[.*\]$/){
					push(@{$region{$chr}}, [$start, $end, $anno]);
				}else{
					next;
				}	
			}elsif($rep eq "none"){
				if($anno =~ /.*:\[.*\]$/){
					next;
				}else{
					push(@{$region{$chr}}, [$start, $end, $anno]);
				}	
			}else{
				die "\n";
			}	
		}else{
			push(@{$region{$chr}}, [$start, $end, $anno]);	
		}	
	}

	if(%cur_trans){
		map { push(@{$region{$cur_chr}}, $_) } @{get_longest($rep, \%cur_trans)};
		%cur_trans = ();
	}

	foreach my $key (keys %region){
        	my @sorted = sort {$a->[0] <=> $b->[0]} @{$region{$key}};
        	$region{$key} = \@sorted;
	}

	close(BED);

	return(\%region);
}

sub get_longest{
	my ($rep, $t) = @_;
	my %trans = %$t;
	my $len = 0;
	my $longest;
	my @res;

	if($rep eq "len"){
		my $start;
		my $end;
		foreach my $key (sort keys %trans){
			my $sum = 0;
			my $min;
			my $max;
			foreach(@{$trans{$key}}){
				if(!defined($min)){
                                        $min = $_->[0];
                                }else{  
                                        $min = $_->[0] if($_->[0] < $min);
                                }
                                if(!defined($max)){
                                        $max = $_->[1];
                                }else{
                                        $max = $_->[1] if($_->[1] > $max);
                                }	
			}
			$sum = $max - $min + 1;
			if($sum > $len){
				$len = $sum;
				$start = $min;
				$end = $max;
				$longest = $key;
			}
		}
		if(defined($longest)){
			push(@res, [$start, $end, $trans{$longest}[0][3]]);
		}
	}else{
		foreach my $key (sort keys %trans){
			my $sum = 0;
			if($rep eq "cdslen"){
				foreach(@{$trans{$key}}){
					if($_->[4] =~ /CDS/){
						$sum += $_->[2];
					}
				}
			}elsif($rep eq "exonlen"){
				foreach(@{$trans{$key}}){
					$sum += $_->[2];
				}
			}
			if($sum > $len){
				$len = $sum;
				$longest = $key;
			}
		}
		if(defined($longest)){
			foreach(@{$trans{$longest}}){
				push(@res, [$_->[0], $_->[1], "$_->[3]:[$longest:$_->[4]]"]);
			}
		}
	}
	return(\@res);
} 

sub covround {
        my ($d) = @_;
        $d = sprintf("%.4g", $d);
        if($d > 0 && $d < 0.1){
                $d = sprintf("%.4f", $d);
        }
        return ($d);
}


1;

