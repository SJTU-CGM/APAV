#!/usr/bin/perl

package APAVsim;

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw/sum/;

sub sim{

	my $usage = "\n\tUsage: apav pavSize --pav <pav_file> [options]

'apav pavSize' is used to estimate the size of pan-genome and core-genome from PAV table.

Necessary input description:
  -i, --pav     <file>	        PAV file.

Options:
  -o, --out     <string>        Output file name.
  -n 		<int>		Specifies the number of random sampling times. 
  				(Default:100)
  --group	<file>		Group file for samples.
  -h, --help                    Print usage page.
  \n";

  	 my $cmdline = $0. " pavSize";
        $cmdline .= " ".join(" ", @ARGV);

        my $stime = `date +"%Y-%m-%d %H:%M:%S"`;
        chomp($stime);

  	my ($pav, $out, $group_file, $help);
	my $number = 100;
	GetOptions(
		'pav|i=s'	=> \$pav,
		'out|o=s'	=> \$out,
		'n=i'		=> \$number,
		'group=s'	=> \$group_file,
		'help|h!'	=> \$help
	) or die $!."\n";

	die $usage if !defined($pav);
	die $usage if $help;
	APAVutils::check_file('--pav/-i', $pav);
	APAVutils::check_file('--group', $group_file) if defined($group_file);
	$out = APAVutils::check_out($out, $pav, ".size");

	open(PAV, "<$pav") or die "Could not open file '$pav'\n";
	open(OUT, ">$out") or die "Could not open file '$out'\n";

	print STDOUT "[".$stime."] [pavSize] Estimate the size of pan-genome and core-genome from PAV table....\n";

        print OUT "##Date: $stime\n";
        print OUT "##CMDline: $cmdline\n";
	if(defined($group_file)){
		print OUT "##Args: --pav $pav --group $group_file --out $out -n $number \n";
	}else{
		print OUT "##Args: --pav $pav --out $out -n $number\n";
	}

	my $sample_n;

	if($group_file){
		open(GROUP, "<$group_file") or die "Could not open file '$group_file'\n";
		my @group;
		my @uniqg;
		my %map;
		my %gtable;
		my %gallp;
		while(<GROUP>){
			$_ =~ s/\r\n$//;
			chomp $_;
			my @arr = split(/\t/, $_);
			die "Error: Can not find phenotype column in $group_file\n." if(scalar(@arr) < 2);
			$map{$arr[0]} = $arr[1]; 
		}
		close(GROUP);
		my %tmp;

		while(<PAV>){
			chomp($_);
			if($_ =~ /^#/){
				next;
			}elsif($_ =~ /^Chr\tStart\tEnd\tLength\tAnnotation/){
				my @header = split(/\t/, $_);
				foreach(@header[(5)..($#header)]){
					if($map{$_}){
						push(@group, $map{$_});
					}else{
						die "Error: Can not find phentype of $_\n";
					}
				}
				@uniqg = grep{++$tmp{$_}<2} @group;
				next;
			}
			my @arr = split(/\t/, $_);
			my @pavs = @arr[(5)..($#arr)];

			$sample_n = scalar(@pavs) if !defined($sample_n);
			die "sample number error\n" if scalar(@pavs) ne $sample_n;
			die "sample number != group number\n" if $sample_n ne scalar(@group);
			
			foreach my $g (@uniqg){
				my @idx = grep{$group[$_] eq $g} 0..$#group;
				my @d = map{$pavs[$_]} @idx;
				if(sum(@d) == 0){
					next;
				}elsif(sum(@d) == scalar(@d)){
					$gallp{$g} ++;	
				}else{
					push(@{$gtable{$g}}, [@d]);
				}
			}
		}

		print OUT "Round\tSampleN\tCore\tPan\tDelta\tGroup\n";

		foreach my $g (@uniqg){
			my $cur_sn = grep{$group[$_] eq $g} 0..$#group;
			my @cur_table = @{$gtable{$g}};
			simprint($number, $cur_sn, \@{$gtable{$g}}, $gallp{$g}, $g);
		}
	}else{
		my @table;
		my $allp = 0;
		while(<PAV>){
			if($_ =~ /^[#|Chr\tStart\tEnd\tLength\tAnnotation\t]/){
				next;
			}
			chomp $_;
			my @arr = split(/\t/, $_);
			my @pavs = @arr[(5)..($#arr)];

			$sample_n = scalar(@pavs) if !defined($sample_n);
			die "sample number error\n" if scalar(@pavs) ne $sample_n;
			if (sum(@pavs) == 0) {
				next;
			}elsif(sum(@pavs) == $sample_n){
				$allp ++;
				next;
			}else{
				push(@table, [@pavs]);
			}
		}
		
		if($#table eq -1){
			die "Can not find dispensable region for estimation.\n";
		}else{
			print OUT "Round\tSampleN\tCore\tPan\tDelta\n";
			simprint($number, $sample_n, \@table, $allp);
		}
	}
	close(PAV);
	close(OUT);

	my $etime = `date +"%Y-%m-%d %H:%M:%S"`;
        chomp($etime);
        print STDOUT "[".$etime."] [pavSize] Finished\n";

}

sub simprint{
	my ($number, $sn, $t, $ap, $g) = @_;
	my @t = @$t;

	foreach my $n(1..$number){
		my @rand = frand($sn);
		my @core = map{$t[$_]->[$rand[0]]} 0..$#t;
		my @pan = @core;
		my $c;
		if(!@core){
			 $c = 0;
		}else{
			$c = sum(@core);
		}
		my $p = $c;
		if($g){
			print OUT join("\t", ($n, "1", $c + $ap, $p + $ap, "", $g))."\n";
		}else{
			print OUT join("\t", ($n, "1", $c + $ap, $p + $ap, ""))."\n";
		}
		my $tmp = $p;
	
		foreach my $i (1..$#rand){
			foreach my $j (0..$#t){
				if($core[$j] == 1 && $t[$j]->[$rand[$i]] != 1){
					$core[$j] = 0;
					$c --;
				}
				if($pan[$j] != 1 && $t[$j]->[$rand[$i]] == 1){
					$pan[$j] = 1;
					$p ++;
				}
			}
			if($g){
				print OUT join("\t", ($n, $i+1, $c + $ap, $p + $ap, $p - $tmp, $g))."\n";
			}else{
				print OUT join("\t", ($n, $i+1, $c + $ap, $p + $ap, $p - $tmp))."\n";
			}
			$tmp = $p;
		}
	}
}


sub frand{
	my @r;
    	my %h;
    	foreach (1..$_[0]){
		$h{$_-1}=1;
    	}
    	my @t=keys(%h);
    	while(@t>0){
		my $i=int(rand(scalar(@t)));
		push @r,$t[$i];
		delete($h{$t[$i]});
		@t=keys(%h);
    	}
    	return @r;
}



1;


