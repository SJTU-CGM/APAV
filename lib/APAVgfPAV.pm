#!/usr/bin/perl

package APAVgfPAV;

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

sub gfpav{

	my $usage = "\n\tUsaege: apav gfamPAV --pav <pav_file> --fam <gene_family_info> [option]

'apav gfamPAV' is used to determine gene family presence-absence based on gene PAV table.

Necessary input description:
  -i, --pav     	<file>	        Gene PAV file.
  --fam			<file>		Gene family table.

Options:
  -o, --out             <string>        Output file name.
  -h, --help                    	Print usage page.
  \n";

  	my $cmdline = $0. " gfamPAV";
        $cmdline .= " ".join(" ", @ARGV);

        my $stime = `date +"%Y-%m-%d %H:%M:%S"`;
        chomp($stime);

        my ($pav, $fam, $out, $help);
        GetOptions(
                'pav|i=s' 	=> \$pav,
		'fam=s'		=> \$fam,
                'out|o=s'       => \$out,
                'help|h!'       => \$help
        ) or die $!."\n";

        die $usage if !defined($pav) & !defined($fam);
        die $usage if $help;
	APAVutils::check_pav_input($pav);
	APAVutils::check_file('--fam/-f', $fam);
	$out = APAVutils::check_out($out, $pav, ".gfpav");
	
	open(PAV, "<$pav") or die "Could not open file '$pav'\n";
	open(FAM, "<$fam") or die "Could not open file '$fam'\n";
	open(OUT, ">$out") or die "Could not open file '$out'\n";
	
	print STDOUT "[".$stime."] [gfamPAV] determine gene family PAV...\n";
	
	print OUT "##Date: $stime\n";
        print OUT "##CMDline: $cmdline\n";
	print OUT "##Args: --pav $pav --fam $fam --out $out \n";

	my %g;
	my @covs;
	my $gid;
	my @arr;
	while(<PAV>){
		next if($_ =~ /^[#\s]+/);
		if($_ =~ /^Chr\tStart\tEnd\tLength\tAnnotation/){
                        print OUT $_;
                        next;
                }
		chomp $_;
		@arr = split(/\t/, $_);
		$gid = $arr[4];
		@covs = @arr[5..$#arr];
		$g{$gid} = [@covs];
	}

	my $f;
	my @res;
	while(<FAM>){
		next if($_ =~ /^[#\s]+/);
		$_ =~ s/\r\n$//;
	        chomp $_;	
		my @arr = split(/\t/, $_);
		my $gf = $arr[0];
		my @genes = split(/,/, $arr[1]);
		@res= ();
		for my $i (0..$#covs){
			$f = 0;
			foreach my $m (@genes){
				if($g{$m}){
					$f = 1 if ($g{$m}[$i] == 1);	
				}else{
					warn "Warn: can not find $m in PAV data.\n";
				}
			}
			push(@res, $f);
		}
		print OUT "-\t-\t-\t-\t".$gf."\t".join("\t", @res)."\n";
		
	}

	close(PAV);
	close(FAM);
	close(OUT);

	my $etime = `date +"%Y-%m-%d %H:%M:%S"`;
        chomp($etime);
        print STDOUT "[".$etime."] [gfamPAV] Finished\n";

}

1;
