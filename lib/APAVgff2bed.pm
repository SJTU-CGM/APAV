#!/usr/bin/perl

package APAVgff2bed;

use strict;
use warnings;
use Getopt::Long;
use APAVutils;

sub gff2bed{

        my $usage = "\n\tUsaege: apav gff2bed --gff <gff_file> [options]

'apav gff2bed' is used to extract the coordinates of gene and gene elements from a GFF format file.

Necessary input description:
  -i, --gff		<file>		Gene annotations in a GFF file.

Options:
  -o, --out		<string>	Output file name.

  --chrl		<file>		Chromosome length file.
  --up_n		<int>		Number of bin(s) in gene upstream region. 
  					(Default:0)
  --up_bin		<int>		The interval width of bin(s) in gene upstream region. 
  					(Default:100)
  --down_n		<int>		Number of bin(s) in gene downstream region. 
  					(Default:0)
  --down_bin		<int>		The interval width of bin(s) in gene downstream region. 
  					(Default:100)
  -h, --help                    	Print usage page.  
 
Warning: '--chrl' is required when setting 'up_n'&'up_bin' or 'down_n'&'down_bin'

\n";

	my $cmdline = $0. " gff2bed";
        $cmdline .= " ".join(" ", @ARGV);

	my $stime = `date +"%Y-%m-%d %H:%M:%S"`;
	chomp($stime);

	my ($gff, $out, $chrl, $help, %chr_len);
	my $up_bin = 100;
	my $up_n = 0;
	my $down_bin = 100;
	my $down_n = 0;
	GetOptions(
		'gff|i=s'	=> \$gff,
		'out|o=s'	=> \$out,
		'chrl=s'	=> \$chrl,
		'up_n=i'	=> \$up_n,
		'up_bin=i'	=> \$up_bin,
		'down_n=i'	=> \$down_n,
		'down_bin=i'	=> \$down_bin,
		'help|h!'	=> \$help
	) or die $!."\n";

	die $usage if !defined($gff);
	die $usage if $help;
	APAVutils::check_file('--gff/-i', $gff);
        $out = APAVutils::check_out($out, $gff, ".bed");
	if(($up_n*$up_bin > 0) || ($down_n*$down_bin > 0)){
		APAVutils::check_file('--chrl', $chrl);
	}

	open(GFF, "<$gff") or die "Could not open file '$gff'\n";
	open(OUT, ">$out") or die "Could not open file '$out'\n";

	print STDOUT "[".$stime."] [gff2bed] Extract the coordinates of gene and gene elements...\n";

        print OUT "##Date: $stime\n";
        print OUT "##CMDline: $cmdline\n";
	if(($up_n*$up_bin > 0) | ($down_n*$down_bin > 0)){
		print OUT "##Args: --gff $gff --out $out --chrl $chrl --up_n $up_n --up_bin $up_bin --down_n $down_n --down_bin $down_bin \n";
	}else{
		print OUT "##Args: --gff $gff --out $out \n";
	}
	
	if($chrl){
		open(CHRL, "<$chrl") or die "Could not open file '$chrl'\n";
	        while(<CHRL>){
			$_ =~ s/\r\n$//;
        	        chomp $_;
                	my @chr_info = split(/\t/, $_);
                	$chr_len{$chr_info[0]} = $chr_info[1];
        	}
	}

	my $cur_gene_id;
	my %cur_tran_map;
	my $cur_tran_count;
	my %cur_exon_info;
	my %cur_lcds_info;
	my $cur_exon_count;
	my $chr;

	my @cur_less_info;
	my @cur_large_info;
	
	while(<GFF>){
		$_ =~ s/\r\n$//;
		if($_ =~ /^[^#\s]+/){
			chomp $_;
			my @info = split(/\t/, $_);
			my $type = $info[2];
			my $start = $info[3];
			my $end = $info[4];
			my $strand = $info[6];
			my $lcds;
			if($type eq "gene"){
				if($info[8] =~ /ID=([^;]+)(;|$)/){
					die "Gene id should not contian ':'\n" if($1 =~ /\:/);
					printRes(\%cur_exon_info, \%cur_lcds_info, $chr, \@cur_less_info, \@cur_large_info, $cur_gene_id) if(defined($cur_gene_id));
					print OUT "$info[0]\t".($start - 1)."\t".$end."\t".$1."\n";
					$chr = $info[0];
					$cur_gene_id = $1;
					%cur_tran_map = ();
					$cur_tran_count = 1;
					%cur_exon_info = ();
					%cur_lcds_info = ();

					if(($up_bin*$up_n  > 0) | ($down_bin*$down_n > 0)){
						die "Can not find length of $chr\n" if(!exists($chr_len{$chr}));
						my ($cur_less_info, $cur_large_info) = getUpDown($cur_gene_id, $chr, $chr_len{$chr}, $strand, $start, $end, $up_bin, $up_n, $down_bin, $down_n);
						@cur_less_info = @$cur_less_info;
						@cur_large_info = @$cur_large_info;

					}else{
						@cur_less_info = ();
						@cur_large_info = ();
					}
				}else{
					die "Can not find gene id";
				}	
			}elsif($type eq "mRNA" | $type eq "transcript"){
				if($info[8] =~ /ID=([^;]+)(;|$)/){
					$cur_tran_map{$1} = "T".$cur_tran_count;
					$cur_tran_count ++;
					$cur_exon_count = 1;
				}else{
					die "Can not find $type id";
				}
			}elsif($type eq "exon"){
				if($info[8] =~ /Parent=([^;]+)(;|$)/){
					if(exists($cur_tran_map{$1})){
						my $cur_tran_id = $cur_tran_map{$1};
						my $cur_exon_id = $cur_tran_id.":exon".$cur_exon_count;
						$cur_exon_count ++;
						if(exists($cur_exon_info{"$start-$end"})){
							push (@{$cur_exon_info{"$start-$end"}}, $cur_exon_id);	
						}else{
							$cur_exon_info{"$start-$end"} = [$cur_exon_id];
						}
					}
				}else{
					die "Can not find exon's parent id";
				}			
			}elsif($type eq "CDS" | $type eq "five_prime_UTR" | $type eq "three_prime_UTR" | $type eq "5UTR" | $type eq "3UTR"){
				if($type eq "CDS"){
					$lcds = "CDS";
				}elsif($type eq "five_prime_UTR" | $type eq "5UTR"){
					$lcds = "5UTR";
				}else{
					$lcds = "3UTR";
				}	
				if($info[8] =~ /Parent=([^;]+)(;|$)/){
					if(exists($cur_tran_map{$1})){
						my $cur_tran_id = $cur_tran_map{$1};
						if(exists($cur_exon_info{"$start-$end"})){
							for(my $i=0; $i<=$#{$cur_exon_info{"$start-$end"}}; $i++){
								if(@{$cur_exon_info{"$start-$end"}}[$i] =~ /^$cur_tran_id:exon/){
									@{$cur_exon_info{"$start-$end"}}[$i] .= "($lcds)";
								}	
							}
						}else{
							for my $key (keys %cur_exon_info){
								if($key =~ /(\d+)-(\d+)/){
									if($start >= $1 && $end <= $2){
										for(my $i=0; $i<=$#{$cur_exon_info{$key}}; $i++){
											if(@{$cur_exon_info{$key}}[$i] =~ /^($cur_tran_id:exon\d+)/){
												my $cur_cds_id = $1.":$lcds";
												@{$cur_exon_info{$key}}[$i] .= "-";
												if(exists($cur_lcds_info{"$start-$end"})){
													push (@{$cur_lcds_info{"$start-$end"}}, $cur_cds_id);
												}else{
													$cur_lcds_info{"$start-$end"} = [$cur_cds_id];
												}
											}
										}
									}
								}
							}
						}
					}else{
						die "Can not find parent";
					}
				}else{
					die "Can not find cds/utr's parent id";
				}
			}
		}
		
	}
	printRes(\%cur_exon_info, \%cur_lcds_info, $chr, \@cur_less_info, \@cur_large_info, $cur_gene_id);

	close(GFF);
	close(CHRL);
	close(OUT);

	my $etime = `date +"%Y-%m-%d %H:%M:%S"`;
        chomp($etime);
        print STDOUT "[".$etime."] [gff2bed] Finished\n";

}

sub printRes{
        my ($exon_info, $lcds_info, $chr, $less_info, $large_info, $gene_id) = @_;
        my @element_info;
        for my $key (keys %$exon_info){
                my @array = grep {$_ =~ /[^-]$/}@{$exon_info->{$key}};
                if(@array){
                        $key =~ /^(\d+)-(\d+)$/;
                        push(@element_info, [$1, $2, join(",",@array)]);
                }
        }
        for my $key (keys %$lcds_info){
                $key =~ /^(\d+)-(\d+)$/;
                push(@element_info, [$1, $2, join(",",@{$lcds_info->{$key}})]);
        }
        @element_info = sort{ $a->[1] <=> $b->[1] }@element_info;
        @element_info = sort{ $a->[0] <=> $b->[0] }@element_info;
        my @less = @$less_info;
        my @large = @$large_info;
        for(my $i=0;$i<=$#less;$i++){
                print OUT $chr."\t".($less[$i][0] - 1)."\t",$less[$i][1]."\t".$less[$i][2]."\n";   ##start
        }
        for(my $i=0;$i<=$#element_info;$i++){
                print OUT $chr."\t".($element_info[$i][0] - 1)."\t",$element_info[$i][1]."\t".$gene_id.":[".$element_info[$i][2]."]"."\n";   ## start
	}
        for(my $i=0;$i<=$#large;$i++){
                print OUT $chr."\t".($large[$i][0] - 1)."\t",$large[$i][1]."\t".$large[$i][2]."\n";   ## start
        }
}

sub getUpDown{
        my ($gene_id, $chr, $chr_len, $strand, $start, $end, $up_bin, $up_n, $down_bin, $down_n) = @_;
        my @res_less;
        my @res_large;
        if($strand eq "+"){     
                for(my $i=$up_n; $i>=1; $i--){
                        if($start - (($i-1)*$up_bin+1) > 1){
                                if($start - $i*$up_bin > 1){
                                        push(@res_less, [$start-$i*$up_bin, $start-(($i-1)*$up_bin+1), $gene_id.":[UP:".(($i-1)*$up_bin+1)."-".$i*$up_bin."]"]);
                                }else{
                                        push(@res_less, [1, $start-(($i-1)*$up_bin+1), $gene_id.":[UP:".(($i-1)*$up_bin+1)."-".($start-1)."]"]);
                                }
                        }
                        
                }
                for(my $i=1; $i<=$down_n; $i++){
                        if($end+(($i-1)*$down_bin+1) < $chr_len){
                                if($end+$i*$down_bin < $chr_len){
                                        push(@res_large, [$end+(($i-1)*$down_bin+1), $end+$i*$down_bin, $gene_id.":[DOWN:".(($i-1)*$down_bin+1)."-".$i*$down_bin."]"]);
                                }else{
                                        push(@res_large, [$end+(($i-1)*$down_bin+1), $chr_len, $gene_id.":[DOWN:".(($i-1)*$down_bin+1)."-".($chr_len-$end)."]"]);            
                                }
                        }
                }
        }elsif($strand eq "-"){
                for(my $i=$down_n; $i>=1; $i--){
                        if($start-(($i-1)*$down_bin+1) > 1){
                                if($start-$i*$down_bin > 1){
                                        push(@res_less, [$start-$i*$down_bin, $start-(($i-1)*$down_bin+1), $gene_id.":[DOWN:".(($i-1)*$down_bin+1)."-".$i*$down_bin."]"]);
                                }else{
                                        push(@res_less, [1, $start-(($i-1)*$down_bin+1), $gene_id.":[DOWN:".(($i-1)*$down_bin+1)."-".($start-1)."]"]);
                                }
                        }
                }
                for(my $i=1; $i<=$up_n; $i++){
                        if($end+(($i-1)*$up_bin+1) < $chr_len){
                                if($end+$i*$up_bin < $chr_len){
                                        push(@res_large, [$end+(($i-1)*$up_bin+1), $end+$i*$up_bin, $gene_id.":[UP:".(($i-1)*$up_bin+1)."-".$i*$up_bin."]"]);
                                }else{
                                        push(@res_large, [$end+(($i-1)*$up_bin+1), $chr_len, $gene_id.":[UP:".(($i-1)*$up_bin+1)."-".($chr_len-$end)."]"]);
                                }
                        }
                }
        }
        return (\@res_less, \@res_large);
}


1;


