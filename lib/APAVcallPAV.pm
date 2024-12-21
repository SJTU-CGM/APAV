#!/usr/bin/perl

package APAVcallPAV;

use strict;
use Getopt::Long;
use Data::Dumper;
use List::Util qw/sum/;
use FindBin qw($Bin);
use Cwd;
use APAVreport;

sub callPAV{

my $usage = "\n\tUsaege: apav callPAV --cov <cov_file> --pheno <phenotype_file> [options]

'apav callPAV' is used to determine PAV based on coverage.

Necessary input description:
  -i, --cov     	<file>		Coverage file.
  --pheno	        <file>		Phenotype file.

Options:
  -o, --out             <string>        Output file prefix.

  --method		<type>		Method of determination: \"fixed\" or \"adaptive\". 
  					(Default:\"fixed\")
  
  fixed method:
  --thre		<float>		Coverage threshold. 
  					(Default:0.5)
  
  adaptive method:
  --mina		<float>		Min absence.  
  					(Default:0.1)
  --iter	        <int>		Maximum number of iterations.  
  					(Default:100)
  
Options for the tracks in genome browser:
  --fa                  <file>		Fasta file of reference genome.
  --gff                 <file>		GFF file is used to add annotation tracks in the genome browser.
  --bamdir          	<dir>		The directory contains mapping results (sorted '.bam' files and index files).
  --slice				Extract the bam file of the target regions, otherwise make symbolic links for the raw bam file.

  -h, --help                    	Print usage page.
  \n";

  	my $cmdline = $0. " callPAV";
        $cmdline .= " ".join(" ", @ARGV);

        my $stime = `date +"%Y-%m-%d %H:%M:%S"`;
        chomp($stime);

        my ($cov, $out_prefix, $help, $phen_file, $fa_file, $gff_file, $bam_dir, $slice);
	my $method = "fixed";
	my $thre = 0.5;
        my $min_absence = 0.1;
	my $max_iter = 100;

        GetOptions(
                'cov|i=s'       => \$cov,
                'out|o=s'       => \$out_prefix,
                'method=s'	=> \$method,
		'thre=f'	=> \$thre,
		'mina=i'	=> \$min_absence,
		'iter=i'	=> \$max_iter,
		'pheno=s'	=> \$phen_file,
		'fa=s'		=> \$fa_file,
		'gff=s'		=> \$gff_file,
		'bamdir=s'	=> \$bam_dir,
		'slice!'	=> \$slice,
                'help|h!'       => \$help
        ) or die $!."\n";

        die $usage if !defined($cov);
        die $usage if $help;
	APAVutils::check_file('--cov/-i', $cov);
	APAVutils::check_file('--pheno/-p', $phen_file);
	$out_prefix = APAVutils::check_out($out_prefix, $cov, "");
	die "Parameter '--method/-m' should be 'fixed' or 'adaptive'\n" if ($method ne "fixed" && $method ne "adaptive");
	APAVutils::check_file('--fa', $fa_file) if defined($fa_file);
	APAVutils::check_file('--gff', $gff_file) if defined($gff_file);
	APAVutils::check_file('--bamdir', $bam_dir) if defined($bam_dir);

	my $out_dir = $out_prefix."_report/";

	mkdir ${out_dir};

	open(COV, "<$cov") or die "Could not open file '$cov'\n";
	open(PAV, ">$out_prefix"."_all.pav") or die "Could not open file '${out_prefix}_all.pav'";
	open(DIST, ">$out_prefix"."_dispensable.pav") or die "Could not open file '${out_prefix}_dispensable.pav'";
	open(DATA, ">$out_dir"."data.json") or die "Could not open file '${out_dir}data.json'";
	open(GREPORT, ">$out_dir"."PAV.html") or die "Could not open file '${out_dir}PAV.html'";
	open(GBED, ">$out_dir"."target.bed") or die "Could not open file '${out_dir}target.bed'";
	open(SREPORT, ">$out_dir"."sample.html") or die "Could not open file '${out_dir}sample.html'";

	print STDOUT "[".$stime."] [callPAV] Determine PAV based on coverage...\n";
        print PAV "##Date: $stime\n";
	print DIST "##Date: $stime\n";
        print PAV "##CMDline: $cmdline\n";
	print DIST "##CMDline: $cmdline\n";
	if($method eq "fixed"){
		print PAV "##Args: --cov $cov --pheno $phen_file --out $out_prefix --method $method -t $thre \n";
		print DIST "##Args: --cov $cov --pheno $phen_file --out $out_prefix --method $method -t $thre \n";
	}else{
		print PAV "##Args: --cov $cov --pheno $phen_file --out $out_prefix --method $method -a $min_absence -n $max_iter \n";
		print DIST "##Args: --cov $cov --pheno $phen_file --out $out_prefix --method $method -a $min_absence -n $max_iter \n";
	}

	my @dtable;
	my @samples;

	print STDOUT "-- Determine PAV\n";
	while(<COV>){
		if($_ =~ /^#/){
			next;
		}
		if($_ =~ /^Chr\tStart\tEnd\tLength\tAnnotation\t/){
			chomp $_;
			my @arr = split(/\t/, $_);
			@samples = @arr[5..$#arr];
			print PAV $_."\n";
			print DIST $_."\n";
			next;
		}
		
		chomp $_;
		my @arr = split(/\t/, $_);
		my $info = join("\t", @arr[0..4]);
		my @covs = @arr[5..$#arr];
		my $pavs = defpav($method, $thre, $min_absence, $max_iter, \@covs);
		foreach (@$pavs){
			if($_ < 1){
				print DIST $info."\t".join("\t", @$pavs)."\n";
				my @mcp = map{@$pavs[$_]."/".$covs[$_]} 0..$#covs;
				my $absence_p = (scalar(@$pavs) - sum(@$pavs)) / scalar(@$pavs);
				$absence_p = sprintf "%.3f",$absence_p;
				$absence_p *= 100;
				push(@dtable, [@arr[4,0..3], $absence_p, @mcp]);
				last;
			}
		}
		print PAV $info."\t".join("\t", @$pavs)."\n";
	}

	system "cp -r ${Bin}/src/js ${out_dir}";
        system "cp -r ${Bin}/src/css ${out_dir}";

	print STDOUT "-- Prepare data for web report\n";
	printGeneReport(\@dtable, \@samples, $out_dir, $fa_file, $gff_file, $bam_dir, $slice);
	print STDOUT "-- Generate web report for PAVs of dispensable regions\n";

	my %p;
	my @phen_name;
	my @phen_flag;
	if($phen_file){
		open(PHEN, "<$phen_file") or die "Could not open file '$phen_file'\n";
		while(<PHEN>){
			$_ =~ s/\r\n$//;
			chomp $_;
			if(!@phen_name){
				@phen_name = split(/\t/, $_);
				shift @phen_name;
				push(@phen_flag, 1) foreach @phen_name;
				next;
			}
			my @arr = split(/\t/, $_);
			foreach (1..$#arr){
				$phen_flag[$_-1] = 0 if(!isNum($arr[$_]));
			}
			$p{$arr[0]} = [@arr];				
		}
	}
	if(%p){
		my @phen_name_str;
		foreach (0..$#phen_name){
			push(@phen_name_str, $phen_name[$_]) if $phen_flag[$_] == 0;
		}
		my @stable;
	
		foreach my $i (0..$#samples){
			my @curphen = @{$p{$samples[$i]}};
			my @curpav = map{my @n = split("/", $_->[$i+6]); $n[0]} @dtable;
			push(@stable, [@curphen, @curpav]);
		}

		my @target_name;
		my @target_info;
		foreach my $i(0..$#dtable){
			push(@target_name, $dtable[$i][0]);
			push(@target_info, [$dtable[$i][0], $dtable[$i][1], $dtable[$i][2], $dtable[$i][3]]);
		}
	
		printSampleReport(\@phen_name, \@phen_name_str, \@stable, \@target_name, \@target_info);
		print STDOUT "-- Generate web report for samples\n";
	}else{
		print STDOUT "Warn: No valid phenotype information is available\n";
	}

	close(COV);
	close(PAV);
	close(DIST);
	close(DATA);
	close(GREPORT);
	close(GBED);
	close(SREPORT);

	system "rm ${out_dir}target.bed" if -e "${out_dir}target.bed";

	my $etime = `date +"%Y-%m-%d %H:%M:%S"`;
        chomp($etime);
        print STDOUT "[".$etime."] [callPAV] Finished\n";
	
}

sub printGeneReport {
	my ($dtable, $samples, $out_dir, $fa_file, $gff_file, $bam_dir, $slice) = @_;
	my @dtable = @$dtable;
	my @lines = map{'"'.join('","', @$_).'"'} @dtable;

	if ($#dtable eq -1){
		system("rm -rf $out_dir");
		die "Can not find dispensable region.\n";
	}
	
	foreach my $i (0..$#dtable){
		my @starts = split(/,/, $dtable[$i][2]);
		my @ends = split(/,/, $dtable[$i][3]);
		my $n = 1;
		if($#starts eq 0){
			print GBED $dtable[$i][1]."\t".($starts[$_]-1)."\t".$ends[$_]."\t".$dtable[$i][0]."\n";
		}else{
			foreach (0..$#starts){
				print GBED $dtable[$i][1]."\t".($starts[$_]-1)."\t".$ends[$_]."\t".$dtable[$i][0]."[".$n."]"."\n";
				$n++;
			}
		}
	}


	my $fa_flag = 0;
	if(defined($fa_file)){
		mkdir "${out_dir}browser";

		system "cp ${fa_file} ${out_dir}browser/reference.fa.gz";
		die "Please install 'samtools' first\n" if(system("command -v samtools > /dev/null 2>&1") != 0);
		print STDOUT "-- Index fasta file\n";
        	system "samtools faidx ${out_dir}browser/reference.fa.gz";
		system "cut -f 1,2 ${out_dir}browser/reference.fa.gz.fai > ref.chrl";
	
		print STDOUT "-- Add annotation track for target regions\n";
		system "sort -k1,1 -k2,2n ${out_dir}target.bed > ${out_dir}target.sorted.bed";
		system "${Bin}/tools/bedToBigBed ${out_dir}target.sorted.bed ref.chrl ${out_dir}browser/target.bb";
		system "rm ${out_dir}target.sorted.bed" if -e "${out_dir}target.sorted.bed"; 
		system "rm ref.chrl";
		$fa_flag = 1;
	}

	my $gff_flag = 0;
	if(defined($gff_file) && defined($fa_file)){
		die "Can not find GFF file\n" if !-e $gff_file;
		$gff_flag = 1;
		print STDOUT "-- Add reference annotation track\n";
		system "cp ${gff_file} ${out_dir}browser/in.gff";
		system '(grep "^#" '.${out_dir}.'browser/in.gff; grep -v "^#" '.${out_dir}.'browser/in.gff | sort -t"`printf \'\t\'`" -k1,1 -k4,4n)  > '.${out_dir}.'browser/reference.gff;';
		system "rm ${out_dir}browser/in.gff";
		system "${Bin}/tools/bgzip ${out_dir}browser/reference.gff";
		system "${Bin}/tools/tabix ${out_dir}browser/reference.gff.gz";
	}

	my $bam_tracks = '';
	if(defined($bam_dir) && defined($fa_file)){
		$bam_dir.="/" unless($bam_dir=~/\/$/);
		my @bams = <$bam_dir*.bam>;
		die "Can not find bam files in the given directory\n" if $#bams == -1;

		my @samples = @bams;
		map{ $_ =~ s/.*\///g }@samples;
		map{ $_ =~ s/.bam$//g }@samples;

		print STDOUT "-- Add alignment tracks\n";
		foreach(@samples){
			if($slice){
				system "samtools view -hb -L ${out_dir}target.bed ${bam_dir}${_}.bam > ${out_dir}browser/${_}.bam";
				system "samtools index ${out_dir}browser/${_}.bam";
			}else{
				my $bam_real = Cwd::realpath("${bam_dir}${_}.bam");
				system "ln -s $bam_real ${out_dir}browser/${_}.bam";
				system "ln -s $bam_real.bai ${out_dir}browser/${_}.bam.bai";
			}
		}

		my @bams_name = @bams;
        	map{ $_ =~ s/.*\///g }@bams_name;
        	map{ $_ =~ s/.bam$//g }@bams_name;
		$bam_tracks = 'let bams = '.'["'.join('","', @bams_name).'"]'.';
                        bams.forEach(function(i){
                                tracks.push(getBamTrack(i, "./browser/"+i+".bam", "./browser/"+i+".bam.bai"));
                        })'
	}

	printf DATA '{"pav":[['.join('],[', @lines).']],';

	my $greport = APAVreport::getPavReport(
		$fa_flag,
		$gff_flag,
		'<th>'.join('</th><th>', @$samples).'</th>', 
		'["'.join('","', @$samples).'"]', 
		$bam_tracks);
	printf GREPORT $greport;
}

sub isNum{
	my ($value) = @_;
	return ($value eq $value+0) ? 1 : 0;
}


sub printSampleReport{

	my ($phen_name, $phen_name_str, $stable, $target_name, $target_info) = @_;
	my @stable = @$stable;
	my @phen_name = @$phen_name;
	my @phen_name_str = @$phen_name_str;

	my @input_phen = map{'<label><input type="checkbox" name="'.$_.'" checked>'.$_.'</label>
		'} @phen_name;
	my @phen_div = map{'<div id="phen-'.$_.'" style="float:left;width:30%;height:300px;"></div>'} @phen_name;
	my @stable_arr = map{'["'.join('","', @{$_}).'"]'} @stable;
	my @target_arr = map{'["","'.join('","', @{$_}).'"]'} @$target_info;

	my @phen_init = map{'let phen'.$_.' = echarts.init(document.getElementById("phen-'.$_.'"));'} @phen_name_str;

	my @phen_draw;
	foreach my $i (0..$#phen_name){
		if(grep {$_ eq $phen_name[$i]} @phen_name_str){
			push(@phen_draw, 'phen'.$phen_name[$i].'.setOption(phenChartOption(getPhenData(sampleTable,'.($i+1).')))');
		}
	}

	printf DATA '"gene":['.join(',', @target_arr).'],'.'"sample":['.join(',', @stable_arr).']}';

	my $sreport = APAVreport::getSampleReport(
		join('<br>', @input_phen), 
		'<th>'.join("</th>\n            <th>", @$target_name).'</th>', 
		join("\n", @phen_div), 
		'["'.join('", "', @phen_name).'"]',  
		join("\n              ", @phen_init), 
		join("\n                              ", @phen_draw));
	print SREPORT $sreport;
}


sub defpav {
	my ($method, $thre, $min_absence, $max_iter, $covs) = @_;
	my @covs = @$covs;
	my @pavs;
	if($method eq "fixed"){
		@pavs = map {($_ >= $thre) ? 1 : 0} @covs;
                        }elsif($method eq "adaptive"){
                                my %tmp;
                                my @uniq = grep {++$tmp{$_} == 1} @covs;
                                my @uniq_s = sort {$a <=>$b} @uniq;

                                my $pav_res;
                                if($#uniq < 2){
                                        @pavs = map {($_ < (1 - $min_absence)) ? 0 : 1} @covs;
                                }else{
                                        my $cov_min = $uniq_s[0];
                                        my $cov_max = $uniq_s[-1];
                                        if($cov_min >= (1 - $min_absence)){
                                                @pavs = map {1} @covs;
                                        }else{
                                                $pav_res = topav($cov_min, $cov_max, \@covs, $max_iter);
                                                @pavs = @$pav_res;
                                        }
                                }
                        }
	return \@pavs;
}


sub topav {
        my ($init_cen1, $init_cen2, $arr, $max_iter) = @_;
        my @idx1;
        my @idx2;
        my $new_cen1;
        my $new_cen2;
        my @res;
        my $n = 1;
        while ($n <= $max_iter){
                for my $i (0..(@$arr-1)){
                        if(abs(@$arr[$i] - $init_cen1) <= abs(@$arr[$i] - $init_cen2)){
                                push(@idx1, $i);
                        }else{
                                push(@idx2, $i);
                        }
                }
                @res = map{(abs($_ - $init_cen1) <= abs($_ - $init_cen2)) ? 0 : 1}@$arr;
                $new_cen1 = sum(@$arr[@idx1]) / scalar(@idx1);
                $new_cen2 = sum(@$arr[@idx2]) / scalar(@idx2);
                if(($new_cen1 == $init_cen1) && ($new_cen2 == $init_cen2)){
                        last;
                }else{
                        $n++;
                        $init_cen1 = $new_cen1;
                        $init_cen2 = $new_cen2;
                        @idx1 = ();
                        @idx2 = ();
                }
        }
        return \@res;
}




1;
