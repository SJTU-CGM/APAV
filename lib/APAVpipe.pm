#!/usr/bin/perl

package APAVpipe;

use strict;
use Getopt::Long;
use Data::Dumper;

sub geneBatch {
	my $usage = "\n\tapav geneBatch --gff <GFF_file> --bamdir <BAM_dir> --pheno <phenotype_file> --out <res_prefix> [options] 

Necessary input description:
  -i, --gff		<file>		GFF file.
  -o, --out 		<string>	Output file prefix.
  -b, --bamdir		<dir>		The directory contains mapping results (sorted '.bam' files).
  --pheno		<file>		Phenotype file.

Options for elements in the upstream and downstream region:
  --chrl		<file>		Chromosome length file.
  --up_n		<int>		Number of bin(s) in gene upstream region.
  					(Default:0)
  --up_bin		<int>		The interval width of bin(s) in gene upstream region.
  					(Default:100)
  --down_n		<int>		Number of bin(s) in gene downstream region.
  					(Default:0)
  --down_bin		<int>		The interval width of bin(s) in gene downstream region.
  					(Default:100)

Options for coverage calculation:
  --rep			<type>		Representative transcript selection methods:
  					\"cdslen\": using the transcript with the longest CDS region;
					\"len\": using the longest transcript;
					\"highcov\": using the transcript with the highest coverage;
					\"none\": using gene body region;
					\"all\": considering all transcript reigon.
                                	(Default:\"cdslen\")
  --mincov		<int>		Depth threshold.
                                	(Default:1)
  --rmele                               Remove analysis on elements.
  --merge				Merge neighboring elements with the same coverage(all 0 or all 1).
  --thread		<int>		Thread number.

Options for PAV Determination:
  --method		<type>		Method of determination: \"fixed\" or \"adaptive\". 
                                        (Default:\"fixed\")
  
  fixed method:
  --thre                <float>         Coverage threshold. 
                                        (Default:0.5)
  
  adaptive method:
  --mina                <float>         Min absence.  
                                        (Default:0.1)
  --iter                <int>           Maximum number of iterations.  
                                        (Default:100)


Options for the tracks in genome browser in PAV report:
  --fa                  <file>          Fasta file of reference genome.
  --cpbam				Copy bam files instead of make symbolic links.


Options for the phenotype association analysis: 
  --p_threshold		<numeric>	The filter threshold for p_value/p_adjusted.
                                        (Default: 0.01)
  --adjust_p				Adjust p_value.


Other Options:
  --fam			<file>		Gene family table.
  --focused_pheno	<string>	The phenotype in focus.
  --example_n		<int>		The number of examples used to show figures drawn by 'pavPlotPhenoVio/pavPlotPhenoBar' commands and 'elePlotCov/elePlotPAV/elePlotDepth' commands.
  					(Default:10)

  -h, --help                            Print usage page.


	\n";

	my ($gff, $bamdir, $out, $pheno);
	my ($chrl, $up_n, $up_bin, $down_n, $down_bin);
	my ($rep, $mincov, $mergecov, $rmele, $thread);
	my ($fa, $gff, $cpbam, $fam);
	my ($method, $thre, $mina, $iter);
	my ($focus, $p_thre, $adjustp);
	my $topn = 10;
	my $help;

	GetOptions(
		'gff=s'		=> \$gff,
		'bamdir|b=s'	=> \$bamdir,
		'out=s'		=> \$out,
		'pheno=s'	=> \$pheno,

		'chrl=s'	=> \$chrl,
		'up_n=i'	=> \$up_n,
		'up_bin=i'	=> \$up_bin,
		'down_n=i'	=> \$down_n,
		'down_bin=i'	=> \$down_bin,

		'rep=s'         => \$rep,
                'mincov=i'      => \$mincov,
		'merge!'	=> \$mergecov,
		'thread=i'	=> \$thread,

		'fa=s'		=> \$fa,
		'gff=s'		=> \$gff,
		'cpbam!'	=> \$cpbam,
		'fam=s'		=> \$fam,
		
		'method=s'	=> \$method,
		'thre=f'	=> \$thre,
		'mina=f'	=> \$mina,
		'iter=i'	=> \$iter,

		'focused_pheno=s'	=> \$focus,
		'rmele!'	=> \$rmele,
		'p_threshold=f'	=> \$p_thre,
		'adjust_p!'	=> \$adjustp,
		'example_n=i'	=> \$topn,

		'help|h!'	=> \$help

	) or die $!."\n";

	die $usage if !defined($gff) & !defined($bamdir) & !defined($out) & !defined($pheno);
        die $usage if $help;

	APAVutils::check_file('--gff/-i', $gff);
	APAVutils::check_file('--pheno', $pheno);
        APAVutils::check_bamdir('--bamdir/-b', $bamdir);
	if((defined($up_n) || defined($up_bin) || defined($down_n) || defined($down_bin))){
                APAVutils::check_file('--chrl', $chrl);
        }

	die "Please install 'samtools' first\n" if(system("command -v samtools > /dev/null 2>&1") != 0);

	die "Error: output directory \"$out\" already exists. To avoid overwriting of existing files, we kindly request that the output directory should not exist.\n" if -e $out;

	system("mkdir $out");

	print STDOUT "\n";

	my @arg_gff = "$0 gff2bed --gff $gff --out $out/$out.bed";
	push @arg_gff, "--chrl $chrl" if defined($chrl);
	push @arg_gff, "--up_n $up_n" if defined($up_n);
	push @arg_gff, "--up_bin $up_bin" if defined($up_bin);
	push @arg_gff, "--down_n $down_n" if defined($down_n);
	push @arg_gff, "--down_bin $down_bin" if defined($down_bin);
	system(join(" ", @arg_gff));

	print STDOUT "\n";

        my @arg_cov = ("$0 staCov --bed $out/$out.bed --bamdir $bamdir --out $out/$out");
        push @arg_cov, "--asgene";
        push @arg_cov, "--rep $rep" if defined($rep);
        push @arg_cov, "--mincov $mincov" if defined($mincov);
	push @arg_cov, "--rmele" if defined($rmele);
	push @arg_cov, "--merge" if defined($mergecov);
	push @arg_cov, "--thread $thread" if defined($thread);	
        system(join(" ", @arg_cov));

	print STDOUT "\n";

        calling($0, '.cov', '', $bamdir, $pheno, $out, $fa, $gff, $method, $thre, $mina, $iter, $cpbam);

	print STDOUT "\n";

	if(!$rmele && $rep ne "len" && $rep ne "highcov" && $rep ne "none"){
		calling($0, '_ele.cov', '_ele', undef, $pheno, $out, undef, undef, $method, $thre, $mina, $iter, $cpbam);
		print STDOUT "\n";
	}

	if(defined($fam)){
		system("$0 gFamPAV --pav ${out}/${out}_all.pav --fam $fam --out $out/$out.fpav");
		print STDOUT "\n";
	}

        sim_analysis($0, $out);

        cov_analysis($0, $pheno, $out);
        common_analysis($0, $pheno, $out, $focus);
        pheno_analysis($0, $pheno, "${out}_dispensable.pav", $out, "phenotype_association", $p_thre, $adjustp, $topn);
	if(!$rmele && $rep ne "len" && $rep ne "highcov" && $rep ne "none"){
		pheno_analysis($0, $pheno, "${out}_ele_dispensable.pav", $out, "phenotype_association_ele", $p_thre, $adjustp, $topn);
		ele_visual($0, $pheno, $bamdir, $gff, $out, $mergecov, $topn);
	}
	
}


sub generalBatch {
	my $usage = "\n\tapav generalBatch --bed <BED_file> --bamdir <BAM_dir> --pheno <phenotype_file> --out <res_prefix> [options]

Necessary input description:
  -i, --gff             <file>          GFF file.
  -o, --out             <string>        Output file prefix.
  -b, --bamdir          <dir>           The directory contains mapping results (sorted '.bam' files).
  --pheno               <file>          Phenotype file.


Options for coverage calculation:
  --mincov              <int>           Depth threshold.
                                        (Default:1)
  --rmele                               Remove analysis on elements.
  --merge				Merge neighboring elements with the same coverage(all 0 or all 1).
  --thread		<int>		Thread number.


Options for PAV Determination:
  --method		<type>          Method of determination: \"fixed\" or \"adaptive\".
                                        (Default:\"fixed\")

  fixed method:
  --thre                <float>         Coverage threshold.
                                        (Default:0.5)

  adaptive method:
  --mina                <float>         Min absence.
                                        (Default:0.1)
  --iter                <int>           Maximum number of iterations.
                                        (Default:100)


Options for the tracks in genome browser:
  --fa                  <file>          Fasta file of reference genome.
  --cpbam                               Copy bam files instead of make symbolic links.


Options for the phenotype association analysis:
  --p_threshold         <numeric>       The filter threshold for p_value/p_adjusted.
                                        (Default: 0.01)
  --adjust_p                            Adjust p_value.


Other Options:
  --focused_pheno       <string>        The phenotype in focus.
  --example_n           <int>           The number of examples used to show figures drawn by 'pavPlotPhenoVio/pavPlotPhenoBar' commands and 'elePlotCov/elePlotPAV/elePlotDepth' commands.

  -h, --help                            Print usage page.

        \n";

	my ($bed, $bamdir, $out, $pheno);
	my ($rep, $mincov, $mergecov, $rmele, $thread);
        my ($fa, $gff, $cpbam);
        my ($method, $thre, $mina, $iter);
        my ($focus, $p_thre, $adjustp);
        my $topn = 10;
        my $help;

        GetOptions(
                'bed=s'         => \$bed,
                'bamdir|b=s'    => \$bamdir,
                'out=s'         => \$out,
                'pheno=s'       => \$pheno,

		'mincov=i'      => \$mincov,
		'merge!'	=> \$mergecov,
		'thread=i'	=> \$thread,

                'fa=s'          => \$fa,
                'gff=s'         => \$gff,
		'cpbam!'	=> \$cpbam,

                'method=s'      => \$method,
                'thre=f'        => \$thre,
                'mina=f'        => \$mina,
                'iter=i'        => \$iter,

                'focused_pheno=s' => \$focus,
		'rmele!'	=> \$rmele,
                'p_threshold=f' => \$p_thre,
                'adjust_p!'     => \$adjustp,
                'example_n=i'   => \$topn,

                'help|h!'       => \$help

        ) or die $!."\n";

        die $usage if !defined($bed) & !defined($bamdir) & !defined($out) & !defined($pheno);
        die $usage if $help;

        APAVutils::check_file('--bed/-i', $bed);
	APAVutils::check_file('--pheno', $pheno);
        APAVutils::check_bamdir('--bamdir/-b', $bamdir);

	die "Please install 'samtools' first\n" if(system("command -v samtools > /dev/null 2>&1") != 0);

	die "Error: output directory \"$out\" already exists. To avoid overwriting of existing files, we kind
ly request that the output directory should not exist.\n" if -e $out;

        system("mkdir $out");

	print STDOUT "\n";

	my @arg_cov = ("$0 staCov --bed $bed --bamdir $bamdir --out $out/$out");
        push @arg_cov, "--mincov $mincov" if defined($mincov);
	push @arg_cov, "--rmele" if defined($rmele);
	push @arg_cov, "--merge" if defined($mergecov);
	push @arg_cov, "--thread $thread" if defined($thread);
        system(join(" ", @arg_cov));

	print STDOUT "\n";

        calling($0, '.cov', '', $bamdir, $pheno, $out, $fa, $gff, $method, $thre, $mina, $iter, $cpbam);

	print STDOUT "\n";

	if(!$rmele){
		calling($0, '_ele.cov', '_ele', undef, $pheno, $out, undef, undef, $method, $thre, $mina, $iter, $cpbam);
		print STDOUT "\n"; 
	}

        sim_analysis($0, $out);
        cov_analysis($0, $pheno, $out);
        common_analysis($0, $pheno, $out, $focus);
        pheno_analysis($0, $pheno, "${out}_dispensable.pav", $out, "phenotype_association", $p_thre, $adjustp, $topn);
        if(!$rmele){
                pheno_analysis($0, $pheno, "${out}_ele_dispensable.pav", $out, "phenotype_association_ele", $p_thre, $adjustp, $topn);
		ele_visual($0, $pheno, $bamdir, 0, $out, $mergecov, $topn);
        }

}

sub ele_visual{
	my ($c, $pheno, $bamdir, $gff, $out, $mergecov, $topn) = @_;
	my $eledir = "element_visualization";

	printLog("[elePlotCov/elePlotPAV/elePlotDepth] Give some examples for showing the coordinate positions and coverage/PAV/depth of elements...");

	system("mkdir $out/$eledir");
	my @ele_list = `grep -v -E "#|Annotation" ${out}/${out}_ele_dispensable.pav |  cut -f 5 | sed 's/:\\[.*\\]//g' | uniq | head -n $topn`;
	
	my $elecov = "${out}/${out}_ele.cov";
	my $elepav = "${out}/${out}_ele_all.pav";
	
	if(defined($mergecov)){
		$elecov = "${out}/${out}_ele.mcov";
	}
	foreach(@ele_list){
                chomp($_);
		my $cur_elecov = "$out/$eledir/${_}.ele.cov";
		my $cur_elepav = "$out/$eledir/${_}.ele.pav";
		if(defined($mergecov)){
			$cur_elecov = "$out/$eledir/${_}.ele.mcov";
		}
		if($gff){
			system("grep -E '${_}' $gff > $out/$eledir/${_}.gff");
			system("grep -E '${_}|Annotation' $elecov > $cur_elecov");
                        system("$c elePlotCov --elecov $cur_elecov --pheno $pheno --gff $out/$eledir/${_}.gff --out $cur_elecov");
                        system("$c elePlotDepth --ele $cur_elecov --bamdir $bamdir --pheno $pheno --gff $out/$eledir/${_}.gff --out $out/$eledir/${_}.ele.depth");
                        system("grep -E '${_}|Annotation' $elepav > $cur_elepav");
                        system("$c elePlotPAV --elepav $cur_elepav --pheno $pheno --gff $out/$eledir/${_}.gff --out $cur_elepav");		
		}else{
			system("grep -E '${_}|Annotation' $elecov > $cur_elecov");
			system("$c elePlotCov --elecov $cur_elecov --pheno $pheno --out $cur_elecov");
			system("$c elePlotDepth --ele $cur_elecov --bamdir $bamdir --pheno $pheno --out $out/$eledir/${_}.ele.depth");
			system("grep -E '${_}|Annotation' $elepav > $cur_elepav");
			system("$c elePlotPAV --elepav $cur_elepav --pheno $pheno --out $cur_elepav");
		}
        }
}


sub calling{

	my ($c, $cov_suffix, $pav_suffix, $bam_dir, $pheno, $out, $fa, $gff, $method, $thre, $mina, $iter, $cpbam) = @_;

	my @arg_callpav = ("$c callPAV --cov $out/${out}${cov_suffix} --out $out/$out${pav_suffix}");
        push @arg_callpav, "--pheno $pheno" if defined($pheno);
        push @arg_callpav, "--fa $fa" if defined($fa);
        push @arg_callpav, "--gff $gff" if defined($gff);
        push @arg_callpav, "--bamdir $bam_dir" if defined($bam_dir);
        push @arg_callpav, "--method $method" if defined($method);
        push @arg_callpav, "--thre $thre" if defined($thre);
        push @arg_callpav, "--mina $mina" if defined($mina);
        push @arg_callpav, "--iter $iter" if defined($iter);
	push @arg_callpav, "--cpbam" if defined($cpbam);
        system(join(" ", @arg_callpav));
}


sub sim_analysis {

	my ($c, $out) = @_;

	system("mkdir ${out}/estimation");

        system("$c pavSize --pav ${out}/${out}_all.pav --out ${out}/estimation/${out}_all.size");

        die "Please install R first\n" if(system("command -v Rscript > /dev/null 2>&1") != 0);
	die "Please install 'APAVplot' R package first\n" if(system("Rscript -e 'library(APAVplot)' > /dev/null 2>&1") != 0);

        printLog("[pavPlotSize] Plot the growth curve of genome estimation...");
        system("$c pavPlotSize --size ${out}/estimation/${out}_all.size --y_title \"Gene Number\" --out ${out}/estimation/${out}_size_curve ");
        system("$c pavPlotSize --size ${out}/estimation/${out}_all.size --y_title \"Gene Number\" --out ${out}/estimation/${out}_size_curve_delta --data_type increasing");
}

sub cov_analysis {
	my ($c, $pheno, $out) = @_;
	system("mkdir ${out}/coverage_visualization");
	printLog("[covPlotHeat] Plot the heatmap of coverage profile...");
        system("$c covPlotHeat --cov ${out}/$out.cov --pheno $pheno --cluster_rows --cluster_columns --out ${out}/coverage_visualization/${out}_cov_heatmap ");
}



sub common_analysis{

	my ($c, $pheno, $out, $focus) = @_;

	system("mkdir ${out}/common_analysis");

        printLog("[pavPlotHist] Plot a histogram to show the classifications and distribution of target regions...");
        system("$c pavPlotHist --pav ${out}/${out}_all.pav --out ${out}/common_analysis/${out}_all_pav_hist");

        printLog("[pavPlotStat] Plot a half-violin to show the number of target regions in samples...");
        if($focus){
                system("$c pavPlotStat --pav ${out}/${out}_all.pav --out ${out}/common_analysis/${out}_all_pav_sta --add_pheno_info $focus" );
        }else{
                system("$c pavPlotStat --pav ${out}/${out}_all.pav --out ${out}/common_analysis/${out}_all_pav_sta" );
        }

	printLog("[pavPlotHeat] Plot a heatmap to show the PAV profile of dispensable regions...");
	system("$c pavPlotHeat --pav ${out}/${out}_dispensable.pav --pheno $pheno --out ${out}/common_analysis/${out}_dispensable_pav_heatmap");

        if($focus){
                printLog("[pavPlotBar] Plot a clustered bar chart to show the classifications of target regions in all samples...");
                system("$c pavPlotBar --pav ${out}/${out}_all.pav --pheno $pheno --add_pheno_info $focus --out ${out}/common_analysis/${out}_all_pav_stackbar");
                printLog("[pavCluster] Cluster samples based on PAV table of dispensable regions...");
                system("$c pavCluster --pav ${out}/${out}_dispensable.pav --pheno $pheno --add_pheno_info $focus --out ${out}/common_analysis/${out}_dispensable_pav_cluster");
                printLog("[pavPCA] Perform PCA analysis for PAV table of dispensable regions...");
                system("$c pavPCA --pav ${out}/${out}_dispensable.pav --pheno $pheno --add_pheno_info $focus --out ${out}/common_analysis/${out}_dispensable_pav_pca");
        }else{
                printLog("[pavPlotBar] Plot a clustered bar chart to show the classifications of target regions in all samples...");
                system("$c pavPlotBar --pav ${out}/${out}_all.pav --out ${out}/common_analysis/${out}_all_pav_stackbar");
                printLog("[pavCluster] Cluster samples based on PAV table of dispensable regions...");
                system("$c pavCluster --pav ${out}/${out}_dispensable.pav --out ${out}/common_analysis/${out}_dispensable_pav_cluster");
                printLog("[pavPCA] Perform PCA analysis for PAV table of dispensable regions...");
                system("$c pavPCA --pav ${out}/${out}_dispensable.pav --out ${out}/common_analysis/${out}_dispensable_pav_pca");
        }

}


sub pheno_analysis {

	my ($c, $pheno, $pav, $out, $dir, $p_thre, $adjustp, $topn) = @_;

	system("mkdir ${out}/$dir");

        printLog("[pavStaPheno] Determine phenotype association...");
        system("$c pavStaPheno --pav ${out}/$pav --pheno $pheno --out ${out}/$dir/${out}_dispensable.phenores");

        printLog("[pavPlotPhenoHeat] Plot the main result of phenotype association analysis in a heatmap...");
        my @arg_phen_heat = ("$c pavPlotPhenoHeat --pav ${out}/$pav --pheno_res ${out}/$dir/${out}_dispensable.phenores --out ${out}/$dir/${out}_dispensable_pheno_heatmap");
        push @arg_phen_heat, "--p_threshold $p_thre" if defined($p_thre);
        push @arg_phen_heat, "--adjust_p" if defined($adjustp);
        system(join(" ", @arg_phen_heat));


        my @phen_list = phen_type($pheno);
        my @phen_str = @{$phen_list[0]};
        my @phen_num = @{$phen_list[1]};

        printLog("[pavPlotPhenoMan] Plot the Manhattan chart for each phenotype...");
        foreach (@phen_str, @phen_num) {
		print STDOUT "-- ${_}\n";
                my @arg_phen_man = ("$c pavPlotPhenoMan --pav ${out}/$pav --pheno $pheno --pheno_res ${out}/$dir/${out}_dispensable.phenores --pheno_name $_ --out ${out}/$dir/${out}_dispensable_pheno_${_}_manhattan");
                push @arg_phen_man, "--adjust_p" if defined($adjustp);
                system(join(" ", @arg_phen_man));
        }

        printLog("[pavPlotPhenoBlock] Plot the 'block chart' for each discrete phenotype...");
        foreach (@phen_str) {
		print STDOUT "-- ${_}\n";
                my @arg_phen_block = ("$c pavPlotPhenoBlock --pav ${out}/$pav --pheno $pheno --pheno_res ${out}/$dir/${out}_dispensable.phenores --pheno_name $_ --out ${out}/$dir/${out}_dispensable_pheno_${_}_block");
                push @arg_phen_block, "--p_threshold $p_thre" if defined($p_thre);
                push @arg_phen_block, "--adjust_p" if defined($adjustp);
                system(join(" ", @arg_phen_block));
        }

        my @phen_res = `head -n $topn ${out}/$dir/${out}_dispensable.phenores`;
        shift @phen_res;
        printLog("[pavPlotPhenoVio/pavPlotPhenoBar] Give some examples for showing the relationship between a genomic region and a certain phenotype in a barplot or violin plot");
        system("mkdir ${out}/$dir/example");
        foreach (@phen_res){
                $_ =~ s/\r\n$//;
                chomp($_);
                my @cur = split(/\t/, $_);
                my $cur_phen = $cur[0];
                my $cur_region = $cur[1];
		$cur_region =~ s/\(/\\\(/g;
		$cur_region =~ s/\)/\\\)/g;

                if(grep {$_ eq $cur_phen} @phen_num){
                        system("$c pavPlotPhenoVio --pav ${out}/$pav --pheno $pheno --region_name $cur_region --pheno_name $cur_phen --out '${out}/$dir/example/${out}_dispensable_pheno_${cur_region}_${cur_phen}_violin'");
                }else{
                        system("$c pavPlotPhenoBar --pav ${out}/$pav --pheno $pheno --region_name $cur_region --pheno_name $cur_phen --out '${out}/$dir/example/${out}_dispensable_pheno_${cur_region}_${cur_phen}_barplot'");
                }
        }

}


sub phen_type {
	my ($file) = @_;
	open(PHEN, "<$file") or die "Could not open file '$file' \n";
	my $header = <PHEN>;
	$header =~ s/\r\n$//;
	chomp($header);
	my @phens = split(/\t/, $header);
	shift @phens;

	my @is_num_col = (1) x scalar(@phens);

	while(<PHEN>){
		$_ =~ s/\r\n$//;
                chomp $_;

		my @fields = split(/\t/, $_);
		shift @fields;
		for (my $col = 0; $col < @fields; $col++) {
        		if ($fields[$col] !~ /^\d+$/) {
            			$is_num_col[$col] = 0;
        		}
    		}
	}

	close(PHEN);

	my @col_str_idx = grep{$is_num_col[$_] == 0} 0..$#phens;
	my @col_str = map{$phens[$_]} @col_str_idx;
	my @col_num_idx = grep{$is_num_col[$_] == 1} 0..$#phens;
	my @col_num = map{$phens[$_]} @col_num_idx;
	
	return (\@col_str, \@col_num);
}



sub printLog{
        my ($info) = @_;
        my $stime = `date +"%Y-%m-%d %H:%M:%S"`;
        chomp($stime);
        print STDOUT "\n[".$stime."] $info\n";
}

1;
