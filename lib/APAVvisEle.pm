#!/usr/bin/perl

package APAVvisEle;

use strict;
use warnings;
use Getopt::Long;
use List::Util qw/min max/;
use Data::Dumper;
use File::Basename;


sub elePlotDepth{

        my $usage = "\n\tUsage: apav elePlotDepth --ele <element_data> --bamdir <bam_dir> [options]

The script will call samtools program, so please make sure samtools is in your PATH.

Necessary input description:
  -i, --ele                     <file>          Element file.
  -b, --bamdir                  <dir>           The directory contains mapping results (sorted '.bam' files).

Option:
  --gff                         <file>          GFF file of ONE target region.
  --pheno			<file>          Phenotype file.
  -o, --out                     <string>        Prefix of the result file.

Visualization options:
  --fig_width           	<numeric>       The width of the figure.
  --fig_height          	<numeric>       The height of the figure.

  --cluster_samples                             Perform clustering on samples.
  --clustering_distance         <type>          Method to measure distance,
                                                choosing from 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.
                                                (Default: 'euclidean')
  --clustering_method           <type>          Method to perform hierarchical clustering,
                                                choosing from 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.
                                                (Default: 'complete')
  --log10                                       Convert depth value to log10(depth+1).

  --ele_region                  <type>          Gene element regions for display,
                                                choosing from 'TN'(transcripts number), 'UP'(upstream) and 'DOWN'(downstream).

  --top_anno_height		<numeric>	The relative height of annotation, ranging from 0 to 1.	
  --left_anno_width		<numeric>  	The relative height of annotation, ranging from 0 to 1.

  --ele_color			<color>		The color of element regions.
  --depth_colors		<colors>	The colors of depth, separated by commas.
  --lowlight_colors		<colors>	The colors for depth of non-element region, separated by commas.

  --gene_colors			<colors>	The colors of gene elements, separated by commas.
                                                The colors are assigned to 'transcript', 'exon', 'CDS', 'five_prime_UTR' and 'three_prime_UTR'.
  --pheno_info_color_list       <key=value>     The colors of phenotype annotations. (eg: 'gender=blue,green')
  --pheno_border		<color>		The color of cell borders in the phenotype annotation.

  --seq_name_size               <numeric>       The size of sequence name.
  --loc_name_size               <numeric>       The size of coordinate numbers.

  --hide_sample_name                            Hide sample names.
  --sample_name_size            <numeric>       The size of sample names.

  --legend_title		<numeric>	The title of legend.
  --legend_title_size           <numeric>       The size of legend title.
  --legend_text_size            <numeric>       The size of legend item labels.

  -h, --help                                    Print usage page.
        
        \n";

        my ($eledata, $bam_dir, $out, $help);

	my $ele_region = "NULL";

        my $gffdata = "NULL";
        my $phenodata = "NULL";

	my $cluster = 0;
        my $clustering_distance = "euclidean";
        my $clustering_method = "complete";

	my $log10 = 0;

	my $fig_width = 10;
        my $fig_height = 8;

        my $top_anno_height = 0.4;
        my $left_anno_width = 0.1;

        my $ele_color = "#A6CEE3";
        my $depth_colors = "#5680ae,white";
        my $lowlight_colors = "gray,white";

        my $gene_colors = "gray70,gray85,#d77478,#ecbf40,#ecbf40";
        my (%pheno_info_color_list, $pheno_info_color_list_str);
	my $pheno_border = "NA";

        my $seq_name_size = "NULL";
        my $loc_name_size = "NULL";

        my $hide_sample_name = 0; ## F
        my $sample_name_size = "NULL";

        my $legend_title = "Depth";
        my $legend_title_size = "NULL";
        my $legend_text_size = "NULL";

	GetOptions(
		'ele|i=s'		=> \$eledata,
		'bamdir|b=s'		=> \$bam_dir,

		'ele_region=s'		=> \$ele_region,

		'out|o=s'		=> \$out,
		'pheno=s'               => \$phenodata,
                'gff=s'                 => \$gffdata,

		'cluster_samples!'	=> \$cluster,
		'clustering_distance=s' => \$clustering_distance,
                'clustering_method=s'   => \$clustering_method,

		'log10!'		=> \$log10,

                'fig_width=f'           => \$fig_width,
                'fig_height=f'          => \$fig_height,

                'top_anno_height=f'     => \$top_anno_height,
                'left_anno_width=f'     => \$left_anno_width,

                'ele_color=s'           => \$ele_color,
                'depth_colors=s'	=> \$depth_colors,
		'lowlight_colors=s'	=> \$lowlight_colors,
                
                'gene_colors=s'         => \$gene_colors,
                'pheno_info_color_list=s'       => \%pheno_info_color_list,
		'pheno_border=s'	=> \$pheno_border,

                'seq_name_size=f'       => \$seq_name_size,
                'loc_name_size=f'       => \$loc_name_size,

                'hide_sample_name!'     => \$hide_sample_name,
                'sample_name_size=f'    => \$sample_name_size,

		'legend_title=s'	=> \$legend_title,		
                'legend_title_size=f'   => \$legend_title_size,
                'legend_text_size=f'    => \$legend_text_size,

                'help|h!'               => \$help
	) or die $!."\n";

	die $usage if !defined($eledata) & !defined($bam_dir);
	die $usage if $help;

	APAVutils::check_file('--ele/-i', $eledata);
	APAVutils::check_file('--pheno', $phenodata) if defined($phenodata);
        APAVutils::check_file('--gff', $gffdata) if defined($gffdata);

        my $head = `grep -v '#' $eledata | head -n 1`;
        if(!($head =~ "^Chr\tStart\tEnd\tLength\tAnnotation")){
                die "Missing the header line\n";
        }
	APAVutils::check_file('--bamdir/-b', $bam_dir);
	$out = APAVutils::check_out($out, "", $eledata.".depth");

        die "Please install 'samtools' first\n" if(system("command -v samtools > /dev/null 2>&1") != 0);
	die "Please install R first\n" if(system("command -v Rscript > /dev/null 2>&1") != 0);
        die "Please install 'APAVplot' R package first\n" if(system("Rscript -e 'library(APAVplot)' > /dev/null 2>&1") != 0);

	if(%pheno_info_color_list){
                for my $key (keys %pheno_info_color_list){
                        $pheno_info_color_list_str.= $key."=".$pheno_info_color_list{$key}.";";
                }
                $pheno_info_color_list_str =~ s/;$//g;
        } else {
                $pheno_info_color_list_str = "NULL";
        }


        my ($chr, $ele_start, $ele_end) = get_ele_range($eledata);
	my ($start, $end);

	if($gffdata ne "NULL"){
		my ($gff_chr, $gff_start, $gff_end) = get_gene_range($gffdata);
		die "Please check the GFF file, the chromosome should be the same as the element data\n" if($gff_chr ne $chr);
		$start = min($ele_start, $gff_start);
		$end = max($ele_end, $gff_end);
	}else{
		$start = $ele_start;
		$end = $ele_end;
	}

        $bam_dir.="/" unless($bam_dir=~/\/$/);
        my @bams = <$bam_dir*.bam>;
	die "Can not find bam files in the given directory\n" if $#bams == -1;

        my @samples = @bams;
        map{ $_ =~ s/.*\///g }@samples;
        map{ $_ =~ s/.bam$//g }@samples;

        my $out_depth = $out.".data";

        system("mkdir $out_depth.tmp");
        system("samtools depth -a -r $chr:$start-$end $bams[0] > $out_depth.tmp/$samples[0]");

        foreach(1..$#bams){
                system("samtools depth -a -r $chr:$start-$end $bams[$_] | awk '{print \$3}' >  $out_depth.tmp/$samples[$_]");
        }

        system("paste ".join(" ", map{ "$out_depth.tmp/".$_} @samples)." > $out_depth");
        system("rm -r $out_depth.tmp");

        system("sed -i '1i\\Chr\tPos\t".join("\t", @samples)."' $out_depth");

	my $dir =  dirname(__FILE__);
        my $exec = $dir."/"."vis_ele.R";

	system("Rscript $exec $dir '$eledata' '$ele_region' '$out' '$phenodata' '$gffdata' depth '$out_depth' '$fig_width' '$fig_height' '$cluster' '$clustering_distance' '$clustering_method' '$log10' '$top_anno_height' '$left_anno_width' '$ele_color' '$depth_colors' '$lowlight_colors' '$gene_colors' '$pheno_info_color_list_str' '$pheno_border' '$seq_name_size' '$loc_name_size' '$hide_sample_name' '$sample_name_size' '$legend_title' '$legend_title_size' '$legend_text_size' 1>/dev/null ");
	

}

sub get_ele_range{

        my ($eledata) = @_;

        my @ccol = `cut -f 1 $eledata | uniq`;
        my @scol = `cut -f 2 $eledata`;
        my @ecol = `cut -f 3 $eledata`;

        shift(@ccol);
        die "Please check the element data, all elements should be located on one sequence\n" if($#ccol > 0);
        my $chr = $ccol[0];
        chomp($chr);

        my @starts = split(/,|;/, join(",", @scol));
        my @ends = split(/,|;/, join(",", @ecol));
        shift(@starts);
        shift(@ends);

        my $start = min(@starts);
        my $end = max(@ends);
        chomp($start);
        chomp($end);

        return ($chr, $start, $end);
}

sub get_gene_range{
	
	my ($gffdata) = @_;

	my @ccol = `cut -f 1 $gffdata | uniq`;
	my @scol = `cut -f 4 $gffdata`;
	my @ecol = `cut -f 5 $gffdata`;
	die "Please check the gff data, all regions should be located on one sequence\n" if($#ccol > 0);
	my $chr = $ccol[0];
	chomp($chr);

	my $start = min(@scol);
	my $end = max(@ecol);
	chomp($start);
	chomp($end);

	return ($chr, $start, $end);
}


sub elePlotPAV {
	my $usage = "\n\tUsage: apav elePlotPAV --elepav <element_PAV_file> [options]

Necessary input desription:
  -i, --elepav                  <file>          Element PAV file of ONE target region.

Options:
  --gff                         <file>          GFF file of gene.
  --pheno                       <file>          Phenotype file.
   -o, --out                    <string>        Figure name.

Visualization options:
  --fig_width                   <numeric>       The width of the figure.
  --fig_height                  <numeric>       The height of the figure.

  --cluster_samples                             Perform clustering on samples.
  --clustering_distance         <type>          Method to measure distance,
                                                choosing from 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.
                                                (Default: 'euclidean')
  --cluster_method              <type>          Method to perform hierarchical clustering,
                                                choosing from 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.
                                                (Default: 'complete')
  --ele_region                  <type>          Gene element regions for display,
                                                choosing from 'TN'(transcripts number), 'UP'(upstream) and 'DOWN'(downstream).

  --top_anno_height             <numeric>       The relative height of annotation, ranging from 0 to 1.
  --left_anno_width             <numeric>       The relative width of annotation, ranging from 0 to 1.

  --ele_color                   <color>         The color of element regions.
  --ele_line_color              <color>         The color for the line pointing to elements.
  --pav_colors                  <colors>        The colors for presence and absence, separated by commas.
  --cell_border			<color>		The color of cell borders in the heat map.
  --gene_colors                 <colors>        The colors for gene elements, separated by commas.
                                                The colors are assigned to 'transcript', 'exon', 'CDS', 'five_prime_UTR' and 'three_prime_UTR'.
  --pheno_info_color_list       <key=value>     The colors for phenotype annotations. (eg: 'gender=blue,green')
  --pheno_border		<color>		The color of cell borders in the phenotype annotation.

  --seq_name_size               <numeric>       The size of sequence name.
  --loc_name_size               <numeric>       The size of coordinate numbers.

  --hide_ele_name                               Hide element names.                             
  --ele_name_size               <numeric>       The size of element names.
  --ele_name_rot              	<numeric>       The rotation of element names.

  --hide_sample_name                            Hide sample names.
  --sample_name_size            <numeric>       The size of sample names.

  --legend_title_size           <numeric>       The size of legend title.
  --legend_text_size            <numeric>       The size of legend item labels.

  -h, --help                                    Print usage page.
	\n";

	my ($eledata, $out, $help);

	my $ele_region = "NULL";

        my $phenodata = "NULL";
        my $gffdata = "NULL";

	my $cluster = 0;
        my $clustering_distance = "euclidean";
        my $clustering_method = "complete";

        my $fig_width = 10;
        my $fig_height = 8;

        my $top_anno_height = 0.4;
        my $left_anno_width = 0.1;

        my $ele_color = "#A6CEE3";
        my $ele_line_color = "gray";
        my $pav_colors = "#5680ae,gray90";
	my $cell_border = "NA";
        my $gene_colors = "gray70,gray85,#d77478,#ecbf40,#ecbf40";
        my (%pheno_info_color_list, $pheno_info_color_list_str);
	my $pheno_border = "NA";

        my $seq_name_size = "NULL";
        my $loc_name_size = "NULL";

        my $hide_ele_name = 0; ## F
        my $ele_name_size = "NULL";
        my $ele_name_rot = 90;

        my $hide_sample_name = 0; ## F
        my $sample_name_size = "NULL";

        my $legend_title_size = "NULL";
        my $legend_text_size = "NULL";


	GetOptions(
                'elepav|i=s'            => \$eledata,
		'ele_region=s'          => \$ele_region,

                'out|o=s'               => \$out,
                'pheno=s'               => \$phenodata,
                'gff=s'                 => \$gffdata,

		'cluster_samples!'	=> \$cluster,
                'clustering_distance=s' => \$clustering_distance,
                'clustering_method=s'   => \$clustering_method,

                'fig_width=f'           => \$fig_width,
                'fig_height=f'          => \$fig_height,

                'top_anno_height=f'     => \$top_anno_height,
                'left_anno_width=f'     => \$left_anno_width,

                'ele_color=s'           => \$ele_color,
                'ele_line_color=s'      => \$ele_line_color,
                'pav_colors=s'          => \$pav_colors,
		'cell_border=s'		=> \$cell_border,
                'gene_colors=s'         => \$gene_colors,
                'pheno_info_color_list=s'       => \%pheno_info_color_list,
		'pheno_border=s'	=> \$pheno_border,		

                'seq_name_size=f'       => \$seq_name_size,
                'loc_name_size=f'       => \$loc_name_size,

                'hide_ele_name!'        => \$hide_ele_name,
                'ele_name_size=f'       => \$ele_name_size,
                'ele_name_rot=f'        => \$ele_name_rot,

                'hide_sample_name!'     => \$hide_sample_name,
                'sample_name_size=f'    => \$sample_name_size,

                'legend_title_size=f'   => \$legend_title_size,
                'legend_text_size=f'    => \$legend_text_size,

                'help|h!'               => \$help
        ) or die $!."\n";

	die $usage if !defined($eledata);
        die $usage if $help;

        APAVutils::check_file('--elepav/-i', $eledata);
	APAVutils::check_file('--pheno', $phenodata) if defined($phenodata);
        APAVutils::check_file('--gff', $gffdata) if defined($gffdata);

        my $head = `grep -v '#' $eledata | head -n 1`;
        if(!($head =~ "^Chr\tStart\tEnd\tLength\tAnnotation")){
                die "Missing the header line\n";
        }
        $out = APAVutils::check_out($out, "", $eledata);

        die "Please install R first\n" if(system("command -v Rscript > /dev/null 2>&1") != 0);
	die "Please install 'APAVplot' R package first\n" if(system("Rscript -e 'library(APAVplot)' > /dev/null 2>&1") != 0);

        if(%pheno_info_color_list){
                for my $key (keys %pheno_info_color_list){
                        $pheno_info_color_list_str.= $key."=".$pheno_info_color_list{$key}.";";
                }
                $pheno_info_color_list_str =~ s/;$//g;
        } else {
                $pheno_info_color_list_str = "NULL";
        }

        my $dir =  dirname(__FILE__);
        my $exec = $dir."/"."vis_ele.R";

        system("Rscript $exec $dir '$eledata' '$ele_region' '$out' '$phenodata' '$gffdata' pav '$fig_width' '$fig_height' '$cluster' '$clustering_distance' '$clustering_method' '$top_anno_height' '$left_anno_width' '$ele_color' '$ele_line_color' '$pav_colors' '$cell_border' '$gene_colors' '$pheno_info_color_list_str' '$pheno_border' '$seq_name_size' '$loc_name_size' '$hide_ele_name' '$ele_name_size' '$ele_name_rot' '$hide_sample_name' '$sample_name_size' '$legend_title_size' '$legend_text_size' 1>/dev/null \n");


}

sub elePlotCov {

	my $usage = "\n\tUsage: apav elePlotCov --elecov <element_coverage_file> [options]

Necessary input desription:
  -i, --elecov	 		<file>		Element coverage file of ONE target region.

Options:
  --gff				<file>		GFF file of gene.
  --pheno			<file>		Phenotype file.
  -o, --out                     <string>        Figure name.

Visualization options:
  --fig_width                   <numeric>       The width of the figure.
  --fig_height                  <numeric>       The height of the figure.

  --cluster_samples                             Perform clustering on samples.
  --clustering_distance         <type>          Method to measure distance,
                                                choosing from 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.
                                                (Default: 'euclidean')
  --cluster_method              <type>          Method to perform hierarchical clustering,
                                                choosing from 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.
                                                (Default: 'complete')
  --ele_region                  <type>          Element region for display,
                                                choosing from 'TN'(transcripts number), 'UP'(upstream) and 'DOWN'(downstream).

  --top_anno_height		<numeric>	The relative height of annotation, ranging from 0 to 1.
  --left_anno_width		<numeric>	The relative width of annotation, ranging from 0 to 1.

  --ele_color			<color>		The color of element regions.
  --ele_line_color		<color>		The color for the line pointing to elements.
  --cov_colors			<colors>	The colors for coverage, separated by commas.
  --cell_border			<color>		The color of cell borders in the heat map.
  --gene_colors			<colors>	The colors for gene elements, separated by commas.
  						The colors are assigned to 'transcript', 'exon', 'CDS', 'five_prime_UTR' and 'three_prime_UTR'.
  --pheno_info_color_list	<key=value>	The colors for phenotype annotations. (eg: 'gender=blue,green')
  --pheno_border		<color>		The color of cell borders in the phenotype annotation.

  --seq_name_size		<numeric>	The size of sequence name.
  --loc_name_size		<numeric>	The size of coordinate numbers.

  --hide_ele_name				Hide element names.				
  --ele_name_size		<numeric>	The size of element names.
  --ele_name_rot		<numeric>	The rotation of element names.

  --hide_sample_name				Hide sample names.
  --sample_name_size		<numeric>	The size of sample names.

  --legend_title_size		<numeric>	The size of legend title.
  --legend_text_size		<numeric>	The size of legend item labels.

  -h, --help					Print usage page.	
	\n";

	my ($eledata, $out, $help);

	my $ele_region = "NULL";

	my $phenodata = "NULL";
	my $gffdata = "NULL";

	my $cluster = 0;
	my $clustering_distance = "euclidean";
        my $clustering_method = "complete";

	my $fig_width = 10;
	my $fig_height = 8;

	my $top_anno_height = 0.4;
	my $left_anno_width = 0.1;

	my $ele_color = "#A6CEE3";
	my $ele_line_color = "gray";
	my $cov_colors = "#5680ae,gray90";
	my $cell_border = "NA";
	my $gene_colors = "gray70,gray85,#d77478,#ecbf40,#ecbf40";
	my (%pheno_info_color_list, $pheno_info_color_list_str);
	my $pheno_border = "NA";

	my $seq_name_size = "NULL";
	my $loc_name_size = "NULL";

	my $hide_ele_name = 0; ## F
	my $ele_name_size = "NULL";
	my $ele_name_rot = 90;

	my $hide_sample_name = 0; ## F
	my $sample_name_size = "NULL";

	my $legend_title_size = "NULL";
	my $legend_text_size = "NULL";

	GetOptions(
		'elecov|i=s'		=> \$eledata,
		'ele_region=s'          => \$ele_region,

		'out|o=s'		=> \$out,
		'pheno=s'		=> \$phenodata,
		'gff=s'			=> \$gffdata,

		'cluster_samples!'	=> \$cluster,
		'clustering_distance=s'	=> \$clustering_distance,
		'clustering_method=s'	=> \$clustering_method,

		'fig_width=f'		=> \$fig_width,
		'fig_height=f'		=> \$fig_height,

		'top_anno_height=f'	=> \$top_anno_height,
		'left_anno_width=f'	=> \$left_anno_width,

		'ele_color=s'		=> \$ele_color,
		'ele_line_color=s'	=> \$ele_line_color,
		'cov_colors=s'		=> \$cov_colors,
		'cell_border=s'		=> \$cell_border,
		'gene_colors=s'		=> \$gene_colors,
		'pheno_info_color_list=s'	=> \%pheno_info_color_list,
		'pheno_border=s'	=> \$pheno_border,
	
		'seq_name_size=f'	=> \$seq_name_size,
		'loc_name_size=f'	=> \$loc_name_size,

		'hide_ele_name!'	=> \$hide_ele_name,
		'ele_name_size=f'	=> \$ele_name_size,
		'ele_name_rot=f'	=> \$ele_name_rot,

		'hide_sample_name!'	=> \$hide_sample_name,
		'sample_name_size=f'	=> \$sample_name_size,

		'legend_title_size=f'	=> \$legend_title_size,
		'legend_text_size=f'	=> \$legend_text_size,

		'help|h!'       	=> \$help
	) or die $!."\n";

	die $usage if !defined($eledata);
	die $usage if $help;

	APAVutils::check_file('--elecov/-i', $eledata);
	APAVutils::check_file('--pheno', $phenodata) if defined($phenodata);
	APAVutils::check_file('--gff', $gffdata) if defined($gffdata);

        my $head = `grep -v '#' $eledata | head -n 1`;
        if(!($head =~ "^Chr\tStart\tEnd\tLength\tAnnotation")){
                die "Missing the header line\n";
        }
	$out = APAVutils::check_out($out, "", $eledata);

	die "Please install R first\n" if(system("command -v Rscript > /dev/null 2>&1") != 0);
	die "Please install 'APAVplot' R package first\n" if(system("Rscript -e 'library(APAVplot)' > /dev/null 2>&1") != 0);

	if(%pheno_info_color_list){
                for my $key (keys %pheno_info_color_list){
                        $pheno_info_color_list_str.= $key."=".$pheno_info_color_list{$key}.";";
                }
                $pheno_info_color_list_str =~ s/;$//g;
        } else {
                $pheno_info_color_list_str = "NULL";
        }

        my $dir =  dirname(__FILE__);
        my $exec = $dir."/"."vis_ele.R";

	system("Rscript $exec $dir '$eledata' '$ele_region' '$out' '$phenodata' '$gffdata' cov '$fig_width' '$fig_height' '$cluster' '$clustering_distance' '$clustering_method' '$top_anno_height' '$left_anno_width' '$ele_color' '$ele_line_color' '$cov_colors' '$cell_border' '$gene_colors' '$pheno_info_color_list_str' '$pheno_border' '$seq_name_size' '$loc_name_size' '$hide_ele_name' '$ele_name_size' '$ele_name_rot' '$hide_sample_name' '$sample_name_size' '$legend_title_size' '$legend_text_size' 1>/dev/null \n");


}


1;
