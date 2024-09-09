#!/usr/bin/perl

package APAVvisPAV;

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use APAVutils;

sub pavPCA {
	my $usage = "\nUsage: apav pavPCA --pav <pav_file> [options]
	
Necessary input description:
  -i, --pav         	<file>		PAV file produced by command 'apav callpav'.
        
Options:
  --pheno       	<file>	        Phenotype file.
  --add_pheno_info      <string>        Add phenotype information.
  -o, --out             <string>        Figure name.
  
Visualization options:
  --fig_width		<numeric>	The width of the figure.
  --fig_height		<numeric>	The height of the figure.

  --pheno_info_colors	<colors>	The colors for the phenotype groups, separated by commas.
  
  --axis_text_size 	<numeric>	The size of tick labels on axis.
  --axis_title_size 	<numeric>	The size of axis title.
  
  --legend_side 	<type>		The position of legend, choosing from 'top', 'bottom', 'right', 'left'.
  --legend_text_size 	<numeric>	The size of legend item labels.
  --legend_title_size 	<numeric>	The size of legend title.

  -h, --help                    	Print usage.
	\n";

	my ($pavdata, $out, $help);

        my $phenodata = "NULL";

	my $fig_width = 8;
        my $fig_height = 6;

	my $add_pheno_info = "NULL";
	my $pheno_info_colors = "NULL";
	my $axis_text_size = "NULL";
	my $axis_title_size = "NULL";
	my $legend_side = "right";
	my $legend_text_size = "NULL";
	my $legend_title_size = "NULL";

	GetOptions(
                'pav|i=s'      		        => \$pavdata,
                'pheno=s'                     	=> \$phenodata,
                'out|o=s'                       => \$out,

		'fig_width=f'			=> \$fig_width,
		'fig_height=f'			=> \$fig_height,
                
		'add_pheno_info=s'		=> \$add_pheno_info,
		'pheno_info_colors=s'		=> \$pheno_info_colors,
		'axis_text_size=f'		=> \$axis_text_size,
		'axis_title_size=f'		=> \$axis_title_size,
		'legend_side=s'			=> \$legend_side,
		'legend_text_size=f'		=> \$legend_text_size,
		'legend_title_size=f'		=> \$legend_title_size,

		'help|h!'                       => \$help
        ) or die $!."\n";
	
	die $usage if !defined($pavdata);

        die $usage if $help;
	APAVutils::check_pav_input($pavdata);
	$out = APAVutils::check_out($out, $pavdata, "_pav_pca");
	APAVutils::check_file("--pheno", $phenodata) if defined($phenodata);

        die "Please install R first\n" if(system("command -v Rscript > /dev/null 2>&1") != 0);
	die "Please install 'APAVplot' R package first\n" if(system("Rscript -e 'library(APAVplot)' > /dev/null 2>&1") != 0);

        my $dir =  dirname(__FILE__);
        my $exec = $dir."/"."vis_pav.R";

	system("Rscript $exec $dir $pavdata $phenodata $out pca '0' '0' '0.1' '0' '0.05' '$fig_width' '$fig_height' '$add_pheno_info' '$pheno_info_colors' '$axis_text_size' '$axis_title_size' '$legend_side' '$legend_text_size' '$legend_title_size' '$fig_width' '$fig_height' 1>/dev/null ");
}

sub pavCluster {
	my $usage = "\nUsage: apav pavCluster --pav <pav_file> [options]
	
Necessary input description:
  -i, --pav         	<file>  	PAV file produced by command 'apav callpav'.
        
Options:
  --pheno		<file>	        Phenotype file.
  --add_pheno_info	<string>	Group the samples by this phenotype.
  -o, --out             <string>        Figure name.

  --clustering_distance	<type>          Method to measure distance,
                                        choosing from 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.
                                        (Default: 'euclidean')
  --clustering_method   <type>          Method to perform hierarchical clustering,
                                        choosing from 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.
                                        (Default: 'complete')

Visualization options:
  --fig_width           <numeric>       The width of the figure.
  --fig_height          <numeric>       The height of the figure.

  --pheno_info_colors	<colors> 	The colors for the phenotype groups, separated by commas.

  --sample_name_size	<numeric>	The size of labels.
  --mult 		<numeric>	A numebr of multiplicative range expansion factors.

  --legend_side 	<type>		The position of legend ('top', 'bottom', 'right', 'left').
  --legend_title_size 	<numeric>	The size of legend title.
  --legend_text_size 	<numeric>	The size of legend item labels.

  -h, --help				Print usage.  
	\n";
	
	my ($pavdata, $out, $help);

        my $phenodata = "NULL";

	my $fig_width = 8;
	my $fig_height = 8;

	my $clustering_distance = "euclidean";
	my $clustering_method = "complete";

	my $add_pheno_info = "NULL";
	my $pheno_info_colors = "NULL";

	my $sample_name_size = 4;
	my $mult = 0.1;
	my $legend_side = "right";
	my $legend_title_size = "NULL";
	my $legend_text_size = "NULL";

	GetOptions(
                'pav|i=s'	               	=> \$pavdata,
                'pheno=s'                     	=> \$phenodata,
                'out|o=s'                       => \$out,

		'fig_width=f'			=> \$fig_width,
		'fig_height=f'			=> \$fig_height,

		'clustering_distance=s'		=> \$clustering_distance,
		'clustering_method=s'		=> \$clustering_method,
		'add_pheno_info=s'		=> \$add_pheno_info,
		'pheno_info_colors=s'		=> \$pheno_info_colors,
		'sample_name_size=f'		=> \$sample_name_size,
		'mult=f'			=> \$mult,
		'legend_side=f'			=> \$legend_side,
		'legend_title_size=f'		=> \$legend_title_size,
		'legend_text_size=f'		=> \$legend_text_size,
                'help|h!'                       => \$help
        ) or die $!."\n";
	
	die $usage if !defined($pavdata);

        die $usage if $help;
	APAVutils::check_pav_input($pavdata);
	$out = APAVutils::check_out($out, $pavdata, "_pav_cluster");
	APAVutils::check_file("--pheno", $phenodata) if defined($phenodata);

        die "Please install R first\n" if(system("command -v Rscript > /dev/null 2>&1") != 0);
	die "Please install 'APAVplot' R package first\n" if(system("Rscript -e 'library(APAVplot)' > /dev/null 2>&1") != 0);

        my $dir =  dirname(__FILE__);
        my $exec = $dir."/"."vis_pav.R";

	system("Rscript $exec $dir $pavdata $phenodata $out cluster '0' '0' '0.1' '0' '0.05' '$fig_width' '$fig_height' '$clustering_distance' '$clustering_method' '$add_pheno_info' '$pheno_info_colors' '$sample_name_size' '$mult' '$legend_side' '$legend_title_size' '$legend_text_size' 1>/dev/null ");
}



sub pavPlotBar {
	my $usage = "\nUsage apav pavPlotBar --pav <pav_file> [options]
	
Necessary input description:
  -i, --pav    	    	<file>		PAV file produced by command 'apav callpav'.   
        
Options:
  --pheno       	<file>		Phenotype file.
  --add_pheno_info      <string>        Group the samples by this phenotype.
  -o, --out             <string>        Figure name.

  --rm_softcore               	 	Exclude `soft-core` when determining the classification of the region.
  --rm_private                 		Exclude `private` when determining the classification of the region.
  --softcore_loss_rate 	<numeric>	Loss rate.
  --use_binomial			Use binomial test when determining `soft-core`.
  --softcore_p_value	<numeric>	P-value (binomial test).

  --clustering_distance <type>          Method to measure distance,
                                        choosing from 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.
                                        (Default: 'euclidean')
  --clustering_method   <type>          Method to perform hierarchical clustering,
                                        choosing from 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.
                                        (Default: 'complete')

Visualization options:
  --show_relative 			Show relative value.
  
  --type_colors 	<colors>	The colors for classifications of the region, separated by commas.
  --phen_info_colors	<colors> 	The colors for the phenotype groups, separated by commas.
  
  --bar_width		<numeric>	The relative width of bars, ranging from 0 to 1.
  --sample_name_size 	<numeric>	The size of sample names.
  
  --legend_side 	<type>		The position of legend ('top', 'bottom', 'right', 'left').
  --legend_title 	<string>	The text for the legend title.
  --legend_title_size 	<numeric>	The size of legend title.
  --legend_text_size 	<numeric>	The size of legend item labels.
  
  --dend_width 		<numeric>	The relative width of dendrogram, ranging from 0 to 1.
  --name_width 		<numeric>	The relative width of sample names, ranging from 0 to 1.

   -h, --help                   	Print usage.
  \n";

	my ($pavdata, $out, $help);

	my $phenodata = "NULL";

	my $fig_width = 10;
	my $fig_height = 10;

        my $rm_softcore = 0; ## 
        my $rm_private = 0; ## 
        my $softcore_loss_rate = 0.1;
	my $use_binomial = 0; ## F
        my $softcore_p_value = 0.05;

	my $show_relative = 0; ## F
	my $add_pheno_info = "NULL";

	my $type_colors = "NULL";
	my $pheno_info_colors = "NULL";

	my $clustering_distance = "euclidean";
	my $clustering_method = "complete";

	my $bar_width = 0.8;
	my $sample_name_size = 3;

	my $legend_side = "right";
	my $legend_title = "Target Region";
	my $legend_title_size = "NULL";
	my $legend_text_size = "NULL";

	my $dend_width = 0.05;
	my $name_width = 0.15;

	GetOptions(
                'pav|i=s'               	=> \$pavdata,
                'pheno=s'                       => \$phenodata,
                'out|o=s'                       => \$out,

		'fig_width=f'           	=> \$fig_width,
                'fig_height=f'          	=> \$fig_height,

                'rm_softcore!'	                => \$rm_softcore,
                'rm_private!'                  	=> \$rm_private,
                'softcore_loss_rate=f'          => \$softcore_loss_rate,
		'use_binomial!'			=> \$use_binomial,
                'softcore_p_value=f'            => \$softcore_p_value,

                'show_relative!'		=> \$show_relative,
		'add_pheno_info=s'		=> \$add_pheno_info,
		'type_colors=s'			=> \$type_colors,
		'pheno_info_colors'		=> \$pheno_info_colors,
		'clustering_distance'		=> \$clustering_distance,
		'clustering_method'		=> \$clustering_method,
		'bar_width=f'			=> \$bar_width,
		'sample_name_size=f'		=> \$sample_name_size,
		'legend_side=s'			=> \$legend_side,
		'legend_title=s'		=> \$legend_title,
		'legend_title_size=f'		=> \$legend_title_size,
		'legend_text_size=f'		=> \$legend_text_size,
		'dend_width=f'			=> \$dend_width,
		'name_width=f'			=> \$name_width,
		'help|h!'                       => \$help
        ) or die $!."\n";

	die $usage if !defined($pavdata);

        die $usage if $help;
	APAVutils::check_pav_input($pavdata);
	$out = APAVutils::check_out($out, $pavdata, "_pav_stackbar");
	APAVutils::check_file("--pheno", $phenodata) if defined($phenodata);

        die "Please install R first\n" if(system("command -v Rscript > /dev/null 2>&1") != 0);
	die "Please install 'APAVplot' R package first\n" if(system("Rscript -e 'library(APAVplot)' > /dev/null 2>&1") != 0);

        my $dir =  dirname(__FILE__);
        my $exec = $dir."/"."vis_pav.R";

	system("Rscript $exec $dir $pavdata $phenodata $out stackbar  '$rm_softcore' '$rm_private' '$softcore_loss_rate' '$use_binomial' '$softcore_p_value' '$fig_width' '$fig_height' '$show_relative' '$add_pheno_info' '$type_colors' '$pheno_info_colors' '$clustering_distance' '$clustering_method' '$bar_width' '$sample_name_size' '$legend_side' '$legend_title' '$legend_title_size' '$legend_text_size' '$dend_width' '$name_width' 1>/dev/null ");

}

sub pavPlotStat {
	my $usage = "\nUsage: apav pavPlotStat --pav <pav_file> [options]

Necessary input description:
  -i, --pav             <file>          PAV file produced by command 'apav callpav'.

Options:
  --pheno		<file>		Phenotype file.
  --add_pheno_info      <string>        Group the samples by this phenotype
  -o, --out             <string>        Figure name.

Visualization options:
  --fig_width           <numeric>       The width of the figure.
  --fig_height          <numeric>       The height of the figure.

  --color		<color>		The color for halfviolin.
  --pheno_info_colors	<colors>	The colors for the phenotype groups, separated by commas.
  
  --x_text_size		<numeric>	The size of tick labels on the x-axis.
  --y_text_size		<numeric>	The size of tick labels on the y-axis.
  --x_title_size	<numeric>	The size of x-axis title.
  --y_title_size	<numeric>	The size of y-axis title.	
	\n";

	my ($pavdata, $out, $help);

	my $phenodata = "NULL";
	my $add_pheno_info = "NULL";

	my $fig_width = 8;
	my $fig_height = 6;

	my $color = "#7e9bc0";
	my $pheno_info_colors = "NULL";

	my $x_text_size = "NULL";
	my $y_text_size = "NULL";
	my $x_title_size = "NULL";
	my $y_title_size = "NULL";

	GetOptions(
		'pav_data|i=s'		=> \$pavdata,
		'pheno|p=s'		=> \$phenodata,
		'out|o=s'		=> \$out,

		'fig_width=f'		=> \$fig_width,
		'fig_height=f'		=> \$fig_height,

		'color=s'		=> \$color,

		'add_pheno_info=s'	=> \$add_pheno_info,
		'pheno_info_colors=s'	=> \$pheno_info_colors,
		
		'x_text_size=f'		=> \$x_text_size,
		'y_text_size=f'		=> \$y_text_size,
		'x_title_size=f'	=> \$x_title_size,
		'y_title_size=f'	=> \$y_title_size,

		'help|h!'		=> \$help
	) or die $!."\n";

	die $usage if !defined($pavdata);

        die $usage if $help;
	APAVutils::check_pav_input($pavdata);
	$out = APAVutils::check_out($out, $pavdata, "_pav_sta");
	APAVutils::check_file("--pheno", $phenodata) if defined($phenodata);

        die "Please install R first\n" if(system("command -v Rscript > /dev/null 2>&1") != 0);
	die "Please install 'APAVplot' R package first\n" if(system("Rscript -e 'library(APAVplot)' > /dev/null 2>&1") != 0);

        my $dir =  dirname(__FILE__);
        my $exec = $dir."/"."vis_pav.R";

	system("Rscript $exec $dir $pavdata $phenodata $out halfviolin '0' '0' '0.1' '0' '0.05' '$fig_width' '$fig_height' '$color' '$add_pheno_info' '$pheno_info_colors' '$x_text_size' '$y_text_size' '$x_title_size' '$y_title_size' 1>/dev/null ");

}


sub pavPlotHist {
	my $usage = "\nUsage: apav pavPlotHist --pav <pav_file> [options]
	
Necessary input description:
  -i, --pav         	<file>        	PAV file produced by command 'apav callpav'.
        
Options:
  -o, --out             <string>        Figure name.
  
  --rm_softcore				Exclude `soft-core` when determining the classification of the region.
  --rm_private				Exclude `private` when determining the classification of the region.
  --softcore_loss_rate	<numeric>	Loss rate.
  --use_binomial                        Use binomial test when determining `soft-core`.
  --softcore_p_value	<numeric>	P-value (binomial test).

Visualization options:
  --fig_width		<numeric>	The width of the figure.
  --fig_height		<numeric>	The height of the figure.

  --hide_ring	 			Hide ring chart.
  --ring_pos_x 		<numeric>	The x-location of the ring chart, ranging from 0 to 1.
  --ring_pos_y 		<numeric>	The y-location of the ring chart, ranging from 0 to 1.
  --ring_r 		<numeric>	The radius of the ring chart, ranging from 0 to 0.5.
  --ring_label_size 	<numeric>	The size of labels on the ring chart.

  --type_colors 	<colors>	Colors for classifications of the region.

  --x_title 		<string>	The text for the x-axis title.
  --x_title_size 	<numeric>	The size of x-axis title.
  --y_title 		<string>	The text for the y-axis title.
  --y_title_size 	<numeric>	The size of y-axis title.
  --x_breaks 		<string>	Break values on the x-axis, separated by commas. (eg: '1,5,10')
  --x_text_size 	<numeric>	The size of tick labels on the x-axis.
  --y_text_size 	<numeric>	The size of tick labels on the y-axis.

  -h, --help		                Print usage.
\n";

	my ($pavdata, $out, $help);

        my $rm_softcore = 0; ##
        my $rm_private = 0; ##
        my $softcore_loss_rate = 0.1;
	my $use_binomial = 0; ## F
        my $softcore_p_value = 0.05;

	my $fig_width = 8;
	my $fig_height = 5;

	my $hide_ring = 0; ## 
	my $ring_pos_x = 0.5;
	my $ring_pos_y = 0.6;
	my $ring_r = 0.3;
	my $ring_label_size = "NA";
	my $type_colors = "NULL";
	my $x_title = "Sample Number";
	my $x_title_size = "NULL";
	my $y_title = "Count";
	my $y_title_size = "NULL";
	my $x_breaks = "NULL";
	my $x_text_size = "NULL";
	my $y_text_size = "NULL";

	GetOptions(
		'pav|i=s'			=> \$pavdata,
                'out|o=s'                       => \$out,

                'rm_softcore!'                 	=> \$rm_softcore,
                'rm_private!'                 	=> \$rm_private,
                'softcore_loss_rate=f'          => \$softcore_loss_rate,
		'use_binomial!'                 => \$use_binomial,
                'softcore_p_value=f'            => \$softcore_p_value,

		'fig_width=f'			=> \$fig_width,
		'fig_height=f'			=> \$fig_height,

		'hide_ring!'			=> \$hide_ring,
		'ring_pos_x=f'			=> \$ring_pos_x,
		'ring_pos_y=f'			=> \$ring_pos_y,
		'ring_r=f'			=> \$ring_r,
		'ring_label_size=f'		=> \$ring_label_size,
		'type_colors=s'			=> \$type_colors,
		'x_title=s'			=> \$x_title,
		'x_title_size=f'		=> \$x_title_size,
		'y_title=s'			=> \$y_title,
		'y_title_size=f'		=> \$y_title_size,
		'x_breaks=s'			=> \$x_breaks,
		'x_text_size=f'			=> \$x_text_size,
		'y_text_size=f'			=> \$y_text_size,
		'help|h!'			=> \$help
	) or die $!."\n";

	die $usage if !defined($pavdata);

        die $usage if $help;
	APAVutils::check_pav_input($pavdata);
	$out = APAVutils::check_out($out, $pavdata, "_pav_hist");

        die "Please install R first\n" if(system("command -v Rscript > /dev/null 2>&1") != 0);
	die "Please install 'APAVplot' R package first\n" if(system("Rscript -e 'library(APAVplot)' > /dev/null 2>&1") != 0);

	my $dir =  dirname(__FILE__);
        my $exec = $dir."/"."vis_pav.R";

	system("Rscript $exec $dir $pavdata 'NULL' $out hist '$rm_softcore' '$rm_private' '$softcore_loss_rate' '$use_binomial' '$softcore_p_value' '$fig_width' '$fig_height' '$hide_ring' '$ring_pos_x' '$ring_pos_y' '$ring_r' '$ring_label_size' '$type_colors' '$x_title' '$x_title_size' '$y_title' '$y_title_size' '$x_breaks' '$x_text_size' '$y_text_size' 1>/dev/null ");

}

sub pavPlotHeat {
	my $usage = "\n Usage: apav pavPlotHeat --pav <pav_file> [options]
	
Necessary input description:
  -i, --pav			<file>		PAV file produced by command 'apav callpav'.
	
Options:
  --pheno			<file>		Phenotype file.
  -o, --out                     <string>        Figure name.

  --rm_softcore                    		Exclude `soft-core` when determining the classification of the region.
  --rm_private                          	Exclude `private` when determining the classification of the region.
  --softcore_loss_rate  	<numeric>       Loss rate.
  --use_binomial                        	Use binomial test when determining `soft-core`.
  --softcore_p_value    	<numeric>       P-value (binomial test).

Visualization options:
  --fig_width           	<numeric>       The width of the figure.
  --fig_height          	<numeric>       The height of the figure.

  --pav_colors			<colors>  	The colors for presence and absence, separated by commas.
  --type_colors 		<colors>	The colors for classifications of the region, separated by commas.
  --region_info_color_list 	<key=value>     The colors for region annotations. (eg: 'chr=black,red')
  --pheno_info_color_list 	<key=value>     The colors for phenotype annotations. (eg: 'gender=blue,green')
  
  --hide_border 				Hide border of blocks in heatmap.
  
  --block_name_size 		<numeric>	The size of block name.
  --block_name_rot		<numeric>	The rotation of block name.
  
  --cluster_rows                                Perform clustering on rows.
  --clustering_distance_rows    <type>          Method of measuring distance when clustring on rows,
                                                choosing from 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski', 'pearson', 'spearman', 'kendall'.
                                                (Default: 'euclidean')
  --clustering_method_rows      <type>          Method to perform hierarchical clustering on rows,
                                                choosing from 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.
                                                (Default: 'complete')
  --row_dend_side               <type>          The position of the row dendrogram, 'left' or 'right'.
  --row_dend_width              <numeric>       The width of the row dendrogram.
  --row_sorted                  <string>        The sorted row names, separated by commas. 
                                                It doesn't work when '--cluster_rows' is used.

  --show_row_names                              Show row names.
  --row_names_side              <type>          The position of row names, 'left' or 'right'.
  --row_names_size              <numeric>       The size of row names.
  --row_names_rot               <numeric>       The rotation of row names.

  --cluster_columns                             Perform clustering on columns.
  --clustering_distance_columns <type>          Method of measuring distance when clustring on columns,
                                                choosing from 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski', 'pearson', 'spearman', 'kendall'.
                                                (Default: 'euclidean')
  --clustering_method_columns   <type>          Method to perform hierarchical clustering on columns,
                                                choosing from 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.
                                                (Default: 'complete')
  --column_dend_side            <type>          The position of the column dendrogram, 'top' or 'bottom'.
  --column_dend_height          <numeric>       The height of the column dendrogram.
  --column_sorted               <string>        The sorted column names, separated by commas. 
                                                It doesn't work when '--cluster_columns'.

  --show_column_names                           Show column names.
  --column_names_side           <type>          The position of column names, 'top' or 'column'.
  --column_names_size           <numeric>       The size of column names.
  --column_names_rot            <numeric>       The rotation of column names.

  --anno_param_row_pheno        <key=value>     Parameters for the phenotype annotation,
                                                including 'show' <type> ('T','F'), 
                                                          'width' <numeric>, 
                                                          'border' <type> ('T', 'F'),
                                                          'name_size' <numeric>, 
                                                          'name_rot <numeric>, 
                                                          'name_side' <type> ('top', 'bottom').
  --anno_param_row_stat         <key=value>     Parameters for the stat annotation of rows,
                                                including 'show' <type> ('T','F'), 
                                                          'width' <numeric>, 
                                                          'border' <type> ('T','F'), 
                                                          'title' <string>, 
                                                          'title_size' <numeric>, 
                                                          'title_rot' <numeric>, 
                                                          'title_side' <type> ('top', 'bottom'), 
                                                          'axis_side' <type> ('top', 'bottom'), 
                                                          'axis_labels_size' <numeric>.
  --anno_param_column_stat      <key=value>     Parameters for the stat annotation of columns,
                                                including 'show' <type> ('T','F'), 
                                                          'height' <numeric>, 
                                                          'border' <type> ('T','F'), 
                                                          'title' <string>, 
                                                          'title_size' <numeric>, 
                                                          'title_rot' <numeric>, 
                                                          'title_side' <type> ('left', 'right'), 
                                                          'axis_side' <type> ('left', 'right'), 
                                                          'axis_labels_size' <numeric>.

  --legend_side                 <type>          The position of legend, choosing from 'top', 'bottom', 'right', 'left'.
  --legend_title                <string>        The text for the legend title.
  --legend_title_size           <numeric>       The size of legend title.
  --legend_text_size            <numeric>       The size of legend item labels.
  --legend_grid_size            <numeric>       The size of legend grid.

  --use_raster                                  Render the heatmap body as a raster image.

  -h, --help                                    Print usage page.

Warning: --region_info_color_list --pheno_info_color_list --anno_param_row_phen --anno_param_column_region --anno_param_row_stat --anno_param_column_stat can be added multiple 
times.
  eg: --anno_param_row_pheno show=T --anno_param_row_pheno width=5 --anno_param_row_pheno name_rot=90
  
  \n";

  	my ($pavdata, $out, $help);

	my $phenodata = "NULL";

	my $rm_softcore = 0; ## 
	my $rm_private = 0; ## 
	my $softcore_loss_rate = 0.1;
	my $use_binomial = 0; ## F
	my $softcore_p_value = 0.05;

	my $fig_width = 10;
        my $fig_height = 6;

	my $region_type = "NULL";
	my $pav_colors = "#5680ae,gray70";
	my $type_colors = "NULL";
        my (%region_info_color_list, $region_info_color_list_str);
        my (%pheno_info_color_list, $pheno_info_color_list_str);

        my $hide_border = 0; ## 
	my $block_name_size = "NULL";
	my $block_name_rot = 0;

        my $cluster_rows = 0; ## False
        my $clustering_distance_rows = "euclidean";
        my $clustering_method_rows = "complete";
        my $row_dend_side = "left";
        my $row_dend_width = 5;
        my $row_sorted = "c()";

        my $show_row_names = 0; ## False
        my $row_names_side = "left";
        my $row_names_size = 10;
        my $row_names_rot = 0;

        my $cluster_columns = 0; ## False
        my $clustering_distance_columns = "euclidean";
        my $clustering_method_columns = "complete";
        my $column_dend_side = "top";
        my $column_dend_height = 5;
        my $column_sorted = "c()";

        my $show_column_names = 0; ## False
        my $column_names_side = "bottom";
        my $column_names_size = 10;
        my $column_names_rot = 90;

        my %anno_param_row_pheno;
        $anno_param_row_pheno{"show"} = "T";
        $anno_param_row_pheno{"width"} = 5;
        $anno_param_row_pheno{"border"} = "F";
        $anno_param_row_pheno{"name_size"} = "NULL";
        $anno_param_row_pheno{"name_rot"} = 90;
        $anno_param_row_pheno{"name_side"} = "top";

	my %anno_param_column_region;
	$anno_param_column_region{"show"} = "T";
	$anno_param_column_region{"height"} = 5;
	$anno_param_column_region{"border"} = "F";
	$anno_param_column_region{"name_size"} = "NULL";
	$anno_param_column_region{"name_rot"} = 0;
	$anno_param_column_region{"name_side"} = "right";

	my %anno_param_row_stat;
        $anno_param_row_stat{"show"} = "T";
        $anno_param_row_stat{"width"} = 10;
        $anno_param_row_stat{"border"} = "F";
        $anno_param_row_stat{"title"} = "Presence\nNumber";
        $anno_param_row_stat{"title_size"} = 10;
        $anno_param_row_stat{"title_side"} = "bottom";
        $anno_param_row_stat{"title_rot"} = 0;
        $anno_param_row_stat{"axis_side"} = "bottom";
        $anno_param_row_stat{"axis_labels_size"} = 8;

        my %anno_param_column_stat;
        $anno_param_column_stat{"show"} = "T";
        $anno_param_column_stat{"height"} = 10;
        $anno_param_column_stat{"border"} = "F";
        $anno_param_column_stat{"title"} = "Presence\nNumber";
        $anno_param_column_stat{"title_size"} = 10;
        $anno_param_column_stat{"title_side"} = "left";
        $anno_param_column_stat{"title_rot"} = 0;
        $anno_param_column_stat{"axis_side"} = "left";
        $anno_param_column_stat{"axis_labels_size"} = 8;

        my $legend_side = "right";
        my $legend_title = "PAV,Region";
        my $legend_title_size = "NULL";
        my $legend_text_size = "NULL";
        my $legend_grid_size = 4;

        my $use_raster = "NULL";
	
	GetOptions(
                'pav|i=s'                       => \$pavdata,
                'pheno=s'                     => \$phenodata,
                'out|o=s'                       => \$out,

		'rm_softcore!'			=> \$rm_softcore,
		'rm_private!'			=> \$rm_private,
		'softcore_loss_rate=f'		=> \$softcore_loss_rate,
		'use_binomial!'                 => \$use_binomial,
		'softcore_p_value=f'		=> \$softcore_p_value,

		'fig_width=f'                   => \$fig_width,
                'fig_height=f'                  => \$fig_height,

		'region_type=s'                 => \$region_type,

                'pav_colors=s'                  => \$pav_colors,
		'type_colors=s'			=> \$type_colors,
                'region_info_color_list=s'      => \%region_info_color_list,
                'pheno_info_color_list'         => \%pheno_info_color_list,

                'hide_border=s'                 => \$hide_border,

		'block_name_size=f'		=> \$block_name_size,
		'block_name_rot=f'		=> \$block_name_rot,

                'cluster_rows!'                 => \$cluster_rows,
                'clustering_distance_rows=s'    => \$clustering_distance_rows,
                'clustering_method_rows=s'      => \$clustering_method_rows,
                'row_dend_side=s'               => \$row_dend_side,
                'row_dend_width=f'              => \$row_dend_width,
                'row_sorted=s'                  => \$row_sorted,

                'show_row_names!'               => \$show_row_names,
                'row_names_side=s'              => \$row_names_side,
                'row_names_size=f'              => \$row_names_size,
                'row_names_rot=f'               => \$row_names_rot,

                'cluster_columns!'              => \$cluster_columns,
                'clustering_distance_columns=s' => \$clustering_distance_columns,
                'clustering_method_columns=s'   => \$clustering_method_columns,
                'column_dend_side=s'            => \$column_dend_side,
                'column_dend_height=f'          => \$column_dend_height,
                'column_sorted=s'               => \$column_sorted,

                'show_column_names!'            => \$show_column_names,
                'column_names_side=s'           => \$column_names_side,
                'column_names_size=f'           => \$column_names_size,
                'column_names_rot=f'            => \$column_names_rot,

                'anno_param_row_pheno=s'        => \%anno_param_row_pheno,
		'anno_param_column_region=s'    => \%anno_param_column_region,
                'anno_param_row_stat=s'         => \%anno_param_row_stat,
                'anno_param_column_stat=s'      => \%anno_param_column_stat,

                'legend_side=s'                 => \$legend_side,
                'legend_title=s'                => \$legend_title,
                'legend_title_size=f'           => \$legend_title_size,
                'legend_text_size=f'            => \$legend_text_size,
                'legend_grid_size=f'            => \$legend_grid_size,

                'use_raster!'                   => \$use_raster,

                'help|h!'       => \$help
        ) or die $!."\n";	

	die $usage if !defined($pavdata);
	die $usage if $help;

	APAVutils::check_pav_input($pavdata);
	$out = APAVutils::check_out($out, $pavdata, "_pav_heatmap");
	APAVutils::check_file("--pheno", $phenodata) if defined($phenodata);

	die "Please install R first\n" if(system("command -v Rscript > /dev/null 2>&1") != 0);
	die "Please install 'APAVplot' R package first\n" if(system("Rscript -e 'library(APAVplot)' > /dev/null 2>&1") != 0);

        if(%region_info_color_list){
                for my $key (keys %region_info_color_list){
                        $region_info_color_list_str.= $key."=".$region_info_color_list{$key}.";";
                }
                $region_info_color_list_str =~ s/;$//g;
        } else {
                $region_info_color_list_str = "NULL";
        }

        if(%pheno_info_color_list){
                for my $key (keys %pheno_info_color_list){
                        $pheno_info_color_list_str.= $key."=".$pheno_info_color_list{$key}.";";
                }
                $pheno_info_color_list_str =~ s/;$//g;
        } else {
                $pheno_info_color_list_str = "NULL";
        }

        my $dir =  dirname(__FILE__);
        my $exec = $dir."/"."vis_pav.R";

	system("Rscript $exec $dir $pavdata $phenodata $out heatmap '$rm_softcore' '$rm_private' '$softcore_loss_rate' '$use_binomial' '$softcore_p_value' '$fig_width' '$fig_height' '$region_type'  '$pav_colors' '$type_colors' '$region_info_color_list_str' '$pheno_info_color_list_str' '$hide_border' '$block_name_size' '$block_name_rot' '$cluster_rows' '$clustering_distance_rows' '$clustering_method_rows' '$row_dend_side' '$row_dend_width' '$row_sorted' '$show_row_names' '$row_names_side' '$row_names_size' '$row_names_rot' '$cluster_columns' '$clustering_distance_columns' '$clustering_method_columns' '$column_dend_side' '$column_dend_height' '$column_sorted' '$show_column_names' '$column_names_side' '$column_names_size' '$column_names_rot' '$anno_param_row_pheno{'show'}' '$anno_param_row_pheno{'width'}' '$anno_param_row_pheno{'border'}' '$anno_param_row_pheno{'name_size'}' '$anno_param_row_pheno{'name_rot'}' '$anno_param_row_pheno{'name_side'}' '$anno_param_column_region{'show'}' '$anno_param_column_region{'height'}' '$anno_param_column_region{'border'}' '$anno_param_column_region{'name_size'}' '$anno_param_column_region{'name_rot'}' '$anno_param_column_region{'name_side'}' '$anno_param_row_stat{'show'}' '$anno_param_row_stat{'width'}' '$anno_param_row_stat{'border'}' '$anno_param_row_stat{'title'}' '$anno_param_row_stat{'title_size'}' '$anno_param_row_stat{'title_side'}' '$anno_param_row_stat{'title_rot'}' '$anno_param_row_stat{'axis_side'}'  '$anno_param_row_stat{'axis_labels_size'}' '$anno_param_column_stat{'show'}' '$anno_param_column_stat{'height'}' '$anno_param_column_stat{'border'}' '$anno_param_column_stat{'title'}' '$anno_param_column_stat{'title_size'}' '$anno_param_column_stat{'title_side'}' '$anno_param_column_stat{'title_rot'}' '$anno_param_column_stat{'axis_side'}' '$anno_param_column_stat{'axis_labels_size'}' '$legend_side' '$legend_title' '$legend_title_size' '$legend_text_size' '$legend_grid_size' '$use_raster'  1>/dev/null  ");


}

1;
