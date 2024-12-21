#!/usr/bin/perl

package APAVvisPheno;

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename;

sub phenosta {
	my $usage = "\nUsage: apav pavStaPheno --pav <pav_file> --pheno <phenotype_file> [options]

Necessary input description:
  -i, --pav         	<file>        	PAV file produced by command 'apav callPAV'.
  --pheno       	<file>        	Phenotype file.
        
Options:
  -o, --out             <string>        Output file.

  --p_adjust_method	<type>		The adjustment methods,
  					choosing from 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'.
					(Default: fdr)
  --parallel_n		<int>		The number of CPU cores used for parallel computing.

  -h, --help                    	Print usage. 
  \n";

   	my ($pavdata, $phenodata, $out, $help);

	my $p_adjust_method = "fdr";
	my $parallel_n = 0;

	GetOptions(
                'pav|i=s'                       => \$pavdata,
                'pheno=s'                       => \$phenodata,
                'out|o=s'                       => \$out,

                'p_adjust_method=s'             => \$p_adjust_method,
                'parallel_n=i'                  => \$parallel_n,

                'help|h!'                       => \$help
        ) or die $!."\n";


	die $usage if !defined($pavdata) & !defined($phenodata);

        die $usage if $help;
	APAVutils::check_pav_input($pavdata);
	APAVutils::check_file('--pheno|-p', $phenodata);
	$out = APAVutils::check_out($out, $pavdata, ".phenores");

        die "Please install R first\n" if(system("command -v Rscript > /dev/null 2>&1") != 0);
	die "Please install 'APAVplot' R package first\n" if(system("Rscript -e 'library(APAVplot)' > /dev/null 2>&1") != 0);

        my $dir =  dirname(__FILE__);
        my $exec = $dir."/"."vis_pheno.R";

        system("Rscript $exec $dir $pavdata $phenodata $out stat \"$p_adjust_method\" \"$parallel_n\"  1>/dev/null ");
	
}

sub phenoPlotHeat{

	my $usage = "\nUsage: apav pavPlotPhenoHeat --pav <pav_file> --pheno_res <phenotype_association_result> [options]

Necessary input description:
  --pav         		<file>		PAV file produced by command 'apav callPAV'.
  --pheno_res   		<file>        	Phenotype association result produced by command 'apav pavStaPheno'.
        
Options:
  -o, --out                     <string>        Figure name. 

  --p_threshold			<numeric>	The threshold of p_value/p_adjusted.
  						(Default: 0.01)
  --adjust_p                    		Adjust p_value.

Visualization options:
  --fig_width                   <numeric>       The width of the figure.
  --fig_height                  <numeric>       The height of the figure.

  --only_show_significant 			Only show p_value/p_adjusted that satisfies the condition.
  --flip 					Flip the cartesian coordinates.
  
  --p_colors 			<colors>	The colors for p_value/p_adjusted, separated by commas.
  --na_col 			<color>		The color for NA values.
  --cell_border_color 		<color>		The color for the border of cells.
  --region_info_color_list 	<colors>	The colors for region annotations. (eg: 'chr=black,red')

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

  --hide_row_names                              Hide row names.
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

  --hide_column_names                           hide column names.
  --column_names_side           <type>          The position of column names, 'top' or 'column'.
  --column_names_size           <numeric>       The size of column names.
  --column_names_rot            <numeric>       The rotation of column names.
  
  --anno_param_region    	<key=value>     Parameters for the region annotation, 
                                                including 'show' <type> ('T', 'F'), 
                                                          'width' <numeric>, 
                                                          'border' <type> ('T', 'F'), 
                                                          'name_size' <numeric>, 
                                                          'name_rot' <numeric>, 
                                                          'name_side' <type> ('left', 'right').
  
  --legend_side                 <type>          The position of legend, choosing from 'top', 'bottom', 'right', 'left'.
  --legend_title                <string>        The text for the legend title.
  --legend_title_size           <numeric>       The size of legend title.
  --legend_text_size            <numeric>       The size of legend item labels.
  --legend_grid_size            <numeric>       The size of legend grid.

  -h, --help                    		Print usage.
  \n";

  	my ($pavdata, $phenores, $out, $help);

	my $fig_width = 4;
        my $fig_height = 10;

	my $p_threshold = 0.01;
        my $adjust_p = 0; ## 
        my $only_show_significant = 0; ## 
        my $flip = 0; ## F

	my $p_colors = "#B95758,#f0d2d0";
	my $na_col = "gray";
	my $cell_border_color = "white";
	my (%region_info_color_list, $region_info_color_list_str);

	my $cluster_rows = 0; ## F
	my $clustering_distance_rows = "euclidean";
	my $clustering_method_rows = "complete";
	my $row_dend_side = "left";
	my $row_dend_width = 5;
	my $row_sorted = "c()";

	my $hide_row_names = 0; ## F
	my $row_names_side = "right";
	my $row_names_size = 10;
	my $row_names_rot = 0;

	my $cluster_columns = 0; ## F
	my $clustering_distance_columns = "euclidean";
	my $clustering_method_columns = "complete";
	my $column_dend_side = "top";
	my $column_dend_height = 5;
	my $column_sorted = "c()";
	
	my $hide_column_names = 0; ## F
	my $column_names_side = "bottom";
	my $column_names_size = 10;
	my $column_names_rot = 90;

	my %anno_param_region;
	$anno_param_region{"show"} = "T"; ## T
	$anno_param_region{"width"} = 5;
	$anno_param_region{"border"} = "F"; ##F
	$anno_param_region{"name_size"} = "NULL";
	$anno_param_region{"name_rot"} = 90;
	$anno_param_region{"name_side"} = "bottom";

	my $legend_side = "right";
	my $legend_title = "PAV,region";
	my $legend_title_size = "NULL";
	my $legend_text_size = "NULL";
	my $legend_grid_size = 4;


        GetOptions(
                'pav|i=s'                       => \$pavdata,
                'pheno_res=s'                   => \$phenores,
                'out|o=s'                       => \$out,

		'fig_width=f'                   => \$fig_width,
                'fig_height=f'                  => \$fig_height,

		'p_threshold=f'			=> \$p_threshold,
		'adjust_p!'			=> \$adjust_p,
		'only_show_significant!'	=> \$only_show_significant,
		'flip!'				=> \$flip,
		
		'p_colors=s'			=> \$p_colors,
		'na_col=s'			=> \$na_col,
		'cell_border_color=s'		=> \$cell_border_color,
		'region_info_color_list=s'	=> \%region_info_color_list,

		'cluster_rows!'			=> \$cluster_rows,
		'clustering_distance_rows=s'	=> \$clustering_distance_rows,
		'clustering_method_rows=s'	=> \$clustering_method_rows,
		'row_dend_side=s'		=> \$row_dend_side,
		'row_dend_width=f'		=> \$row_dend_width,
		'row_sorted=s'			=> \$row_sorted,

		'hide_row_names!'		=> \$hide_row_names,
		'row_names_side=s'		=> \$row_names_side,
		'row_names_size=f'		=> \$row_names_size,
		'row_names_rot=f'		=> \$row_names_rot,

		'cluster_columns!'		=> \$cluster_columns,
		'clustering_distance_columns=s'	=> \$clustering_distance_columns,
		'clustering_method_columns=s'	=> \$clustering_method_columns,
		'column_dend_side=s'		=> \$column_dend_side,
		'column_dend_height=f'		=> \$column_dend_height,
		'column_sorted=s'		=> \$column_sorted,

		'hide_column_names!'		=> \$hide_column_names,
		'column_names_side=s'		=> \$column_names_side,
		'column_names_size=f'		=> \$column_names_size,
		'column_names_rot=f'		=> \$column_names_rot,

		'anno_param_region=s'		=> \%anno_param_region,

		'legend_side=s'			=> \$legend_side,
		'legend_title=s'		=> \$legend_title,
		'legend_title_size=f'		=> \$legend_title_size,
		'legend_text_size=f'		=> \$legend_text_size,
		'legend_grid_size=f'		=> \$legend_grid_size,

                'help|h!'                       => \$help
        ) or die $!."\n";

	die $usage if !defined($pavdata) & !defined($phenores) & !defined($out);

        die $usage if $help;

	APAVutils::check_pav_input($pavdata);
	APAVutils::check_file('--pheno_res', $phenores);
	$out = APAVutils::check_out($out, $phenores, "_pheno_heatmap");

	my @pheno_filter;
	if($adjust_p){
		@pheno_filter = `awk '\$4 < $p_threshold' $phenores`;
	}else{
		@pheno_filter = `awk '\$3 < $p_threshold' $phenores`;
	}
	if($#pheno_filter < 0){
		die "No results after filtering at the current threshold.\n";
	}

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

        my $dir =  dirname(__FILE__);
        my $exec = $dir."/"."vis_pheno.R";

	system("Rscript $exec $dir $pavdata 'NULL' '$out' heatmap '$phenores' '$fig_width' '$fig_height' '$p_threshold' '$adjust_p' '$only_show_significant' '$flip' '$p_colors' '$na_col' '$cell_border_color' '$region_info_color_list_str' '$cluster_rows' '$clustering_distance_rows' '$clustering_method_rows' '$row_dend_side' '$row_dend_width' '$row_sorted' '$hide_row_names' '$row_names_side' '$row_names_size' '$row_names_rot' '$cluster_columns' '$clustering_distance_columns' '$clustering_method_columns' '$column_dend_side' '$column_dend_height' '$column_sorted' '$hide_column_names' '$column_names_side' '$column_names_size' '$column_names_rot' '$anno_param_region{'show'}' '$anno_param_region{'width'}' '$anno_param_region{'border'}' '$anno_param_region{'name_size'}' '$anno_param_region{'name_rot'}' '$anno_param_region{'name_side'}' '$legend_side' '$legend_title' '$legend_title_size' '$legend_text_size' '$legend_grid_size' 1>/dev/null ");

}


sub phenoPlotBlock{
	my $usage = "\nUsage: apav pavPlotPhenoBlock --pav <pav_file> --pheno <phenotype file> --pheno_res <phenotype_association_result> --phen_name <phenotype name> [options]

Necessary input description:
  -i, --pav         		<string>        PAV file produced by command 'apav callPAV'.
  --pheno           		<file>          Phenotype file.
  --pheno_res   		<file>	        Phenotype association result produced by command 'apav pavStaPheno'.
  --pheno_name                  <string>	Phenotype name.	

Options:
  -o, --out                     <string>        Figure name. 

  --p_threshold                 <numeric>       The threshold of p_value/p_adjusted.
                                                (Default: 0.01)
  --adjust_p                                    Adjust p_value.

Visualization options:
  --fig_width                   <numeric>       The width of the figure.
  --fig_height                  <numeric>       The height of the figure.

  --only_show_significant                       Only show p_value/p_adjusted that satisfies the condition.
  --flip                                        Flip the cartesian coordinates.

  --per_colors 			<colors>	The colors for absence percentage, separated by commas.
  --na_col                      <color>         The color for NA values.
  --cell_border_color           <color>         The color for the border of cells.
  --region_info_color_list      <colors>        The colors for region annotations. (eg: 'chr=black,red')

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

  --hide_row_names                              Hide row names.
  --row_names_side              <type>          The position of row names, 'left' or 'right'.
  --row_names_size              <numeric>       The size of row names.
  --row_names_rot               <numeric>       The rotation of row names.

  --cluster_columns                             Perform clustering on columns.
  --clustering_distance_columns <type>          Method of measuring distance when clustring on columns.
                                                choosing from 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski', 'pearson', 'spearman', 'kendall'.
                                                (Default: 'euclidean')
  --clustering_method_columns   <type>          Method to perform hierarchical clustering on columns.
                                                choosing from 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.
                                                (Default: 'complete')
  --column_dend_side            <type>          The position of the column dendrogram, 'top' or 'bottom'.
  --column_dend_height          <numeric>       The height of the column dendrogram.
  --column_sorted               <string>        The sorted column names, separated by commas.
                                                It doesn't work when '--cluster_columns'.

  --hide_column_names                           Hide column names.
  --column_names_side           <type>          The position of column names, 'top' or 'column'.
  --column_names_size           <numeric>       The size of column names.
  --column_names_rot            <numeric>       The rotation of column names.

  --anno_param_region           <key=value>     Parameters for the region annotation, 
                                                including 'show' <type> ('T', 'F'), 
                                                          'width' <numeric>, 
                                                          'border' <type> ('T', 'F'), 
                                                          'name_size' <numeric>, 
                                                          'name_rot' <numeric>, 
                                                          'name_side' <type> ('left', 'right').
  
  --legend_side                 <type>          The position of legend, choosing from 'top', 'bottom', 'right', 'left'.
  --legend_title                <string>        The text for the legend title.
  --legend_title_size           <numeric>       The size of legend title.
  --legend_text_size            <numeric>       The size of legend item labels.
  --legend_grid_size            <numeric>       The size of legend grid.

  -h, --help                                    Print usage.	
	\n";

	my ($pavdata, $phenores, $phenodata, $pheno_name, $out, $help);

	my $fig_width = 10;
        my $fig_height = 6;

	my $p_threshold = 0.01;
        my $adjust_p = 0; ## 
        my $only_show_significant = 0; ## T
        my $flip = 0; ## F

        my $per_colors = "#d6deeb,#376da1";
        my $na_col = "gray";
        my $cell_border_color = "white";
        my (%region_info_color_list, $region_info_color_list_str);

        my $cluster_rows = 0; ## F
        my $clustering_distance_rows = "euclidean";
        my $clustering_method_rows = "complete";
        my $row_dend_side = "left";
        my $row_dend_width = 5;
        my $row_sorted = "c()";

        my $hide_row_names = 0; ## F
        my $row_names_side = "right";
        my $row_names_size = 10;
        my $row_names_rot = 0;

        my $cluster_columns = 0; ## F
        my $clustering_distance_columns = "euclidean";
        my $clustering_method_columns = "complete";
        my $column_dend_side = "top";
        my $column_dend_height = 5;
        my $column_sorted = "c()";
        
        my $hide_column_names = 0; ## F
        my $column_names_side = "bottom";
        my $column_names_size = 10;
        my $column_names_rot = 90;

        my %anno_param_region;
        $anno_param_region{"show"} = "T";
        $anno_param_region{"width"} = 5;
        $anno_param_region{"border"} = "F";
        $anno_param_region{"name_size"} = "NULL";
        $anno_param_region{"name_rot"} = 90;
        $anno_param_region{"name_side"} = "bottom";

	my $legend_side = "right";
        my $legend_title = "PAV,region";
        my $legend_title_size = "NULL";
        my $legend_text_size = "NULL";
        my $legend_grid_size = 4;	


	GetOptions(
                'pav|i=s'                       => \$pavdata,
		'pheno=s'			=> \$phenodata,
                'pheno_res=s'                   => \$phenores,
		'pheno_name=s'			=> \$pheno_name,
                'out|o=s'                       => \$out,

		'fig_width=f'                   => \$fig_width,
                'fig_height=f'                  => \$fig_height,

                'p_threshold=f'                 => \$p_threshold,
                'adjust_p!'                     => \$adjust_p,
                'only_show_significant!'        => \$only_show_significant,
                'flip!'                         => \$flip,

		'per_colors=s'			=> \$per_colors,
		'na_col=s'			=> \$na_col,
		'cell_border_color=s'		=> \$cell_border_color,
		'region_info_color_list=s'	=> \%region_info_color_list,

		'cluster_rows!'                 => \$cluster_rows,
                'clustering_distance_rows=s'    => \$clustering_distance_rows,
                'clustering_method_rows=s'      => \$clustering_method_rows,
                'row_dend_side=s'               => \$row_dend_side,
                'row_dend_width=f'              => \$row_dend_width,
                'row_sorted=s'                  => \$row_sorted,

                'hide_row_names!'               => \$hide_row_names,
                'row_names_side=s'              => \$row_names_side,
                'row_names_size=f'              => \$row_names_size,
                'row_names_rot=f'               => \$row_names_rot,

                'cluster_columns!'              => \$cluster_columns,
                'clustering_distance_columns=s' => \$clustering_distance_columns,
                'clustering_method_columns=s'   => \$clustering_method_columns,
                'column_dend_side=s'            => \$column_dend_side,
                'column_dend_height=f'          => \$column_dend_height,
                'column_sorted=s'               => \$column_sorted,

                'hide_column_names!'            => \$hide_column_names,
                'column_names_side=s'           => \$column_names_side,
                'column_names_size=f'           => \$column_names_size,
                'column_names_rot=f'            => \$column_names_rot,

                'anno_param_region=s'           => \%anno_param_region,

                'legend_side=s'                 => \$legend_side,
                'legend_title=s'                => \$legend_title,
                'legend_title_size=f'           => \$legend_title_size,
                'legend_text_size=f'            => \$legend_text_size,
                'legend_grid_size=f'            => \$legend_grid_size,

                'help|h!'                       => \$help
	) or die $!."\n";

	die $usage if !defined($pavdata) & !($phenodata) & !defined($phenores)  & !defined($pheno_name);

        die $usage if $help;
	APAVutils::check_pav_input($pavdata);
	APAVutils::check_file('--pheno', $phenodata);
        APAVutils::check_file('--pheno_res', $phenores);
	APAVutils::check_arg('--pheno_name', $pheno_name);
	$out = APAVutils::check_out($out, $phenores, "_pheno_".$pheno_name."_block");

	my @pheno_filter;
	if($adjust_p){
                @pheno_filter = `awk '\$4 < $p_threshold' $phenores`;
        }else{
                @pheno_filter = `awk '\$3 < $p_threshold' $phenores`;
        }
        if($#pheno_filter < 0){
                die "No results after filtering at the current threshold.\n";
        }

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

        my $dir =  dirname(__FILE__);
        my $exec = $dir."/"."vis_pheno.R";

	system("Rscript $exec $dir $pavdata $phenodata '$out' block '$phenores' '$fig_width' '$fig_height' '$pheno_name' '$p_threshold' '$adjust_p' '$only_show_significant' '$flip' '$per_colors' '$na_col' '$cell_border_color' '$region_info_color_list_str' '$cluster_rows' '$clustering_distance_rows' '$clustering_method_rows' '$row_dend_side' '$row_dend_width' '$row_sorted' '$hide_row_names' '$row_names_side' '$row_names_size' '$row_names_rot' '$cluster_columns' '$clustering_distance_columns' '$clustering_method_columns' '$column_dend_side' '$column_dend_height' '$column_sorted' '$hide_column_names' '$column_names_side' '$column_names_size' '$column_names_rot' '$anno_param_region{'show'}' '$anno_param_region{'width'}' '$anno_param_region{'border'}' '$anno_param_region{'name_size'}' '$anno_param_region{'name_rot'}' '$anno_param_region{'name_side'}' '$legend_side' '$legend_title' '$legend_title_size' '$legend_text_size' '$legend_grid_size' 1>/dev/null ");

}



sub phenPlotManhattan{
	
	my $usage = "\nUsage: apav pavPlotPhenoMan --pav <pav_file> --pheno_res <phenotype_association_result> --pheno_name <phenotype> [options]
	
Necessary input description:
  -i, --pav         	<file>        	PAV file produced by command 'apav callPAV'.
  --pheno               <file>          Phenotype file.
  
  --pheno_res   	<file>        	Phenotype association result produced by command 'apav pavStaPheno'.
  --pheno_name          <string>	Phenotype name.

Options:
  -o, --out             <string>        Figure name.
  --adjust_p                            Adjust p_value.
  --highlight_top_n     <int>           The top `n` points will be highlighted.

Visualization Settings:
  --fig_width           <numeric>       The width of the figure.
  --fig_height          <numeric>       The height of the figure.

  --highlight_text_size <numeric>	The size of labels on highlight points.

  --point_size 		<numeric>	The size of points.
  --x_text_size 	<numeric>	The size of tick labels on x-axis.
  --x_text_angle	<numeric>	The angle of tick labels.
  --x_title_size 	<numeric>	The size of x-axis title.
  --y_text_size 	<numeric>	The size of tick labels on y-axis.
  --y_title_size 	<numeric>	The size of y-axis title.

  -h, --help				Print usage.  
	\n";


	my ($pavdata, $phenodata, $phenores, $pheno_name, $out, $help);

	my $fig_width = 10;
        my $fig_height = 6;

	my $adjust_p = 0; ## 
	my $highlight_top_n = 5;
	my $highlight_text_size = 4;
	my $point_size = 1.5;
	my $x_text_size = "NULL";
	my $x_text_angle = 0;
	my $x_title_size = "NULL";
	my $y_text_size = "NULL";
	my $y_title_size = "NULL";


	GetOptions(
		'pav|i=s'                       => \$pavdata,
		'pheno=s'			=> \$phenodata,
		'pheno_name=s'                  => \$pheno_name,
                'pheno_res=s'                   => \$phenores,
                'out|o=s'                       => \$out,

		'fig_width=f'           	=> \$fig_width,
                'fig_height=f'          	=> \$fig_height,

		'adjust_p!'			=> \$adjust_p,
		'highlight_top_n=i'		=> \$highlight_top_n,
		'highlight_text_size=f'		=> \$highlight_text_size,

		'point_size=f'			=> \$point_size,
		'x_text_size=f'			=> \$x_text_size,
		'x_text_angle=f'		=> \$x_text_angle,
		'x_title_size=f'		=> \$x_title_size,
		'y_text_size=f'			=> \$y_text_size,
		'y_title_size=f'		=> \$y_title_size,
		'help|h!'			=> \$help
	) or die $!."\n";

	die $usage if !defined($pavdata) & !defined($phenores) & !defined($pheno_name);
        die $usage if $help;

	APAVutils::check_pav_input($pavdata);
	APAVutils::check_file('--pheno', $phenodata);
        APAVutils::check_file('--pheno_res', $phenores);
	APAVutils::check_arg('--pheno_name', $pheno_name);
	$out = APAVutils::check_out($out, $phenores, "_pheno_".$pheno_name."_manhattan");

        die "Please install R first\n" if(system("command -v Rscript > /dev/null 2>&1") != 0);
	die "Please install 'APAVplot' R package first\n" if(system("Rscript -e 'library(APAVplot)' > /dev/null 2>&1") != 0);

        my $dir =  dirname(__FILE__);
        my $exec = $dir."/"."vis_pheno.R";

        system("Rscript $exec $dir $pavdata $phenodata '$out' manhattan $phenores \"$pheno_name\" '$fig_width' '$fig_height' \"$adjust_p\" \"$highlight_top_n\" \"$highlight_text_size\" \"$point_size\" \"$x_text_size\" \"$x_text_angle\" \"$x_title_size\" \"$y_text_size\" \"$y_title_size\" 1>/dev/null ");


}


sub phenoPlotBar {
	my $usage = "\nUsage: apav pavPlotPhenoBar --pav <pav_file> --pheno <phenotype_file> --region_name <region> --pheno_name <phenotype> [options]

Necessary input description:
  -i, --pav         	<file>        	PAV file produced by command 'apav callPAV'.
  --pheno		<file>		Phenotype file.

  --pheno_name 		<string>	The name of phenotype.
  --region_name 	<string>	The name of target region.

Option:
  -o, --out		<string>	Figure name.
 
Visualization options:
  --fig_width           <numeric>       The width of the figure.
  --fig_height          <numeric>       The height of the figure.

  --pav_colors 		<colors>	The colors for presence and absence.
  --bar_width 		<numeric>	The relative width of bars, ranging from 0 to 1.

  --x_text_size 	<numeric>	The size of tick labels on x-axis.
  --x_title_size 	<numeric>	The size of x-axis title.
  --y_text_size 	<numeric>	The size of tick labels on y-axis.
  --y_title_size 	<numeric>	The size of y-axis title.

  --legend_side 	<type>		The position of legend, choosing from 'top', 'bottom', 'right', 'left'.
  --legend_title_size 	<numeric>	The size of legend title.
  --legend_text_size 	<numeric>	The size of legend item labels.
	
  -h, --help				Print usage.
	\n";

	my ($pavdata, $phenodata, $out, $help);
	my ($pheno_name, $region_name);

	my $fig_width = 8;
        my $fig_height = 6;

	my $pav_colors = "gray70,steelblue";
	my $bar_width = 0.8;

	my $x_text_size = "NULL";
	my $x_title_size = "NULL";
	my $y_text_size = "NULL";
	my $y_title_size = "NULL";

	my $legend_side = "top";
	my $legend_title_size = "NULL";
	my $legend_text_size = "NULL";


	GetOptions(
		'pav|i=s'		=> \$pavdata,
		'pheno=s'		=> \$phenodata,
		'out|o=s'		=> \$out,

		'fig_width=f'           => \$fig_width,
                'fig_height=f'          => \$fig_height,

		'pheno_name=s'		=> \$pheno_name,
		'region_name=s'		=> \$region_name,

		'pav_colors=s'		=> \$pav_colors,
		'bar_width=f'		=> \$bar_width,

		'x_text_size=f'		=> \$x_text_size,
		'x_title_size=f'	=> \$x_title_size,
		'y_text_size=f'		=> \$y_text_size,
		'y_title_size=f'	=> \$y_title_size,

		'legend_side=s'		=> \$legend_side,
		'legend_title_size=f'	=> \$legend_title_size,
		'legend_text_size=f'	=> \$legend_text_size,
		'help|h!'		=> \$help
	) or die $!."\n";


	die $usage if !defined($pavdata) & !defined($phenodata) & !defined($pheno_name) & !defined($region_name);
        die $usage if $help;

        APAVutils::check_pav_input($pavdata);
	APAVutils::check_file('--pheno', $phenodata);
	APAVutils::check_arg('--pheno_name', $pheno_name);
	APAVutils::check_arg('--region_name', $region_name);
	$out = APAVutils::check_out($out, $pavdata, "_pheno_".$pheno_name."_".$region_name."_bar");

	die "Please install R first\n" if(system("command -v Rscript > /dev/null 2>&1") != 0);
	die "Please install 'APAVplot' R package first\n" if(system("Rscript -e 'library(APAVplot)' > /dev/null 2>&1") != 0);
        
	my $dir =  dirname(__FILE__);
        my $exec = $dir."/"."vis_pheno.R";

	system("Rscript $exec $dir $pavdata $phenodata '$out' bar '$pheno_name' '$region_name' '$fig_width' '$fig_height' '$pav_colors' '$bar_width' '$x_text_size' '$x_title_size' '$y_text_size' '$y_title_size' '$legend_side' '$legend_title_size' '$legend_text_size' 1>/dev/null ");

}

sub phenoPlotViolin {
	
	my $usage = "\nUsage: apav pavPlotPhenoVio --pav <pav_file> --pheno <phenotype_file> --region_name <region> --pheno_name <phenotype> [options]

Necessary input description:
  -i, --pav             <file>          PAV file produced by command 'apav callPAV'.
  --pheno               <file>          Phenotype file.

  --pheno_name          <string>        The name of phenotype.
  --region_name         <string>        The name of target region.

Option:
  -o, --out             <string>        Figure name.

Visualization Settings:
  --fig_width           <numeric>       The width of the figure.
  --fig_height          <numeric>       The height of the figure.
  
  --pav_colors          <colors>        The colors for presence and absence.

  --x_text_size         <numeric>       The size of tick labels on x-axis.
  --x_title_size        <numeric>       The size of x-axis title.
  --y_text_size         <numeric>       The size of tick labels on y-axis.
  --y_title_size        <numeric>       The size of y-axis title.

  --legend_side         <type>          The position of legend, choosing from 'top', 'bottom', 'right', 'left'.
  --legend_title_size   <numeric>       The size of legend title.
  --legend_text_size    <numeric>       The size of legend item labels.

  -h, --help				Print usage.		
	\n";

	my ($pavdata, $phenodata, $pheno_name, $region_name, $out, $help);

	my $fig_width = 8;
        my $fig_height = 6;

	my $pav_colors = "gray70,steelblue";
	my $x_text_size = "NULL";
	my $x_title_size = "NULL";
	my $y_text_size = "NULL";
	my $y_title_size = "NULL";
	my $legend_side = "top";
	my $legend_title_size = "NULL";
	my $legend_text_size = "NULL";


	GetOptions(
		'pav|i=s'               => \$pavdata,
                'pheno=s'               => \$phenodata,
                'out|o=s'               => \$out,

		'fig_width=f'           => \$fig_width,
                'fig_height=f'          => \$fig_height,
	
		'pheno_name=s'		=> \$pheno_name,
		'region_name=s'		=> \$region_name,

		'pav_colors=s'		=> \$pav_colors,
		'x_text_size=f'		=> \$x_text_size,
		'x_title_size=f'	=> \$x_title_size,
		'y_text_size=f'		=> \$y_text_size,
		'y_title_size=f'	=> \$y_title_size,
		
		'legend_side=s'		=> \$legend_side,
		'legend_title_size=f'	=> \$legend_title_size,
		'legend_text_size=f'	=> \$legend_text_size,

		'help|h!'		=> \$help
	) or die $!."\n";

	die $usage if !defined($pavdata) & !defined($phenodata) & !defined($pheno_name) & !defined($region_name);
        die $usage if $help;

        APAVutils::check_pav_input($pavdata);
        APAVutils::check_file('--pheno', $phenodata);
        APAVutils::check_arg('--pheno_name', $pheno_name);
        APAVutils::check_arg('--region_name', $region_name);
	$out = APAVutils::check_out($out, $pavdata, "_pheno_".$pheno_name."_".$region_name."_violin");

	die "Please install R first\n" if(system("command -v Rscript > /dev/null 2>&1") != 0);
        die "Please install 'APAVplot' R package first\n" if(system("Rscript -e 'library(APAVplot)' > /dev/null 2>&1") != 0);

        my $dir =  dirname(__FILE__);
        my $exec = $dir."/"."vis_pheno.R";

	system("Rscript $exec $dir $pavdata $phenodata '$out' violin '$pheno_name' '$region_name' '$fig_width' '$fig_height' '$pav_colors' '$x_text_size' '$x_title_size' '$y_text_size' '$y_title_size' '$legend_side' '$legend_title_size' '$legend_text_size' 1>/dev/null ");


}


1;
