#!/usr/bin/perl

package APAVvisCov;

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename;

sub covPlotHeat {

	my $usage = "\nUsage: apav covPlotHeat --cov <cov_file> -o <figure_name> [options]
	
Necessary input description:
  -i, --cov			<file>		Coverage file.

Options:
  -p, --pheno			<file>		Phenotype file.
  -o, --out                     <string>        Figure name.

Visualization options:
  --fig_width           	<numeric>       The width of the figure.
  --fig_height          	<numeric>       The height of the figure.

  --cov_colors			<colors>	The colors for coverage, separated by commas.
  						(Default: 'white,#F8766D')
  --region_info_color_list	<key=value>	The colors for region annotations. (eg: 'chr=black,red')
  --pheno_info_color_list	<key=value>	The colors for phenotype annotations. (eg: 'gender=blue,green')

  --border					Draw border for main heatmap. 

  --cluster_rows				Perform clustering on rows.
  --clustering_distance_rows	<type>		Method of measuring distance when clustring on rows,
  						choosing from 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski', 'pearson', 'spearman', 'kendall'.
						(Default: 'euclidean')
  --clustering_method_rows	<type>		Method to perform hierarchical clustering on rows,
  						choosing from 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.
						(Default: 'complete')
  --row_dend_side		<type>		The position of the row dendrogram, 'left' or 'right'.
  --row_dend_width		<numeric>	The width of the row dendrogram.
  --row_sorted			<string>	The sorted row names, separated by commas. 
  						It doesn't work when '--cluster_rows' is used.

  --show_row_names				Show row names.
  --row_names_side		<type>		The position of row names, 'left' or 'right'.
  --row_names_size		<numeric>	The size of row names.
  --row_names_rot		<numeric>	The rotation of row names.

  --cluster_columns				Perform clustering on columns.
  --clustering_distance_columns	<type>		Method of measuring distance when clustring on columns.
  						choosing from 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski', 'pearson', 'spearman', 'kendall'.
                                                (Default: 'euclidean')
  --clustering_method_columns	<type>		Method to perform hierarchical clustering on columns.
  						choosing from 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.
                                                (Default: 'complete')
  --column_dend_side		<type>		The position of the column dendrogram, 'top' or 'bottom'.
  --column_dend_height		<numeric>	The height of the column dendrogram.
  --column_sorted		<string>	The sorted column names, separated by commas. 
  						It doesn't work when '--cluster_columns'.

  --show_column_names				Show column names.
  --column_names_side		<type>		The position of column names, 'top' or 'column'.
  --column_names_size		<numeric>	The size of column names.
  --column_names_rot		<numeric>	The rotation of column names.

  --anno_param_row_pheno	<key=value>	Parameters for the phenotype annotation,
  						including 'show' <type> ('T','F'), 
							  'width' <numeric>, 
							  'border' <type> ('T', 'F'),
					       		  'name_size' <numeric>, 
							  'name_rot <numeric>, 
							  'name_side' <type> ('top', 'bottom').
  --anno_param_row_stat		<key=value>	Parameters for the stat annotation of rows,
  						including 'show' <type> ('T','F'), 
							  'width' <numeric>, 
							  'border' <type> ('T','F'), 
							  'title' <string>, 
							  'title_size' <numeric>, 
							  'title_rot' <numeric>, 
							  'title_side' <type> ('top', 'bottom'), 
							  'axis_side' <type> ('top', 'bottom'), 
							  'axis_labels_size' <numeric>.
  --anno_param_column_stat	<key=value>	Parameters for the stat annotation of columns,
 						including 'show' <type> ('T','F'), 
							  'height' <numeric>, 
							  'border' <type> ('T','F'), 
							  'title' <string>, 
							  'title_size' <numeric>, 
							  'title_rot' <numeric>, 
							  'title_side' <type> ('left', 'right'), 
							  'axis_side' <type> ('left', 'right'), 
							  'axis_labels_size' <numeric>.

  --legend_side			<type>		The position of legend, choosing from 'top', 'bottom', 'right', 'left'.
  --legend_title		<string>	The text for the legend title.
  --legend_title_size		<numeric>	The size of legend title.
  --legend_text_size		<numeric>	The size of legend item labels.
  --legend_grid_size		<numeric>	The size of legend grid.

  --use_raster					Render the heatmap body as a raster image.

  -h, --help					Print usage page.

Warning: --region_info_color_list --pheno_info_color_list --anno_param_row_phen --anno_param_column_region --anno_param_row_stat --anno_param_column_stat can be added multiple times.
  eg: --anno_param_row_pheno show=T --anno_param_row_pheno width=5 --anno_param_row_pheno name_rot=90
  \n";

	my ($covdata, $out, $help);

	my $phenodata = "NULL";

	my $fig_width = 10;
	my $fig_height = 6;

	my $cov_colors = "white,steelblue";
	my (%region_info_color_list, $region_info_color_list_str);
	my (%pheno_info_color_list, $pheno_info_color_list_str);

	my $border = 0; ## T->F
	my $cluster_rows = 0; ## F
	my $clustering_distance_rows = "euclidean";
	my $clustering_method_rows = "complete";
	my $row_dend_side = "left";
	my $row_dend_width = 5;
	my $row_sorted = "c()";

	my $show_row_names = 0; ## F
	my $row_names_side = "left";
	my $row_names_size = 10;
	my $row_names_rot = 0;

	my $cluster_columns = 0; ## F
	my $clustering_distance_columns = "euclidean";
	my $clustering_method_columns = "complete";
	my $column_dend_side = "top";
	my $column_dend_height = 5;
	my $column_sorted = "c()";

	my $show_column_names = 0; ## F
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
	$anno_param_row_stat{"show"} = "F";
	$anno_param_row_stat{"width"} = 10;
	$anno_param_row_stat{"border"} = "F";
	$anno_param_row_stat{"title"} = "Coverage";
	$anno_param_row_stat{"title_size"} = 10;
	$anno_param_row_stat{"title_side"} = "bottom";
	$anno_param_row_stat{"title_rot"} = 0;
	$anno_param_row_stat{"axis_side"} = "bottom";
	$anno_param_row_stat{"axis_labels_size"} = 8;

	my %anno_param_column_stat;
	$anno_param_column_stat{"show"} = "F";
	$anno_param_column_stat{"height"} = 10;
	$anno_param_column_stat{"border"} = "F";
	$anno_param_column_stat{"title"} = "Coverage";
	$anno_param_column_stat{"title_size"} = 10;
	$anno_param_column_stat{"title_side"} = "left";
	$anno_param_column_stat{"title_rot"} = 0;
	$anno_param_column_stat{"axis_side"} = "left";
	$anno_param_column_stat{"axis_labels_size"} = 8;

	my $legend_side = "right";
	my $legend_title = "Coverage";
	my $legend_title_size = "NULL";
	my $legend_text_size = "NULL";
	my $legend_grid_size = 4;

	my $use_raster = "NULL";	

	
	GetOptions(
                'cov|i=s'			=> \$covdata,
		'pheno|p=s'			=> \$phenodata,
                'out|o=s'			=> \$out,

		'fig_width=f'             	=> \$fig_width,
                'fig_height=f'            	=> \$fig_height,
		
		'cov_colors=s'			=> \$cov_colors,
		'region_info_color_list=s'	=> \%region_info_color_list,
		'pheno_info_color_list=s'	=> \%pheno_info_color_list,

		'border=s'			=> \$border,
		'cluster_rows!'			=> \$cluster_rows,
		'clustering_distance_rows=s'	=> \$clustering_distance_rows,
		'clustering_method_rows=s'	=> \$clustering_method_rows,
		'row_dend_side=s'		=> \$row_dend_side,
		'row_dend_width=f'		=> \$row_dend_width,
		'row_sorted=s'			=> \$row_sorted,

		'show_row_names!'		=> \$show_row_names,
		'row_names_side=s'		=> \$row_names_side,
		'row_names_size=f'		=> \$row_names_size,
		'row_names_rot=f'		=> \$row_names_rot,
		'cluster_columns!'		=> \$cluster_columns,
		'clustering_distance_columns=s'	=> \$clustering_distance_columns,
		'clustering_method_columns=s'	=> \$clustering_method_columns,
		'column_dend_side=s'		=> \$column_dend_side,
		'column_dend_height=f'		=> \$column_dend_height,
		'column_sorted=s'		=> \$column_sorted,

		'show_column_names!'		=> \$show_column_names,
		'column_names_side=s'		=> \$column_names_side,
		'column_names_size=f'		=> \$column_names_size,
		'column_names_rot=f'		=> \$column_names_rot,

		'anno_param_row_pheno=s'	=> \%anno_param_row_pheno,
		'anno_param_column_region=s'	=> \%anno_param_column_region,
		'anno_param_row_stat=s'		=> \%anno_param_row_stat,
		'anno_param_column_stat=s'	=> \%anno_param_column_stat,

		'legend_side=s'			=> \$legend_side,
		'legend_title=s'		=> \$legend_title,
		'legend_title_size=f'		=> \$legend_title_size,
		'legend_text_size=f'		=> \$legend_text_size,
		'legend_grid_size=f'		=> \$legend_grid_size,

		'use_raster!'			=> \$use_raster,

		'help|h!'       => \$help
	) or die $!."\n";
	
	die $usage if !defined($covdata);

	die $usage if $help;
	APAVutils::check_file('--cov/-i', $covdata);
	my $head = `grep -v '#' $covdata | head -n 1`;
	if(!($head =~ "^Chr\tStart\tEnd\tLength\tAnnotation")){
		die "Please make sure the input file is produced by command 'apav staCov', 'apv staElecov' or 'apav mergeElecov'\n";
	}
	$out = APAVutils::check_out($out, $covdata, "_cov_heatmap");
	APAVutils::check_file('--pheno/-p', $phenodata) if defined($phenodata);

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
        my $exec = $dir."/"."vis_cov.R";

	system("Rscript $exec $dir $covdata $phenodata $out \"$cov_colors\" \"$region_info_color_list_str\" \"$pheno_info_color_list_str\" \"$border\" \"$cluster_rows\" \"$clustering_distance_rows\" \"$clustering_method_rows\" \"$row_dend_side\" \"$row_dend_width\" \"$row_sorted\" \"$show_row_names\" \"$row_names_side\" \"$row_names_size\" \"$row_names_rot\" \"$cluster_columns\" \"$clustering_distance_columns\" \"$clustering_method_columns\" \"$column_dend_side\" \"$column_dend_height\" \"$column_sorted\" \"$show_column_names\" \"$column_names_side\" \"$column_names_size\" \"$column_names_rot\" \"$anno_param_row_pheno{\"show\"}\" \"$anno_param_row_pheno{\"width\"}\" \"$anno_param_row_pheno{\"border\"}\" \"$anno_param_row_pheno{\"name_size\"}\" \"$anno_param_row_pheno{\"name_rot\"}\" \"$anno_param_row_pheno{\"name_side\"}\" \"$anno_param_column_region{\"show\"}\" \"$anno_param_column_region{\"height\"}\" \"$anno_param_column_region{\"border\"}\" \"$anno_param_column_region{\"name_size\"}\" \"$anno_param_column_region{\"name_rot\"}\" \"$anno_param_column_region{\"name_side\"}\" \"$anno_param_row_stat{\"show\"}\" \"$anno_param_row_stat{\"width\"}\" \"$anno_param_row_stat{\"border\"}\" \"$anno_param_row_stat{\"title\"}\" \"$anno_param_row_stat{\"title_size\"}\" \"$anno_param_row_stat{\"title_side\"}\" \"$anno_param_row_stat{\"title_rot\"}\" \"$anno_param_row_stat{\"axis_side\"}\" \"$anno_param_row_stat{\"axis_labels_size\"}\" \"$anno_param_column_stat{\"show\"}\" \"$anno_param_column_stat{\"height\"}\" \"$anno_param_column_stat{\"border\"}\" \"$anno_param_column_stat{\"title\"}\" \"$anno_param_column_stat{\"title_size\"}\" \"$anno_param_column_stat{\"title_side\"}\" \"$anno_param_column_stat{\"title_rot\"}\" \"$anno_param_column_stat{\"axis_side\"}\" \"$anno_param_column_stat{\"axis_labels_size\"}\" \"$legend_side\" \"$legend_title\" \"$legend_title_size\" \"$legend_text_size\" \"$legend_grid_size\" \"$use_raster\"  \"$fig_width\" \"$fig_height\" 1>/dev/null \n");


}




1;



