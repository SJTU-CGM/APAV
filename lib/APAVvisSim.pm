#!/usr/bin/perl

package APAVvisSim;

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename;


sub simPlotCurve{

	my $usage = "\n\tUsage: apav pavPlotSim -i <sim_output_file> [options]
	
Necessary input desription:
  -i, --simout		<file>		Result file produced by command 'apav sim'.

Option:
  -o, --out		<string>	Figure name.

Visualization Options:
  
  --chart_type		<type>		Chart type. It should be one of 'errorbar', 'jitter' and 'ribbon'.
  					(Default:'ribbon')
  --data_type		<type>		The type of data displayed. It should be 'count' or 'change'.
  					(Default: 'count')

  --fig_width		<numeric>	The width of the figure.
  --fig_height		<numeric>	The height of the figure.

  --x_title		<string>	The text for the x-axis title.
  --y_title		<string>	The text for the y-axis title.
  --x_title_size	<numeric>	The size of x-axis title.
  --y_title_size	<numeric>	The size of y-axis title.
  --x_breaks		<string>	The break values on the x-axis, separated by commas. (eg: '1,5,10')
  --y_breaks		<string>	The break values on the y-axis, separated by commas. (eg: '10,20,30')
  --x_text_size		<numeric>	The size of tick labels on x-axis.
  --y_text_size		<numeric>	The size of tick labels on x-axis.

  --hide_legend				Hide the legend.
  --legend_side		<type>		The position of legend, choosing from 'top', 'bottom', 'left', 'right'.(Default: right)
  --legend_title	<string>	The text for the legend title.
  --legend_title_size	<numeric>	The size of legend title.
  --legend_text_size	<numeric>	The size of legend item labels.
   
  --errorbar_width	<numeric>	The relative width of errorbar, ranging from 0 to 1.
  --errorbar_size	<numeric>	The size of errorbar.
  --errorbar_color	<colors>	The colors for errorbar, separated by commas.
  --errorbar_alpha	<numeric>	The opacity of errorbar, ranging from 0 to 1.
  --errorbar_point_size	<numeric>	The size of point representing the mean value.
  --errorbar_point_color<colors>	The colors for point representing the mean value, separated by commas.

  --jitter_width	<numeric>	The numeric vector giving the relative width of jittered points, ranging from 0 to 1.
  --jitter_size		<numeric>	The size of jittered points.
  --jitter_color	<colors>	The colors for jittered points, separated by commas.
  --jitter_alpha	<numeric>	The opacity of jittered points, ranging from 0 to 1.
  --jitter_point_size	<numeric>	The size of point representing the mean value.
  --jitter_point_color	<colors>	The colors for point representing the mean value, separated by commas.

  --path_size		<numeric>	The size of path.
  --path_color		<colors>	The colors for path, separated by commas.
  --ribbon_fill		<colors>	The colors for ribbon, separated by commas.
  --ribbon_alpha	<numeric>	The opacity of ribbon, ranging from 0 to 1.

  -h, --help				Print usage page.
  \n";

	my ($simdata, $out, $help);
	
	my $chart_type = "ribbon";
	my $data_type = "count";

	my $fig_width = 6;
	my $fig_height = 4;

	my $x_title = "Sample Number";
	my $y_title = "NULL";
	my $x_title_size = "NULL";
    	my $y_title_size = "NULL";
    	my $x_breaks = "NULL";
    	my $y_breaks = "NULL";
    	my $x_text_size = "NULL";
    	my $y_text_size = "NULL";

    	my $hide_legend = 0;  ## F
    	my $legend_side = "right";
    	my $legend_title = "NULL";
    	my $legend_title_size = "NULL";
    	my $legend_text_size = "NULL";

    	my $errorbar_width = 0.8;
    	my $errorbar_size = 1;
    	my $errorbar_color = "NULL";
    	my $errorbar_alpha = 0.8;
    	my $errorbar_point_size = 2;
    	my $errorbar_point_color = "NULL";

    	my $jitter_width = 0.8;
    	my $jitter_size = 1;
    	my $jitter_color = "NULL";
    	my $jitter_alpha = 0.1;
    	my $jitter_point_size = 2;
    	my $jitter_point_color = "NULL";

    	my $path_size = 1;
    	my $path_color = "NULL";
    	my $ribbon_fill = "NULL";
    	my $ribbon_alpha = 0.5;

	GetOptions(
		'simout|i=s'	=> \$simdata,
		'out|o=s'	=> \$out,

		'chart_type=s'	=> \$chart_type,
		'data_type=s'	=> \$data_type,

		'fig_width=f'		=> \$fig_width,
		'fig_height=f'		=> \$fig_height,

		'x_title=s'		=> \$x_title,
		'y_title=s'		=> \$y_title,
		'x_title_size=f'	=> \$x_title_size,
		'y_title_size=f'	=> \$y_title_size,
		'x_breaks=s'		=> \$x_breaks,
		'y_breaks=s'		=> \$y_breaks,
		'x_text_size=f'		=> \$x_text_size,
		'y_text_size=f'		=> \$y_text_size,

		'hide_legend!'		=> \$hide_legend,
		'legend_side=s'		=> \$legend_side,
		'legend_title=s'	=> \$legend_title,
		'legend_title_size=f'	=> \$legend_title_size,
		'legend_text_size=f'	=> \$legend_text_size,

		'legend_text_size=f'	=> \$legend_text_size,
		'errorbar_size=f'	=> \$errorbar_size,
		'errorbar_color=s'	=> \$errorbar_color,
		'errorbar_alpha=f'	=> \$errorbar_alpha,
		'errorbar_point_size=f'	=> \$errorbar_point_size,
		'errorbar_point_color=s'=> \$errorbar_point_color,

		'jitter_width=f'	=> \$jitter_width,
		'jitter_size=f'		=> \$jitter_size,
		'jitter_color=s'	=> \$jitter_color,
		'jitter_alpha=f'	=> \$jitter_alpha,
		'jitter_point_size=f'	=> \$jitter_point_size,
		'jitter_point_color=s'	=> \$jitter_point_color,

		'path_size=f'		=> \$path_size,
		'path_color=s'		=> \$path_color,
		'ribbon_fill=s'		=> \$ribbon_fill,
		'ribbon_alpha=f'	=> \$ribbon_alpha,

		'help|h!'	=> \$help
	) or die $!."\n";

	die $usage if !defined($simdata);

	die $usage if $help;
	APAVutils::check_file('--simout/-i', $simdata);
	my $head = `grep -v '#' $simdata | head -n 1`;
	if($head ne "Round\tSampleN\tCore\tPan\tDelta\n" && $head ne "Round\tSampleN\tCore\tPan\tDelta\tGroup\n"){
		die "Please make sure the input file is produced by command 'apav pavSim'\n";
	}
	$out = APAVutils::check_out($out, $simdata, "_sim_curve");

	die "Please install R first\n" if(system("command -v Rscript > /dev/null 2>&1") != 0);
	die "Please install 'APAVplot' R package first\n" if(system("Rscript -e 'library(APAVplot)' > /dev/null 2>&1") != 0);
	
	my $dir =  dirname(__FILE__);
	my $exec = $dir."/"."vis_sim.R";

	system("Rscript $exec $dir \"$simdata\" \"$out\" \"$chart_type\" \"$data_type\" \"$x_title\" \"$y_title\" \"$x_title_size\" \"$y_title_size\" \"$x_breaks\" \"$y_breaks\" \"$x_text_size\" \"$y_text_size\" \"$hide_legend\" \"$legend_side\" \"$legend_title\" \"$legend_title_size\" \"$legend_text_size\" \"$errorbar_width\" \"$errorbar_size\" \"$errorbar_color\" \"$errorbar_alpha\" \"$errorbar_point_size\" \"$errorbar_point_color\" \"$jitter_width\" \"$jitter_size\" \"$jitter_color\" \"$jitter_alpha\" \"$jitter_point_size\" \"$jitter_point_color\" \"$path_size\" \"$path_color\" \"$ribbon_fill\" \"$ribbon_alpha\" \"$fig_width\" \"$fig_height\" 1>/dev/null ");


}

1;


