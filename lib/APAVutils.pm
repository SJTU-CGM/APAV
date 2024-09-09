#!/usr/bin/perl

package APAVutils;

use strict;
use warnings;
use Getopt::Long;

sub check_arg {
	my ($param, $in) = @_;
	die "Parameter '$param' not found. Please check your input and try again.\n" if !defined($in);
}

sub check_file {
	my ($param, $in) = @_;
	check_arg($param, $in);
	if($in ne "NULL"){
		die "\"$in\": No such file or directory.\n" if !(-e $in);
	}
}

sub check_out {
	my ($out, $in, $suffix) = @_;
	my $res;
	if(!defined($out)){
		$res = $in;
	       	$res =~ s/.*\///; 
		$res =~ s/\.[^\.]*$//;
		$res .= $suffix;
	}else{
		$res = $out;
	}
	return $res;
}

sub check_pav_input {
	my ($pavdata) = @_;
	check_file('--pav|-i', $pavdata);
	my $head = `grep -v '#' $pavdata | head -n 1`;
        if(!($head =~ "^Chr\tStart\tEnd\tLength\tAnnotation")){
                die "Please make sure the input file is produced by command 'apav callPAV'.\n";
        }
}

sub check_bamdir {
	my ($param, $bamdir) = @_;
	APAVutils::check_file($param, $bamdir);
        $bamdir.="/" unless($bamdir=~/\/$/);
        my @bams = <$bamdir*.bam>;
        die "Can not find bam files in the given directory\n" if $#bams == -1;
}

1;
