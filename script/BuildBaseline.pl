#! /urs/bin/perl -w
use strict;
use Getopt::Long;

my ($input, $outdir, $prefix, $number,$dupratio,$help);

GetOptions(
                        "i:s" => \$input,
                        "o:s" => \$outdir,
                        "s:s" => \$prefix,
                        "n:i" => \$number,
                        "d:f" => \$dupratio,
                        "help|?" => \$help,
);

if(!$input || $help){
        die &USAGE;
}
sub USAGE{
        my $usage =<<"USAGE"

usage: perl $0 [options]
options:
      -i <s> input file
      -o <s> output dir
      -s <i> prefix
      -n <i> min sample number
      -d <f> min duplicate ratio
example: perl $0 -i  cutoff.txt -n 25 -d 0.1  -o ./ -s Panel
author : zhao.lili\@genecast.com.cn
USAGE
}

$number||=25;
$dupratio||=0.1;

my %hash;
open FILE,"$input" or die "[ERROR] file $input open failed!\n";;
<FILE>;
while(<FILE>){
	chomp;
	my @tmp=split /\t/;
	$tmp[1]=sprintf ("%.2f",$tmp[1]);
	if ($tmp[1]>=$dupratio && $tmp[2]>=$number){
		$hash{$tmp[0]}{$tmp[1]}="$tmp[6]|$tmp[5]|$tmp[4]|$tmp[2]|$tmp[1]";
	}
}	
close FILE;

open OUT,">$outdir/$prefix\_bMSI_BASELINE.txt";
foreach my $gene(keys %hash){
	print OUT  "$gene";
	foreach my $dup(sort {$a<=>$b} keys %{$hash{$gene}}){
		print OUT "\t$hash{$gene}{$dup}";
	}
	print OUT "\n";
}
close OUT;
