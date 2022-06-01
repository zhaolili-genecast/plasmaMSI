#! /urs/bin/perl -w 
use strict;
use Statistics::Basic qw(:all);
use Getopt::Long;
my ($bed,$model,$indir, $outdir, $prefix, $help);

GetOptions(
                        "p:s" => \$bed,
                        "m:s" => \$model,
			"i:s" => \$indir,
                        "o:s" => \$outdir,
                        "s:s" => \$prefix,
                        "help|?" => \$help,
);

if(!$indir || $help){
        die &USAGE;
}

sub USAGE{
        my $usage =<<"USAGE"

usage: perl $0 [options]
options:
      -p <s> bed file
      -i <s> input dir
      -o <s> output dir
      -s <i> prefix

example: perl $0 -p Panel.MSI.bed -m Panel_bMSI_Model.txt -i ./stat  -o ./ -s Panel
author : zhao.lili\@genecast.com.cn
USAGE
}


my %fra;
my %indel;
open BED,"$bed" or die "[ERROR] file $bed open failed!\n";;
while(<BED>){
	chomp;
	next if(/^#/ || /^$/);
	my @tmp=split /\t/;
	$fra{$tmp[3]}=$tmp[5];
	$indel{$tmp[3]}=$tmp[4];
}
close BED;

my @file=glob ("$indir/*.kld.txt");

my %hash;
my %name;
my %samplenum;


foreach my $file(@file){
	my $name1=(split /\//,$file)[-1];
	my ($name,$frac)=split /_down/,$name1;
	$frac=~s/.stat.txt$//;
	open FILE,"$file";
	<FILE>;
	while(<FILE>){
		chomp;
		my @tmp=split /\t/;
		next if (!exists $indel{$tmp[1]});
		my $indel=0;
		my @type=split/\s+/,$tmp[-1];
                foreach my $type(@type){
                        my ($a,$p,$d)=split /\:/,$type;
                        if ($a!=0){
                                $indel=$indel+$d;
                        }else{
				$indel=$indel+$fra{$tmp[1]}*$d;
			}
                }
                next if ($indel<$indel{$tmp[1]});
		
		my $info="$tmp[1]\t$tmp[4]";
		my $sampleinfo="$name\t$tmp[1]\t$tmp[4]";
		if (! exists $name{$sampleinfo} && $tmp[3]!=0){
			$samplenum{$tmp[1]}{$tmp[4]}++;
			$name{$sampleinfo}=1;
			$hash{$info}.=";$tmp[3]";
		}
	}
	close FILE;
}
open OUT,">$outdir/$prefix.kld.dis.txt";
foreach my $marker0(sort keys %hash){
	my ($gene,$dup)=split /\t/,$marker0;
	$hash{$marker0}=~s/^;//;
	my @tmp=split /\;/,$hash{$marker0};
	my $num=@tmp;
	next if ($num<20);
	my $kld;
	print OUT "$gene\t$dup\t$num\t";
	foreach my $info(sort{$a<=>$b} @tmp){
		$kld.= "$info,";

	}
	$kld=~s/,$//;
	print OUT "$kld\n";
}
close OUT;

