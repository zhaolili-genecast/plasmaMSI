#! /urs/bin/perl -w
use strict;
use Getopt::Long;

my ($baseline, $bed, $infile, $outfile, $help);
GetOptions(
                        "p:s" => \$bed,
                        "b:s" => \$baseline,
			"i:s" => \$infile,
                        "f:s" => \$outfile,
                        "h|?" => \$help,
);

if(!$bed || !$baseline ||$help){
        die &USAGE;
}

sub USAGE{
        my $usage =<<"USAGE"

usage: perl $0 [options]
options:
      -p <s> bed file
      -b <s> baseline file
      -i <s> inputfile 
      -f <s> output file
example: perl $0 -p Panel.bMSI.bed -b Panel_bMSI_BASELINE.txt -i sample.kld.txt  -f sample_msi.txt 
author : zhao.lili\@genecast.com.cn
USAGE
}
open BASELINE,"$baseline" or die "[ERROR] file $baseline  open failed!\n";
my %baseline;
my %basemarker;
while(<BASELINE>){
	chomp;
	my @tmp=split /\t/;
        my $marker=shift @tmp;
        my $min=(split /\|/,$tmp[0])[3];
        my $max=(split /\|/,$tmp[-1])[3];
	foreach my $info(@tmp){
		my ($sd,$average,$cutoff,$num,$dup)=split /\|/,$info;
		$baseline{$marker}{$dup}=$cutoff;
	}
	$basemarker{$marker}=1;
}
close BASELINE;

my %indeldp;
my %fra;
open BED,"$bed" or die "[ERROR] file $bed open failed!\n";
while(<BED>){
        chomp;
        next if (/Start/);
        my @tmp=split /\t/;
	$indeldp{$tmp[3]}=$tmp[4];
	$fra{$tmp[3]}=$tmp[5];
}
close BED;

open OUT,">$outfile" or die "[ERROR] file $outfile open failed!\n";
my $head;
foreach my $gene(sort keys %basemarker){
	$head.="\t$gene";
}
print OUT "sample_id\tunstable_loci\tpassing_loci\tmsi_score\tmsi_status$head\n";

open FILE,"$infile";
my $name=(split /\//,$infile)[-1];
$name=~s/\.kld\.txt//;
my $passloci=0;
my $positive=0;
my %marker=();
while(<FILE>){
	chomp;
	my @tmp=split /\t/;
	next if (! exists $baseline{$tmp[1]});
	next if (! exists $indeldp{$tmp[1]});
	my $gene=$tmp[1];
	my $depth=$tmp[2]; 
	my $kld=$tmp[3];
	my $dup_ratio=$tmp[4];
	next if (! exists $baseline{$tmp[1]}{$dup_ratio});
	my @infos=split /\s+/,$tmp[5];
	my $indeldepth=0;
	foreach my $info(@infos){
		my ($a,$p,$d)=split /\:/,$info;
		if ($a==0){
			$indeldepth+=$d*$fra{$gene};
		}else{
			$indeldepth+=$d;
		}
	}
	if ($indeldepth >=$indeldp{$gene}){
        	$passloci++;
		my $kldcut=$baseline{$gene}{$dup_ratio};
		if ($kld>$kldcut){
			$marker{$gene}=1;
			$positive++;
		}else{
			$marker{$gene}=0;
		}
	}
}
close FILE;
	
my $status="QNS";
my $pratio="NA";
if ($passloci>0){
	$pratio=sprintf("%.2f",$positive/$passloci);
}

if ($passloci>=15){
        if ($positive>=4 && $pratio >=0.1){
                $status="MSI-H";
        }else{$status="MSS";}
}else{
        if ($positive>=4){
               	$status="MSI-H";
        }else{
                $status="QNS";
        }
}
print OUT "$name\t$positive\t$passloci\t$pratio\t$status";

foreach my $gene(sort keys %basemarker){
        if (!exists $marker{$gene}){
                print OUT "\t";
        }else{
                print OUT "\t$marker{$gene}"
        }
}
print OUT "\n";	
close OUT;
