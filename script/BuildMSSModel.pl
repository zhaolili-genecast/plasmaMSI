#! /urs/bin/perl -w
use strict;
use Statistics::Basic qw(:all);
use List::Util qw/sum max/;
use Getopt::Long;

my ($bed,$indir, $outdir, $prefix, $number,$dupratio,$help);

GetOptions(
                        "p:s" => \$bed,
                        "i:s" => \$indir,
                        "o:s" => \$outdir,
                        "s:s" => \$prefix,
			"n:i" => \$number,
			"d:f" => \$dupratio,
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
      -i <s> input dir containing allele files
      -o <s> output dir
      -s <i> prefix
      -n <i> min sample number
      -d <f> min duplicate ratio
example: perl $0 -p Panel.bMSI.bed -i ./stat -n 25 -d 0.1  -o ./ -s Panel 
author : zhao.lili\@genecast.com.cn
USAGE
}

$number||=25;
$dupratio||=0.1; 
open BED,"$bed" or die "[ERROR] file $bed open failed!\n";
my %fra;
my %indel;
while(<BED>){
	chomp;
	my @tmp=split /\t/;
	$fra{$tmp[3]}=$tmp[5];
	$indel{$tmp[3]}=$tmp[4];
}
close BED;


open OUT,"> $outdir/$prefix\_bMSI_Model.txt" or die "[ERROR] file $outdir/$prefix\_bMSI_Model.txt open failed!\n";
print OUT "Gene\tIndelLength:AlleleFraction:AlleleFractionSD|Indel number|sample number|dupratio\n";

my (%pro,%dep,%index,%number);

my @file=glob "$indir/*.allele.txt";
foreach my $file(@file){
	my $n=(split /\//,$file)[-1];
        my ($name,$frac)=split /_down/,$n;
	open FILE,"$file";
	
	<FILE>;
	while(<FILE>){
		chomp;
		my @tmp=split /\t/;
		my $indel=0;
		next if (!exists $fra{$tmp[1]});
		my @type=split/\s+/,$tmp[-1];
		foreach my $type(@type){
                        my $d=(split /\:/,$type)[2];
                        my $allele=(split /\:/,$type)[0];
			if ($allele!=0){
				$indel=$indel+$d;
			}else{$indel=$indel+$fra{$tmp[1]}*$d;}
                }
		next if ($indel<$indel{$tmp[1]});
		my $dup=$tmp[3];
		my $sample_dup="$name\t$tmp[1]\t$dup";
		next if (exists $index{$sample_dup});
		$index{$sample_dup}=1;
		$number{$tmp[1]}{$dup}++;
		foreach my $type(@type){
			my ($a,$p,$d)=split /\:/,$type;
			$pro{$tmp[1]}{$dup}{$a}{$name}=$d;
		}	 
	}
	close FILE;
}

foreach my $gene(sort keys %pro){
	my $line="";
	foreach my $dup(sort{$a<=>$b} keys %{$pro{$gene}}){
		my @finalp;
		my @finalsd;
		my @finalalle;
		foreach my $allele(sort{$a<=>$b} keys %{$pro{$gene}{$dup}}){
			my @p;
			foreach my $sample (keys %{$pro{$gene}{$dup}{$allele}}){
				push @p,$pro{$gene}{$dup}{$allele}{$sample};
			}
			my $numtmp=@p;
			if ($numtmp<$number{$gene}{$dup}){
				my $distance=$number{$gene}{$dup}-$numtmp;
				for (my $i=1;$i<=$distance;$i++){
					push @p,0;
				}
			}
					
			if ($numtmp/$number{$gene}{$dup} <0.1){
				next;
			}
			my $meanp=sprintf("%0.2f",mean(@p));
			my $sdp=sprintf("%0.2f",stddev(@p));
			if($meanp >1){				
				push @finalp,$meanp;
				push @finalsd,$sdp;
				push @finalalle,$allele;
			}
		}

		my $indeldepth=0;
		if ($dup>=$dupratio && $number{$gene}{$dup}>= $number){
			for(my $i=0;$i<@finalalle;$i++){
				$line.="$finalalle[$i]:$finalp[$i]:$finalsd[$i]|";
				if ($finalalle[$i] !=0){
					$indeldepth=$indeldepth+$finalp[$i];
				}else{
					$indeldepth=$indeldepth+$fra{$gene}*$finalp[$i];
				}
			}
		$line.="$indeldepth|$number{$gene}{$dup}|$dup\t";
		}
	}
	if ($line ne ""){
		$line=~s/\t$//;
        	print OUT "$gene\t$line\n";
        }
}
close OUT;
