#! /ure/bin/perl -w
use strict;
use Statistics::Basic qw(:all);
use Getopt::Long;

my ($file, $window, $help);

GetOptions(
                        "i:s" => \$file,
                        "s:i" => \$window,
                        "help|?" => \$help,
);

if(!$file || $help){
        die &USAGE;
}

sub USAGE{
        my $usage =<<"USAGE"

usage: perl $0 [options]
options:
      -i <s> input file [intermediate.file]
      -s <i> window size
example: perl $0 -i Panel.kld.dis.txt  -s 5  >  Panel.kld.txt
author : zhao.lili\@genecast.com.cn
USAGE
}


my $N=$window;
my $cut=($window-1)/2;

open FILE,"$file" or die "[ERROR] file $file open failed!\n";
my %line;
#print "#Marker\tsd|average|cutoff|number|duplicateratio\n";
while(<FILE>){
	chomp;
	next if (/^#/);
	my @tmp=split /\t/;
	$line{$tmp[0]}.="$tmp[3]|$tmp[2]|$tmp[1]\t";
}
close FILE;

foreach my $marker(sort keys %line){
	$line{$marker}=~s/\t$//;
	my @tmp=split /\t/,$line{$marker};
	my (@kld,@num,@dup);
	foreach my $info(@tmp){
		my($kld,$num,$dup)=(split /\|/,$info)[0,1,2];
		push @kld,$kld;
		push @num,$num;
		push @dup,$dup;
	}
	my $number=@kld;
	for (my $i=0;$i<$number;$i++){
		my $tmpnum;
		my $tmpkld;
		if ($i<$cut){
			print "$marker\t$dup[$i]\t$num[$i]\t$kld[$i]\n";
		}elsif($i>=$cut and $i<$number-$cut){
			for (my $j=$i-$cut;$j<=$i+$cut;$j++){
				$tmpnum+=$num[$j];
				$tmpkld.=",$kld[$j]";
                             }
			$tmpkld=~s/^,//;
                	print "$marker\t$dup[$i]\t$tmpnum\t$tmpkld\n";
                }elsif($i>=$number-$cut){
			print "$marker\t$dup[$i]\t$num[$i]\t$kld[$i]\n";
			
		}
	}
}
close FILE;
