#!/urs/bin/perl -w
use strict;
use List::Util qw/min max/;
use Getopt::Long;

my ($bed, $model, $indir, $outdir, $help);

GetOptions(
                        "p:s" => \$bed,
                     	"m:s" => \$model,
			"i:s" => \$indir,
                        "o:s" => \$outdir,
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
      -m <s> MSS model file
      -i <s> input dir containing allele file
      -o <s> output dir
example: perl $0 -p Panel.bMSI.bed -m Panel_bMSI_Model.txt -i ./stat   -o ./ 
author : zhao.lili\@genecast.com.cn
USAGE
}

open BED,"$bed" or die "[ERROR] file $bed open failed!\n";
my %fra;
while(<BED>){
	chomp;
	next if (/^#/);
	my @tmp=split /\t/;
	$fra{$tmp[3]}=$tmp[5];
}
close BED;

open MODE,"$model" or die "[ERROR] file $model open failed!\n";
my %maxdup;
my %baseline;
my $e=1;

while(<MODE>){
	chomp;
	next if (/^#/|| /^Gene/);
	my @tmp=split /\t/;
        my $marker=shift @tmp;	
	$maxdup{$marker}=(split /\|/,$tmp[-1])[-1];
	foreach my $info(@tmp){
		my @info=split /\|/,$info;
		my $dup= pop @info;
		$dup=sprintf("%0.2f",$dup);
		pop @info;pop @info;
		my $infos=join "|",@info;
		$baseline{$marker}{$dup}=$infos;
	}
}
close MODE;

my @file=glob "$indir/*.allele.txt";
foreach my $file(@file){
	open STAT,"$file";
	my $name=(split /\//,$file)[-1];
	$name=~s/\.allele\.txt/\.kld\.txt/;
	open OUT, "> $outdir/$name";
	<STAT>;
	print OUT "Position\tName\tDepth\tKLD\tDup_Ratio\tIndelLength:AlleleFraction:SupportingCalls\n";
	while(<STAT>){
        	chomp;
        	my @tmp=split /\t/;
		my $dup=$tmp[3];
        	next if (!exists $baseline{$tmp[1]}{$tmp[3]});
		my %allele=();
	
		my %p=();
		my @ptmp=split /\s+/,$tmp[-1];
		my $pd=0;
		my $pnum=0;
        	foreach my $p(@ptmp){
			my ($allele,$num)=(split /\:/,$p)[0,2];
                        if ($allele!=0){
				$pd=$pd+$num;   
        	                $p{$allele}=$num;
                	        $allele{$allele}=1;
				$pnum++;
			}else{
				if($fra{$tmp[1]}>0){ 
					$pd=$pd+$fra{$tmp[1]}*$num;
					$p{$allele}=$fra{$tmp[1]}*$num;
					$allele{$allele}=1;
					$pnum++;
				}
			}
		}
		
		my %q=();
		my @qtmp=split /\|/,$baseline{$tmp[1]}{$tmp[3]};
		my $qd=0;
		my $qnum=0;
		foreach my $q(@qtmp){
			my ($allele,$num)=(split /\:/,$q)[0,1];
                        if ($allele!=0){
				$qd=$qd+$num;
                        	$q{$allele}=$num;
                        	$allele{$allele}=1;
				$qnum++;
			}else{
				if ($fra{$tmp[1]}>0){
                              		$qd=$qd+$fra{$tmp[1]}*$num;
                                	$q{$allele}=$fra{$tmp[1]}*$num;
                                	$allele{$allele}=1;
                                	$qnum++;  
				}
                        }
                }
		
		my $all=keys %allele;
		my $qnon=$all-$qnum;
		my $pnon=$all-$pnum;
			
		my $kld=0;
		foreach my $allele(sort {$a<=>$b} keys %allele){
                        if ($qnon>0){
				if (!exists $q{$allele}){
					$q{$allele}=$e;
					$qd=$qd+$e;
				}
                         }
                        if ($pnon>0){
                                if (!exists $p{$allele}){
					$p{$allele}=$e;
					$pd=$pd+$e;
                                }
                        }
                }
		
                foreach my $allele(sort {$a<=>$b} keys %allele){
			my $pp=$p{$allele}/$pd;
                        my $qp=$q{$allele}/$qd;
			if ($p{$allele}==0 || $q{$allele}==0){
			print "$tmp[1]\t$allele\n";
			}
                        $kld+=$pp*((log($pp)/log(2))-(log($qp)/log(2)));
		}
		my $line="$tmp[0]\t$tmp[1]\t$tmp[2]\t$kld\t$dup\t$tmp[4]";
		print OUT "$line\n";
	}
	close STAT;
	close OUT;
}
