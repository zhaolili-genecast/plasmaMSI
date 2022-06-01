#! /urs/bin/perl -w
use strict;
use List::Util;
use Getopt::Long;
use Statistics::Basic qw(:all);

my ($bed, $bam, $flank, $num, $outdir, $sample,$samtools, $help);
GetOptions(
			"p:s" => \$bed,
                        "b:s" => \$bam,
			"f:i"=> \$flank,
                        "n:i" => \$num,
			"o:s" => \$outdir,
			"s:s" => \$sample,
			"t:s" => \$samtools,
                        "h|?" => \$help,
);

$flank ||= 2;
$num||= 2;
$samtools||="samtools";
if(!$bed || !$bam || $help){
        die &USAGE;
}

my %hash;
my %chr;
my %pos;
my %markerlen;
open BED,"$bed" or die "[ERROR] file $bed open failed!\n";
open BEDT,">$outdir/$sample\.bed.tmp" or die "[ERROR] file $outdir/$sample\.bed.tmp open failed!\n";
while(<BED>){
	chomp;
	next if (/start/ || /^#/);
	my @tmp=split /\t/;
	$hash{$tmp[3]}=$_;
	$chr{$tmp[0]}=1;
	my $pos="$tmp[0]:$tmp[1]-$tmp[2]";
	$pos{$tmp[3]}=$pos;
	my $extends=$tmp[1];
	my $extende=$tmp[2];
	print BEDT "$tmp[0]\t$extends\t$extende\n";
	$markerlen{$tmp[3]}=$tmp[2]-$tmp[1];
}
close BED;
close BEDT;

my %ref2reads;
my %read;
my %dupread;
my %insert;
my %dup;

open SAM,"$samtools view -F 0xB00 -L $outdir/$sample\.bed.tmp -q 1 $bam|" or die "[ERROR] file $bam open failed!\n";### remove multi alignment
open OUT1,">$outdir/$sample\.reads.txt" or die "[ERROR] file $outdir/$sample\.result.txt open failed!\n";
print OUT1 "ReadID\tMsequence\tReadSequence\n";
open OUT2,">$outdir/$sample\.allele.txt"or die "[ERROR] file $outdir/$sample\.stat.txt open failed!\n";
print OUT2 "Position\tName\tDepth\tDup_Ratio\tIndelLength:AlleleFraction:SupportingCalls\n";
open DUP, "> $outdir/$sample\.DUP.txt" or die "[ERROR] file $outdir/$sample\.DUP.txt open failed!\n";

my %all;
my %dupall;
while(<SAM>){
	chomp;
	my @tmp=split /\t/;
	my $info="$tmp[0]\t$tmp[9]";
	my $StartAndDirection;
	my $StartAndEnd;
        my $flagall=$tmp[1];
	if(exists $chr{$tmp[2]} && !($flagall & 4)){
		%insert=();
		my $cigar=$tmp[5];$cigar=~s/S/S;/g;$cigar=~s/M/M;/g;$cigar=~s/D/D;/g;$cigar=~s/I/I;/g;$cigar=~s/H/H;/g;
                my @cigar=split /;/,$cigar;
		
		my $lenforend=0;
		my $endflag=0;
		foreach my $c(@cigar){
			if($c=~/(\d+)M/||$c=~/(\d+)D/){
				$lenforend+=$1;
			}
		}
		
                if ($tmp[6] eq "=" && $tmp[8]!=0){
                        if ($tmp[8]>0 ){
                                $StartAndDirection="$tmp[2]:$tmp[3]:$tmp[8]:+";
				my $End=$tmp[3]+$tmp[8]-1;
				$StartAndEnd="$tmp[2]:$tmp[3]:$End";
				
                        }elsif($tmp[8]<0){
                                $tmp[8]=0-$tmp[8];
                                $StartAndDirection="$tmp[2]:$tmp[7]:$tmp[8]:-";
				my $End=$tmp[7]+$tmp[8]-1;
				$StartAndEnd="$tmp[2]:$tmp[7]:$End";
                        }
                }else{  
			next;
                        my $cigar=$tmp[5];$cigar=~s/S/S;/g;$cigar=~s/M/M;/g;$cigar=~s/D/D;/g;$cigar=~s/I/I;/g;$cigar=~s/H/H;/g;
                        my @cigar=split /;/,$cigar;
                        my $lenforend=0;
                        foreach my $c(@cigar){
                                if($c=~/(\d+)M/||$c=~/(\d+)D/){
                                        $lenforend+=$1;
                                }
                        }
                	$StartAndDirection="$tmp[2]:$tmp[3]:$lenforend:$tmp[6]:0";
                }

		foreach my $gene(sort keys %hash){
                        my @str=split /\t/,$hash{$gene};
			my ($dupchr,$dupstart,$dupend)=split /\:/,$StartAndEnd;
                        if ($tmp[2] eq $str[0] && ($str[1]+1-$dupstart>=$flank) && ($dupend -$str[2]>=$flank)) {
                                $endflag=1
                        }
                }
                next if ($endflag==0);

		my $readsstart=0; my $lengthref=0;my $lengthreads=0;
                if ($cigar[0]=~/(\d+)S/){## first is soft clip;
                	$readsstart=$1;
		}elsif($cigar[0]=~/(\d+)M/){## first is map;
			$lengthref=$1;
			$lengthreads=$1;
			for(my $j=1;$j<=$lengthref;$j++){
	 			$ref2reads{$j}=$j;
 			}
		}

                for (my $i=1;$i<@cigar;$i++){
               	        if ($cigar[$i]=~/(\d+)M/){
                       	        my $dis=$lengthreads-$lengthref;
				$lengthref+=$1;
				$lengthreads+=$1;
				my $lengthreftmp=$1;
				my $start=$lengthref-$lengthreftmp+1;
				for(my $j=$start;$j<=$lengthref;$j++){
					$ref2reads{$j}=$readsstart+$j+$dis;
				}
                        }
                        if ($cigar[$i]=~/(\d+)I/){
                        	my $dis=$lengthreads-$lengthref;
				$lengthreads+=$1;
				my $lengthreadstmp=$1;
				my $start=$lengthreads-$lengthreadstmp+1;
				for(my $j=$start;$j<=$lengthreads;$j++){
                                	$ref2reads{$lengthref}=$readsstart+$j+$dis;
					my $pos=$lengthref+$tmp[3]-1;
					$insert{$pos}=1
                                }
			}
			if ($cigar[$i]=~/(\d+)D/){
				my $dis=$lengthreads-$lengthref;
				$lengthref+=$1;
				my $lengthreftmp=$1;
				my $start=$lengthref-$lengthreftmp+1;
				for(my $j=$start;$j<=$lengthref;$j++){
                                	$ref2reads{$j}=$readsstart+$lengthreads+$dis;
					my $pos=$lengthref+$tmp[3]-1;
                                }
			}
                }
		foreach my $gene(sort keys %hash){
			next if (exists $read{$gene}{$tmp[0]});###### remove mate reads 
			
			my @dupinfopri=split /\:/,$StartAndDirection;
			my $dupposition;
                        if ($dupinfopri[-1] eq "0"){
				$dupposition="$gene\t$dupinfopri[0]\t$dupinfopri[1]\t$dupinfopri[2]:$dupinfopri[3]:s";
			}else{
                                $dupposition="$gene\t$dupinfopri[0]\t$dupinfopri[1]\t$dupinfopri[2]";
			}
			
			my ($dupchr,$dupstart,$dupend)=split /\:/,$StartAndEnd;
			my @str=split /\t/,$hash{$gene};
			if ($tmp[2] eq $str[0]){
			    if ($str[1]+1-$dupstart>=$flank && $dupend-$str[2]>=$flank){
				if (!exists $dupread{$gene}{$tmp[0]}){
					$dupread{$gene}{$tmp[0]}=1;
					$dupall{$dupposition}++;
				}
				
				if (($str[1]+1-$tmp[3]>=$flank) && (($tmp[3]+$lengthref-1)-$str[2]>=$flank)){
					if (!(%insert && !exists $insert{$str[1]})){
						my $sp=($str[1]-1)-$tmp[3]+1;## ($str[1]+1-2)-$tmp[3]+1: 1-base before number bp of $str[1]
						my $ep=$str[2]-$tmp[3]+1; ## $str[2]-$tmp[3]+1: 1-base of $str[2]
						my $rl=$ref2reads{$ep}-$ref2reads{$sp}-1; ## $ref2reads{$ep} - $ref2reads{$sp} +1-2.
						if ($str[1]+1-$tmp[3]==$flank){
							my $pre=substr($tmp[9],0,$flank+1);
							my @pre=split //,$pre;
							my $base=$pre[-1];
							my $flag=0;
							for (my $i=0;$i<@pre-1;$i++){
								if ($pre[$i] ne $base ){
									$flag=1;
								}
							}
							next if ($flag==0);
						}
						if (($tmp[3]+$lengthref-1)-$str[2]==$flank){
							my $post=substr($tmp[9],-($flank+1));
							my @post=split //,$post;
							my $base=$post[0];
							my $flag=0;
							for (my $i=1;$i<@post;$i++){
								if ($post[$i] ne $base ){
									$flag=1;
								}
							}
							next if ($flag==0);
						}
						my $substart=$ref2reads{$sp}+1; ## $str[1] psoation  0-base
						my $subreads=substr($tmp[9],$substart,$rl);
						my $Len=length($subreads);
						print OUT1 "$tmp[0]\t$subreads\t$tmp[9]\t$tmp[5]\t$gene\n";	
						$read{$gene}{$tmp[0]}=1;
						
						my @dupinfo=split /\:/,$StartAndDirection;
                                                my $position;
                                                if ($dupinfo[-1] eq "0"){
                                                        $position="$gene\t$dupinfo[0]\t$dupinfo[1]\t$dupinfo[2]:$dupinfo[3]:s";
                                                }else{
                                                        $position="$gene\t$dupinfo[0]\t$dupinfo[1]\t$dupinfo[2]";
                                                }
						$dup{$position}.="$Len ";		
					}
				}
			}
			}	
		}
	}
}
close SAM;
close OUT1;
unlink "$outdir/$sample\.bed.tmp";

my %mode;
my %mult_mode;
foreach my $p( sort keys %dup){
	my @tmp=split /\t/,$p;
	my $gene=$tmp[0];
	my $pos="$tmp[1]\t$tmp[2]\t$tmp[3]";
	$dup{$p}=~s/\s+$//;
	my @len=split /\s+/,$dup{$p};
	my $nanumber=$dupall{$p}-@len;
	if ($nanumber>0){
		for (my $i=0;$i<$nanumber;$i++){
			push @len,-1;
		}
	}
	next if ($tmp[3]=~/\:s$/);

	my $num4mode=mode(@len);
	$num4mode=~s/\[//;$num4mode=~s/\]$//;$num4mode=~s/\,//g;
	my @num4mode=split /\s+/,$num4mode;
	if (@num4mode ==1 && $num4mode>0){
		$mode{$gene}{$num4mode}++;
		$all{$gene}+=@len;		
		my $dupp=join " ",@len;
		print DUP "$p\t$dupp\n";
	}elsif(@num4mode ==1 && $num4mode<0){
		my $dupp=join " ",@len;
		print DUP "*$p\t$dupp\n";
	
	}else{
		$mult_mode{$gene}{$pos}=$num4mode;
	}
}
close DUP;

foreach my $gene(sort keys %mode){
	my $depth=0;
	my $allele;	
	foreach my $length(sort keys %{$mode{$gene}}){
                if ($mode{$gene}{$length} >= $num ){
                        $depth=$depth+$mode{$gene}{$length};
                }
        }

        my $dup_ratio=sprintf("%.2f",($all{$gene}-$depth)/$all{$gene});
        if ($depth>=$num){
                foreach my $length(sort {$a<=>$b} keys %{$mode{$gene}}){
                        my $ratio=$mode{$gene}{$length}/$depth;
                        my $len=$length-$markerlen{$gene};
                        if ($mode{$gene}{$length}>=$num){
                                $allele.="$len:$ratio:$mode{$gene}{$length} ";
                        }
                }
		$allele=~s/\s+$//;
		print OUT2 "$pos{$gene}\t$gene\t$depth\t$dup_ratio\t$allele\n";
	}
}
close OUT2;

sub USAGE{
        my $usage =<<"USAGE"

usage: perl $0 [options]
options:
      -p <s> bed format file
      -b <s> alignment bam file
      -f <i> flank size, default [2]
      -n <i> consistency numbern,default[2]
      -s <s> prefix of result files
      -t <s> samtools PATH 
      -o <s> output
example: perl $0 -p Panel.bMSI.bed  -b sample.sorted.mkdup.realign.bam  -f 2 -n 2 -o MSI -s sample20171207-B_I25
author : zhao.lili\@genecast.com.cn

USAGE
}
