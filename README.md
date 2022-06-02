# bMSI-CAST

bMSI-CAST (blood MSI Caller Adjusted with Sequence duplicaTes), an NGS-based  MSI detection method is compatible plasma samples and does not require matched normal samples.

# Installation
bMSI-CAST is written in Perl. Directly download the entire directory and run bmsicast.sh to perform MSI dectection on a single test sample. No extra compilation is required.

# Usage

Version 0.10
## baseline construction

User can run the script ConstructeBaseline.sh to build the MSS Model and Baseline.

Usage: ConstructeBaseline.sh <bed file> <model inputdir> <baseline inputdir> <outputdir> <ci> <prefix>

##Options
-`bed file` \<string\> the list of markers to those presented in your capture design,format:"chr\tstart\tend\tname\tdepthcutoff\tfraction".
-`model inputdir` \<string\> the path of directory containg alleledistribution files for building MSS Model.User generates there files with Dup_AlleleCaller.pl
-`baseline inputdir` \<string\> the path of directory containg alleledistribution files for building BASELINE.User generates there files with Dup_AlleleCaller.pl
-`outputdir` \<string\> output directory. 
-`ci` \<string\>  Confidence interval to determine MSI-H/MSS.
-`prefix` \<string\> prefix for output files.


## run MSI detection

Usage: bmsicast.sh <database> <bam>  <outputdir> <sample> <prefix>

## Options
-`database` \<string\> the directory that contains baseline databases.  
-`bam` \<string\> the path of input BAM file.  
-`outputdir` \<string\> output directory.
-`sample` \<string\> sample name.
-`prefix` \<string\> prefix for baseline databases.  
 
# Example
An example command to run MSI detection on a single pre-deduplication BAM is:
```
perl $0 -database <database> -bam `bam` -flank 2 -num 2 -minsize 1 -quality 1  -outdir MSI -prefix `prefix`
```
The program takes a pre-deduplication BAM as input, perform duplication removal, call allele distribution, call loci-level MSI status and sample-level MSI status.

Output files include:
1. `sample`.reads.txt
Each entry contains spanning reads information (readID, microsatellite sequence, read sequence, CIGAR and loci name) from the input BAM file.
 
2. `sample`.DUP.txt
Each row contains the loci name, reference name, aligned starting position, insert size, lengths of each members pertaining to a family.

3. `sample`.allele.txt
Each row represents info. of all fragments covering a locus, including dup-ratios, and allele length offsets compared to the germline allele and corresponding frequencies.

4. `sample`.kld.txt

Each row represents info. of all fragments covering a locus, including KLD value, dup-ratios, and allele length offsets compared to the germline allele and corresponding frequencies.

5. `prefix`.msi.txt
One line sample-level MSI calling results.


# Contact
We can be reached by: Lili Zhao, zhao.lili@genecast.com.cn; Hongyu Xie, xie.hongyu@genecast.com.cn
