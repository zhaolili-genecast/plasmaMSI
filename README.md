# plasmaMSI

plasmaMSI is an algorithm to detect microsatellite instability in plasma cell-free DNA from NGS data.

# Installation
plasmaMSI is written in Perl. Directly download the entire directory and run plasmaMSI.sh to perform MSI dectection on a single test sample. No extra compilation is required.

# Usage

## ConstructeBaseline.sh: Baseline construction

User can run the script ConstructeBaseline.sh to build the Baseline MSS Model.

Usage: 
```bash
ConstructeBaseline.sh <bed file> <model inputdir> <baseline inputdir> <outputdir> <ci> <prefix>
```

### Options
-`bed file` \<string\> the list of markers to those presented in your capture design,format:"chr\tstart\tend\tname\tdepthcutoff\tfraction".\
-`model inputdir` \<string\> the path of directory containg alleledistribution files for building MSS Model.User generates there files with Dup_AlleleCaller.pl \
-`baseline inputdir` \<string\> the path of directory containg alleledistribution files for building BASELINE.User generates there files with Dup_AlleleCaller.pl \
-`outputdir` \<string\> output directory. \
-`ci` \<string\>  Confidence interval to determine MSI-H/MSS.\
-`prefix` \<string\> prefix for output files.


## plasmaMSI.sh: Call if a sample is MSI or MSS

Usage: 
```bash
plasmaMSI.sh \<database\> \<bam\>  \<outputdir\> \<sample\> \<prefix\>
```

### Options
-`database` \<string\> the directory that contains baseline databases. \
-`bam` \<string\> the path of input BAM file. \
-`outputdir` \<string\> output directory. \
-`sample` \<string\> sample name. \
-`prefix` \<string\> prefix for baseline databases.
 
### Other examples
An example command to run MSI detection on a single pre-deduplication BAM is:
```bash
plasmaMSI.sh  ../database  sample.sorted.mkdup.realign.bam  ../examples sample test
```
The program takes a pre-deduplication BAM as input, perform duplication removal, call allele distribution, call loci-level MSI status and sample-level MSI status.

Output files include:

1. `sample`.reads.txt: Each line contains spanning reads information (readID, microsatellite sequence, read sequence, CIGAR and loci name) from the input BAM file.
 
2. `sample`.DUP.txt: Each line contains the loci name, reference name, aligned starting position, insert size, lengths of each members pertaining to a family.

3. `sample`.allele.txt: Each line represents information of all fragments covering a locus, including dup-ratios, and allele length offsets compared to the germline allele and corresponding frequencies.

4. `sample`.kld.txt: Each line represents information of all fragments covering a locus, including KLD value, dup-ratios, and allele length offsets compared to the germline allele and corresponding frequencies.

5. `prefix`.msi.txt: Sample-level MSI calling results.

# Publication
- Fengchang Huang, Lili Zhao, Hongyu Xie, Tiancheng Han, Jian Huang, Xiaoqing Wang, Jun Yang, Yuanyuan Hong, Jingchao Shu, Jianing Yu, Qingyun Li, Yu S. Huang, Weizhi Chen, Ji He, Wenliang Li. plasmaMSI: a systematic method to detect next-generation sequencing-based microsatellite instability in plasma cell-free DNA (2024 Under review)



# Contact
- Lili Zhao (zhaolili_607@126.com)
- Yu S. Huang (polyactis@gmail.com)
