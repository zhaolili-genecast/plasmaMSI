#!/bin/bash
if [ -z "$5" ]; then
        echo
        echo "Usage: $0 <baseline> <bam>  <outputdir> <sample> <prefix>"
        echo
        exit 1
fi
baseline=$1
bam=$2
outputdir=$3
sample=$4
prefix=$5
set -e

WORKFLOW_DIR="$(dirname $(realpath $0))"
if [ ! -d ${outputdir} ]; then
	mkdir -p ${outputdir}
fi

##### caculate dupratio && call alleles  #####
perl ${WORKFLOW_DIR}/Dup_AlleleCaller.pl  -p ${baseline}/${prefix}.bMSI.bed  -b ${bam}  -f 2 -n 2 -o ${outputdir}  -s ${sample}

##### caculate kld #####  
perl ${WORKFLOW_DIR}/CalculateKLD.pl -p ${baseline}/${prefix}.bMSI.bed  -m ${baseline}/${prefix}_bMSI_Model.txt -i ${outputdir}  -o ${outputdir} 
#### clacssify msi status ###
perl ${WORKFLOW_DIR}/Classify.pl -p ${baseline}/${prefix}.bMSI.bed -b ${baseline}/${prefix}_bMSI_BASELINE.txt  -i ${outputdir}/${sample}.kld.txt  -f ${outputdir}/${sample}.msi.txt  
