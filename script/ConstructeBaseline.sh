#!/bin/bash
if [ -z "$6" ]; then
        echo
        echo "Usage: $0 <bed file> <model inputdir> <baseline inputdir> <outputdir> <ci> <prefix>"
        echo
        exit 1
fi
bed=$1
minputdir=$2
binputdir=$3
outputdir=$4
ci=$5
prefix=$6
set -e

kld=${outputdir}/kld
WORKFLOW_DIR="$(dirname $(realpath $0))"
if [ ! -d ${kld} ]; then
        mkdir -p ${kld}
fi

##### Build MSS Model #####
perl ${WORKFLOW_DIR}/BuildMSSModel.pl -p ${bed} -i ${minputdir} -n 25 -d 0.1 -o ${outputdir} -s /${prefix}
##### Build BASELINE ##### 
perl ${WORKFLOW_DIR}/CalculateKLD.pl -p ${bed}  -m ${outputdir}/${prefix}_bMSI_Model.txt -i ${binputdir} -o ${kld}
perl ${WORKFLOW_DIR}/CountKLD.pl -p ${bed} -i ${kld} -o ${outputdir} -s ${prefix}
perl ${WORKFLOW_DIR}/Smooth.pl -i ${outputdir}/${prefix}.kld.dis.txt -s 5 > ${outputdir}/${prefix}.kld.txt
Rscript ${WORKFLOW_DIR}/KLDCutoff.R ${outputdir}/${prefix}.kld.txt ${ci} ${outputdir}/${prefix}.forbaseline.txt
perl ${WORKFLOW_DIR}/BuildBaseline.pl -i ${outputdir}/${prefix}.forbaseline.txt  -n 25 -d 0.1  -o ${outputdir} -s ${prefix}
