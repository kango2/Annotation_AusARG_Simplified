#!/bin/bash
#PBS -N test
#PBS -l ncpus=1,mem=1gb,walltime=00:01:00,storage=gdata/if89+gdata/xl04,jobfs=1GB
#PBS -j oe


cd /g/data/xl04/jc4878/github/Annotation_AusARG_Simplified

module use /g/data/if89/apps/modulefiles
module load Augustus/3.4.0 perllib/v5.26.3 blat/37 RepeatMasker/4.1.2-p1 scipio/1.4 pblat/2.5 pslCDnaFilter/0 parallel/20191022 seqkit/2.5.1 seqtk/1.3 diamond/2.1.9 python3-as-python

echo -e "${AUGUSTUS_CONFIG}"
mkdir -p /g/data/xl04/jc4878/github/Annotation_AusARG_Simplified/PBS
