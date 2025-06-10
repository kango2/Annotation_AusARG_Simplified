#!/bin/bash
#PBS -N minimap2
#PBS -l ncpus=16,walltime=1:00:00,storage=gdata/if89+gdata/xl04,mem=100GB
#PBS -j oe


set -e
module use /g/data/if89/apps/modulefiles
module unload minimap2 samtools
module load minimap2/2.26 samtools/1.19.2

# Check if transcriptome ends with .fa or .fasta
if [[ $transcriptome == *.fa ]]; then
    export base=$(basename "$transcriptome" .fa)
elif [[ $transcriptome == *.fasta ]]; then
    export base=$(basename "$transcriptome" .fasta)
else
    echo "Unsupported transcriptome file extension."
    exit 1
fi


# These alignments are super fast, takes about 1 minute
minimap2 -t ${PBS_NCPUS} -ax splice:hq \
${genome} \
${transcriptome} | samtools sort -@ ${PBS_NCPUS} -O BAM -o ${output}/${base}_2genome.bam
