#!/bin/bash
#PBS -N 2TrainingAugustus
#PBS -l ncpus=8,walltime=48:00:00,storage=gdata/if89+gdata/xl04,mem=80GB,jobfs=80GB
#PBS -j oe



### Species name, dont use whitespace
export species="Species_name"
### working directory, keep same as in GenerateTrainingGene
export workingdir="/path/to/workingdir"
### Path to genome, keep same as in 1GenerateTrainingGene
export genome="/g/data/xl04/jc4878/Bassiana_publication/YAHS/genome/BASDU_HifiASM_YAHS_SUP_CurV1.1.fasta.masked"


cd ${workingdir}
mkdir -p ${workingdir}/Augustus
mkdir -p ${workingdir}/Augustus/training
mkdir -p ${workingdir}/Augustus/config

### Loading modules
module use /g/data/if89/apps/modulefiles
module load Augustus/3.4.0 perllib/v5.26.3 blat/37 RepeatMasker/4.1.2-p1 scipio/1.4 pblat/2.5 pslCDnaFilter/0 python3-as-python

rsync -a $AUGUSTUS_CONFIG_PATH/ Augustus/config/
export AUGUSTUS_CONFIG_PATH=${workingdir}/Augustus/config

cd ${workingdir}/Augustus/training
gff2gbSmallDNA.pl ${workingdir}/TrainingGene/Training_gene_NoOverlap.gff3 ${genome} 2000 ${workingdir}/Augustus/training/training_unfiltered.gb
new_species.pl --species=${species}
grep -h -vE '^(#|$)' ${workingdir}/Augustus/config/species/${species}/${species}_metapars.cfg \
${workingdir}/Augustus/config/species/${species}/${species}_metapars.utr.cfg \
> ${workingdir}/Augustus/config/species/${species}/${species}_metapars_and_utr.cfg

### Filter out unsatified genes in training gene file
etraining --species=${species} --UTR=on --print_utr=on --stopCodonExcludedFromCDS=f ${workingdir}/Augustus/training/training_unfiltered.gb 2> training_unfiltered.err
cat training_unfiltered.err | perl -pe 's/.*in sequence (\S+): .*/$1/' > badgenes.lst
filterGenes.pl badgenes.lst ${workingdir}/Augustus/training/training_unfiltered.gb  > ${workingdir}/Augustus/training/training_filtered.gb 

### randomly split the filtered training files into 100 + the rest, 100 will be used for evaluating accuracy
### This generates two additional file
### training_filtered.gb.train (original number - 100)
### training_filtered.gb.test (100)
randomSplit.pl ${workingdir}/Augustus/training/training_filtered.gb 100

### Then split the .train into 500 for optimising augustus parametersm, the 500 further splits into 200 for intermediate evaluation and 300 rest for only training (which is fast)
### This generates two additional file
### training_filtered.gb.train.train.train (300)
### training_filtered.gb.train.train.test (200)
randomSplit.pl training_filtered.gb.train 500
###            training_filtered.gb.train.test 500
###            training_filtered.gb.train.train (original number - 100 -500)
randomSplit.pl training_filtered.gb.train.test 200
###            training_filtered.gb.train.test.test 200
###            training_filtered.gb.train.test.train 300



### an initial training with genes not used in evaluation
etraining --species=${species} --UTR=on --print_utr=on --stopCodonExcludedFromCDS=f ${workingdir}/Augustus/training/training_filtered.gb.train.test
### Use 100 for first evaluation
augustus --species=${species} --UTR=on --print_utr=on --stopCodonExcludedFromCDS=f ${workingdir}/Augustus/training/training_filtered.gb.test | tee first_evaluation.out
grep -A 36 Evaluation first_evaluation.out > ${workingdir}/log/first_evaluation.report


# Now we optimize with 500 genes, 200 for evaluation and all 500 for training, the max 5 rounds has been chosen but it will finish earlier if no improvement are found
optimize_augustus.pl \
--species=${species} \
--cpus=${PBS_NCPUS} \
--rounds=5 \
${workingdir}/Augustus/training/training_filtered.gb.train.test.test \
--onlytrain=${workingdir}/Augustus/training/training_filtered.gb.train.test.train \
--UTR=on \
--stopCodonExcludedFromCDS=f \
--metapars=${workingdir}/Augustus/config/species/${species}/${species}_metapars_and_utr.cfg \
--cleanup=1 \
> ${workingdir}/OptimiseAugustus.running

# Retrain after optimization
etraining --species=${species} --UTR=on --print_utr=on --stopCodonExcludedFromCDS=f ${workingdir}/Augustus/training/training_filtered.gb.train.test

# final evaluation
augustus --species=${species} --UTR=on --print_utr=on --stopCodonExcludedFromCDS=f ${workingdir}/Augustus/training/training_filtered.gb.test | tee final_evaluation.out
grep -A 36 Evaluation final_evaluation.out > ${workingdir}/log/final_evaluation.report



# rename log file from .running to .done
mv ${workingdir}/OptimiseAugustus.running ${workingdir}/log/OptimiseAugustus.done
