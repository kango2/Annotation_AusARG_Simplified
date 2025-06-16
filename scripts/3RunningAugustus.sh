#!/bin/bash
#PBS -N 3RunningAugustus
#PBS -l ncpus=32,walltime=24:00:00,storage=gdata/if89+gdata/xl04,mem=100GB
#PBS -j oe



# ### Species name, keep same as in 2TrainingAugustus.sh
# export species="Species_name"
# ### working directory, keep same as in 2TrainingAugustus & 1GenerateTrainingGene
# export workingdir="/path/to/workingdir"
# ### Path to genome, keep same as in 2GenerateTrainingGene & 1GenerateTrainingGene
# export genome="/g/data/xl04/jc4878/Bassiana_publication/YAHS/genome/BASDU_HifiASM_YAHS_SUP_CurV1.1.fasta.masked"
# ### Path to augustus config folder
# # keep below unchanged if $workingdir is same as in 2TrainingAugustus
# # Point this to your previously trained config folder if coming from 1GenerateTrainingGene
# export AUGUSTUS_CONFIG="${workingdir}/Augustus/config"
# ### Path to custom extrinsic.cfg file
# export CUSTOM_EXTRINSIC="/g/data/xl04/jc4878/github/Annotation_AusARG/extrinsic.MPE_modded.cfg"
# ### Path to uniprot_swissprot diamond database
# export uniprot_sprot="/g/data/xl04/hrp561/basdurnaseq/uniprot_sprot.diamond.db.dmnd"
# ### Path to uniprot_trembl diamond database
# export uniprot_trembl="/g/data/xl04/hrp561/basdurnaseq/uniprot_trembl.diamond.db.dmnd"





cd ${workingdir}
mkdir -p ${workingdir}/Augustus
mkdir -p ${workingdir}/Augustus/annotation
mkdir -p ${workingdir}/Augustus/annotation/parallel
mkdir -p ${workingdir}/Augustus/annotation_functional
mkdir -p ${workingdir}/Augustus/Split_genome
mkdir -p ${workingdir}/RunningAugustus
mkdir -p ${workingdir}/RunningAugustus/script
mkdir -p ${workingdir}/log

###################################################################################
module use /g/data/if89/apps/modulefiles
module load Augustus/3.4.0 perllib/v5.26.3 blat/37 RepeatMasker/4.1.2-p1 scipio/1.4 pblat/2.5 pslCDnaFilter/0 parallel/20191022 seqkit/2.5.1 seqtk/1.3 diamond/2.1.9 python3-as-python

### Splitting each fasta file in the multi-fasta genomes to a single fasta file
## Drop comment in fasta header and split 
seqtk seq -C ${genome} | seqkit split -i -O ${workingdir}/Augustus/Split_genome
## Rename each individual file into the scaffold name they contain
for f in ${workingdir}/Augustus/Split_genome/*part*; do
NAME=$(grep ">" $f); mv $f ${workingdir}/Augustus/Split_genome/${NAME#>}.fa
done

seqtk comp ${genome} | \
cut -f1,2 | awk -v var="$workingdir" '{print var"/Augustus/Split_genome/"$1".fa\t"var"/Hints/Hints.gff3\t1\t"$2}' > ${workingdir}/RunningAugustus/GenomeChrIDs.lst

## the temp output where augustus spit its output to
export augDir=${workingdir}/Augustus/annotation/parallel
export augCall="augustus --species=${species} --softmasking=1 --UTR=on --print_utr=on --extrinsicCfgFile=${CUSTOM_EXTRINSIC} --exonnames=off --stopCodonExcludedFromCDS=f --AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG}"

createAugustusJoblist.pl --sequences=${workingdir}/RunningAugustus/GenomeChrIDs.lst --wrap="#" --overlap=2000000 \
--chunksize=20000000 --outputdir=$augDir/ --joblist=${workingdir}/RunningAugustus/jobs.lst \
--jobprefix=${workingdir}/RunningAugustus/script/AugustusScript --command "$augCall"
chmod +x ${workingdir}/RunningAugustus/script/AugustusScript*


## Running gene prediction in parallel
parallel --jobs ${PBS_NCPUS} --no-notice "nice {}" < ${workingdir}/RunningAugustus/jobs.lst

for x in $(cat ${workingdir}/RunningAugustus/jobs.lst)
do
    cat $(cat $x | perl -ne 'if(m/--outfile=(\S*) --errfile/){print $1}')
done | join_aug_pred.pl > ${workingdir}/Augustus/annotation/augustus.gff

## Process the output
cd ${workingdir}/Augustus/annotation
getAnnoFasta.pl augustus.gff
grep -v "#" augustus.gff > augustus.gtf
gtf2gff.pl < augustus.gtf --out=augustus.gff3 --gff3
cp augustus.aa ${workingdir}/Augustus_peptide.fasta

## cat error files into log, this will likely contain warnings about fully soft-masked sequences in the genome. This is normal
cat ${workingdir}/Augustus/annotation/parallel/*.err > ${workingdir}/log/augustus.err


cd ${workingdir}/Augustus/annotation_functional
## Running diamond blastp to uniprot_sprot and uniprot_trembl
diamond blastp --db ${uniprot_sprot} \
--out blast_uniprot_sprot.out \
--query ${workingdir}/Augustus/annotation/augustus.aa \
--un uniprot_sprot_unaligned.fa --unfmt fasta --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles \
--max-target-seqs 1 --max-hsps 1 --evalue 1e-3 --threads ${PBS_NCPUS}

diamond blastp --db ${uniprot_trembl} \
--out blast_uniprot_trembl.out \
--query ${workingdir}/Augustus/annotation_functional/uniprot_sprot_unaligned.fa \
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles \
--max-target-seqs 1 --max-hsps 1 --evalue 1e-3 --threads ${PBS_NCPUS}

cat blast_uniprot_sprot.out blast_uniprot_trembl.out > blast_uniprot.out

## combine gene ID with uniprot ID and gene name from blast results. Prints "unknown" if no uniprot hit
join -1 1 -2 1 -t $'\t' -a 1 -e unknown -o1.1,2.2,2.3 <(awk 'BEGIN {FS="\t"; OFS="\t"} $3 == "gene" {sub(/;/, "", $9); sub(/^ID=/, "", $9); {print $9".t1"}}' ${workingdir}/Augustus/annotation/augustus.gff3 | sort -k1,1) \
<(awk -F'\t' '{if ($13 ~ /GN=/) {split($13,a,"GN="); split(a[2],b," "); split($2, c, "|"); $2 = c[2]; print $1"\t"$2"\t" b[1]} else {split($2, c, "|"); print $1"\t"c[2]"\tunknown"}}' ${workingdir}/Augustus/annotation_functional/blast_uniprot.out | sort -k1,1) | \
sort -k2 -V | awk -F'\t' -v OFS='\t' '{sub(/\.t1$/,"",$1); print "ID="$0}' | awk '{print $1";\t"$2"\t"$3}' > geneID_to_uniprotID.tabular

## add uniprot ID to gff3 file as a new attribute in 9th column
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2";gene_name="$3; next} $9 in a {$9=$9"uniprot_ID="a[$9]";"}1' geneID_to_uniprotID.tabular ${workingdir}/Augustus/annotation/augustus.gff3 | \
sed $'1i##gff-version 3' > ${workingdir}/Augustus_annotation.gff3
cd ${workingdir}

## generate a tabular file with gene ID, midpoint, gene length, strand, uniprot ID and gene name
awk -F'\t' '$3 == "gene" {midpoint = int(($4 + $5) / 2); split($9,a,"ID="); split(a[2],b,";"); split($9,c,"uniprot_ID="); split(c[2],d,";"); split($9,e,"gene_name="); split(e[2],f,";"); print $1"\t"midpoint"\t"$5-$4"\t"$7"\t"b[1]"\t"d[1]"\t"f[1]}' ${workingdir}/Augustus_annotation.gff3 | \
sed $'1iContig\tMiddle position\tGene length\tStrand\tGene ID\tUniprot accession number\tGene name' > ${workingdir}/Augustus_gene_table.tabular
