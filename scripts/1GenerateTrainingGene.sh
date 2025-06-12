#!/bin/bash
#PBS -N 2TrainingAugustus
#PBS -l ncpus=12,walltime=48:00:00,storage=gdata/if89+gdata/xl04,mem=190GB,jobfs=10GB
#PBS -j oe



module use /g/data/if89/apps/modulefiles
module unload minimap2 samtools
module load minimap2/2.26 
module load samtools/1.19.2
module load perllib/v5.26.3 genometools/1.6.2
module load seqtk/1.3
module load AGAT/1.4.0
module load pythonlib/3.9.2
module load cdhit/4.8.1
module load bedtools/2.31.0
module load python3-as-python parallel/20191022

export LANG=""
export scriptdir=${repository_path}/scripts
cd ${workingdir}
mkdir -p ${workingdir}/TrainingGene
mkdir -p ${workingdir}/TrainingGene/Workingdir
mkdir -p ${workingdir}/TrainingGene/Workingdir/highQ_cDNA
mkdir -p ${workingdir}/TrainingGene/Workingdir/highQ_CDS
mkdir -p ${workingdir}/TrainingGene/Workingdir/Alignment_cDNA
mkdir -p ${workingdir}/TrainingGene/Workingdir/Alignment_CDS
mkdir -p ${workingdir}/TrainingGene/Workingdir/tmp
mkdir -p ${workingdir}/TrainingGene/Workingdir/Merged
mkdir -p ${workingdir}/Hints
mkdir -p ${workingdir}/Hints/Workingdir
mkdir -p ${workingdir}/Hints/Workingdir/weakerQ_cDNA
mkdir -p ${workingdir}/Hints/Workingdir/weakerQ_CDS
mkdir -p ${workingdir}/Hints/Workingdir/Alignment_cDNA
mkdir -p ${workingdir}/Hints/Workingdir/Alignment_CDS
mkdir -p ${workingdir}/Hints/Workingdir/tmp
mkdir -p ${workingdir}/Hints/Workingdir/Merged
mkdir -p ${workingdir}/log
mkdir -p ${workingdir}/log/PBS

grep '>' ${TRANSLATION_OUTPUT}/*.cds.all.fa | \
perl -lne '$_=~/>(\S+) len=(\d+) .* (\S+)/; $id = $1; $transcriptLen = $2; @info = split(":",$3); if ($info[3] eq "yes" and $info[4] eq "yes" and $info[5] >= 0.95 and $info[5] <= 1.05) { $id =~ /(\S+)_i\d+/; $gid=$1; print "$gid\t$id\t".join("\t",@info)."\t".($info[1]-$info[0])."\t$transcriptLen"; }' | \
sort -k1,1 -k9,9nr -k10,10nr | perl -lne '@a=split("\t",$_); unless (exists $sel{$a[0]}) { print join("\t",@a); $sel{$a[0]}="" }' > ${workingdir}/TrainingGene/Workingdir/highQ.filtered.transcripts.tsv

for i in $(ls ${TRANSLATION_OUTPUT}/*.cds.all.fa); do
    seqtk subseq -l 60 ${i} <(cut -f2 ${workingdir}/TrainingGene/Workingdir/highQ.filtered.transcripts.tsv) > ${workingdir}/TrainingGene/Workingdir/highQ_CDS/$(basename ${i} .cds.all.fa).cds.highQ.fa
done 2>&1 | tee > ${workingdir}/log/filter.log
for i in $(ls ${TRANSLATION_OUTPUT}/*.cdna.all.fa); do
    seqtk subseq -l 60 ${i} <(cut -f2 ${workingdir}/TrainingGene/Workingdir/highQ.filtered.transcripts.tsv) > ${workingdir}/TrainingGene/Workingdir/highQ_cDNA/$(basename ${i} .cdna.all.fa).cdna.highQ.fa
done 2>&1 | tee >> ${workingdir}/log/filter.log

grep '>' ${TRANSLATION_OUTPUT}/*.cds.all.fa | \
perl -lne '$_=~/>(\S+) len=(\d+) .* (\S+)/; $id = $1; $transcriptLen = $2; @info = split(":",$3); if ($info[3] eq "yes" and $info[4] eq "yes" and ($info[5] < 0.95 or $info[5] > 1.05)) { $id =~ /(\S+)_i\d+/; $gid=$1; print "$gid\t$id\t".join("\t",@info)."\t".($info[1]-$info[0])."\t$transcriptLen"; }' | \
sort -k1,1 -k9,9nr -k10,10nr | perl -lne '@a=split("\t",$_); unless (exists $sel{$a[0]}) { print join("\t",@a); $sel{$a[0]}="" }' > ${workingdir}/Hints/Workingdir/weakerQ.filtered.transcripts.tsv

for i in $(ls ${TRANSLATION_OUTPUT}/*.cds.all.fa); do
    seqtk subseq -l 60 ${i} <(cut -f2 ${workingdir}/Hints/Workingdir/weakerQ.filtered.transcripts.tsv) > ${workingdir}/Hints/Workingdir/weakerQ_CDS/$(basename ${i} .cds.all.fa).cds.weakerQ.fa
done 2>&1 | tee > ${workingdir}/log/filterWeaker.log
for i in $(ls ${TRANSLATION_OUTPUT}/*.cdna.all.fa); do
    seqtk subseq -l 60 ${i} <(cut -f2 ${workingdir}/Hints/Workingdir/weakerQ.filtered.transcripts.tsv) > ${workingdir}/Hints/Workingdir/weakerQ_cDNA/$(basename ${i} .cdna.all.fa).cdna.weakerQ.fa
done 2>&1 | tee >> ${workingdir}/log/filterWeaker.log

export output=${workingdir}/TrainingGene/Workingdir/Alignment_CDS
for i in $(ls ${workingdir}/TrainingGene/Workingdir/highQ_CDS/*.cds.highQ.fa); do
    export base=$(basename "$i" .fa)
    echo -e "minimap2 -t 1 -ax splice:hq ${genome} ${i} | samtools sort -@ 1 -O BAM -o ${output}/${base}_2genome.bam" >> ${workingdir}/minimap2_command.list
done
export output=${workingdir}/TrainingGene/Workingdir/Alignment_cDNA
for i in $(ls ${workingdir}/TrainingGene/Workingdir/highQ_cDNA/*.cdna.highQ.fa); do
    export base=$(basename "$i" .fa)
    echo -e "minimap2 -t 1 -ax splice:hq ${genome} ${i} | samtools sort -@ 1 -O BAM -o ${output}/${base}_2genome.bam" >> ${workingdir}/minimap2_command.list
done
##
export output=${workingdir}/Hints/Workingdir/Alignment_CDS
for i in $(ls ${workingdir}/Hints/Workingdir/weakerQ_CDS/*.cds.weakerQ.fa); do
    export base=$(basename "$i" .fa)
    echo -e "minimap2 -t 1 -ax splice:hq ${genome} ${i} | samtools sort -@ 1 -O BAM -o ${output}/${base}_2genome.bam" >> ${workingdir}/minimap2_command.list
done
export output=${workingdir}/Hints/Workingdir/Alignment_cDNA
for i in $(ls ${workingdir}/Hints/Workingdir/weakerQ_cDNA/*.cdna.weakerQ.fa); do
    export base=$(basename "$i" .fa)
    echo -e "minimap2 -t 1 -ax splice:hq ${genome} ${i} | samtools sort -@ 1 -O BAM -o ${output}/${base}_2genome.bam" >> ${workingdir}/minimap2_command.list
done
cat ${workingdir}/minimap2_command.list | parallel --jobs ${PBS_NCPUS} {}
rm ${workingdir}/minimap2_command.list


# Convert bam to gff
for i in $(ls ${workingdir}/TrainingGene/Workingdir/Alignment_cDNA/*.bam); do
  agat_convert_minimap2_bam2gff.pl -i ${i} -o ${workingdir}/TrainingGene/Workingdir/Alignment_cDNA/$(basename ${i} .bam).gff3
done 2>&1 | tee > ${workingdir}/log/PBS/convertBAMtoGFF3.log
for i in $(ls ${workingdir}/TrainingGene/Workingdir/Alignment_CDS/*.bam); do
  agat_convert_minimap2_bam2gff.pl -i ${i} -o ${workingdir}/TrainingGene/Workingdir/Alignment_CDS/$(basename ${i} .bam).gff3
done 2>&1 | tee >> ${workingdir}/log/convertBAMtoGFF3.log

## Same but weaker hints
for i in $(ls ${workingdir}/Hints/Workingdir/Alignment_cDNA/*.bam); do
  agat_convert_minimap2_bam2gff.pl -i ${i} -o ${workingdir}/Hints/Workingdir/Alignment_cDNA/$(basename ${i} .bam).gff3
done 2>&1 | tee > ${workingdir}/log/convertBAMtoGFF3Weaker.log
for i in $(ls ${workingdir}/Hints/Workingdir/Alignment_CDS/*.bam); do
  agat_convert_minimap2_bam2gff.pl -i ${i} -o ${workingdir}/Hints/Workingdir/Alignment_CDS/$(basename ${i} .bam).gff3
done 2>&1 | tee >> ${workingdir}/log/convertBAMtoGFF3Weaker.log



for i in $(ls ${workingdir}/TrainingGene/Workingdir/Alignment_cDNA/*_2genome.gff3); do
  python ${scriptdir}/processminimap2gff3.py ${i} | \
  # Pipe output to awk command to replace ID= or Parent= with the string in column 10, get rid of column 10, then replace cDNA_match_part with exon
  awk 'BEGIN {FS="\t"; OFS="\t"} {
    if ($0 ~ /^#/) {
        print
        next
    }
    if ($3 == "cDNA_match") {
        $3="gene"
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
            if (kv[1] == "ID") {
                parts[i] = kv[1]"=" $10".gene"
            }
        }
        $9 = parts[1]
        sub(/;$/, "", $9)
        print $1"\t"$2"\tmRNA\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tID="$10";Parent="$10".gene"
    } else if ($3 == "cDNA_match_part") {
        $3="exon"
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
            if (kv[1] == "Parent") {
                parts[i] = kv[1]"=" $10
            }
            else if (kv[1] == "ID") {
                parts[i] = kv[1]"=" $10".exonOnLine"NR
            }
        }
        $9 = parts[1] ";" parts[2]
        sub(/;$/, "", $9)
    }
    print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9
}' > ${workingdir}/TrainingGene/Workingdir/Alignment_cDNA/$(basename ${i} .gff3)_processed.gff3
done 2>&1 | tee > ${workingdir}/log/processGFF3.log
for i in $(ls ${workingdir}/TrainingGene/Workingdir/Alignment_CDS/*_2genome.gff3); do
  python ${scriptdir}/processminimap2gff3.py ${i} | \
  # Pipe output to awk command to replace ID= or Parent= with the string in column 10, get rid of column 10, then replace cDNA_match_part with exon
  awk 'BEGIN {FS="\t"; OFS="\t"} {
    if ($0 ~ /^#/) {
        print
        next
    }
    if ($3 == "cDNA_match") {
        $3="gene"
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
            if (kv[1] == "ID") {
                parts[i] = kv[1]"=" $10".gene"
            }
        }
        $9 = parts[1]
        sub(/;$/, "", $9)
        print $1"\t"$2"\tmRNA\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tID="$10";Parent="$10".gene"
    } else if ($3 == "cDNA_match_part") {
        $3="CDS"
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
            if (kv[1] == "Parent") {
                parts[i] = kv[1]"=" $10
            }
            else if (kv[1] == "ID") {
                parts[i] = kv[1]"=" $10".CDSOnLine"NR
            }
        }
        $9 = parts[1] ";" parts[2]
        sub(/;$/, "", $9)
    }
    print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9
}' > ${workingdir}/TrainingGene/Workingdir/Alignment_CDS/$(basename ${i} .gff3)_processed.gff3
done 2>&1 | tee >> ${workingdir}/log/processGFF3.log

## Same but for weaker hints
for i in $(ls ${workingdir}/Hints/Workingdir/Alignment_cDNA/*_2genome.gff3); do
  python ${scriptdir}/processminimap2gff3.py ${i} | \
  # Pipe output to awk command to replace ID= or Parent= with the string in column 10, get rid of column 10, then replace cDNA_match_part with exon
  awk 'BEGIN {FS="\t"; OFS="\t"} {
    if ($0 ~ /^#/) {
        print
        next
    }
    if ($3 == "cDNA_match") {
        $3="gene"
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
            if (kv[1] == "ID") {
                parts[i] = kv[1]"=" $10".gene"
            }
        }
        $9 = parts[1]
        sub(/;$/, "", $9)
        print $1"\t"$2"\tmRNA\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tID="$10";Parent="$10".gene"
    } else if ($3 == "cDNA_match_part") {
        $3="exon"
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
            if (kv[1] == "Parent") {
                parts[i] = kv[1]"=" $10
            }
            else if (kv[1] == "ID") {
                parts[i] = kv[1]"=" $10".exonOnLine"NR
            }
        }
        $9 = parts[1] ";" parts[2]
        sub(/;$/, "", $9)
    }
    print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9
}' > ${workingdir}/Hints/Workingdir/Alignment_cDNA/$(basename ${i} .gff3)_processed.gff3
done 2>&1 | tee > ${workingdir}/log/processGFF3Weaker.log
for i in $(ls ${workingdir}/Hints/Workingdir/Alignment_CDS/*_2genome.gff3); do
  python ${scriptdir}/processminimap2gff3.py ${i} | \
  # Pipe output to awk command to replace ID= or Parent= with the string in column 10, get rid of column 10, then replace cDNA_match_part with exon
  awk 'BEGIN {FS="\t"; OFS="\t"} {
    if ($0 ~ /^#/) {
        print
        next
    }
    if ($3 == "cDNA_match") {
        $3="gene"
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
            if (kv[1] == "ID") {
                parts[i] = kv[1]"=" $10".gene"
            }
        }
        $9 = parts[1]
        sub(/;$/, "", $9)
        print $1"\t"$2"\tmRNA\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tID="$10";Parent="$10".gene"
    } else if ($3 == "cDNA_match_part") {
        $3="CDS"
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
            if (kv[1] == "Parent") {
                parts[i] = kv[1]"=" $10
            }
            else if (kv[1] == "ID") {
                parts[i] = kv[1]"=" $10".CDSOnLine"NR
            }
        }
        $9 = parts[1] ";" parts[2]
        sub(/;$/, "", $9)
    }
    print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9
}' > ${workingdir}/Hints/Workingdir/Alignment_CDS/$(basename ${i} .gff3)_processed.gff3
done 2>&1 | tee >> ${workingdir}/log/processGFF3Weaker.log









# [ Section 2]
## Find transcript where the cDNA alignment completely cover 100% of the CDS alignment in bedtools intersect
## Calculate UTR from exon/CDS difference
## Infer introns for hints
for i in $(ls ${workingdir}/TrainingGene/Workingdir/Alignment_cDNA/*_processed.gff3); do
    awk '$3 == "gene"' ${i} | awk 'BEGIN {FS="\t"; OFS="\t"} 
    {
        if ($3 == "gene") {
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            }
            $9 = parts[1]
            sub(/\.gene$/, "", $9)
            sub(/^ID=/, "", $9)
        }
        print $1"\t"$4"\t"$5"\t"$9"\t"$6"\t"$7
    }' > ${workingdir}/TrainingGene/Workingdir/tmp/cDNA.bed
    awk '$3 == "gene"' ${workingdir}/TrainingGene/Workingdir/Alignment_CDS/$(basename ${i} cdna.highQ_2genome_processed.gff3)cds.highQ_2genome_processed.gff3 | awk 'BEGIN {FS="\t"; OFS="\t"} 
    {
        if ($3 == "gene") {
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            }
            $9 = parts[1]
            sub(/\.gene$/, "", $9)
            sub(/^ID=/, "", $9)
        }
        print $1"\t"$4"\t"$5"\t"$9"\t"$6"\t"$7
    }' > ${workingdir}/TrainingGene/Workingdir/tmp/CDS.bed
    bedtools intersect -F 1.0 -a ${workingdir}/TrainingGene/Workingdir/tmp/cDNA.bed -b ${workingdir}/TrainingGene/Workingdir/tmp/CDS.bed -wa -wb -s | awk '$4 == $10 {print $4}' > ${workingdir}/TrainingGene/Workingdir/tmp/both.txt
    echo '##gff-version 3' > ${workingdir}/TrainingGene/Workingdir/tmp/merge.gff3
    grep -w -f ${workingdir}/TrainingGene/Workingdir/tmp/both.txt ${i} >> ${workingdir}/TrainingGene/Workingdir/tmp/merge.gff3
    awk '$3 == "CDS"' ${workingdir}/TrainingGene/Workingdir/Alignment_CDS/$(basename ${i} cdna.highQ_2genome_processed.gff3)cds.highQ_2genome_processed.gff3 | grep -w -f ${workingdir}/TrainingGene/Workingdir/tmp/both.txt >> ${workingdir}/TrainingGene/Workingdir/tmp/merge.gff3
    gt gff3 -sortlines -tidy -retainids ${workingdir}/TrainingGene/Workingdir/tmp/merge.gff3 | gt gff3 -retainids -sort > ${workingdir}/TrainingGene/Workingdir/tmp/merge2.gff3
    python3 ${scriptdir}/addUTRs.py ${workingdir}/TrainingGene/Workingdir/tmp/merge2.gff3 ${workingdir}/TrainingGene/Workingdir/tmp/merge2_with_UTRs.gff3
    gt gff3 -retainids -addintrons -sort ${workingdir}/TrainingGene/Workingdir/tmp/merge2_with_UTRs.gff3 | \
    awk 'BEGIN {FS="\t"; OFS="\t"} {
        if ($3 == "intron") {
            $2="genometools"
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            $9 = parts[1]".intronOn"NR";"parts[1]
            sub(/^Parent=/, "ID=", $9)
            print
            }
        }
        else {
            print
        }
    }' > ${workingdir}/TrainingGene/Workingdir/tmp/merge3_with_UTRs.gff3
    mv ${workingdir}/TrainingGene/Workingdir/tmp/merge3_with_UTRs.gff3 ${workingdir}/TrainingGene/Workingdir/Merged/$(basename ${i} cdna.highQ_2genome_processed.gff3)highQ_with_UTRs.gff3
done 2>&1 | tee > ${workingdir}/log/mergeGFF3.log
awk '$1 != "warning:"' ${workingdir}/log/mergeGFF3.log > ${workingdir}/log/mergeGFF3.log2

## Same but for weaker hints
for i in $(ls ${workingdir}/Hints/Workingdir/Alignment_cDNA/*_processed.gff3); do
    awk '$3 == "gene"' ${i} | awk 'BEGIN {FS="\t"; OFS="\t"} 
    {
        if ($3 == "gene") {
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            }
            $9 = parts[1]
            sub(/\.gene$/, "", $9)
            sub(/^ID=/, "", $9)
        }
        print $1"\t"$4"\t"$5"\t"$9"\t"$6"\t"$7
    }' > ${workingdir}/Hints/Workingdir/tmp/cDNA.bed
    awk '$3 == "gene"' ${workingdir}/Hints/Workingdir/Alignment_CDS/$(basename ${i} cdna.weakerQ_2genome_processed.gff3)cds.weakerQ_2genome_processed.gff3 | awk 'BEGIN {FS="\t"; OFS="\t"} 
    {
        if ($3 == "gene") {
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            }
            $9 = parts[1]
            sub(/\.gene$/, "", $9)
            sub(/^ID=/, "", $9)
        }
        print $1"\t"$4"\t"$5"\t"$9"\t"$6"\t"$7
    }' > ${workingdir}/Hints/Workingdir/tmp/CDS.bed
    bedtools intersect -F 1.0 -a ${workingdir}/Hints/Workingdir/tmp/cDNA.bed -b ${workingdir}/Hints/Workingdir/tmp/CDS.bed -wa -wb -s | awk '$4 == $10 {print $4}' > ${workingdir}/Hints/Workingdir/tmp/both.txt
    echo '##gff-version 3' > ${workingdir}/Hints/Workingdir/tmp/merge.gff3
    grep -w -f ${workingdir}/Hints/Workingdir/tmp/both.txt ${i} >> ${workingdir}/Hints/Workingdir/tmp/merge.gff3
    awk '$3 == "CDS"' ${workingdir}/Hints/Workingdir/Alignment_CDS/$(basename ${i} cdna.weakerQ_2genome_processed.gff3)cds.weakerQ_2genome_processed.gff3 | grep -w -f ${workingdir}/Hints/Workingdir/tmp/both.txt >> ${workingdir}/Hints/Workingdir/tmp/merge.gff3
    gt gff3 -sortlines -tidy -retainids ${workingdir}/Hints/Workingdir/tmp/merge.gff3 | gt gff3 -retainids -sort > ${workingdir}/Hints/Workingdir/tmp/merge2.gff3
    python3 ${scriptdir}/addUTRs.py ${workingdir}/Hints/Workingdir/tmp/merge2.gff3 ${workingdir}/Hints/Workingdir/tmp/merge2_with_UTRs.gff3
    gt gff3 -retainids -addintrons -sort ${workingdir}/Hints/Workingdir/tmp/merge2_with_UTRs.gff3 | \
    awk 'BEGIN {FS="\t"; OFS="\t"} {
        if ($3 == "intron") {
            $2="genometools"
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            $9 = parts[1]".intronOn"NR";"parts[1]
            sub(/^Parent=/, "ID=", $9)
            print
            }
        }
        else {
            print
        }
    }' > ${workingdir}/Hints/Workingdir/tmp/merge3_with_UTRs.gff3
    mv ${workingdir}/Hints/Workingdir/tmp/merge3_with_UTRs.gff3 ${workingdir}/Hints/Workingdir/Merged/$(basename ${i} cdna.weakerQ_2genome_processed.gff3)weakerQ_with_UTRs.gff3
done 2>&1 | tee > ${workingdir}/log/mergeGFF3weaker.log
awk '$1 != "warning:"' ${workingdir}/log/mergeGFF3weaker.log > ${workingdir}/log/mergeGFF3weaker.log2








# [ Section 3 ]
## Getting rid of redundant protein sequences in different transcriptome using cd-hit, protein with sequence identity >= 80% are clustered and a representative sequence is given
## Note: non-redundancy is important in optimising augustus, the criteria given in the link below is >= 80% identical on aa level
### https://vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/training.html
rm ${workingdir}/TrainingGene/Workingdir/tmp/*
for i in $(ls ${TRANSLATION_OUTPUT}/*.pep.all.fa); do
    seqtk subseq -l 60 ${i} <(cut -f2 ${workingdir}/TrainingGene/Workingdir/highQ.filtered.transcripts.tsv) >> ${workingdir}/TrainingGene/Workingdir/tmp/transcriptome.pep.highQ.fa
done

### Change -T to ${PBS_NCPUS} if running on PBS
cd-hit -i ${workingdir}/TrainingGene/Workingdir/tmp/transcriptome.pep.highQ.fa -o ${workingdir}/TrainingGene/Workingdir/tmp/transcriptome.pep.highQ.representative.fa -c 0.8 -aS 0.9 -g 1 -d 0 -n 3 -M 75000 -T ${PBS_NCPUS}

## same but weaker hints
rm ${workingdir}/Hints/Workingdir/tmp/*
for i in $(ls ${TRANSLATION_OUTPUT}/*.pep.all.fa); do
    seqtk subseq -l 60 ${i} <(cut -f2 ${workingdir}/Hints/Workingdir/weakerQ.filtered.transcripts.tsv) >> ${workingdir}/Hints/Workingdir/tmp/transcriptome.pep.weakerQ.fa
done
cd-hit -i ${workingdir}/Hints/Workingdir/tmp/transcriptome.pep.weakerQ.fa -o ${workingdir}/Hints/Workingdir/tmp/transcriptome.pep.weakerQ.representative.fa -c 0.8 -aS 0.9 -g 1 -d 0 -n 3 -M 75000 -T ${PBS_NCPUS}



seqtk seq -C ${workingdir}/TrainingGene/Workingdir/tmp/transcriptome.pep.highQ.representative.fa | grep ">" | sed 's/>//g' | sort > ${workingdir}/TrainingGene/Workingdir/tmp/transcriptome.pep.highQ.representative.list
cp ${workingdir}/TrainingGene/Workingdir/tmp/transcriptome.pep.highQ.representative.list ${workingdir}/TrainingGene

## Same but weaker hints
seqtk seq -C ${workingdir}/Hints/Workingdir/tmp/transcriptome.pep.weakerQ.representative.fa | grep ">" | sed 's/>//g' | sort > ${workingdir}/Hints/Workingdir/tmp/transcriptome.pep.weakerQ.representative.list
cp ${workingdir}/Hints/Workingdir/tmp/transcriptome.pep.weakerQ.representative.list ${workingdir}/Hints



### grep these out from gff3
grep --no-filename -w -f ${workingdir}/TrainingGene/Workingdir/tmp/transcriptome.pep.highQ.representative.list ${workingdir}/TrainingGene/Workingdir/Merged/*highQ_with_UTRs.gff3 > ${workingdir}/TrainingGene/Training_gene.gff3

rm ${workingdir}/TrainingGene/Workingdir/tmp/*
### Find overlapping genes, get rid of all of them. 1kbp are added to either side of the gene so the 2kp flanking region generated from these gene later on will have no known coding region
awk '$3 == "gene"' ${workingdir}/TrainingGene/Training_gene.gff3 | awk 'BEGIN {FS="\t"; OFS="\t"} 
{
    if ($3 == "gene") {
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
        }
        $9 = parts[1]
        sub(/\.gene$/, "", $9)
        sub(/^ID=/, "", $9)
    }
    print $1"\t"$4-1000"\t"$5+1000"\t"$9"\t"$6"\t"$7
}' | awk 'BEGIN {FS="\t"; OFS="\t"} {if ($2 <= 0) {$2 = "1"} {print}}' > ${workingdir}/TrainingGene/Workingdir/tmp/Training_gene.bed

bedtools intersect -a ${workingdir}/TrainingGene/Workingdir/tmp/Training_gene.bed -b ${workingdir}/TrainingGene/Workingdir/tmp/Training_gene.bed -wa -wb | \
awk '$4 != $10 {print $10}' | sort | uniq > ${workingdir}/TrainingGene/Workingdir/tmp/overlapping.list
grep -w -v -f ${workingdir}/TrainingGene/Workingdir/tmp/overlapping.list ${workingdir}/TrainingGene/Training_gene.gff3 > ${workingdir}/TrainingGene/Training_gene_NoOverlap.gff3

### Now Process the general gff3 for into gff3 hints format, We will not be getting rid of overlap and instead will provide these and let augustus decide which is better
awk 'BEGIN {FS="\t"; OFS="\t"} {
        if ($3 == "intron") {
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            }
            $9 = parts[2]";pri=4;src=E"
            sub(/^Parent=/, "grp=", $9)
            $6 = "10"
            print
        }
        else if ($3 == "exon") {
            $3="exonpart"
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            }
            $9 = parts[2]";pri=4;src=E"
            sub(/^Parent=/, "grp=", $9)
            $6 = "10"
            print
            }
        else if ($3 == "five_prime_UTR" || $3 == "three_prime_UTR") {
            $3="UTRpart"
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            }
            $9 = parts[2]";pri=4;src=E"
            sub(/^Parent=/, "grp=", $9)
            $6 = "10"
            print
        }
    }' ${workingdir}/TrainingGene/Training_gene.gff3 > ${workingdir}/Hints/Hints.gff3

## same for weaker hints
grep --no-filename -w -f ${workingdir}/Hints/transcriptome.pep.weakerQ.representative.list ${workingdir}/Hints/Workingdir/Merged/*weakerQ_with_UTRs.gff3 | \
awk 'BEGIN {FS="\t"; OFS="\t"} {
        if ($3 == "intron") {
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            }
            $9 = parts[2]";pri=3;src=E"
            sub(/^Parent=/, "grp=", $9)
            $6 = "5"
            print
        }
        else if ($3 == "exon") {
            $3="exonpart"
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            }
            $9 = parts[2]";pri=3;src=E"
            sub(/^Parent=/, "grp=", $9)
            $6 = "5"
            print
            }
        else if ($3 == "five_prime_UTR" || $3 == "three_prime_UTR") {
            $3="UTRpart"
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            }
            $9 = parts[2]";pri=3;src=E"
            sub(/^Parent=/, "grp=", $9)
            $6 = "5"
            print
        }
    }' >> ${workingdir}/Hints/Hints.gff3



# Continue with TrainingAugustus.sh
## Shared variables are:
### ${workingdir}
### ${genome}
