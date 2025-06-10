diamond makedb --threads ${PBS_NCPUS} --db uniprot_sprot.diamond.db --in uniprot_sprot.fasta.gz

diamond makedb --threads ${PBS_NCPUS} --db uniprot_trembl.diamond.db --in uniprot_trembl.fasta.gz

# Database in Hardip's directory here
## /g/data/xl04/hrp561/basdurnaseq/uniprot_sprot.diamond.db.dmnd
## /g/data/xl04/hrp561/basdurnaseq/uniprot_trembl.diamond.db.dmnd
