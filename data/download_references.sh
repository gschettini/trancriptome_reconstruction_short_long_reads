#!/bin/bash

# Download reference files - Bos taurus (OPTIONAL)

#Ensembl fasta
cd /reference/fasta && ftp https://ftp.ensembl.org/pub/current_fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa.gz && gunzip *.gz

#Ensembl annotation
cd /reference/annotation && ftp https://ftp.ensembl.org/pub/current_gtf/bos_taurus/Bos_taurus.ARS-UCD1.3.111.gtf.gz && gunzip *.gz

#NCBI fasta
cd /reference/fasta && ftp https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/all_assembly_versions/GCF_002263795.2_ARS-UCD1.3/GCF_002263795.2_ARS-UCD1.3_genomic.fna.gz && gunzip *.gz

#NCBI annotation
cd /reference/annotation && ftp https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/all_assembly_versions/GCF_002263795.2_ARS-UCD1.3/GCF_002263795.2_ARS-UCD1.3_genomic.gff.gz && gunzip *.gz

#Formatting to have the same chromosome name (ENSEMBL and NCBI)
cd /reference/annotation && wget https://hgdownload.soe.ucsc.edu/hubs/GCF/002/263/795/GCF_002263795.2/GCF_002263795.2.chromAlias.txt && awk '{print $1, $2}' GCF_002263795.2.chromAlias.txt | sed '1d' > GCF_002263795.2.chromAlias_2.txt && awk 'FNR==NR { patterns[$1]=$2; next } $1 in patterns { sub($1, patterns[$1]) } 1' GCF_002263795.2.chromAlias_2.txt GCF_002263795.2_ARS-UCD1.3_genomic.gff > NCBI_GCF_002263795.2_ARS-UCD1.3_genomic.gff

cd /reference/fasta && awk 'FNR==NR { patterns[$1]=$2; next } /^>/ { for (pattern in patterns) gsub(pattern, patterns[pattern]) } 1' /reference/annotation/GCF_002263795.2.chromAlias_2.txt GCF_002263795.2_ARS-UCD1.3_genomic.fna > NCBI_GCF_002263795.2_ARS-UCD1.3_genomic.fna

#NONCODE v5 (lncRNA database)
cd /reference/annotation && wget http://v5.noncode.org/datadownload/NONCODEv5_bosTau6.lncAndGene.bed.gz && gunzip *.gz && bedToGenePred NONCODEv5_bosTau6.lncAndGene.bed genepred.gp && genePredToGtf file genepred.gp NONCODEv5_bosTau6.lncAndGene.gtf && sed 's/^chr//' NONCODEv5_bosTau6.lncAndGene.gtf > a_NONCODEv5_bosTau6.lncAndGene.gtf && rm genepred.gp NONCODEv5_bosTau6.lncAndGene.gtf && mv a_NONCODEv5_bosTau6.lncAndGene.gtf NONCODEv5_bosTau6.lncAndGene.gtf

#non-redudant protein NCBI database
cd /reference/others/nr_db && wget https://ftp.ncbi.nlm.nih.gov/blast/db/v5/FASTA/nr.gz && wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && tar -xzf taxdump.tar.gz && rm taxdump.tar.gz && wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz

#RNAsamba training file
cd /reference/others && curl -O https://raw.githubusercontent.com/apcamargo/RNAsamba/master/data/partial_length_weights.hdf5