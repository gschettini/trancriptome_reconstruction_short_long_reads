#!/bin/bash
folder_a=/test/set_a
folder_b=/test/set_b
n_cores=30
memory_ram=150
genome_size=2.8g
path_to_index=/reference/fasta/hisat/index/Bos_taurus.ARS-UCD1.3.111
splice_sites=/reference/fasta/hisat/splice_sites/splice_sites.txt
reference_genome_ensembl=/reference/fasta/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa
gtf_ensembl=/reference/annotation/Bos_taurus.ARS-UCD1.3.111.gtf
reference_genome_ncbi=/reference/fasta/NCBI_GCF_002263795.2_ARS-UCD1.3_genomic.fna
gff_ncbi=/reference/annotation/NCBI_GCF_002263795.2_ARS-UCD1.3_genomic.gff
gtf_lncRNA=/reference/annotation/NONCODEv5_bosTau6.lncAndGene.gtf
reference_GMAP_ensembl=/reference/GMAP/ensembl
reference_GMAP_ncbi=/reference/GMAP/ncbi
adapters_SR_a=/Trimmomatic-0.39/adapters/NexteraPE-PE.fa
adapters_SR_b=/Trimmomatic-0.39/adapters/NexteraPE-PE.fa
organism_1=9913,9915,30522
organism_1_e_value=0.000001
organism_1_perc_id=90
organism_1_cov=90
organism_2=9606,10090
organism_2_e_value=0.000001
organism_2_perc_id=70
organism_2_cov=90
organism_3=40674
organism_3_e_value=0.000001
organism_3_perc_id=70
organism_3_cov=90

#Make all necessary directories##
if [ $folder_a != "/test/set_a" ]; then
  echo "Making directories in the new folder_a"
  mkdir $folder_a
  cd $folder_a
  mkdir SR_fastq && mkdir SR_fastq/trimming && mkdir SR_fastq/alignment && mkdir SR_fastq/preprocessed
  mkdir LR_nanopore && mkdir LR_nanopore/canu
  mkdir de_novo_assembly
  mkdir annotation_GMAP && mkdir annotation_GMAP/ensembl && mkdir annotation_GMAP/ncbi
  mkdir counting
  mkdir comparisons && mkdir comparisons/ensembl && mkdir comparisons/ncbi && mkdir comparisons/lncRNA
  mkdir functional_characterization
  mkdir R_outputs
  cp /test/set_a/files_a.yaml $folder_a/
else
  echo "Default folder_a set."
fi

if [ -z "$folder_b" ]; then
  echo "Folder_b was not set"
fi

if [[ ( "$folder_b" != "/test/set_b" ) && ( "$folder_b" != "" ) ]]; then
  echo "Making directories in the new folder_b"
  mkdir $folder_b
  cd $folder_b || exit 1
  mkdir SR_fastq && mkdir SR_fastq/trimming && mkdir SR_fastq/alignment && mkdir SR_fastq/preprocessed
  mkdir LR_nanopore && mkdir LR_nanopore/canu
  mkdir de_novo_assembly
  mkdir annotation_GMAP && mkdir annotation_GMAP/ensembl && mkdir annotation_GMAP/ncbi
  mkdir counting
  mkdir comparisons && mkdir comparisons/ensembl && mkdir comparisons/ncbi && mkdir comparisons/lncRNA
  mkdir functional_characterization
  mkdir R_outputs
  cp /test/set_b/files_b.yaml $folder_b/
fi

if [[ ("$folder_b" == "/test/set_b" ) ]]; then
  echo "Default folder_b set."
fi

sed -i "s|/test/set_a|$folder_a|g" $folder_a/files_a.yaml
sed -i "s|/test/set_b|$folder_b|g" $folder_b/files_b.yaml
