#!/bin/bash
source ./settings.sh

##Build Hisat2 splice_sites file##
if test -f $splice_sites; then
  echo "Hisat2 splice_sites file already exist."
  echo "Skipping the Hisat2 splice_sites building."
else
  echo "Hisat2 splice_sites file does not exist."
  echo "Building Hisat2 splice_sites now..."
  python3.7 /hisat2-2.2.1/hisat2_extract_splice_sites.py $gtf_ensembl > $splice_sites
fi

##Build Hisat2 index##
if [ -n "$(ls -A /reference/fasta/hisat/index/*.ht2 2>/dev/null)" ]; then
  echo "Hisat2 index file already exist."
  echo "Skipping Hisat2 index building."
else
  echo "Hisat2 index file does not exist."
  echo "Building Hisat2 index now..."
  /hisat2-2.2.1/hisat2-build -p $(echo $n_cores) --ss $splice_sites $reference_genome_ensembl $path_to_index
fi


##Trimming - set A##
for sample in $(ls $folder_a/SR_fastq/*.fastq* | awk -F '_R[12]' '{print $1}' | sort | uniq); do

read1=$(ls $folder_a/SR_fastq/*.fastq* | grep $sample | grep "_R1")
read2=$(ls $folder_a/SR_fastq/*.fastq* | grep $sample | grep "_R2")

read1_file=$(basename "$read1" | cut -f 1 -d '.')
read2_file=$(basename "$read2" | cut -f 1 -d '.')
extension="${read1#*.}"

echo 'Trimming from set A now'
echo ${read1_file}
echo ${read2_file}

output_R1=$folder_a/SR_fastq/trimming/$(echo $read1_file).trimmed.$(echo $extension)
output_R1_un=$folder_a/SR_fastq/trimming/$(echo $read1_file).untrimmed.$(echo $extension)
output_R2=$folder_a/SR_fastq/trimming/$(echo $read2_file).trimmed.$(echo $extension)
output_R2_un=$folder_a/SR_fastq/trimming/$(echo $read2_file).untrimmed.$(echo $extension)

java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads $(echo $n_cores) $read1 $read2 $output_R1 $output_R1_un $output_R2 $output_R2_un ILLUMINACLIP:$adapters_SR_a:5:10:10 SLIDINGWINDOW:3:25 MINLEN:40 AVGQUAL:35
done;

rm $folder_a/SR_fastq/trimming/*.untrimmed.fastq*

##Alignment - set A##

cd $folder_a/SR_fastq/alignment

for sample in $(ls $folder_a/SR_fastq/trimming/*.trimmed.fastq* | awk -F '_R[12]' '{print $1}' | sort | uniq); do

read1=$(ls $folder_a/SR_fastq/trimming/*.trimmed.fastq*| grep $sample | grep "_R1")
read2=$(ls $folder_a/SR_fastq/trimming/*.trimmed.fastq* | grep $sample | grep "_R2")

read1_file=$(basename "$read1" | cut -f 1 -d '.' | sed "s/_R1//g" )

echo $read1
echo $read2
echo 'Aligning from set A now'

/hisat2-2.2.1/hisat2  -p $(echo $n_cores) -k 1 -x $path_to_index -1 $read1 -2 $read2 --bowtie2-dp 2 --score-min L,0,-1 --dta --known-splicesite-infile $splice_sites  -S $(echo $read1_file)_alignment.sam --no-unal --summary-file $(echo $read1_file)_alignment_summary.txt
done;

##Preprocessing - set A##

echo "Preprocessing SR from set A now"

cd $folder_a/SR_fastq/alignment

for sample in $(ls $folder_a/SR_fastq/alignment/*.sam | awk -F '_alignment' '{print $1}' | sort); do

sample_file=$(basename "$sample" | cut -f 1 -d '.')

echo ${sample_file}

samtools sort -@ $(echo $n_cores) -O BAM -l 9 -o $(echo ${sample_file})_alignment.bam $(echo ${sample_file})_alignment.sam
samtools index $(echo ${sample_file})_alignment.bam
samtools view --threads $(echo $n_cores) -b -h -F 1796 -o $(echo ${sample_file})_alignment.filtered.bam $(echo ${sample_file})_alignment.bam
samtools sort -@ $(echo $n_cores) -O BAM -l 9 -o $(echo ${sample_file})_alignment.filtered.sorted.bam $(echo ${sample_file})_alignment.filtered.bam
/biobambam2/bin/bammarkduplicates I=$(echo ${sample_file})_alignment.filtered.sorted.bam O=$(echo ${sample_file})_alignment.filtered.sorted.undup.bam markthreads=$(echo $n_cores) level=9 rmdup=1 index=1 dupindex=0 verbose=0
/biobambam2/bin/bamtofastq filename=$(echo ${sample_file})_alignment.filtered.sorted.undup.bam F=$folder_a/SR_fastq/preprocessed/$(echo ${sample_file})_R1.fastq F2=$folder_a/SR_fastq/preprocessed/$(echo ${sample_file})_R2.fastq

done;

cd $folder_a/SR_fastq/preprocessed

for sample in $(ls $folder_a/SR_fastq/preprocessed/*R*.fastq | awk -F '_R[12]' '{print $1}' | sort | uniq); do


read1=$(ls $folder_a/SR_fastq/preprocessed/*R*.fastq | grep $sample | grep "_R1.fastq")
read2=$(ls $folder_a/SR_fastq/preprocessed/*R*.fastq | grep $sample | grep "_R2.fastq")

read1_file=$(basename "$read1" | cut -f 1 -d '.')
read2_file=$(basename "$read2" | cut -f 1 -d '.')
sample_name=$(echo $read1_file | awk -F '_R1' '{print $1}' )

echo $sample_name

#Error correction / Removing contamination
/bbmap/bbduk.sh in1=$(echo ${read1_file}).fastq in2=$(echo ${read2_file}).fastq2 out1=$(echo ${read1_file}).trimmed.unmatched.fastq.gz out2=$(echo ${read2_file}).trimmed.unmatched.fastq.gz ref=/bbmap/resources/phix174_ill.ref.fa.gz,/bbmap/resources/sequencing_artifacts.fa.gz t=$(echo $n_cores) k=21 hdist=1 stats=stats_$sample_name.txt

#quality filtering
/bbmap/bbduk.sh in1=$(echo ${read1_file}).trimmed.unmatched.fastq.gz in2=$(echo ${read2_file}).trimmed.unmatched.fastq.gz out1=$(echo ${read1_file}).filtered.trimmed.fastq.gz out2=$(echo ${read2_file}).filtered.trimmed.fastq.gz qtrim=rl trimq=30 minlength=100 t=$(echo $n_cores) minbasequality=25

#normalization and error correction
/bbmap/bbnorm.sh in1=$(echo ${read1_file}).filtered.trimmed.fastq.gz in2=$(echo ${read2_file}).filtered.trimmed.fastq.gz out1=$(echo ${read1_file}).normalized.filtered.trimmed.fastq out2=$(echo ${read2_file}).normalized.filtered.trimmed.fastq target=30 maxdepth=30 mindepth=10 prefilter=t bits=16 prefilter ecc=t t=$(echo $n_cores)

done;

rm *_R1.fastq *_R2.fastq *.trimmed.unmatched.fastq.gz *_R1.filtered.trimmed.fastq *_R2.filtered.trimmed.fastq

##ONT long-read correction##
echo "Correcting LR - set A"
/canu/build/bin/canu -correct -p set_a -d $folder_a/LR_nanopore/canu/correction genomeSize=$(echo $genome_size) hapMemory=$(echo $memory_ram) hapThreads=$(echo $n_cores) cormhapMemory=$(echo $memory_ram) corThreads=$(echo $n_cores) batMemory=$(echo $memory_ram) batThreads=$(echo $n_cores) stopOnLowCoverage=0 minReadLength=200 minOverlapLength=200 corMhapSensitivity=high corOutCoverage=all corOverlapper=minimap minInputCoverage=0 -nanopore-raw $folder_a/LR_nanopore/*.fastq

##ONT long-read trimming##
echo "Trimming LR - set A"
/canu/build/bin/canu -trim -p set_a -d $folder_a/LR_nanopore/canu/trimming genomeSize=$(echo $genome_size) mmapMemory=$(echo $memory_ram) mmapThreads=$(echo $n_cores) obtmmapMemory=$(echo $memory_ram) batMemory=$(echo $memory_ram) batThreads=$(echo $n_cores) stopOnLowCoverage=0  minReadLength=200 minOverlapLength=200 utgOverlapper=minimap corOutCoverage=all utgMhapSensitivity=high minInputCoverage=0 -nanopore-corrected $folder_a/LR_nanopore/canu/correction/set_a.correctedReads.fasta.gz 

##de novo assembly - set A ##
echo "de novo assembly - set A"
/SPAdes-3.14.1-Linux/bin/spades.py --rna --dataset $folder_a/files_a.yaml -o $folder_a/de_novo_assembly/ -t $(echo $n_cores) -k 21, 33, 55, 77
seqkit stats $folder_a/de_novo_assembly/hard_filtered_transcripts.fasta > $folder_a/de_novo_assembly/hard_filtered_transcripts_summary.txt

##Comparisons ##
#Build index ensembl
if [ -n "$(ls -A /reference/GMAP/ensembl/*ensembl_index.coords 2>/dev/null)" ];  then
  echo "GMAP ensembl index file already exist."
  echo "Skipping GMAP ensembl index building."
else
  echo "GMAP ensembl index file does not exist."
  echo "Building GMAP ensembl index now..."
  gmap_build -D $reference_GMAP_ensembl $reference_genome_ensembl --genomedb=ensembl_index
fi

#Build index ncbi
if [ -n "$(ls -A /reference/GMAP/ncbi/*ncbi_index.coords 2>/dev/null)" ]; then
  echo "GMAP ncbi index file already exist."
  echo "Skipping GMAP ncbi index building."
else
  echo "GMAP ncbi index file does not exist."
  echo "Building GMAP ncbi index now..."
gmap_build -D $reference_GMAP_ncbi $reference_genome_ncbi --genomedb=ncbi_index
fi



##GMAP alignment - ENSEMBL (SET A)
cd $folder_a/annotation_GMAP/ensembl
gmap -D $reference_GMAP_ensembl/ensembl_index -d ensembl_index -t $(echo $n_cores) $folder_a/de_novo_assembly/hard_filtered_transcripts.fasta -f gff3_gene --microexon-spliceprob=1 --nofails --quality-protocol=illumina --suboptimal-score=0.99 --min-identity=0.90 --npaths=1 > set_a_ensembl.gff3
gmap -D $reference_GMAP_ensembl/ensembl_index -d ensembl_index -t $(echo $n_cores) $folder_a/de_novo_assembly/hard_filtered_transcripts.fasta -f samse --microexon-spliceprob=1 --nofails --quality-protocol=illumina --suboptimal-score=0.99 --min-identity=0.90 --npaths=1 > set_a_ensembl.sam

#   Filtering sequences with more than 10 mismatches (SET A)
awk ' $0 ~ /mismatches=[0-9];/ ' set_a_ensembl.gff3 > set_a_ensembl_filtered.gff3

##Filter - Rscript 1
Rscript /Rscripts/script1_a.R

rm $folder_a/annotation_GMAP/ensembl/set_a_ensembl.gff3 $folder_a/annotation_GMAP/ensembl/set_a_ensembl_filtered.gff3

##Compare our new gff file with ENSEMBL gff file (SET A)

cd $folder_a/comparisons/ensembl
/gffcompare/gffcompare -r $gtf_ensembl -T $folder_a/annotation_GMAP/ensembl/set_a_ensembl_filtered_2.gff3 -o "set_a_ensembl"

awk '/class_code "u";/{print}' set_a_ensembl.annotated.gtf > unknown_set_a_ensembl.annotated.gtf

#   total unknown
echo "Total unknown transcripts - ENSEMBL - SET A"
wc -l unknown_set_a_ensembl.annotated.gtf

##Filter unknown transcripts from assembly file (SET A)

cd $folder_a/comparisons/ensembl
awk '{print $10}' unknown_set_a_ensembl.annotated.gtf > flag_u_set_a_ensembl.txt
sed -e 's/\"//g' flag_u_set_a_ensembl.txt > flag_u_set_a_ensembl_2.txt
sed -e 's/\;//g' flag_u_set_a_ensembl_2.txt > flag_u_set_a_ensembl_3.txt
sed -e 's/\.mrna*//g' flag_u_set_a_ensembl_3.txt > flag_u_set_a_ensembl_4.txt
sed -e 's/\NODE/>NODE/' flag_u_set_a_ensembl_4.txt | uniq > flag_u_set_a_ensembl_5.txt
sed -e 's/..$//g' flag_u_set_a_ensembl_5.txt > flag_u_set_a_ensembl_6.txt
mv flag_u_set_a_ensembl_6.txt flag_u_set_a_ensembl2.txt
rm flag_u_set_a_ensembl_*.txt

cd $folder_a/de_novo_assembly
file="$folder_a/comparisons/ensembl/flag_u_set_a_ensembl2.txt"

parallel --jobs $(echo $n_cores) --eta '
  count={#}
  line={1}

  temp_file1="1_$count.fasta"
  temp_file2="2_$count.fasta"

  sed -n -e "/$line/,/>NODE/ p" hard_filtered_transcripts.fasta > "$temp_file1"
  if [ -s "$temp_file1" ]; then
    head -n -1 "$temp_file1" > "$temp_file2"
  fi' :::: "$file"

cat 2_*.fasta > flag_u_set_a_ensembl.fasta

rm 1_*.fasta
rm 2_*.fasta

##GMAP alignment - NCBI (SET A)

cd $folder_a/annotation_GMAP/ncbi
gmap -D $reference_GMAP_ncbi/ncbi_index -d ncbi_index -t $(echo $n_cores) $folder_a/de_novo_assembly/flag_u_set_a_ensembl.fasta -f gff3_gene --microexon-spliceprob=1 --nofails --quality-protocol=illumina --suboptimal-score=0.99 --min-identity=0.90 --npaths=1 > set_a_ncbi.gff3
gmap -D $reference_GMAP_ncbi/ncbi_index -d ncbi_index -t $(echo $n_cores) $folder_a/de_novo_assembly/flag_u_set_a_ensembl.fasta -f samse --microexon-spliceprob=1 --nofails --quality-protocol=illumina --suboptimal-score=0.99 --min-identity=0.90 --npaths=1 > set_a_ncbi.sam

#   Filtering sequences with more than 10 mismatches (SET A)
awk ' $0 ~ /mismatches=[0-9];/ ' set_a_ncbi.gff3 > set_a_ncbi_filtered.gff3

##Filter - Rscript 2
Rscript /Rscripts/script2_a.R

##Compare our new gff file with NCBI gff file (SET A)

cd $folder_a/comparisons/ncbi
/gffcompare/gffcompare -r $gff_ncbi -T $folder_a/annotation_GMAP/ncbi/set_a_ncbi_filtered_2.gff3 -o "set_a_ncbi"

awk '/class_code "u";/{print}' set_a_ncbi.annotated.gtf > unknown_set_a_ncbi.annotated.gtf

#   total unknown
echo "Total unknown transcripts - NCBI - SET A"
wc -l unknown_set_a_ncbi.annotated.gtf

##Filter unknown transcripts from assembly file (SET A)

cd $folder_a/comparisons/ncbi
awk '{print $10}' unknown_set_a_ncbi.annotated.gtf > flag_u_set_a_ncbi.txt
sed -e 's/\"//g' flag_u_set_a_ncbi.txt > flag_u_set_a_ncbi_2.txt
sed -e 's/\;//g' flag_u_set_a_ncbi_2.txt > flag_u_set_a_ncbi_3.txt
sed -e 's/\.mrna*//g' flag_u_set_a_ncbi_3.txt > flag_u_set_a_ncbi_4.txt
sed -e 's/\NODE/>NODE/' flag_u_set_a_ncbi_4.txt | uniq > flag_u_set_a_ncbi_5.txt
sed -e 's/..$//g' flag_u_set_a_ncbi_5.txt > flag_u_set_a_ncbi_6.txt
mv flag_u_set_a_ncbi_6.txt flag_u_set_a_ncbi2.txt
rm flag_u_set_a_ncbi_*.txt

cd $folder_a/de_novo_assembly
file="$folder_a/comparisons/ncbi/flag_u_set_a_ncbi2.txt"

parallel --jobs $(echo $n_cores) --eta '
  count={#}
  line={1}

  temp_file1="1_$count.fasta"
  temp_file2="2_$count.fasta"

  sed -n -e "/$line/,/>NODE/ p" flag_u_set_a_ensembl.fasta > "$temp_file1"
  if [ -s "$temp_file1" ]; then
    head -n -1 "$temp_file1" > "$temp_file2"
  fi' :::: "$file"

cat 2_*.fasta > flag_u_set_a_ncbi.fasta

rm 1_*.fasta
rm 2_*.fasta

##Filter - Rscript 3
Rscript /Rscripts/script3_a.R

##Formating new unknown gtf file before compare with lncRNA database (NONCODEv5) - Set A##

cd $folder_a/comparisons/ncbi
sed -e 's/\NODE/"NODE/g' set_a_ncbi_filtered_3.gff3 > set_a_ncbi_filtered_3_2.gtf
sed -e 's/\; xloc/"; xloc/g' set_a_ncbi_filtered_3_2.gtf > set_a_ncbi_filtered_3_3.gtf
sed -e 's/\; gene_id/"; gene_id/g' set_a_ncbi_filtered_3_3.gtf > set_a_ncbi_filtered_3_4.gtf
sed -e 's/\; exon_number/"; exon_number/g' set_a_ncbi_filtered_3_4.gtf > set_a_ncbi_filtered_3_5.gtf
sed -e 's/\; class_code/"; class_code/g' set_a_ncbi_filtered_3_5.gtf > set_a_ncbi_filtered_3_6.gtf
mv set_a_ncbi_filtered_3_6.gtf set_a_ncbi_filtered_4.gtf
rm set_a_ncbi_filtered_3_*.gtf

##Compare our gff file with lncRNA database (NONCODEv5) - Set A##

cd $folder_a/comparisons/lncRNA

/gffcompare/gffcompare -r $gtf_lncRNA -T $folder_a/comparisons/ncbi/set_a_ncbi_filtered_4.gtf -o "set_a_lncRNA"

awk '/class_code "u";/{print}' set_a_lncRNA.annotated.gtf > unknown_set_a_lncRNA.annotated.gtf
#   total unknown
wc -l unknown_set_a_lncRNA.annotated.gtf

##Filter - Rscript 4
Rscript /Rscripts/script4_a.R

##Putting quotes before filter - Set A##

cd $folder_a/comparisons/lncRNA
sed -e 's/\NODE/"NODE/g' set_a_lncRNA_filtered.gff3 > set_a_lncRNA_filtered_2.gff3
sed -e 's/\; xloc /"; xloc "/g' set_a_lncRNA_filtered_2.gff3 > set_a_lncRNA_filtered_3.gff3
sed -e 's/\; gene_id/"; gene_id/g' set_a_lncRNA_filtered_3.gff3 > set_a_lncRNA_filtered_4.gff3
sed -E 's/; exon_number ([[:alnum:]]+)/"; exon_number "\1"/g' set_a_lncRNA_filtered_4.gff3 > set_a_lncRNA_filtered_5.gff3
sed -e 's/\; class_code u/"; class_code "u"/g' set_a_lncRNA_filtered_5.gff3 > set_a_lncRNA_filtered_6.gff3
sed -E 's/; tss_id TSS([[:alnum:]]+)/; tss_id "TSS\1"/g' set_a_lncRNA_filtered_6.gff3 > set_a_lncRNA_filtered_7.gff3

mv set_a_lncRNA_filtered_7.gff3 set_a_all_flags_u.gff3
rm set_a_lncRNA_filtered_*.gff3

##Filtering GTF file (genes with >300bp & >1 exon) - Set A##
cd $folder_a/comparisons/lncRNA
#remove transcripts less than 300bp
/gffread/gffread -l 300 set_a_all_flags_u.gff3 -o set_a_all_flags_u_2.gff3
#discard single exons
/gffread/gffread -U -T set_a_all_flags_u_2.gff3 -o set_a_all_flags_u_3.gtf

##Filter (Upstream and Downstream < 5kb known genes)- Rscript 5
Rscript /Rscripts/script5_a.R

## Putting quotes before concatenate gtf - Set A##
cd $folder_a/comparisons/lncRNA

sed -e 's/\NODE/"NODE/g' set_a_all_flags_u_clean.gtf > set_a_all_flags_u_clean_2.gtf
sed -e 's/\; gene_id/"; gene_id/g' set_a_all_flags_u_clean_2.gtf > set_a_all_flags_u_clean_3.gtf
sed -E 's/.path([[:alnum:]]+)/.path\1"/g' set_a_all_flags_u_clean_3.gtf > set_a_all_flags_u_clean_4.gtf

mv set_a_all_flags_u_clean_4.gtf set_a_all_flags_u_clean2.gtf
rm set_a_all_flags_u_clean_*.gtf

##Concatenating gtf - Set A##
cd $folder_a/comparisons/lncRNA

#concatenating gtf overlapping exon boundaries (SET A)
/gffread/gffread -MKY --cset -T set_a_all_flags_u_clean2.gtf -o set_a_all_flags_u_clean_concatenated.gtf -d info_dup_set_a 
#discard single exons (SET A)
/gffread/gffread -U -T set_a_all_flags_u_clean_concatenated.gtf -o set_a_all_flags_u_clean_concatenated2.gtf 

##Quality filtering (remove transcripts without count) - Set A##

cd $folder_a/SR_fastq/alignment

for sample in $(ls $folder_a/SR_fastq/alignment/*_alignment.filtered.sorted.undup.bam | sort); do
echo $sample
sample_name=$(basename "$sample" | cut -f 1 -d '.')
sample_name=$(echo $sample_name | awk -F '_alignment' '{print $1}' )
echo $sample_name
/subread-2.0.6-Linux-x86_64/bin/featureCounts -s 1 -a $folder_a/comparisons/lncRNA/set_a_all_flags_u_clean_concatenated2.gtf -o $folder_a/counting/filter_$(echo $sample_name).count -F 'GTF' -t 'exon' -g 'gene_id' --ignoreDup -p -T $(echo $n_cores) $sample
done;

##Filter (transcripts not counted)- Rscript 6
Rscript /Rscripts/script6_a.R

#Set A
/gffread/gffread $folder_a/comparisons/lncRNA/set_a_all_flags_u_clean_concatenated2.gtf --nids $folder_a/R_outputs/not_counted_transcripts_set_a.txt -T -o $folder_a/comparisons/lncRNA/set_a_all_flags_u_final.gtf

##Merging close transcripts##

#Set A
cd $folder_a/comparisons/lncRNA/

/gffread/gffread set_a_all_flags_u_final.gtf -o set_a_all_flags_u_final2.gff
awk '{ if ($3 == "transcript") print $0 }' set_a_all_flags_u_final2.gff > set_a_all_flags_u_final2_transcripts.gff
/gffread/gffread set_a_all_flags_u_final2_transcripts.gff --bed -o set_a_all_flags_u_final2_transcripts.bed
sort -k1,1 -k2,2n set_a_all_flags_u_final2_transcripts.bed > set_a_all_flags_u_final2_transcripts2.bed
bedtools merge -i set_a_all_flags_u_final2_transcripts2.bed -header -s -d 300 -c 4,5,6 -o collapse,sum,distinct -delim "+" > merged_set_a_all_flags_u_final2_transcripts2.bed
bedToGenePred merged_set_a_all_flags_u_final2_transcripts2.bed genepred.gp
genePredToGtf file genepred.gp merged_set_a_all_flags_u_final2_transcripts2.gtf
/gffread/gffread merged_set_a_all_flags_u_final2_transcripts2.gtf -MQK --cset -T -o  set_a_all_flags_u_final3.gtf 
awk '{ if ($3 == "transcript") print $0 }' set_a_all_flags_u_final3.gtf  > only_transcripts_set_a_all_flags_u_final3.gtf

rm set_a_all_flags_u_final2.gff set_a_all_flags_u_final2_transcripts.gff set_a_all_flags_u_final2_transcripts.bed set_a_all_flags_u_final2_transcripts2.bed merged_set_a_all_flags_u_final2_transcripts2.bed genepred.gp merged_set_a_all_flags_u_final2_transcripts2.gtf

##Get FASTA sequences from final gtf coordinates##

#Set A
cd $folder_a/comparisons/lncRNA/
/gffread/gffread -w set_a_all_flags_u_final3.fasta -g $reference_genome_ensembl set_a_all_flags_u_final3.gtf

##Identifying orthologs - (Humans + Mice + Mammals) Possible family gene member (Bos taurus + Bos indicus + Bos taurusXBos Indicus)

if test -f /reference/others/nr_db/nr_db.dmnd; then
  echo "Diamond database already exist."
  echo "Skipping diamond database building."
else
  echo "Diamond database does not exist."
  echo "Building Diamond database now..."
  cd /reference/others/nr_db/
  diamond makedb --in nr.gz --db /reference/others/nr_db/nr_db --taxonmap prot.accession2taxid.FULL.gz --taxonnodes nodes.dmp --taxonnames names.dmp
fi

#Set A
cd $folder_a/functional_characterization
 
#Cattle or Organism 1
if [[ -v organism_1 ]];  then  
  echo "Organism 1 was set:"
  echo $organism_1
  diamond blastx -d /reference/others/nr_db/nr_db -q $folder_a/comparisons/lncRNA/set_a_all_flags_u_final3.fasta --taxonlist $organism_1 -e $organism_1_e_value --subject-cover $organism_1_cov -o set_a_nr_1.tsv -p $n_cores -k 1 -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids qcovhsp stitle --unal 1 --id $organism_1_perc_id
  awk '{print $1}' set_a_nr_1.tsv > set_a_nr_1_list_no_match.tsv
  awk '$2 ~ /^[[:alnum:]]+/' set_a_nr_1.tsv > set_a_nr_1_match.tsv
  /seqtk/seqtk subseq $folder_a/comparisons/lncRNA/set_a_all_flags_u_final3.fasta set_a_nr_1_list_no_match.tsv > set_a_all_flags_u_no_match.fasta 
else  
  echo "Organism 1 was not set. Skipping..."  
fi

#Human and Mouse | or Organism 2
if [[ -v organism_2 ]];  then  
  echo "Organism 2 was set:"
  echo $organism_2
  diamond blastx -d /reference/others/nr_db/nr_db -q set_a_all_flags_u_no_match.fasta --taxonlist $organism_2 -e $organism_2_e_value --subject-cover $organism_2_cov -o set_a_nr_2.tsv -p $n_cores -k 1 -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids qcovhsp stitle --unal 1 --id $organism_2_perc_id
  awk '{print $1}' set_a_nr_2.tsv > set_a_nr_2_list_no_match.tsv
  awk '$2 ~ /^[[:alnum:]]+/' set_a_nr_2.tsv > set_a_nr_2_match.tsv
  /seqtk/seqtk subseq set_a_all_flags_u_no_match.fasta set_a_nr_2_list_no_match.tsv > set_a_all_flags_u_no_match2.fasta
else  
  echo "Organism 2 was not set. Skipping..."  
fi

#Mammals or Organism 3
if [[ -v organism_3 ]];  then  
  echo "Organism 3 was set:"
  echo $organism_3
  diamond blastx -d /reference/others/nr_db/nr_db -q set_a_all_flags_u_no_match2.fasta --taxonlist $organism_3 -e $organism_3_e_value --subject-cover $organism_3_cov -o set_a_nr_3.tsv -p $n_cores -k 1 -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids qcovhsp stitle --unal 1 --id $organism_3_perc_id
  awk '{print $1}' set_a_nr_3.tsv > set_a_nr_3_list_no_match.tsv
  awk '$2 ~ /^[[:alnum:]]+/' set_a_nr_3.tsv > set_a_nr_3_match.tsv
  /seqtk/seqtk subseq set_a_all_flags_u_no_match2.fasta set_a_nr_3_list_no_match.tsv > set_a_all_flags_u_no_match3.fasta
else  
  echo "Organism 3 was not set. Skipping..."  
fi
 
##Filtering ortologous proteins## - Rscript7
Rscript /Rscripts/script7_a.R

##Identifying coding vs noncoding potential - RNAsamba##

#Set A

cd $folder_a/functional_characterization

awk '{print $1}' all_hits_set_a.txt | sed -e 's/\"//g' > protein_list.csv
if test -f $folder_a/functional_characterization/set_a_all_flags_u_no_match3.fasta; then
  seqkit grep -iv -f protein_list.csv $folder_a/functional_characterization/set_a_all_flags_u_no_match3.fasta --quiet > set_a_all_flags_u_no_orthologous.fasta
fi
if [[ ! -e $folder_a/functional_characterization/set_a_all_flags_u_no_match3.fasta && -e $folder_a/functional_characterization/set_a_all_flags_u_no_match2.fasta ]]; then
  seqkit grep -iv -f protein_list.csv $folder_a/functional_characterization/set_a_all_flags_u_no_match2.fasta --quiet > set_a_all_flags_u_no_orthologous.fasta
fi
if ! test -f $folder_a/functional_characterization/set_a_all_flags_u_no_match2.fasta; then
  seqkit grep -iv -f protein_list.csv $folder_a/functional_characterization/set_a_all_flags_u_no_match1.fasta --quiet > set_a_all_flags_u_no_orthologous.fasta
fi

rnasamba classify output.tsv set_a_all_flags_u_no_orthologous.fasta /reference/others/partial_length_weights.hdf5

awk '$4 ~ /^coding/' output.tsv > set_a_all_flags_u_no_orthologous_coding.txt
awk '$4 ~ /^noncoding/' output.tsv > set_a_all_flags_u_no_orthologous_noncoding.txt
echo "The remaining SET A transcript loci without orthologous proteins similiarity was classified in:"
echo "coding_potential:"
wc -l set_a_all_flags_u_no_orthologous_coding.txt
echo "non_coding_potential:"
wc -l set_a_all_flags_u_no_orthologous_noncoding.txt

##Summary functional characterization## - Rscript 8
Rscript /Rscripts/script8_a.R
