#!/bin/bash

source ./settings.sh

#Set A - SR
sratoolkit.3.0.10-ubuntu64/bin/fastq-dump --split-files -O $folder_a/SR_fastq/ SRR5645291
sratoolkit.3.0.10-ubuntu64/bin/fastq-dump --split-files -O $folder_a/SR_fastq/ SRR5645294

#Set A - LR
#sratoolkit.3.0.7-ubuntu64/bin/fastq-dump -O $folder_a/LR_fastq/
#sratoolkit.3.0.7-ubuntu64/bin/fastq-dump -O $folder_a/LR_fastq/

#Set B - SR
sratoolkit.3.0.10-ubuntu64/bin/fastq-dump --split-files -O $folder_b/SR_fastq/ SRR23565574
sratoolkit.3.0.10-ubuntu64/bin/fastq-dump --split-files -O $folder_b/SR_fastq/ SRR23565576

#Set B - LR
#sratoolkit.3.0.7-ubuntu64/bin/fastq-dump -O $folder_b/LR_fastq/
#sratoolkit.3.0.7-ubuntu64/bin/fastq-dump -O $folder_b/LR_fastq/


#Rename filenames - Set A
cd $folder_a/SR_fastq

rename 's/_1/_R1/g' *.fastq*
rename 's/_2/_R2/g' *.fastq*

#Rename filenames - Set B
cd $folder_b/SR_fastq

rename 's/_1/_R1/g' *.fastq*
rename 's/_2/_R2/g' *.fastq*
