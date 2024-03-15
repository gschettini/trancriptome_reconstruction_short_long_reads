[![GitHub release (latest by date)](https://img.shields.io/github/v/release/gschettini/trancriptome_reconstruction_short_long_reads)](https://github.com/gschettini/trancriptome_reconstruction_short_long_reads/releases/)
[![GitHub Downloads](https://img.shields.io/github/downloads/gschettini/trancriptome_reconstruction_short_long_reads/total.svg?style=social&logo=github&label=Download)](https://github.com/gschettini/trancriptome_reconstruction_short_long_reads/releases)
[![Docker pulls](https://img.shields.io/docker/pulls/gschettini/integrating_sr_lr_workflow_transcriptome_reconstruction.svg?label=Docker%20pulls&color=blue)](https://hub.docker.com/r/gschettini/integrating_sr_lr_workflow_transcriptome_reconstruction)

<font size=20>__Integrating Illumina short-reads and ONT long-reads - Transcriptome reconstruction__</font> 

1. [About](#sec1)</br>
2. [Installation](#sec2)</br> 
3. [Pre-requirements](#sec3)</br>
4. [Usage](#sec4)</br>
    4.1. [Two-set pipeline](#sec4.1)</br>
    4.2. [One-set pipeline](#sec4.2)</br>
5. [Outputs](#sec5)</br>
    5.1. [Pre-processing Illumina Short-reads](#sec5.1)</br>
    5.2. [Pre-processing ONT Long-reads](#sec5.2)</br>
    5.3. [<i>de novo</i> Assembly](#sec5.3)</br>
    5.4. [Annotation Comparisons](#sec5.4)</br>
    5.5. [Functional Characterization](#sec5.5)</br>
6. [Citations](#sec6)</br>

<a name="sec1"></a>
# About

The "Integrating Illumina short-reads and ONT long-reads workflow" is a comprehensive pipeline, which uses RNA sequencing data from both 150bp paired-end Illumina and Oxford Nanopore Technologies reads. This pipeline aims to reconstruct the transcriptome, identify potential unknown transcripts, and perform functional classification.

> [!NOTE]  
> This work was partly supported by funding from the Agricultural Genome to Phenome Initiative (AG2PI), which is funded by USDA-NIFA award 2022-70412-38454. However, USDA-NIFA and AG2PI funding bodies had no role in the design and development of the workflow.<br><br>
> <img src="https://github.com/gschettini/trancriptome_reconstruction_short_long_reads/assets/107274474/847664c0-32ed-401c-89e0-7d5c593ffd92" width="100" height="55"><img src="https://github.com/gschettini/trancriptome_reconstruction_short_long_reads/assets/107274474/6629b405-5fcc-4862-8f1e-cd78f83d4967" width="100" height="55">

<a name="sec2"></a>
# Installation
```sh
docker pull gschettini/integrating_sr_lr_workflow_transcriptome_reconstruction
```

<b>OR</b>

```sh
mkdir ~/integrating_SR_LR
cd ~/integrating_SR_LR
git clone https://github.com/gschettini/trancriptome_reconstruction_short_long_reads.git
docker build -t integrating_SR_LR .
```

<a name="sec3"></a>
# Pre-requirements

* Download all necessary files/scripts for this workflow <i>(Step required only if this GitHub repository was cloned)</i>:
```sh
cd /
wget https://github.com/gschettini/trancriptome_reconstruction_short_long_reads/blob/main/integrating_sr_lr_files.tar.gz
tar -xzvf integrating_sr_lr_files.tar.gz
```

Also, it is required to download/copy into the Docker container sample sequences (Illumina PE short-reads and ONT long-reads).

* Download sample sequences <br><i>(It is an optional step since the user can copy files into the docker container)</i>:

```sh
./download_sequences.sh
```

* Download reference files <br><i>(It is an optional script since the user can copy files into the docker container)</i>:

```sh
./download_references.sh
```

`download_sequences.sh` and `download_references` can be modified for other organisms and samples.


> [!WARNING]  
> Depending on the organism and amount of sample data, there will be steps requiring a large amount of memory RAM.
>  - Hisat2 indexes: It is possible that this step cannot be completed due to insufficient memory RAM.<br><br>
>        - Option 1: Download pre-indexed genomes from [Hisat2 website](https://daehwankimlab.github.io/hisat2/download/) into index folder (`/reference/fasta/hisat/index/`) and change `path_to_index` variable in `setting.sh`. <br><br>
>        - Option 2a: Adapt the workflow. It can be done by removing `--ss $splice_sites` option from the line 38 (`/hisat2-2.2.1/hisat2-build -p $(echo $n_cores) --ss $splice_sites $reference_genome_ensembl $path_to_index`) in `/data/integrating_SR_LR.rmd` or from the line 21 in (`master/docker_script.sh` or `/master/docker_script_a.sh`). <br><br>
>        - Option 2b: Further adapt the workflow. Manually set Hisat2 build options related to memory RAM allocation (e.g. `--bmaxdivn` and `--dcv`). For more details, access [Hisat2 manual](https://daehwankimlab.github.io/hisat2/manual/).<br><br><br>
>  - <i>de novo</i> Assemblies: It is possible that the step cannot be completed due to insufficient memory RAM.<br><br>
>        - Option 1: Adapt the workflow. Set more strict thresholds during the Pre-processing step by modifying `/data/integrating_SR_LR.rmd` (L199-L206 and L258-L265) or `master/docker_script.sh` (L149-L156 and L197-L204) or `master/docker_script_a.sh` (L103-L109). For more information related to parameters, check [BBtools User guide](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/).<br>


<a name="sec4"></a>
# Usage

There are two ways to use the workflow. 

-  Two-set pipeline with two sets of samples (e.g. Oocytes vs. Blastocysts), which includes a final comparison and determination of unknown loci shared and/or exclusive for each set.

-  One-set pipeline with only one set of samples (e.g. Oocytes), which does not include a final comparison and determination of unknown loci shared and/or exclusive for each set.

> [!TIP]
> An optional [Rmarkdown code](https://github.com/gschettini/trancriptome_reconstruction_short_long_reads/blob/main/data/integrating_SR_LR.rmd) is provided `/data/integrating_SR_LR.rmd` and can be easily adapted in case of limitations.


The file `settings.sh` is a shell script that can be modified according to the user's needs.


```sh
#!/bin/bash
folder_a=/test/set_a
folder_b=/test/set_b
n_cores=30
memory_ram=150
genome_size=2.8g
path_to_index=/reference/fasta/hisat/index/Bos_taurus.ARS-UCD1.2.110
splice_sites=/reference/fasta/hisat/splice_sites/splice_sites.txt
reference_genome_ensembl=/reference/fasta/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa
gtf_ensembl=/reference/annotation/Bos_taurus.ARS-UCD1.2.110.gtf
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

...

```

`folder_a` and `folder_b` can be changed to any folder you want to use to conduct the workflow.
If you only need a <b>one-way</b> pipeline, you can set `folder_b` as blank by leaving it empty (e.g. `folder_b= `).

<i> In case of any modification in `folder_a` and/or `folder_b`, `settings.sh` needs to be executed to create all necessary (sub-)directories into the new paths. </i>

```sh
./settings.sh
```

`n_cores` and `memory_ram` allow the user to set the number of threads and allocate a specific amount of memory RAM in Gb.

`genome_size` determines the genome size of the selected organism in Gbases. This value will be used by [Canu](#sec6) during the ONT long-reads correction step.

`path_to_index` and `splice_sites` are files that will be produced by [Hisat2](#sec6) and be used in the workflow. Both can be modified to a different name/organism.

`reference_genome_ensembl`, `gtf_ensembl`, `reference_genome_ncbi`, `gff_ncbi` and `gtf_lncRNA` are necessary files for alignment ([Hisat2](#sec6) and [GMAP](#sec6)) and comparisons ([GFFread](#sec6) and [GFFcompare](#sec6)) steps.

`reference_GMAP_ensembl` and `reference_GMAP_ncbi` are folders where [GMAP](#sec6) will generate the index files to align and generate new annotation files.

`adapters_SR_a` and `adapters_SR_b` set adapters used for paired-end 150bp Illumina short-reads which will be used for adapter trimming ([Trimmomatic](#sec6)). 
The list of adapters can be found in `/Trimmomatic-0.39/adapters`.


`Organism_1`, `Organism_2`, and `Organism_3` set the taxonomic IDs for the functional characterization through protein sequence similarity against NCBI non-redundant protein database ([NCBI](#sec6)) by [DIAMOND](#sec6). Leaving any of these `Organism` variables empty will disable the search for this specific organism.<br>
<i>Also, it is possible to search for multiple taxonomic IDs by separating them with commas. (e.g. organism_1=9913,9915,30522)</i>

`Organism_1_e_value`, `Organism_2_e_value`, and `Organism_3_e_value` sets an e-value threshold for each organism set.

`Organism_1_perc_id`, `Organism_2_perc_id`, and `Organism_3_perc_id` sets a percentage of identity threshold between the query and subject for each organism set.

`Organism_1_cov`, `Organism_2_cov`, and `Organism_3_cov` set a subject coverage threshold for each organism set.

<a name="sec4.1"></a>
## Two-set pipeline

First, `folder_a` and `folder_b` are required to be set in `settings.sh` to activate a two-set pipeline.

Second, it is necessary to modify/add the path to Illumina short-read samples in the `files_a.yaml` and/or `files_b.yaml`, which will be used in the [RNAspades](#sec6)

-  /test/set_a/files_a.yaml
-  /test/set_b/files_b.yaml

<i>In case of modifications/addition of samples,<b>ONLY</b> change the SRA accession IDs:

(e.g. SRR5645291 -> sample_001).</i>

The first part of the path to the files will be changed automatically if `folder_a` and/or `folder_b` are modified (<i>/test/set_a</i>).

In case of more samples, copy the path/to/files and modify the SRA accession/sample name.

<i> Those files still do not exist if you did not run the workflow. However, they will be generated and used along with the de novo assembly. </i>

```yaml
    [
      {
        orientation: "fr",
        type: "paired-end",
        right reads: [
          "/test/set_a/SR_fastq/preprocessed/SRR5645291_R1.normalized.filtered.trimmed.fastq",
          "/test/set_a/SR_fastq/preprocessed/SRR5645294_R1.normalized.filtered.trimmed.fastq"
        ],
        left reads: [
          "/test/set_a/SR_fastq/preprocessed/SRR5645291_R2.normalized.filtered.trimmed.fastq",
          "/test/set_a/SR_fastq/preprocessed/SRR5645294_R2.normalized.filtered.trimmed.fastq"
        ]
      },
      {
        type: "nanopore",
        single reads: [
          "/test/set_a/LR_nanopore/canu/trimming/set_a.trimmedReads.fasta.gz"
        ]
      }
    ]
```

Then, with all the settings and required files (references and samples) set, it is possible to start the workflow by `integrating_SR_LR.sh`.

```sh
./integrating_SR_LR.sh
```


<a name="sec4.2"></a>
## One-set pipeline

First, it is necessary to leave `folder_b` in `settings.sh` as blank to activate a one-set pipeline.

Second, It is required to modify/add the path to samples in the `files_a.yaml` similar in the [two-set pipeline](#sec4.1) .

-  /test/set_a/files_a.yaml

Then, with all the settings and required files (references and samples) set, it is possible to start the workflow by `integrating_SR_LR.sh`.

```sh
./integrating_SR_LR.sh
```

<a name="sec5"></a>
# Outputs

<a name="sec5.1"></a>
## Pre-processing Illumina Short-reads
<b>Trimming</b><br>
`$folder_a/SR_fastq/trimming/*.trimmed.fastq*` and `$folder_b/SR_fastq/trimming/*.trimmed.fastq*`: Adapter-trimmed FASTQ files from your Illumina PE short-reads by [Trimmomatic](#sec6).

<b>Filtering</b><br>
`$folder_a/SR_fastq/alignment/*._alignment.filtered.sorted.undup.bam` and `$folder_b/SR_fastq/alignment/*_alignment.filtered.sorted.undup.bam `: BAM files from alignment using [Hisat2](#sec6), sorted, filtered (-F [1796](https://broadinstitute.github.io/picard/explain-flags.html)), and with duplicated reads removed by [SAMtools](#sec6) and [Biombambam2](#sec6).

<b>Pre-processing FASTQ files</b><br>
`$folder_a/SR_fastq/preprocessed/*R*.normalized.filtered.trimmed.fastq` and `$folder_b/SR_fastq/preprocessed/*R*.normalized.filtered.trimmed.fastq`: Pre-processed FASTQ files from your Illumina PE short-reads by [BBtools](#sec6) (Contaminant removed (sequencing artifacts), Error corrected, Quality filtered (`minbasequality` = 25 ) and normalized (read coverage maximum of 30x)).

<a name="sec5.2"></a>
## Pre-processing ONT Long-reads
<b>Correction</b><br>
`$folder_a/LR_nanopore/canu/correction/set_a.correctedReads.fasta.gz` and `$folder_b/LR_nanopore/canu/correction/set_b.correctedReads.fasta.gz`: Corrected ONT Long-reads by [Canu](#sec6).

<b>Trimming</b><br>
`$folder_a/LR_nanopore/canu/trimming/set_a.trimmedReads.fasta.gz` and `$folder_b/LR_nanopore/canu/trimming/set_b.trimmedReads.fasta.gz`: Corrected and Trimmed ONT Long-reads by [Canu](#sec6).

<a name="sec5.3"></a>
## <i>de novo</i> Assembly
<b>FASTA file</b><br>
`$folder_a/de_novo_assembly/hard_filtered_transcripts.fasta` and `$folder_b/de_novo_assembly/hard_filtered_transcripts.fasta`: Long and reliable sequence reads from <i>de novo</i> assembly by [rnaSPAdes](#sec6).<br>
The remaining outputs present in these folders can be found in [rnaSPAdes](https://cab.spbu.ru/files/release3.14.1/rnaspades_manual.html) and [SPAdes](https://cab.spbu.ru/files/release3.14.1/manual.html) manual.<b> However, they will not be used in our workflow.</b>

<b>Summary</b><br>
`$folder_a/de_novo_assembly/hard_filtered_transcripts_summary.txt` and  `$folder_b/de_novo_assembly/hard_filtered_transcripts_summary.txt`: <i>de novo</i> assembly summary by [seqkit](#sec6).


<a name="sec5.4"></a>
## Annotation Comparisons

## ENSEMBL

`$folder_a/comparisons/ensembl/set_a_ensembl.annotated.gtf` and `$folder_b/comparisons/ensembl/set_b_ensembl.annotated.gtf` are comparisons with ENSEMBL annotation file outputs, generated by [GFFcompare](#sec6).<br>

`$folder_a/comparisons/ensembl/unknown_set_a_ensembl.annotated.gtf` and `$folder_b/comparisons/ensembl/unknown_set_b_ensembl.annotated.gtf` are a list of unknown transcripts (flag "-u"), which were used for the following comparisons.<br>

<i>More details about flags can be found at [GFFcompare website](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml#transfrag-class-codes).</i>

## NCBI/RefSeq

`$folder_a/comparisons/ncbi/set_a_ncbi.annotated.gtf` and `$folder_b/comparisons/ncbi/set_b_ncbi.annotated.gtf` are comparisons with NCBI annotation file outputs, generated by [GFFcompare](#sec6).<br>

`$folder_a/comparisons/ncbi/unknown_set_a_ncbi.annotated.gtf` and `$folder_b/comparisons/ncbi/unknown_set_b_ncbi.annotated.gtf` are a list of unknown transcripts (flag "-u"), which were used for the following comparisons.<br>

## NONCODEv5 (lncRNA)

`$folder_a/comparisons/lncRNA/set_a_lncRNA.annotated.gtf` and `$folder_b/comparisons/lncRNA/set_b_lncRNA.annotated.gtf` are comparisons with NONCODEv5 lncRNA annotation file outputs, generated by [GFFcompare](#sec6).<br>

`$folder_a/comparisons/lncRNA/unknown_set_a_lncRNA.annotated.gtf` and `$folder_b/comparisons/lncRNA/unknown_set_b_lncRNA.annotated.gtf` are a list of unknown transcripts (flag "-u").<br>

## Quality filtering after comparisons

`$folder_a/comparisons/lncRNA/set_a_all_flags_u_final3.gtf` and `$folder_b/comparisons/lncRNA/set_b_all_flags_u_final3.gtf` are a final annotation file after filtered out any transcript with less than 300bp of length and only one exon, close to known genes (<5kb upstream/downstream) and reduced redundancy (merged immediate transcripts <300bp) by [GFFread](#sec6), [BEDtools](#sec6), and [GenomicRanges](#sec6). Also, any transcript with zero raw counts in more than 90% of samples (done by [FeatureCounts](#sec6) with Illumina SR samples) was filtered out.

`$folder_a/comparisons/lncRNA/only_transcripts_set_a_all_flags_u_final3.gtf` and `$folder_b/comparisons/lncRNA/only_transcripts_set_b_all_flags_u_final3.gtf` are a list of unknown transcripts from the last annotation file .<br>

<a name="sec5.5"></a>
## Functional Characterization

`$folder_a/functional_characterization/summary_transcript_loci_set_a` and `$folder_b/functional_characterization/summary_transcript_loci_set_b` are the final summary of unknown transcript, which includes hits found by [DIAMOND](#sec6) based on the thresholds (percentage of identity and e-score) and `organism` set in `settings.sh`. Also, it includes potential coding/non-coding classification done by [RNAsamba](#sec6).<br>

<i>If Two-set was enabled, it will generate additional files, `$folder_a/R_outputs/set_a_transcripts_in_common_with_set_b.txt`, `$folder_a/R_outputs/set_a_transcripts_exclusive.txt` and `$folder_b/R_outputs/set_b_transcripts_exclusive.txt`, which will have a list of transcripts in common between the two sets of samples and a list of exclusive ones, respectively.</i>

<a name="sec6"></a>
# Citations

-   BBtools - Bushnell, B. (2018) https://jgi.doe.gov/data-and-tools/software-tools/bbtools/

-   BEDtools - Quinlan, A. R., & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841–842. https://doi.org/10.1093/bioinformatics/btq033

-   biobambam2 - Tischler, G., & Leonard, S. (2014). biobambam: tools for read pair collation based algorithms on BAM files. Source Code for Biology and Medicine, 9(1), 13. https://doi.org/10.1186/1751-0473-9-13

-   DIAMOND - Buchfink, B., Reuter, K., & Drost, H.-G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods, 18(4), 366–368. https://doi.org/10.1038/s41592-021-01101-x

-   Ensembl database - Flicek, P., Amode, M. R., Barrell, D., Beal, K., Billis, K., Brent, S., Carvalho-Silva, D., Clapham, P., Coates, G., Fitzgerald, S., Gil, L., Girón, C. G., Gordon, L., Hourlier, T., Hunt, S., Johnson, N., Juettemann, T., Kähäri, A. K., Keenan, S., … Searle, S. M. J. (2014). Ensembl 2014. Nucleic Acids Research, 42(D1), D749–D755. https://doi.org/10.1093/nar/gkt1196

-   Featurecounts - Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), 923–930. https://doi.org/10.1093/bioinformatics/btt656

-   GenomicRanges - Lawrence, M., Huber, W., Pagès, H., Aboyoun, P., Carlson, M., Gentleman, R., Morgan, M. T., & Carey, V. J. (2013). Software for Computing and Annotating Genomic Ranges. PLoS Computational Biology, 9(8), e1003118. https://doi.org/10.1371/journal.pcbi.1003118

-   GFFcompare - Pertea, G., & Pertea, M. (2020). GFF Utilities: GffRead and GffCompare. F1000Research, 9, 304. https://doi.org/10.12688/f1000research.23297.2

-   GFFread - Pertea, G., & Pertea, M. (2020). GFF Utilities: GffRead and GffCompare. F1000Research, 9, 304. https://doi.org/10.12688/f1000research.23297.2

-   GMAP - Wu, T. D., & Watanabe, C. K. (2005). GMAP: a genomic mapping and alignment program for mRNA and EST sequences. Bioinformatics, 21(9), 1859–1875. https://doi.org/10.1093/bioinformatics/bti310

-   Hisat2 - Kim, D., Paggi, J. M., Park, C., Bennett, C., & Salzberg, S. L. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nature Biotechnology, 37(8), 907–915. https://doi.org/10.1038/s41587-019-0201-4

-   htmlwidgets - Vaidyanathan R., Xie Y., Allaire JJ., Cheng J. https://CRAN.R-project.org/package=htmlwidgets

-   NCBI database - Sayers, E. W., Bolton, E. E., Brister, J. R., Canese, K., Chan, J., Comeau, D. C., Connor, R., Funk, K., Kelly, C., Kim, S., Madej, T., Marchler-Bauer, A., Lanczycki, C., Lathrop, S., Lu, Z., Thibaud-Nissen, F., Murphy, T., Phan, L., Skripchenko, Y., … Sherry, S. T. (2022). Database resources of the national center for biotechnology information. Nucleic Acids Research, 50(D1), D20–D26. https://doi.org/10.1093/nar/gkab1112

-   NONCODEv5 - Fang, S., Zhang, L., Guo, J., Niu, Y., Wu, Y., Li, H., Zhao, L., Li, X., Teng, X., Sun, X., Sun, L., Zhang, M. Q., Chen, R., & Zhao, Y. (2018). NONCODEV5: a comprehensive annotation database for long non-coding RNAs. Nucleic Acids Research, 46(D1), D308–D314. https://doi.org/10.1093/nar/gkx1107

-   Parallel - Tange, O. (2023). GNU Parallel 20230522 ('Charles') [stable]. Zenodo. https://doi.org/10.5281/zenodo.7958356

-   RNAsamba - Camargo, A. P., Sourkov, V., Pereira, G. A. G., & Carazzolle, M. F. (2020). RNAsamba: neural network-based assessment of the protein-coding potential of RNA sequences. NAR Genomics and Bioinformatics, 2(1). https://doi.org/10.1093/nargab/lqz024

-   RNAspades - Prjibelski, A. D., Puglia, G. D., Antipov, D., Bushmanova, E., Giordano, D., Mikheenko, A., Vitale, D., & Lapidus, A. (2020). Extending rnaSPAdes functionality for hybrid transcriptome assembly. BMC Bioinformatics, 21(S12), 302. https://doi.org/10.1186/s12859-020-03614-2

-   rtracklayer - Lawrence M, Gentleman R, Carey V (2009). “rtracklayer: an R package for interfacing with genome browsers.” Bioinformatics, 25, 1841-1842. doi:10.1093/bioinformatics/btp328, http://bioinformatics.oxfordjournals.org/content/25/14/1841.abstract.

-   rvest - Wickham H. https://CRAN.R-project.org/package=rvest

-   SAMtools - Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352

-   Seqkit - Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962

-   Seqtk -  Seqtk: Toolkit for processing sequences in FASTA/Q formats. GitHub. https://github.com/lh3/seqtk

-   SRAtoolkit - SRAtoolkit: The SRA Toolkit and SDK from NCBI is a collection of tools and libraries for using data in the INSDC Sequence Read Archives. GitHub. https://github.com/ncbi/sra-tools

-   stringr - Wickham H. https://CRAN.R-project.org/package=stringr 

-   Trimmomatic - Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114–2120. https://doi.org/10.1093/bioinformatics/btu170

-   UCSC Genome browser and Tools - Karolchik, D., Hinrichs, A. S., Furey, T. S., Roskin, K. M., Sugnet, C. W., Haussler, D., & James Kent, W. (2004). The UCSC Table Browser data retrieval tool. Nucleic Acids Research, 32(90001), 493D – 496. https://doi.org/10.1093/nar/gkh103
