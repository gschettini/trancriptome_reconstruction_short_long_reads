# Download base image ubuntu 22.04
FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive

# LABEL about the custom image
LABEL maintainer="gschettini@vt.edu"
LABEL version="0.1"
LABEL description="This is a custom Docker Image for Integrating Illumina and ONT sequence workflow."

# Install base utilities
RUN apt-get -y update \
    && apt-get install -y build-essential \
    && apt-get install -y curl apt-transport-https ca-certificates software-properties-common \
    && apt-get install -y vim nano wget curl ftp bzip2 gzip screen \
    && apt-get install -y r-base r-base-dev \
    && apt-get install -y default-jre \
    && apt-get install -y parallel rename \
    && apt-get install -y libmaus2-2 libcurl4-openssl-dev autotools-dev autoconf automake libtool m4 git cmake libssl-dev libfontconfig1-dev libxml2-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    && rm -rf /var/lib/apt/lists/* \
    && add-apt-repository universe \
    && add-apt-repository -y ppa:deadsnakes/ppa \
    && apt-get -y update \
    && apt-get install -y curl python3.7 python3.7-dev python3.7-distutils \
    && update-alternatives --install /usr/bin/python python /usr/bin/python3.7 1 \
    && update-alternatives --set python /usr/bin/python3.7 \
    && curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
    && python3.7 get-pip.py \
    && pip install --no-cache-dir nvidia-cudnn-cu11 \
    && pip install --no-cache-dir tensorflow \
    && git clone https://github.com/gt1/libmaus2.git && mv /libmaus2 /lib && cd /lib/libmaus2 && libtoolize && aclocal && autoheader && automake --force-missing --add-missing && autoconf && ./configure --prefix=/lib/libmaus2 && make && make install && cd /

# Download/Install packages
RUN wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download && mv download hisat2-2.2.1-Linux_x86_64.zip && unzip hisat2-2.2.1-Linux_x86_64.zip \
    && wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && unzip Trimmomatic-0.39.zip \
    && wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && tar xf samtools-1.18.tar.bz2 && cd /samtools-1.18 && ./configure && make && make install && cd / \
    && git clone https://github.com/gt1/biobambam2.git && cd /biobambam2 && autoreconf -i -f && ./configure --with-libmaus2=/lib/libmaus2 --prefix=/biobambam2 && make install && cd / \
    && git clone https://github.com/arq5x/bedtools2.git && cd /bedtools2 && make && make install && cd / \
    && git clone https://github.com/marbl/canu.git && cd /canu/src && make && cd / \
    && wget https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-Linux-x86_64.tar.gz/download && mv download subread-2.0.6-Linux-x86_64.tar.gz && tar xzf subread-2.0.6-Linux-x86_64.tar.gz \
    && git clone https://github.com/gpertea/gffcompare.git && cd /gffcompare && make release && cd / \ 
    && git clone https://github.com/gpertea/gffread.git && cd /gffread && make release && cd / \
    && wget https://sourceforge.net/projects/bbmap/files/BBMap_39.06.tar.gz/download && mv download BBMap_39.06.tar.gz && tar xzf BBMap_39.06.tar.gz \
    && wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2023-07-20.tar.gz && tar xvzf gmap-gsnap-2023-07-20.tar.gz && cd /gmap-2023-07-20 && ./configure && make && make install && cd / \
    && pip install --no-cache-dir rnasamba \
    && wget http://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz && tar xzf diamond-linux64.tar.gz -C /bin/ \
    && apt-get install -y ncbi-blast+ \
    && wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.10/sratoolkit.3.0.10-ubuntu64.tar.gz && tar -xzf sratoolkit.3.0.10-ubuntu64.tar.gz && rm sratoolkit.3.0.10-ubuntu64.tar.gz \
    && wget https://github.com/ablab/spades/releases/download/v3.14.1/SPAdes-3.14.1-Linux.tar.gz && tar -xzf SPAdes-3.14.1-Linux.tar.gz \
    && git clone https://github.com/lh3/minimap2.git && cd minimap2 && make && cp /minimap2/minimap2 /bin && cd / \
    && git clone https://github.com/lh3/seqtk.git && cd seqtk && make && cd / \
    && wget https://github.com/shenwei356/seqkit/releases/download/v2.7.0/seqkit_linux_amd64.tar.gz && tar -xzf seqkit_linux_amd64.tar.gz && rm seqkit_linux_amd64.tar.gz && mv seqkit /bin/ && cd / \
    && wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v385/bedToGenePred && chmod +x bedToGenePred && mv bedToGenePred /bin/ \
    && wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v385/genePredToGtf && chmod +x genePredToGtf && mv genePredToGtf /bin/

# Install R-packages
RUN Rscript -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org")'
RUN Rscript -e 'BiocManager::install(pkgs = c("stringr", "rtracklayer", "GenomicRanges"), site_repository ="https://cloud.r-project.org", force = T, update = T)'

#Remove all unnecessary files
RUN rm *.gz \
    && rm *.zip \
    && rm *.bz2

# Set the default command to run when the container starts
CMD ["/bin/bash"]
