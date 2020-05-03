FROM centos:7

ENV PATH=/apps/exseek/bin:/apps/ucsc_tools:/opt/microsoft/ropen/3.5.3/lib64/R/bin:/apps/anaconda3/bin:/apps/anaconda2/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

COPY CentOS-Base.repo /etc/yum.repos.d
RUN yum group install -y 'Development Tools' \
    && yum install -y epel-release \
    && yum install -y curl screen \
         yum-utils \
        coreutils which binutils vim nano openssl-devel Xorg-devel mesa-libGL-devel \
        libxml-devel libxml2-devel \
        java-1.8.0-openjdk wget pigz gawk \
    && yum-builddep -y R \
    && mkdir -p /apps/archive /apps/ucsc_tools


# install anaconda3
COPY Anaconda3-2019.03-Linux-x86_64.sh /apps/archive
RUN bash /apps/archive/Anaconda3-2019.03-Linux-x86_64.sh -b -p /apps/anaconda3 

# install pip packages
RUN /apps/anaconda3/bin/pip install skrebate mlxtend tqdm jinja2 snakemake pyBigWig multiqc

# UCSC kent tools
RUN mkdir -p /apps/ucsc_tools \
    && wget -q -P /apps/ucsc_tools/ 'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToGenePred' \
        'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf' \
        'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig' \
        'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph' \
    && chmod +x /apps/ucsc_tools/*

# install R
COPY microsoft-r-open-3.5.3.tar.gz /apps/archive
RUN tar -C /apps/archive -zxf /apps/archive/microsoft-r-open-3.5.3.tar.gz \
    && bash /apps/archive/microsoft-r-open/install.sh -s -a -u

# install R packages
COPY install-R-packages.R /apps/archive/
RUN cd /apps/archive; Rscript /apps/archive/install-R-packages.R

# install conda packages
RUN /apps/anaconda3/bin/conda install -y 'python<3.7' numpy scipy scikit-learn 'openssl<1.1' \
    && /apps/anaconda3/bin/conda install -y pandas matplotlib seaborn h5py \
    && /apps/anaconda3/bin/conda install -c bioconda -y bedtools samtools star subread bowtie2 rsem bamtools cutadapt picard gffread gffcompare \
    && /apps/anaconda3/bin/conda install -c bioconda -y fastx_toolkit biopython

# cleanup 
RUN /apps/anaconda3/bin/conda clean -y --all \
    && yum clean all \
    && rm -rf /apps/archive

# Install exseek-pipeline
COPY exseek_pipeline-1.0.0-py3-none-any.whl /apps
RUN /apps/anaconda3/bin/pip install /apps/exseek_pipeline-1.0.0-py3-none-any.whl \
    && rm -f /apps/exseek_pipeline-1.0.0-py3-none-any.whl

