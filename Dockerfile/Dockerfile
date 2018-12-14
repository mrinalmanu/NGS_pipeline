# Base Image
FROM biocontainers/biocontainers:latest
USER root
RUN conda install bcftools=1.3.1
RUN apt-get update
RUN apt-get install docker.io -y
RUN apt-get install python-psutil -y
RUN apt-get install python-numpy -y
RUN apt-get install python-pandas -y
RUN apt-get install tabix
RUN docker pull broadinstitute/gatk
RUN docker pull quay.io/biocontainers/vt:0.57721--hf74b74d_1
RUN docker pull opengenomics/varscan:latest
RUN docker pull wasa000/bwa_samtools131
RUN docker pull bcftools:1.9
RUN docker pull willmclaren/ensembl-vep:latest
RUN docker pull broadinstitute/gatk
RUN git clone https://github.com/ptriska/NGS_pipeline/ /home/NGS_pipeline/


