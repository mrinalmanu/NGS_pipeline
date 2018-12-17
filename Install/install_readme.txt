1 Installation
All scripts are in https://github.com/ptriska/NGS_pipeline/
1.1 Build docker base container
docker image build -t [container name] --no-cache .
1.2 Pull docker images with tools
run pull_images.sh. This will pull following images with tools:
vt genome tools: 	quay.io/biocontainers/vt:0.57721--hf74b74d_1
varscan: 		opengenomics/varscan:latest
bwa+samtools1.3.1:	wasa000/bwa_samtools131
bcftools1.9: 		wasa000/bcftools:1.9
vep annotator: 	willmclaren/ensembl-vep:latest
GATK4: 		broadinstitute/gatk
1.3 Run base container
Base container is started by bash script run_ngs_container.sh. Container mounts several locations:
    • docker socker (-v /var/run/docker.sock:/var/run/docker.sock)
    • location where input data will be stored (container will mount this location as /home/input)
    • location where output data will be stored (container will mount this location as /home/output)
    • location with reference files (container will mount this location as /home/Reference)
mandatory flags:
-r: Reference dir 
-i: Input docker mount point folder
-o: Output docker mount point folder
-I: Base container image
-n: Container name (username)*

*make sure that container name is not being used, try docker ps -a |grep [container name].
2 Prepare run
2.1 Create sample list
Create a comma-separated text file with mandatory header. Next rows contain samples, one sample per row.
NOTE: input and output location is mounted at /home/input (output). Keep this in mind when specifying location of FASTQ files and specify only subfolders! Example: If your FASTQ is in /home/analysis/FASTQ/sample1/sample1.fastq, and your mounted input location is /home/analysis, then set your input location in sample file as /home/input/FASTQ/sample1/sample1.fastq
Header fields are:
SAMPLE_ID,			- any text string (no white spaces)
PIPELINE,			- pipeline_ID [GATK4_germline_v1, GATK4_somatic_TN_v1,GATK_varscan_TN_v1]
GENOME_BUILD,		- [hg19, GRCh38]
FASTQ1,			- location of normal fastq1, if tumor-only run: leave empty
FASTQ2,			- location of normal fastq2, if tumor-only run: leave empty
FASTQ1_TUMOR,		- location of tumor fastq1, if germline run: leave empty
FASTQ2_TUMOR,		- location of tumor fastq2, if germline run: leave empty
CLEANUP,			- set YES is wish to remove intermediate BAM
 LIB_ID,			- library_ID any text string (not important)
PL_ID,			- PL string (not important)
PU_ID				- PU string (not important)
2.2 Review docker images versions (optional)
Docker images versions file is located in /home/NGS_pipeline/docker_images_versions.py
For each pipeline, you can specify exact docker image to be used for each tool. For example, if you wish to use different version of samtools in somatic and germline pipeline, you can specify for each pipeline which image of samtools should be used. This can improve customization and reproducibility of results. Docker images are specified in docker_images_versions.txt in /home/NGS_pipeline.
2.3 Review config_file (optional)
Config_file is located in /home/NGS_pipeline/pipelines/config_file.py
It contains names of indexes and reference genome fasta files in Reference folder. Make sure these point to existing files in your mounted Reference directory. Maximum threads available can be set by changing of variable "max_nr_threads".
3 Run samples
Run is invoked by script run_pipeline.py. Script takes only one argument: sample list.
4 Documentation of scripts
Pipelines are organized into modules in /pipelines. Each pipeline module contains a class named after the pipeline (eg. GATK4_somatic_TN_v1). This class contains function run_pipeline() which initiates processing of sample.
4.1 Management script: run_pipeline.py
Script has three parts: first parses docker_images_versions.txt, second parses sample list and creates an instance of class of respective pipeline for every sample in sample list. Third part checks available cores and submits samples for processing by invoking run_pipeline() function.
4.2 Pipeline modules
Pipeline modules are in /home/NGS_pipeline/pipelines
All pipeline modules contain self-named class Pipeline class contains 
    • init function
    • run-in-docker function (runs tools in containers)
    • run_pipeline function
Function run_pipeline() contains
    1. functions, which make building blocks of the pipeline
    2. PIPELINE WORKFLOW part, which invokes functions in correct order. In tumor-normal somatic pipelines, pipeline workflow consists of pre-process, which is done in parallel, and then join variant calling and annotation.
