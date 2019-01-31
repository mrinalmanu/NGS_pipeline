### GATK4 RNA pipeline version 1.0
### Pipeline created: 01/23/2019

import sys
import subprocess
import logging
import os
from gatk4_docker import gatk_docker
from multiprocessing import Process
from config_file import *
from lock_module import *
import time


def check_path(path):                		## Check if file exists
			if os.path.exists(path) == True:
				return True
			else:
				return False
				
def checkContainer(container_name):
		containers = subprocess.check_output(['docker','ps','-a']).decode(encoding="437")
		if container_name in containers.split():
			subprocess.call(["docker","rm",container_name])
		else:
			return				

class GATK_RNA_STAR():
	def __init__(self, sample_name, genome_build, f1, f2,tf1,tf2, cleanup,lib_ID, pl_ID, pu_ID, docker_images_dict):
		self.sample_name = sample_name
		self.fastq1_tumor = tf1
		self.fastq2_tumor = tf2
		self.fastq1 = f1
		self.fastq2 = f2
		if (f1 =="") or (f2 ==""):
			print (self.sample_name+": ERROR: FASTQ file is missing, aborting.")
			quit()
		if (tf1 !="") or (tf2 !=""):
			print (self.sample_name+": ERROR: tumor FASTQ file specified in germline pipeline, aborting.")
			quit()
		self.input_folder=input_folder
		self.threads = max_threads_process
		self.ram = "50000"
		self.open_files = []
		self.genome_build = genome_build
		self.cleanup = cleanup
		self.docker_images_dict = docker_images_dict
		self.lib_ID = lib_ID
		self.pl_ID = pl_ID
		self.pu_ID = pu_ID
		
		self.output_folder = output_folder+"/"+self.sample_name+"/output/"
		if check_path(self.output_folder)==False:
			print ("MKDIR", self.output_folder)			
			os.system("mkdir -p "+self.output_folder)		
		os.system("chmod 777 "+self.output_folder)	
		
		# logger setup
		LOG_FORMAT = "%(levelname)s\t%(asctime)s\t%(module)s\t%(funcName)s\t%(message)s"
		logging.basicConfig(filename = None,
												level = logging.DEBUG,
												format = LOG_FORMAT)
		self.logger = logging.getLogger()
	
	
	def run_in_docker(self, cmd, image, threads_needed, stdout=None, stderr=None):
		""" Run a command inside docker container"""
		
		if len(cmd)<2:
			container_name = self.sample_name+"_vt"
		else:
			container_name = self.sample_name+"_"+cmd[0]
		checkContainer(container_name)
		
		if len(cmd[0].split(" "))>1:
			if cmd[0].split(" ")[1]=="STAR":
				dcmd = ["docker", "run","--name",container_name,
								"-v", "{}:{}".format(self.input_folder, self.input_folder),
								"-v", "{}:{}".format(output_folder, output_folder),
								"-v", "{}:{}".format(reference_folder, reference_folder),
								'--ipc="host"',
								self.docker_images_dict[image],
								"/bin/bash -c "
								]
		else:
			dcmd = ["docker", "run","--name",container_name,
								"-v", "{}:{}".format(self.input_folder, self.input_folder),
								"-v", "{}:{}".format(output_folder, output_folder),
								"-v", "{}:{}".format(reference_folder, reference_folder),
								self.docker_images_dict[image]
								]
											
		dcmd += cmd
		if stdout is not None:
			stdout = open(stdout,"w")
			self.open_files.append(stdout)
		if stderr is not None:
			stderr = open(stderr,"w")
			self.open_files.append(stderr)
			
		############ LOCK DOWN NEEDED THREADS #############
		threads_needed = int(threads_needed)
		wait_go = "wait"
		while wait_go !='go':
			time.sleep(10)
			wait_go = check_threads(batch_ID, container_name,"start",(-1)*threads_needed, max_nr_threads)	
		
		########### RUN COMMAND IN DOCKER #################
		self.logger.info("Running: "+" ".join(dcmd))
		os.system(" ".join(dcmd))
		checkContainer(container_name)
		
		########### RELEASE LOCK DOWN CORES ###############
		wait_go = "wait"
		while wait_go !='go':
			time.sleep(10)
			wait_go = check_threads(batch_ID, container_name,"finish", threads_needed, max_nr_threads)	
	###################### END OF RUN COMMANDS IN DOCKER CONTAINERS ###################
			
	def run_pipeline(self):
		
		open("pipelines/"+batch_ID+".CPUlog",'w') ###open empty log file
		### Genome build selection ###
		if self.genome_build == "GRCh38":
			bwa_index = bwa_index_GRCh38
			reference_fasta = reference_fasta_GRCh38
			dbsnp_vcf = dbsnp_vcf_GRCh38
		if self.genome_build == "hg19":
			bwa_index = bwa_index_hg19
			reference_fasta = reference_fasta_hg19
			dbsnp_vcf = dbsnp_vcf_hg19
		
		### Pipeline auxiliary files ###
		unsorted_bam = self.output_folder+self.sample_name+".STAR.Aligned.sortedByCoord.out.bam"
		sorted_bam = self.output_folder+self.sample_name+".STAR.Aligned.sortedByCoord.out.bam"
		sort_err = self.output_folder+self.sample_name+".sort.err"
		markDuplicates_bam = self.output_folder+self.sample_name+".MarkDuplicates.bam"
		markDuplicates_metrics = self.output_folder+self.sample_name+".MarkDuplicates-metrics.txt"
		ReadGroups_bam = self.output_folder+self.sample_name+".ReadGroups.bam"
		BaseRecalibrator_metrics = self.output_folder+self.sample_name+".BaseRecalibrator-metrics.txt"
		BaseRecalibrator_bam = self.output_folder+self.sample_name+".BaseRecalibrator.bam"
		vcf = self.output_folder+self.sample_name+".vcf.gz"
		########################################
		
		def align_STAR():
			# Align to reference genome (bwa) and sort (samtools)
			
			# 1st pass
			first_pass_command = [" STAR --genomeDir "
								+bwa_index+" "
								+"--readFilesIn "
								+self.fastq1+" "
								+self.fastq2+" "
								+"--genomeLoad LoadAndKeep "
							#	+"--twopassMode Basic "
								+"--runThreadN "
								+self.threads
								]
			
		
			#2nd pass
			second_pass_command = [" STAR --runMode genomeGenerate --genomeDir 2ndpass_"+self.sample_name+"/ " 
								+"--genomeFastaFiles " 
								+reference_fasta+" "
								+"--sjdbFileChrStartEnd SJ.out.tab "
								+"--sjdbOverhang 75 "
							#	+"--genomeLoad LoadAndKeep "
								+"--runThreadN "+self.threads
								]
			#final pass					
			final_alignment_command = [" STAR --genomeDir 2ndpass_"+self.sample_name+"/ "
								+"--readFilesIn "
								+self.fastq1+" "
								+self.fastq2+" "
								+"--runThreadN "
								+self.threads+" "
								+"--chimJunctionOverhangMin 15 "
								+"--chimSegmentMin 15 "
								+"--outSAMtype BAM SortedByCoordinate "
								+"--outFileNamePrefix "
								+self.output_folder+self.sample_name+".STAR."
								]
			
			star_command = '"'+"".join(first_pass_command+[" ; rm -r 2ndpass_"+self.sample_name+" ; mkdir 2ndpass_"+self.sample_name+" ; "]+second_pass_command+[" ; "]+final_alignment_command)+'"'
			self.run_in_docker([star_command],"star", self.threads)							
			
			
			
				
		
		def addReadGroups():
			# Add read groups (GATK)
			parameters_dict = {
								"input":unsorted_bam,
								"output":ReadGroups_bam,
								"RGLB":self.lib_ID,
								"RGPL":self.pl_ID,
								"RGPU":self.pu_ID,
								"RGSM":self.sample_name
								}
			ReadGroups_log = self.output_folder+self.sample_name+".ReadGroups.log"
			gatk_docker("gatk_add_read_groups", parameters_dict,
						ReadGroups_log, self.ram,self.docker_images_dict["gatk"], )
	
		def buildRecalibrator():
			# Base recalibrator (GATK)
			parameters_dict = {
								"input":sorted_bam,
								"output":BaseRecalibrator_metrics,
								"reference":reference_fasta,
								"known-sites":dbsnp_vcf,
								"use-original-qualities":"true"
								}
			BaseRecalibrator_log = output_folder+self.sample_name+".BaseRecalibrator.log"
			gatk_docker("gatk_build_recalibrator", parameters_dict,
						BaseRecalibrator_log, self.ram,self.docker_images_dict["gatk"])

			
		def applyRecalibrator():
			# Base recalibrator - applying model (GATK)
			parameters_dict = {
								"input":sorted_bam,
								"output":BaseRecalibrator_bam,
								"bqsr":BaseRecalibrator_metrics,
								"use-original-qualities":"true",
								"static-quantized-quals":"10",
								"static-quantized-quals":"20",
								"static-quantized-quals":"30"
								}
			ApplyBQSR_log = output_folder+self.sample_name+".ApplyBQSR.log"
			gatk_docker("gatk_apply_recalibrator", parameters_dict,
						ApplyBQSR_log, self.ram,self.docker_images_dict["gatk"])
		
		def variantCalling_HaplotypeCaller(input_bam):
			# Haplotype caller
			parameters_dict = {
								"input":input_bam,
								"output":input_bam.replace(".bam",".vcf.gz"),
								"dont-use-soft-clipped-bases":"true",
								"stand-call-conf":"20",
								"reference":reference_fasta
							#	"emit-reference-confidence":"true",
							#	"max-alternate-alleles":"3"
								}
			HaplotypeCaller_log = input_bam+".HaplotypeCaller.log"
			gatk_docker("gatk_haplotype_caller", parameters_dict,
						HaplotypeCaller_log, self.ram,self.docker_images_dict["gatk"])
						
		

			
		def validateSam():
			# Validate alignment file integrity
			parameters_dict = {
								"input":output_folder+self.sample_name+".bwa.bam",
								"MODE":"SUMMARY"
								}
			validate_log = output_folder+self.sample_name+".validateSamFile.log"
			gatk_docker(gatk_validate_sam, parameters_dict,
						HaplotypeCaller_log, self.ram,self.docker_images_dict["gatk"])

		def indexBam(path):
			if check_path(path) == True:
				index_cmd = ['samtools','index',path]
				index_log= output_folder+self.sample_name+".samtoolsindex.log"
				#self.run_in_docker(index_cmd, stderr=index_log)
				self.run_in_docker(index_cmd,"samtools", self.threads)
			else:
				raise RuntimeError("BaseRecalibrator bam file not found")
				quit()
		def sortBam(bam_to_sort, sorted_bam, t_n):
			if check_path(bam_to_sort) == True:
				sort_cmd = ['samtools','sort', bam_to_sort]
				sort_cmd += ['-o', sorted_bam]
				sort_cmd += ['-@', self.threads]
				sort_log= self.output_folder+self.sample_name+".samtoolssort.log"
				#self.run_in_docker(index_cmd, stderr=index_log)
				self.run_in_docker(sort_cmd, "samtools", self.threads)
			else:
				raise RuntimeError("unsorted bam file not found")
				quit()		
		
		def processVcf(input_vcf, filename):
			cmd= ['/bin/bash -c "zcat '+input_vcf+' | vt decompose -s - | vt normalize -q - -n -r '+reference_fasta+' 2> /dev/null " | bgzip -c > '+filename+'normalized.vcf.gz']
			self.run_in_docker(cmd, "vt", self.threads)
		
		def extractChromosome(bam_file, chromosome):
			index_cmd = ['samtools','view', bam_file, "chr"+chromosome, "-b","-o", bam_file.replace("bam","chr"+chromosome+".bam")]
			index_cmd += ['-@', self.threads]
			index_log= self.output_folder+self.sample_name+".samtools_split.log"
			#self.run_in_docker(index_cmd, stderr=index_log)
			self.run_in_docker(index_cmd,"samtools", self.threads)
		
		def joinChromosomeVcfs():
			joint_vcf = self.output_folder+self.sample_name+".STAR.vcf.gz"
			partial_vcf_list = []
			for chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']:
				partial_vcf_list.append(BaseRecalibrator_bam.replace("bam","chr"+chromosome+".vcf.gz"))
			os.system("bcftools concat -o "+joint_vcf+" -O z "+" ".join(partial_vcf_list))	
					
		#***********************************
		# Pipeline workflow ###
		#***********************************
		
		# Check if alignment file exists
		if check_path(unsorted_bam)==False:
			self.logger.info(self.sample_name+" sample: align")
			align_STAR()
			
		else:
			self.logger.info("alignment file present, skipping")
#			validateSam()
		
		# Check if readgroups were added
		if check_path(self.output_folder+self.sample_name+".ReadGroups.bam")==False:
			self.logger.info(self.sample_name+" sample: Add readgroups")
			addReadGroups()
		else:	
			self.logger.info("readgroups ok, skipping...")
	
	#	#Sort Bam
	#	self.logger.info("sort")
	#	sortBam(ReadGroups_bam, sorted_bam, "foo")
	#	print ("")
		
		# Check if recalibrator model was built
		if check_path(self.output_folder+self.sample_name+".BaseRecalibrator-metrics.txt")==False:
			buildRecalibrator()
		else:	
			self.logger.info("recalibrator model already build, skipping...")	

		# Check if recalibrator model was applied
		if check_path(self.output_folder+self.sample_name+".BaseRecalibrator.bam")==False:
			self.logger.info(self.sample_name+" sample: Base recalibrator apply")
			applyRecalibrator()
		else:	
			self.logger.info("recalibrator model already applied, skipping...")		
		 
		# Index bam before variant calling
		indexBam(BaseRecalibrator_bam)
		
		# Split BAM to chromosomes for parallel variant calling
	#	if chromosome_parallel==True:
		for chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']:
			extractChromosome(BaseRecalibrator_bam, chromosome)
		
		for chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']:
			indexBam(BaseRecalibrator_bam.replace("bam","chr"+chromosome+".bam"))	
		
		proc = []
		for chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']:
			p = Process(target = variantCalling_HaplotypeCaller, args=[BaseRecalibrator_bam.replace("bam","chr"+chromosome+".bam")])
			p.start()
			proc.append(p)

		for p in proc:
			p.join()
		
		# Join chromosome vcf together
		
		joinChromosomeVcfs()
			
		# PROCESS VCF
		vcf_file = self.output_folder+self.sample_name+".STAR.vcf.gz"
		if check_path(vcf_file)==True:
			self.logger.info(self.sample_name+" sample: vcf process")		
			processVcf(vcf_file, vcf_file.replace("vcf.gz",""))
		
		# Remove intermediate GATK-produced bams
		if self.cleanup=="YES":
			try:
				subprocess.Popen("rm "+self.output_folder+self.sample_name+".MarkDuplicates.bam",shell=True,stderr=subprocess.PIPE)
			except:
				pass
			try:
				subprocess.Popen("rm "+self.output_folder+self.sample_name+".ReadGroup.bam",shell=True,stderr=subprocess.PIPE)
			except:
				pass
			try:
				subprocess.Popen("rm "+self.output_folder+self.sample_name+".BaseRecalibrator.bam",shell=True,stderr=subprocess.PIPE)
			except:
				pass
				
		self.logger.info("sample "+self.sample_name+" pipeline finished.")	
