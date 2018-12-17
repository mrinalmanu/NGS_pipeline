### GATK4 somatic pipeline for tumor-normal pairs version 1.0
### Pipeline created: 10/18/2018

import sys
import subprocess
import logging
import os
from gatk4_docker import gatk_docker
from multiprocessing import Process
from config_file import *


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
			
class GATK4_somatic_TN_v1():
	def __init__(self, sample_name, genome_build, f1, f2,tf1,tf2, cleanup,lib_ID, pl_ID, pu_ID, docker_images_dict):

		self.sample_name= sample_name
		self.fastq1_tumor = tf1
		self.fastq2_tumor = tf2
		self.fastq1 = f1
		self.fastq2 = f2
		if (f1 =="") or (f2 =="") or (tf1 =="") or (tf2 ==""):
			print >>sys.stderr, t_n_sample_name+": FASTQ file is missing, aborting."
			quit()
		self.input_folder= input_folder
		self.threads = max_nr_threads
		self.ram = "50000"
		self.open_files = []
		self.genome_build = genome_build
		self.cleanup = cleanup
		self.docker_images_dict = docker_images_dict
		self.lib_ID = lib_ID
		self.pl_ID = pl_ID
		self.pu_ID = pu_ID
		self.output_folder = output_folder+"/"+self.sample_name+"/GATK4_somatic/output/"
		
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
	
	def run_in_docker(self, cmd,t_n, image, stdout=None, stderr=None):
		""" Run a command inside docker container"""
		
		if len(cmd)<2:
			container_name = self.sample_name+"_"+t_n+"_vt"
		else:
			container_name = self.sample_name+"_"+t_n+"_"+cmd[0]
	
		checkContainer(container_name)
		dcmd = ["docker", "run","--name",container_name,
						"-v", "{}:{}".format(self.input_folder, self.input_folder),
						"-v", "{}:{}".format(output_folder, output_folder),
						"-v", "{}:{}".format(reference_folder, reference_folder),
						self.docker_images_dict[image]]
		dcmd += cmd
	
		if stdout is not None:
			stdout = open(stdout,"w")#

			self.open_files.append(stdout)
		if stderr is not None:
			stderr = open(stderr,"w")
			self.open_files.append(stderr)
		self.logger.info("Running: "+" ".join(dcmd))
		os.system(" ".join(dcmd))
		checkContainer(container_name)
	
	
	def run_pipeline(self):
		### Genome build selection ###
		if self.genome_build == "GRCh38":
			bwa_index = bwa_index_GRCh38
			reference_fasta = reference_fasta_GRCh38
			dbsnp_vcf = dbsnp_vcf_GRCh38
			print "VEP REFERENCE FOR GRCH38 not present, exiting."
			quit()
			
		if self.genome_build == "hg19":
			bwa_index = bwa_index_hg19
			reference_fasta = reference_fasta_hg19
			dbsnp_vcf = dbsnp_vcf_hg19
			vep_reference_folder = vep_reference_folder_hg19
		
		
		def variantCalling_mutect2():
			# variant caller
			vcf = self.output_folder+self.sample_name+".mutect2.vcf.gz"
			parameters_dict = {
								"input_tumor":self.output_folder+self.sample_name+"_r_tumor.bam",
								"input_normal":self.output_folder+self.sample_name+"_r_normal.bam",
								"normal-sample":self.sample_name+"_r_normal",
								"tumor-sample":self.sample_name+"_r_tumor",
								"output":vcf,
								"reference":reference_fasta,
							#	"--disable-read-filter":""
							#	"dbsnp":
							#	"emit-reference-confidence":"true",
							#	"max-alternate-alleles":"3"
								}
			mutect2_log = self.output_folder+self.sample_name+".Mutect2.log"
			gatk_docker("gatk_mutect2", parameters_dict,
						mutect2_log, self.ram,self.docker_images_dict["gatk"])
		
		def variantCalling_mutect2_partial(chromosome):
			# variant caller
			vcf = self.output_folder+self.sample_name+".chr"+chromosome+".mutect2.vcf.gz"
			parameters_dict = {
								"input_tumor":self.output_folder+self.sample_name+"_r_tumor.chr"+chromosome+"_RG.bam",
								"input_normal":self.output_folder+self.sample_name+"_r_normal.chr"+chromosome+"_RG.bam",
								"normal-sample":self.sample_name+"_r_normal.chr"+chromosome+"_RG",
								"tumor-sample":self.sample_name+"_r_tumor.chr"+chromosome+"_RG",
								"output":vcf,
								"reference":reference_fasta,
							#	"--disable-read-filter":""
							#	"dbsnp":
							#	"emit-reference-confidence":"true",
							#	"max-alternate-alleles":"3"
								}
			mutect2_log = self.output_folder+self.sample_name+".chr"+chromosome+".Mutect2.log"
			gatk_docker("gatk_mutect2", parameters_dict,
						mutect2_log, self.ram,self.docker_images_dict["gatk"])				

		
	
		def validateSam(t_n):
			# Validate alignment file integrity
			parameters_dict = {
								"input":self.output_folder+t_n_sample_name+".bwa."+t_n+".bam",
								"MODE":"SUMMARY"
								}
			validate_log = self.output_folder+t_n_sample_name+".validateSamFile."+t_n+".log"
			gatk_docker(gatk_validate_sam, parameters_dict,
						HaplotypeCaller_log, self.ram,self.docker_images_dict["gatk"])
		
		def processVcf(input_vcf, filename):
			# decompose and normalize variants in vt 
			cmd= ['/bin/bash -c "zcat '+input_vcf+' | vt decompose -s - | vt normalize -q - -n -r '+reference_fasta+' 2> /dev/null " | bgzip -c > '+filename+'normalized.vcf.gz']
			self.run_in_docker(cmd, "asa","vt")
		
		def filterMutectCalls(vcf_file):
			# perform variant filtering
			parameters_dict = {
								"variant":vcf_file,
								"output":vcf_file.replace(".vcf.gz",".filtered.vcf.gz")
								}
			filter_log = self.output_folder+self.sample_name+".filter.log"
			gatk_docker("gatk_filter_mutect", parameters_dict,
						filter_log, self.ram,self.docker_images_dict["gatk"])
		
		def annotateVep(vcf_file):
			# annotation in Variant Effect Predictor
			annotate_cmd = ["vep",
							"--af", 
							"--af_1kg", 
							"--af_esp",
							"--af_gnomad", 
							"--appris", 
							"--biotype", 
							"--check_existing", 
							"--pick", 
							"--distance 5000", 
							"--merged", 
							"--polyphen b", 
							"--pubmed", 
							"--regulatory", 
							"--sift b", 
							"--species homo_sapiens", 
							"--symbol", 
							"--tsl", 
							"--cache", 
							"--offline", 
							"--dir_cache", vep_reference_folder, 
							"--fork", self.threads, 
							"--input_file",vcf_file, 
							"--force_overwrite", 
							"--fasta", reference_fasta, 
							"-o", vcf_file.replace(".vcf.gz",".annotated.vcf.gz"), 
							"--vcf"]
								
			self.run_in_docker(annotate_cmd,"asa", "vep")
		
		def indexVcf(vcf_file):
			# index vcf file
			index_cmd=["bcftools","index","-t","-f",vcf_file]
			self.run_in_docker(index_cmd, "asa","bcftools")
		
		def sortVcf(vcf_file):
			# index vcf file
			index_cmd=["bcftools","sort","-O z","-o",vcf_file.replace(".vcf.gz",".sort.vcf.gz"),vcf_file]
			self.run_in_docker(index_cmd, "asa","bcftools")	
	
		def joinChromosomeVcfs():
			joint_vcf = self.output_folder+self.sample_name+".mutect2.vcf.gz"
			partial_vcf_list = []
			for chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']:
				
				if check_path(self.sample_name+"_samplenames_txt")==True:
					os.system("rm "+self.sample_name+"_samplenames_txt")
				os.system("echo "+self.sample_name+"_normal >> "+self.sample_name+"_samplenames_txt")
				os.system("echo "+self.sample_name+"_tumor >> "+self.sample_name+"_samplenames_txt")
				os.system("bcftools reheader -s "+self.sample_name+"_samplenames_txt "+self.output_folder+self.sample_name+".chr"+chromosome+".mutect2.vcf.gz -o "+self.output_folder+self.sample_name+".chr"+chromosome+".reheader.vcf.gz")
				os.system("rm "+self.sample_name+"_samplenames_txt")
				
			for chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']:
				partial_vcf_list.append(self.output_folder+self.sample_name+".chr"+chromosome+".reheader.vcf.gz")
			self.logger.info("Joining vcf: Running: bcftools concat -o "+joint_vcf+" -O z "+" ".join(partial_vcf_list))
			os.system("bcftools concat -o "+joint_vcf+" -O z "+" ".join(partial_vcf_list))		
				
		def preprocess(t_n):
		#### pre-processing for tumor and normal separately ###
			t_n_sample_name= self.sample_name+"_"+t_n
			### Pre-process auxiliary files ###
			unsorted_bam = self.output_folder+self.sample_name+".bwa."+t_n+".bam"
			sorted_bam = self.output_folder+self.sample_name+".bwa.sorted."+t_n+".bam"
			sort_err = self.output_folder+self.sample_name+".sort."+t_n+".err"
			markDuplicates_bam = self.output_folder+self.sample_name+".MarkDuplicates."+t_n+".bam"
			markDuplicates_metrics = self.output_folder+self.sample_name+".MarkDuplicates-metrics."+t_n+".txt"
			ReadGroups_bam = self.output_folder+self.sample_name+".ReadGroups."+t_n+".bam"
			BaseRecalibrator_metrics = self.output_folder+self.sample_name+".BaseRecalibrator-metrics."+t_n+".txt"
			BaseRecalibrator_bam = self.output_folder+self.sample_name+"_r_"+t_n+".bam"
			########################################
			
			def align_BWA(t_n):
				#Align to reference genome (bwa) and sort (samtools)
				if t_n =="tumor":
					fastq1 = self.fastq1_tumor
					fastq2 = self.fastq2_tumor
				if t_n =="normal":
					fastq1 = self.fastq1
					fastq2 = self.fastq2	
			
				self.run_in_docker(["python align_bwa_sort_samtools.py "
					+reference_fasta+" "
					+fastq1+" "
					+fastq2+" "
					+self.threads+" "
					+self.output_folder+t_n_sample_name+" "
					+unsorted_bam],
					t_n,
					"bwa_samtools131")
			
			def markDuplicates(t_n):
				# Mark duplicates (GATK)
				parameters_dict = {
									"input":unsorted_bam,
									"output":markDuplicates_bam,
									"mark_dupl_metrics":markDuplicates_metrics,
									"optical_duplicate_pixel_dist":"2500",
									"assume_sort_order":"queryname",
									"clear_DT":"false",
									"add_pg_tag_to_reads":"false"
										}
				markDuplicates_log = self.output_folder+t_n_sample_name+".MarkDuplicates.log"
				gatk_docker("gatk_mark_duplicates",parameters_dict,
							markDuplicates_log, self.ram,self.docker_images_dict["gatk"])
			
			def addReadGroups(t_n):
				# Add read groups (GATK)
				parameters_dict = {
									"input":markDuplicates_bam,
									"output":ReadGroups_bam,
									"RGLB":self.lib_ID,
									"RGPL":self.pl_ID,
									"RGPU":self.pu_ID,
									"RGSM":self.sample_name+"_r_"+t_n
									}
				ReadGroups_log = self.output_folder+t_n_sample_name+".ReadGroups."+t_n+".log"
				gatk_docker("gatk_add_read_groups", parameters_dict,
							ReadGroups_log, self.ram,self.docker_images_dict["gatk"])
			
			def addReadGroups_partial_bam(chromosome, input_bam):
				# Add read groups (GATK)
				parameters_dict = {
									"input":input_bam,
									"output":input_bam.replace(".bam","_RG.bam"),
									"RGLB":self.lib_ID,
									"RGPL":self.pl_ID,
									"RGPU":self.pu_ID,
									"RGSM":input_bam.replace(".bam","_RG").split("/")[-1]
									}
				ReadGroups_log = self.output_folder+t_n_sample_name+".ReadGroups.log"
				gatk_docker("gatk_add_read_groups", parameters_dict,
							ReadGroups_log, self.ram,self.docker_images_dict["gatk"])				
		
			def buildRecalibrator(t_n):
				# Base recalibrator (GATK)
				parameters_dict = {
									"input":sorted_bam,
									"output":BaseRecalibrator_metrics,
									"reference":reference_fasta,
									"known-sites":dbsnp_vcf,
									"use-original-qualities":"true"
									}
				BaseRecalibrator_log = self.output_folder+t_n_sample_name+".BaseRecalibrator."+t_n+".log"
				gatk_docker("gatk_build_recalibrator", parameters_dict,
							BaseRecalibrator_log, self.ram,self.docker_images_dict["gatk"])
	
				
			def applyRecalibrator(t_n):
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
				ApplyBQSR_log = self.output_folder+t_n_sample_name+".ApplyBQSR."+t_n+".log"
				gatk_docker("gatk_apply_recalibrator", parameters_dict,
							ApplyBQSR_log, self.ram,self.docker_images_dict["gatk"])
			
			def sortBam(unsorted_bam, sorted_bam, t_n):
				if check_path(unsorted_bam) == True:
					sort_cmd = ['samtools','sort','-@',self.threads, unsorted_bam]
					sort_cmd += ['-o', sorted_bam]
					sort_log= self.output_folder+self.sample_name+".samtoolssort."+t_n+".log"
					#self.run_in_docker(index_cmd, stderr=index_log)
					self.run_in_docker(sort_cmd,t_n, self.docker_images_dict["samtools"])
				else:
					raise RuntimeError("unsorted bam file not found")
					quit()
		
		
			def indexBam(path, t_n):
				if check_path(path) == True:
					index_cmd = ['samtools','index',path]
					index_log= self.output_folder+self.sample_name+".samtoolsindex."+t_n+".log"
					#self.run_in_docker(index_cmd, stderr=index_log)
					self.run_in_docker(index_cmd,t_n,self.docker_images_dict["samtools"])
				else:
					print path
					raise RuntimeError("BaseRecalibrator bam file not found")
					quit()

			def extractChromosome(bam_file, chromosome):
				index_cmd = ['samtools','view', bam_file, "chr"+chromosome, "-b","-o", bam_file.replace("bam","chr"+chromosome+".bam")]
				index_cmd += ['-@', '40']
				index_log= self.output_folder+self.sample_name+".samtools_split.log"
				#self.run_in_docker(index_cmd, stderr=index_log)
				self.run_in_docker(index_cmd,t_n,self.docker_images_dict["samtools"])
									
			
			#***********************************
			# pre-process workflow ###
			#***********************************
			# ALIGNMENT BWA (Check if alignment file exists)
			if check_path(self.output_folder+self.sample_name+".bwa."+t_n+".bam")==False:
				self.logger.info(self.sample_name+"align")
				align_BWA(t_n)
			
			else:
				self.logger.info(self.sample_name+t_n+" sample:alignment done, skipping...")
			print ("")
			
			# MARK DUPLICATES (Check if duplicated reads were marked)
			if check_path(self.output_folder+self.sample_name+".MarkDuplicates."+t_n+".bam")==False:
				self.logger.info("mark_dupl")
				markDuplicates(t_n)
			else:	
				self.logger.info(t_n+" sample:duplicates marked, skipping...")
			print ("")
			
			# ADD READGROUPS (Check if readgroups were added)
			if check_path(self.output_folder+self.sample_name+".ReadGroups."+t_n+".bam")==False:
				self.logger.info(self.sample_name+"add_readgroups")
				addReadGroups(t_n)
			else:	
				self.logger.info(self.sample_name+t_n+" sample: readgroups ok, skipping...")
			print ("")
			
			#Sort Bam
			if check_path(sorted_bam)==False:
				self.logger.info(self.sample_name+"sort")
				sortBam(ReadGroups_bam, sorted_bam, t_n)
			print ("")
			
			# BASE RECALIBRATOR (Check if recalibrator model was built)
			if check_path(self.output_folder+self.sample_name+".BaseRecalibrator-metrics."+t_n+".txt")==False:
				self.logger.info(self.sample_name+"build_recalib")
				buildRecalibrator(t_n)
			else:	
				self.logger.info(t_n+" sample: recalibrator model already build, skipping...")	
			print ("")
			
			# APPLY BQSR (Check if recalibrator model was applied)
			if check_path(BaseRecalibrator_bam)==False:
				self.logger.info("apply_recalib")		
				applyRecalibrator(t_n)
			else:	
				self.logger.info(self.sample_name+t_n+" sample: recalibrator model already applied, skipping...")		
			print ("")
			 
			# Index bam before variant calling
			if check_path(BaseRecalibrator_bam+".bai")==False:
				
				self.logger.info(self.sample_name+t_n+" sample: indexing")		
				indexBam(BaseRecalibrator_bam, t_n)
			else:
				self.logger.info(self.sample_name+t_n+" sample: indexed, skipping...")		
			print ("")
			
			# # Split BAM to chromosomes for parallel variant calling
			# #	if chromosome_parallel==True:
			# for chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']:
				# extractChromosome(BaseRecalibrator_bam, chromosome)
		
			# for chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']:
				# addReadGroups_partial_bam(chromosome, BaseRecalibrator_bam.replace("bam","chr"+chromosome+".bam"))
		
			# for chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']:
				# indexBam(BaseRecalibrator_bam.replace(".bam",".chr"+chromosome+"_RG.bam"),t_n)	
			
			
	
			#### END OF pre-processing function for tumor and normal separately ###
	
		
		# Initiate pre-process for tumor and normal in parallel
		proc = []
		for t_n in ['tumor','normal']:
			p = Process(target = preprocess, args=[t_n])
			p.start()
			proc.append(p)

		for p in proc:
			p.join()
		
		### PARALLEL TUMOR- NORMAL VARIANT CALLING ###########################################################
		
		# proc = []
		# for chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']:
			# p = Process(target = variantCalling_mutect2_partial, args=[chromosome])
			# p.start()
			# proc.append(p)

		# for p in proc:
			# p.join()

#		# VARIANT CALLING MUTECT2 WITHOUT PARALLELISATION (Check if variant calling was performed)
#		if check_path(self.output_folder+self.sample_name+".mutect2.vcf.gz")==False:
#			self.logger.info(self.sample_name+" Mutect2 variant calling")
#			variantCalling_mutect2()
#		else:	
#			self.logger.info(self.sample_name+"variant calling already done, skipping...")
		
		
		# join partial chromosome vcfs
#		joinChromosomeVcfs()

		
		
		# PROCESS VCF in VT (check if already processed)
		if check_path(self.output_folder+self.sample_name+".mutect2.normalized.vcf.gz")==False:
			vcf_file = self.output_folder+self.sample_name+".mutect2.vcf.gz"
			if check_path(vcf_file)==True:
				processVcf(vcf_file, vcf_file.replace("vcf.gz",""))
		else:
			self.logger.info(self.sample_name+": vcf processing already done, skipping...")
		
		sortVcf(self.output_folder+self.sample_name+".mutect2.normalized.vcf.gz")
		
		# Index normalized vcf (bcftools index vcf)
		if check_path(self.output_folder+self.sample_name+".mutect2.normalized.sort.vcf.gz.tbi")==False:
			self.logger.info(self.sample_name+" sample: indexing")		
			indexVcf(self.output_folder+self.sample_name+".mutect2.normalized.sort.vcf.gz")
		else:
			self.logger.info(self.sample_name+" sample: normalized vcf indexed, skipping...")		
			print ("")
		
		
		## FILTER MUTECT CALLS (Check if already filtered)
		if check_path(self.output_folder+self.sample_name+".mutect2.normalized.sort.filtered.vcf.gz")==False:		
			self.logger.info(self.sample_name+" Mutect2 calls filtering")
			filterMutectCalls(self.output_folder+self.sample_name+".mutect2.normalized.sort.vcf.gz")
		else:
			self.logger.info(self.sample_name+": calls filtered, skipping...")
			
		## ANNOTATE CALLS WITH VEP (Check if already annotated)
		if check_path(self.output_folder+self.sample_name+".mutect2.normalized.sort.filtered.vep.vcf.gz")==False:		
			self.logger.info(self.sample_name+" annotating with VEP")
			annotateVep(self.output_folder+self.sample_name+".mutect2.normalized.sort.filtered.vcf.gz")
		else:
			self.logger.info(self.sample_name+": annotated, skipping...")	
				
		# Remove intermediate GATK-produced bams
		
		
			
		
		if self.cleanup=="YES":
			try:
				subprocess.Popen("rm "+self.output_folder+t_n_sample_name+".MarkDuplicates."+t_n+".bam",shell=True,stderr=subprocess.PIPE)
			except:
				pass
			try:
				subprocess.Popen("rm "+self.output_folder+t_n_sample_name+".ReadGroup."+t_n+".bam",shell=True,stderr=subprocess.PIPE)
			except:
				pass
			try:
				subprocess.Popen("rm "+self.output_folder+t_n_sample_name+".BaseRecalibrator."+t_n+".bam",shell=True,stderr=subprocess.PIPE)
			except:
				pass
				
		self.logger.info("sample "+self.sample_name+"_"+t_n+" pipeline finished.")	
	
