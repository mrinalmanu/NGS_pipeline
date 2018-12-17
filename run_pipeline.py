import sys
sys.path.append('/home/NGS_pipeline/pipelines')
from pipelines import *
from GATK4_germline_v1_all_docker import GATK4_germline_v1
from GATK4_somatic_TN_v1 import GATK4_somatic_TN_v1
#from GATK4_somatic_TO_v1 import GATK4_somatic_TO_v1
from GATK_varscan_TN_v1 import GATK_varscan_TN_v1
from multiprocessing import Process
import time
import multiprocessing

################ INPUT SAMPLE FILE ######################
sample_file=open(sys.argv[1]).read().splitlines()
docker_images_versions=open("/home/NGS_pipeline/docker_images_versions.txt").read().split('#')
minimum_cores = 10

################ PARSE DOCKER IMAGE VERSIONS FROM FILE ##################
latest_dict={}
pipelines_tools_dict={}
for paragraph in docker_images_versions:
	if paragraph =="":
		continue
	else:
		paragraph=paragraph.split('\n')
		if paragraph[0]=="MOST_RECENT_IMAGES":
			for item in paragraph[1:]:
				item=item.split()
				if len(item)!=2:
					continue
				else:
					latest_dict.update({item[0]:item[1]})
		else:
			pipeline_dict={}
			for item in paragraph[1:]:
				item=item.split()
				if len(item)!=2:
					continue
				else:
					if item[1]=="latest":
						if item[0] in latest_dict:
							item[1]=latest_dict[item[0]]
						else:
							#print >> sys.stderr,"ERROR: tool "+item[0]+" defined in pipeline "+paragraph[0]+" does not have specified latest version in "+sys.argv[2]	
							quit()
					pipeline_dict.update({item[0]:item[1]})
			pipelines_tools_dict.update({paragraph[0]:pipeline_dict})
###################### END IMAGES PARSING #############################
###################### PARSE SAMPLE LIST ##############################
instances_list=[]
for line in sample_file[1:]:
	if list(line)[0]=="#":
		continue
	lline=line.split(',')
	sample_ID = lline[0]
	pipeline_ID = lline[1]
	genome_build = lline[2]
	fastq1 = lline[3]
	fastq2 = lline[4]
	fastq1_tumor = lline[5]
	fastq2_tumor = lline[6]
	cleanup = lline[7]
	lib = lline[8]
	pl = lline[9]
	pu = lline[10]
	
	print (fastq1)
	print (fastq2)
	print (fastq1_tumor)
	print (fastq2_tumor)
##################### END SAMPLE LIST PARSING #########################

############### CREATE INSTANCES OF PIPELINE CLASS FOR SAMPLES ########
	if pipeline_ID == "GATK4_germline_v1":
		instances_list.append(GATK4_germline_v1(sample_ID, genome_build, fastq1, fastq2, fastq1_tumor, fastq2_tumor, cleanup,lib,pl,pu, pipelines_tools_dict[pipeline_ID]))
	elif pipeline_ID == "GATK4_somatic_TN_v1":
		instances_list.append(GATK4_somatic_TN_v1(sample_ID, genome_build, fastq1, fastq2, fastq1_tumor, fastq2_tumor, cleanup,lib,pl,pu, pipelines_tools_dict[pipeline_ID]))
	elif pipeline_ID == "GATK4_somatic_TO_v1":
		instances_list.append(GATK4_somatic_TO_v1(sample_ID, genome_build, fastq1, fastq2, fastq1_tumor, fastq2_tumor, cleanup,lib,pl,pu, pipelines_tools_dict[pipeline_ID]))
	elif pipeline_ID == "GATK_varscan_TN_v1":
		instances_list.append(GATK_varscan_TN_v1(sample_ID, genome_build, fastq1, fastq2, fastq1_tumor, fastq2_tumor, cleanup,lib,pl,pu, pipelines_tools_dict[pipeline_ID]))
	
	else:
		print ("Sample "+sample_ID+": pipeline "+pipeline+" not available")


##################### SEND SAMPLES TO PIPELINE #########################
def machine_busy(nr_threads): ##Check load on machine before submiting job
	import psutil
	import numpy as np
	cores_available= multiprocessing.cpu_count()
		
	busycores_matrix=[]	
	for nr in range(5):
		values = psutil.cpu_percent(interval=3, percpu=True)
		busycores = 0	
		for core in values:
			if core >= 50:
				busycores += 1
		busycores_matrix.append(busycores)

	
	print "all available cores: ",cores_available	
	print "nr of busy cores: ",np.mean(busycores_matrix)
	if cores_available - np.mean(busycores_matrix)< nr_threads:
		print "machine busy"
		return True
	else:
		return False


def run(instance):
	instance.run_pipeline()
################## RUN SAMPLE IF ENOUGH CORES ###########
for instance in instances_list:
	wait_for_cores = machine_busy(minimum_cores)
	while wait_for_cores==True:
		time.sleep(5)
		wait_for_cores = machine_busy(minimum_cores)
	print "CPU resources sufficient, submitting sample", instance.sample_name
	proc = []
	p = Process(target=run, args=[instance])
	p.start()
	proc.append(p)

for p in proc:
	p.join()
