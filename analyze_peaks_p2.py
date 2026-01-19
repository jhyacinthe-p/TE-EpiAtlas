import pickle#cPickle as pickle
import sys
print(sys.version)
import pandas as pd
import numpy
import argparse
import os
import time
import itertools
from pybedtools import BedTool as bedtool
import pybedtools
import multiprocessing

#####
# Author:   Jeffrey Hyacinthe
# Email:    jeffrey.hyacinthe@mail.mcgill.ca
# Date:     2020/10/30
# Note:     Python version == 3.6
#####
"""
Intersects a given peak file with TE and returns the overlap count in comparison to associated random simulated samples.
returns a csv table, a csv of peak distribution relative to tss and a bed file of intersection with repeatmaker

Usage:
analyze_peaks.py input_file -g genome -l length -n size -r iteration_count -o output_file -p processes

input_file: [str] standard bed file
genome: default=38 [int] 19 or 38
lenght: default=200 [int] lenght of region directly around peak to consider
size: default=-1 [int] number of simulated peaks per trial. -1 sets to input peak count
iteration_count: default=10 [int] random iteration count, recommended 1000 (Slow).
output_file: default=analyze_peaks_output.csv [str] csv table output file
processes: parralelize the work using multiple processes
"""

if __name__ == "__main__":
	#=======================================================
	#               Class & Functions
	#=======================================================
	def memory_usage_psutil():
		# return the memory usage in MB
		import psutil
		process = psutil.Process(os.getpid())
		mem = process.memory_info()[0] / float(2 ** 20)
		return mem

	def timestamp(text):
		this_time=time.time()
		global last_time
		print(text+" | Since last: "+str(this_time-last_time)+"s | Total : "+str(this_time-start_time)+"s" )
		last_time = this_time
		return this_time

	def linear_mean(x, _mean, n):
		if n==1:
			return x
		return (x + (n-1)*_mean)/n

	def linear_S(x, _S, n, _mean, _mean_prev):
		return _S + (x - _mean_prev)*(x - _mean)

	def linear_sd(x, _sd, n, _mean, _mean_prev):
		_S = (_sd**2)*n
		S = linear_S(x, _S, n, _mean, _mean_prev)
		return (S/n)**(1/2)

	def linear_var(x, _var, n, _mean):
		return ((n-2)/(n-1))*_var+((x-_mean)**2)/n


	def add_to_joint_table(te_table, _joint_table):
		"""
		Add a TE table to a common, cummulative table
		Calculates the mean and sd for instance count and bases of TE linearly.
		For the actual sample, it reports the actual count since there is only one sample.
		For the random simulation, it will report the mean and sd of all iterations
		"""
		global counter
		with counter.get_lock():
			counter.value += 1
		n = counter.value
		i=0
		joint_table={}
		with(open(te_table.fn, "r")) as f:
			for line in f:
				te, count, bases = line.strip().split()
				count = int(count)
				bases = int(bases)
				#	#Calculate Variances firsts to use the "previous" (n-1) mean value
				joint_table[te] = _joint_table[te]
				mean_prev = joint_table[te]["count"]
				joint_table[te]["count"] = linear_mean(count, joint_table[te]["count"], n)
				joint_table[te]["count_sd"] = linear_sd(count, joint_table[te]["count_sd"], n, joint_table[te]["count"], mean_prev)
				mean_prev = joint_table[te]["bases"]
				joint_table[te]["bases"] = linear_mean(bases, joint_table[te]["bases"], n)
				joint_table[te]["bases_sd"] = linear_sd(bases, joint_table[te]["bases_sd"], n, joint_table[te]["bases"], mean_prev)
				joint_table[te]["more_in_sample"] += (sample_te_table[te]["count"]>count)
				_joint_table[te] = joint_table[te]
				i+=1


	def tableit(in_bed):
		#Turn a bed file into a list of TE with instance count and cummulative base count
		#Start at 1 indexing
		field_count = len(in_bed[0].fields)
		te_table = in_bed.groupby(g=[field_count-3], c=field_count, o=['count','sum'])
		return te_table

	def merge_keeplongest(file):
		# Passes through bed file and keeps only longest overlap of duplicate A entries
		# Assumes sorted bedfile
		_tmp = []
		with open(file, 'r') as f:
			line = f.readline()
			line_next = line
			line_held = line
			while True:
				line = line_next
				line_next = f.readline()
				if not line_next:
					_tmp.append(line_held.strip().split())
					break  # EOF
				line1 = line.strip().split()
				line2 = line_next.strip().split()
				#If next line not same peak, add held line and move to next
				if line1[0:3]!=line2[0:3]:
					_tmp.append(line_held.strip().split())
					line_held = line_next
					continue
				#if next  line has larger overlap, set it as held line
				if int(line2[-1])>int(line1[-1]):
					line_held = line_next
		_tmp = sorted(_tmp, key=lambda x:(x[-4], x[0], int(x[1]), int(x[2])))
		return bedtool(_tmp)

	def add_trial(intersected_bed, joint_table):
		"""
		Add a random simulation to the joint cummulative table
		Cleans up the intersection, format it into a count table and delete all associated file
		"""
		print("add trial start")
		cleaned_intersect = merge_keeplongest(intersected_bed.fn)
		te_table = tableit(cleaned_intersect)
		add_to_joint_table(te_table, joint_table)
		#clean up
		for b in [intersected_bed, cleaned_intersect, te_table]:
			if os.path.exists(b.fn):
				os.remove(b.fn)
			else:
				print("The file does not exist: "+str(b.fn))
		print("add trial end")

	def create_rmsk_info(file):
		"""
		From repeat masker file, generate a mapping from repeat name to its class and family
		Count and distribute the repeats (bases and count) in their family and class
		"""
		file_obj = open(file,"r")
		linecount = 0
		rep_to_class = {}
		rep_to_fam = {}
		fam_to_class = {}
		rmsk_data = {"repeat":{},"repeat_instance":{},"total":{"bases":0,"count":0}}
		for line in file_obj:
			cols = line.strip().split("	")
			rep = cols[10]
			rclass = cols[11]
			rfam = cols[12]
			if linecount>0:#skip header
				rstart = int(cols[6])
				rend = int(cols[7])
				base_count = rend-rstart
				repi = rep+"_"+cols[5]+"_"+cols[6]
				#============== Repeat
				if rep not in rmsk_data["repeat"]:
					rmsk_data["repeat"][rep]={"count":0,"bases":0}
				rmsk_data["repeat"][rep]["count"] += 1
				rmsk_data["repeat"][rep]["bases"] += base_count
				#============== Repeat instance
				if repi not in rmsk_data["repeat_instance"]:
					rmsk_data["repeat_instance"][repi]={"count":0,"bases":0}
				rmsk_data["repeat_instance"][repi]["count"] += 1
				rmsk_data["repeat_instance"][repi]["bases"] += base_count
				#=============== Total
				rmsk_data["total"]["count"] += 1
				rmsk_data["total"]["bases"] += base_count
				#Repeat to Family and Class Mapping
				if rep not in rep_to_class:
					rep_to_class[rep] = rclass
				if rep not in rep_to_fam:
					rep_to_fam[rep] = rfam
				if rfam not in fam_to_class:
					fam_to_class[rfam] = rclass
			linecount+= 1
		file_obj.close()
		return rmsk_data, rep_to_class, rep_to_fam, fam_to_class

	def get_annotation_distributions(_inbed):
		"""
		Get peak counts contained in defined region annotations to generate equivalent random samples
		returns peaks within the annotation(annot_bed), the entire annotation regions (annot_zone_bed)
		and counts within annotation (annot_count). 
		"""
		annot_bed = {}
		annot_zone_bed = {}
		annot_count = {}
		annot_count["all"]={}
		annot_count["all"]["count"]=0
		prefix = "lib/annot_zones/annot_zone"
		#Get peak distribution within the interations
		remain_bed=bedtool(_inbed)
		for annot in annot_table:
			fname=prefix+"_"+annot+".sorted.bed"
			annot_zone_bed[annot] = bedtool(fname)
			annot_bed[annot] = remain_bed.intersect(annot_zone_bed[annot], u=True ,sorted=True, stream=True)
			remain_bed = remain_bed.intersect(annot_zone_bed[annot], v=True ,sorted=True)
			annot_count[annot] = {"count":len(annot_bed[annot])}
			print(str(annot)+" has "+str(annot_count[annot]["count"])+" peaks")
			annot_count["all"]["count"] += annot_count[annot]["count"]
		#Turn count distribution to proportions
		for annot in annot_table:
			annot_count[annot]["ratio"] = annot_count[annot]["count"]/(1.0*annot_count["all"]["count"])
		annot_count["all"]["ratio"] = read_count
		#Save the distribution
		df_distribution = pd.DataFrame(annot_count).transpose()
		df_distribution.to_csv(out_csv_pre_ext+"_distribution.csv")
		return annot_bed, annot_zone_bed, annot_count

	def bed_resize_intervals(peak_file, length):
		"""
		resize all intervals in bed to given lenght
		"""
		peak_file_obj = bedtool(peak_file)
		new_line = []
		for interval in peak_file_obj:
			mid = int((interval.start+interval.stop)/2)
			new_start =  max(mid - int(length/2), 0)
			new_end =  mid+ 1 + int(length/2) #should be chromosome max 
			modified_line = '	'.join([interval[0]]+[str(new_start)]+[str(new_end)]+interval[3:])
			new_line.append(modified_line)
		new_bed = bedtool(new_line)#,  from_string=True
		return new_bed

	def get_lenght(in_bed):
		"""
		Get total cummulative lenght of all peaks
		"""
		tsize=0
		peak_file_obj=in_bed
		for interval in peak_file_obj:
			start=interval.start
			end=interval.stop
			tsize+=end-start
		return tsize

	def get_lenght_f(in_bed):
		"""
		Get total cummulative lenght of all peaks. file variant.
		"""
		tsize=0
		peak_file_obj=open(in_bed.fn,"r")
		for interval in peak_file_obj:
			cols = interval.split("	")
			start = int(cols[1])
			end = int(cols[2])
			tsize+=end-start
		peak_file_obj.close()
		return tsize


	def differential_results(sample_te_table, trials_te_table):
		"""
		Obtain counts and other metrics of all repeats through a dict of dict.[repeat][feature]
		returns CSV ready dataframe object
		"""
		results = {}			
		for te in all_repeats:
			results[te] = {}
			#Common Variables
			results[te]["name"] = te
			results[te]["superfamily"] = repeat_to_superfamily[te]
			results[te]["class"] = repeat_to_class[te]
			results[te]["read_count"] = read_count
			results[te]["read_bases"] = read_bases
			results[te]["random_read_count"] = annot_count["all"]["count"]#random_read_count
			results[te]["random_read_bases"] = annot_count["all"]["count"]*length#random_read_bases
			results[te]["iterations"] = iter_count
			#Sample Variables
			results[te]["count"] = sample_te_table[te]["count"]
			results[te]["bases"]  = sample_te_table[te]["bases"]
			#Random simulation variables
			results[te]["time_over"] = trials_te_table[te]["more_in_sample"]
			results[te]["random_count"] = trials_te_table[te]["count"]
			results[te]["random_count_sd"] = trials_te_table[te]["count_sd"]#**(1/2)
			results[te]["random_bases"]  = trials_te_table[te]["bases"]
			results[te]["random_bases_sd"] = trials_te_table[te]["bases_sd"]#**(1/2)
			#Repeat Masker variables
			results[te]["rmsk_count"] = sample_te_table[te]["count"]
			results[te]["rmsk_bases"]  = sample_te_table[te]["bases"]
		#Convert to dataframe and reorder columns
		df_results = pd.DataFrame(results).transpose()
		df_results = df_results[["name","superfamily","class","count","random_count","random_count_sd","time_over","rmsk_count","bases",
		"random_bases","random_bases_sd","rmsk_bases","read_count","read_bases","random_read_count","random_read_bases",
		"iterations"]]
		df_results = df_results.sort_values(by=["class","superfamily","name"])
		return df_results

	def build_preshuffle():
		"""
		Generate a template bed file for each annotation region with the right lenght and peak count
		will be shuffled in each iteration in the "build_random" bed generation
		"""
		preshuffle_beds = {}
		for annot in annot_table:
			random_str=[iii for iii in ["chrX 1 "+str(length+1)+" "+str(0)]*annot_count[annot]["count"]] 
			b = bedtool("\n".join(random_str), from_string=True)
			preshuffle_beds[annot] = b
		return preshuffle_beds

	def build_random():
		"""
		Generate a random sample with matching annotation distribution to input sample
		"""
		rand_annot_bed_list = []
		for annot in preshuffle_beds:
			a = preshuffle_beds[annot].shuffle(g=ref_genome_file, incl=annot_zone_bed[annot].fn)
			rand_annot_bed_list.append(a)#.fn
			print(str(annot)+" effective size: "+str(annot_count[annot]["count"]))
		print("Full effective size: "+str(annot_count["all"]["count"]*iter_count))
		random_bed=rand_annot_bed_list[0].cat(*rand_annot_bed_list[1:], postmerge=False, stream=True)

		for annot_file in rand_annot_bed_list:
			if os.path.exists(annot_file.fn):
				os.remove(annot_file.fn)
			else:
				print("The file does not exist")
		del rand_annot_bed_list[:]
		random_bed = random_bed.sort(stream=True)
		return random_bed

	def iterate_random(infile):
		"""
		Performs one iteration of creating an annotation matching random sample,
		tabulating its repeat overlap and adding it to a cummulative total.
		Usually parrallelized
		"""
		#print("in iterate"+str(infile))
		it_time=time.time()
		random_bed = build_random()
		intersect_random_bed = random_bed.intersect(infile, wo=True, sorted=True)
		add_trial(intersect_random_bed, trials_te_table)
		print("Iteration Time: "+str(time.time()-it_time) )
		return 1
		
	def givename(feature):
		feature.name = str(i)
		return feature

	def getScriptPath():
		###Returns path to the script
		return os.getcwd()#path.dirname(__file__)


	#=======================================================
	#                      BEGIN
	#=======================================================
	version = 0.5
	#=======================================================
	#                  OS & DIR SETUP
	#=======================================================
	if os.name == 'posix':
		current_dir = '.'
		NL="\n"
	else :
		current_dir = getScriptPath()
		NL="\r\n"

	path = current_dir
	#=======================================================
	#                   Arguments
	#=======================================================
	print(">>> Analyze TE Peaks v."+str(version))
	print(' '.join(os.sys.argv[0:]))
	print()
	parser = argparse.ArgumentParser(description='Optional app description')
	parser.add_argument('input_file',nargs='?', type=str, default="", help='Input bed')
	parser.add_argument('-g',nargs='?', dest='ref_genome', type=int, default=38, help='int number of human Reference genome used')
	parser.add_argument('-l',nargs='?', dest='length', type=int, default=200, help='length to set peaks')
	parser.add_argument('-n',nargs='?', dest='size', type=int, default=-1, help='desired number of peaks in random samples. -1 same as input')
	parser.add_argument('-r',nargs='?', dest='iter_count', type=int, default=10, help='number of random sample iteration')
	parser.add_argument('-o',nargs='?', dest='output_file', type=str, default="analyze_peaks_output.csv", help='output file for generated random set')
	parser.add_argument('-p',nargs='?', dest='processes', type=int, default=1, help='Number of processes to use in parallel')
	args = parser.parse_args()
	peak_file = args.input_file
	genome_id = args.ref_genome
	iter_count = args.iter_count
	length = args.length
	size = args.size
	out_csv = args.output_file
	print("Input file: "+str(peak_file))
	#create folder if needed
	if "/" in out_csv:
		print("Output file: "+str(out_csv))
		out_csv_dir = "/".join(out_csv.split("/")[:-1])
		print("Output dir: "+str(out_csv_dir))
		if not os.path.exists(out_csv_dir):
			os.makedirs(out_csv_dir)

	out_csv_pre_ext = ".".join(out_csv.split(".")[:-1])

	#Annotations (to be file derived)
	annot_table={"tss":{"min":-1000,"max":1000},
	"promoter":{"min":-5000,"max":1000},
	"intragenic":{"min":0,"max":0},
	"proximal":{"min":-10000,"max":10000},
	"distal":{"min":-100000,"max":100000},
	"desert":{"min":-9999999999,"max":9999999999}}
	annot_table_ordered = ["intragenic", "proximal", "distal", "desert"]
	annot_table_tss_ordered = ["tss","promoter"]
	print("-----------------------------------------------")
	#=======================================================
	#            Genome & Associated files
	#=======================================================
	#genome & associated files
	refseq = "lib/refseq"+str(genome_id)+".sorted.bed"
	refseq_tss = "lib/refseq"+str(genome_id)+"_tss.sorted.bed"
	repeat_bed = "lib/hg"+str(genome_id)+"_repeats.sorted.bed"
	ref_genome = "hg"+str(genome_id)
	ref_genome_file = "lib/human."+ref_genome+".genome"
	repfile = "lib/rmsk_hg"+str(genome_id)+".txt"
	rmsk_data_file = "lib/rmsk_data_hg"+str(genome_id)+".cpickle"
	#Load repeat masker data
	if os.path.isfile(rmsk_data_file):
		file_obj=open(path+"/"+rmsk_data_file,"rb")
		rmsk_data, repeat_to_class, repeat_to_superfamily, fam_to_class= pickle.load(file_obj)
		file_obj.close()
		print("Loaded repeat masker data")
	else:
		#repfile = "rmsk_hg38.txt"#"hg19_repeats.sorted.bed"
		rmsk_data, repeat_to_class, repeat_to_superfamily, fam_to_class = create_rmsk_info(repfile)
		file_obj=open(path+"/"+rmsk_data_file ,"wb")
		pickle.dump((rmsk_data, repeat_to_class, repeat_to_superfamily, fam_to_class) , file_obj)
		file_obj.close()
		print( "Created repeat masker data")
	all_repeats = rmsk_data["repeat"].keys()


	#=======================================================
	#                   CODE START
	#=======================================================
	"""
	analyze_peaks.py input_file -g genome -l length -n size -r iteration_count -o output_file -p processes
	"""
	save_main_intersect=True
	start_time = time.time()
	last_time = start_time
	timestamp("\nStart")
	counter = multiprocessing.Value('i', 0)
	main_bed = bed_resize_intervals(peak_file, length)#bedtool(peak_file)#
	timestamp("Loaded and resized sample bed")
	read_count = main_bed.count()
	read_bases = get_lenght_f(main_bed)
	print("Base count: "+str(read_bases))
	#if size is -1 (default) set size equal to input bed
	if size==-1:
		size = main_bed.count()
	print("Peak count: "+str(size))
	main_bed = main_bed.sort()
	timestamp("Sorted sample")
	repeat_bed = bedtool(repeat_bed)
	timestamp("Loaded Repeat bed file")
	intersect_main_bed = main_bed.intersect(repeat_bed, wo=True, sorted=True)#, stream=True
	timestamp("Intersect sample and repeats ")

	sample_te_table = {}
	for te in all_repeats:
		sample_te_table[te] = {"count":0, "count_sd":0, "bases":0, "bases_sd":0, "more_in_sample":0}
	main_cleaned_intersect = merge_keeplongest(intersect_main_bed.fn)
	if save_main_intersect:
		main_cleaned_intersect.saveas(out_csv_pre_ext+"_cleaned_intersection.bed")
	main_te_table = tableit(main_cleaned_intersect)
	add_to_joint_table(main_te_table, sample_te_table)
	print("++++++++++++++++++++++++++++++++++++++++++++")
	#Get peak distribution within the annotations. Also output result in file
	annot_bed, annot_zone_bed, annot_count = get_annotation_distributions(main_bed)
	timestamp("Done peak annotation distributions")
	preshuffle_beds = build_preshuffle()
	timestamp("Preshuffled beds generated")
	trials_te_table = multiprocessing.Manager().dict()
	#initiate the trial joint table
	counter = multiprocessing.Value('i', 0)
	for te in all_repeats:
		trials_te_table[te] = {"count":0, "count_sd":0, "bases":0, "bases_sd":0, "more_in_sample":0}
	print("Process: "+str(args.processes))
	print("++++++++++++++++++++++++++++++++++++++++++++")
	try:
		paradata = list(repeat_bed.parallel_apply(iter_count, iterate_random, (repeat_bed, ), {}, processes=args.processes))
		#paradata = list(bedtools.parallel_apply(repeat_bed, method = iterate_random, method_args=(repeat_bed, ), {},iterations = iter_count, processes=args.processes))
	except KeyboardInterrupt:
		print("You cancelled the program!")
		os.sys.exit(1)
	timestamp("Done "+str(iter_count)+" trials")
	print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	#Build result table from the sample TE table and the Random trials TE table
	results = differential_results(sample_te_table, trials_te_table)#tabulate_results(sample_trial, random_trials)
	timestamp("Done Tabulating ")
	results.to_csv(out_csv)

	print("Finished.")
	print("total time: "+str(time.time()-start_time) )
	#print("Total Mem used: "+str(memory_usage_psutil()))
