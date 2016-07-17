
#############################
## Import required modules ##
#############################
""" 
- argparse, collections, csv, re, and sys are part of the standard Python library
- biopython, numpy, sympy, and cyvcf are easily installed via conda and pip
"""

import argparse
import collections
import csv
import re
import sys
from Bio import SeqIO
import numpy as np
import sympy
import cyvcf


########################################
## Import files from the command line ##
########################################
parser = argparse.ArgumentParser(description="This program calculates nucleotide"\
								" diversity separately for autosomes and sex-linked"\
								" chromosomes.  It requires, as input, a minimum of:"\
								" a VCF of polymorphic sites, a list of sex-linked contigs, and"\
								" bed files containing callable regions for each population.  It is assumed"\
								" that all filtering (both VCF and bed file) has been done elsewhere.")

#Print help/usage if no arguments are supplied
if len(sys.argv)==1:
	parser.print_usage()
	sys.exit(1)
	
# A function to check if bootstrap input is non-negative 
def check_negative(boot_int):
    value = int(boot_int)
    if value < 0:
         raise argparse.ArgumentTypeError("%s is an invalid value.  Input for --bootstrap must be nonnegative integer" % boot_int)
    return value
	
# Parse the command line
parser.add_argument("--vcf", required=True, 
					help="REQUIRED. Input VCF file.  Can be gzipped.")
parser.add_argument("--callable_regions", nargs="*", required=True, 
					help="REQUIRED. Bed files (no header) containing callable regions for each population (to accurately calculate pi). Order of files must correspond exactly to the order of the population lists")
parser.add_argument("--bootstrap", default=0, type=check_negative, 
					help="Default is 0. If n > 0, returns 95% confidence interval of distribution of n replicates with replacement")
parser.add_argument("--male_list", nargs="*", default=None, 
					help="Default is None (i.e., all females).  Provide file listing males. List males using exact names (capitalization, etc.) as found in the vcf file.  Allows script to properly handle X chromosomes")
parser.add_argument("--population_lists", nargs="*", default=None,
					help="Default is None. If None, script will treat all individuals as coming from the same population.  Else, will read files and treat all individuals in each file as coming from the same population.")
parser.add_argument("--x_linked_scaffolds", nargs="*", default=None,
					help="Default is None. Provide a file listing (either in a single column, or horizontally in a single row separated by tabs). If None, script will treat all scaffolds as autosomal.  Else, will treat all scaffolds in file as X-linked.")
parser.add_argument("--scaffold_sites_filter", type=int, default=0,
					help="Default is 0.  Minimum number of callable sites for a scaffold/chromosome to be included in analyses.  Scaffold is removed if this minimum is not met in any individual.")
parser.add_argument("--min_cov", type=int, default=0,
					help="Default is 0.  Minimum read depth for a sample at a site for that genotype to be included.")
parser.add_argument("--QD", type=float, default=0,
					help="Default is 0. Minimum quality by depth (QD) value.")
parser.add_argument("--FS", type=float, default=10000.0,
					help="Default is 10000. Maximum Fisher strand (FS) score.")
parser.add_argument("--QUAL", type=float, default=0.0,
					help="Default is 0.  Minimum site quality score.")
parser.add_argument("--MAPQ", type=float, default=0.0,
					help="Default is 0. Minimum mean mapping quality score.")
args = parser.parse_args()

if len(args.population_lists) != len(args.callable_regions):
	print "Order and number of --callable_region files must match --population_lists files exactly"
	sys.exit(1)
	
else:
	print "Found %d populations (based on --callable_region and --population_lists files)" % len(args.population_lists)
	print "Matched the following files (if incorrect, be sure --callable_region and --population_lists files are ordered in the same way - i.e., --callable_regions pop1_callable.bed pop2_callable.bed --population_lists pop1_ids.txt pop2_ids.txt):"
	for idx,val in enumerate(args.callable_regions):
		print "%s     %s" % (args.callable_regions[idx], args.population_lists[idx])

print "Reading input files..."
print ""

# Processes X-linked scaffolds
with open(args.x_linked_scaffolds[0],"r") as f:
	x_linked = [item for sublist in list(csv.reader(f,delimiter="\t")) for item in sublist]
	x_linked = [x.strip() for x in x_linked]
	while "" in x_linked:
		x_linked.remove("")

# Process input list of males
# Includes a number of cleaning and filtering steps to clean up accidental whitespace, and
#		read from both vertical and horizonal lists (horizontal have to be tab separated)
# If no input lists are provided, it treats all samples in the VCF as female
if args.male_list != None:
	with open(args.male_list[0],"r") as f:
		males = [item for sublist in list(csv.reader(f,delimiter="\t")) for item in sublist]
		males = [x.strip() for x in males]
		while "" in males:
			males.remove("")		
else:
	males = []

# Open vcf file and initiate the CyVCF parser
vcf_file = open(args.vcf,"r")
vcf_reader = cyvcf.Reader(vcf_file)

# Process input population lists
# Includes a number of cleaning and filtering steps to clean up accidental whitespace, and
#		read from both vertical and horizonal lists (horizontal have to be tab separated)
# If no input lists are provided, it treats all samples in the VCF as coming from the same
#		population

if args.population_lists != None:
	populations = []
	no_populations = False
	for i in args.population_lists:
		with open(i,"r") as f:
			temp_pop = [item for sublist in list(csv.reader(f,delimiter="\t")) for item in sublist]
			temp_pop = [x.strip() for x in temp_pop]
			while "" in temp_pop:
				temp_pop.remove("")
		populations.append([temp_pop,i,[[],[]],[[],[]]])
else:
	populations = [[vcf_reader.samples,args.vcf,[[],[]],[[],[]]]]
	


#########################################################
## Calculate total callable sequence for each DNA type ##
#########################################################

print ""
print "Calculating total callable sequence and filtering low-coverage scaffolds..."
print ""

auto_seq = 0
x_seq = 0

auto_contigs = {}
x_contigs = {}

# Process input bed files of callable regions

for idx,i in enumerate(args.callable_regions):
	callable_dict = {}
	with open(i,"r") as f:
		for j in csv.reader(f,delimiter="\t"):
			if j[0] in callable_dict:
				callable_dict[j[0]] += int(j[2]) - int(j[1])
			else:
				callable_dict[j[0]] = int(j[2]) - int(j[1])
	initial_scaffolds = len(callable_dict)
	filtered_dict = {k:v for (k,v) in callable_dict.iteritems() if v >= args.scaffold_sites_filter}
	passing_scaffolds = len(filtered_dict)
	populations[idx].append(filtered_dict)
	
	print "For population corresponding to %s - %s:" % (args.callable_regions[idx], args.population_lists[idx])
	print "Initial scaffolds: %d" % initial_scaffolds
	print "Removed scaffolds (not enough sequence coverage): %d" % (initial_scaffolds - passing_scaffolds)
	print "Remaining scaffolds: %d" % passing_scaffolds
	print ""
					

########################################################################
## Functions to calculated diversity and conduct bootstrap resampling ##
########################################################################

def pi_overall(tot_diff, k, sequence_length):
	""" Calculates mean nucleotide diversity, pi, using equation from Box 1.3 in 
	Charlesworth and Charleswoth (2010):
	(tot_diff / (k choose 2)) / sequence_length
	where:
		tot_diff = the total number of differences observed in all pairwise comparisons
		k = the number of chromosomes sampled (e.g. k = 6 for 3 diploid individuals)
		sequence_length = the number of bases analyzed per sequence (should be same for
										all sequences)
	"""
	if k == 0:
		return 0
	elif k == 1:
		return 0
	else:
		numerator = float(tot_diff) / ((float(k) * (float(k) - 1)) / 2.0)
		return numerator / float(sequence_length)
	

def pi_site(allele_count_list):
	"""Function calculates pi from Hohenlohe et al. (2010)
	
	pi = 1 - sum((n_i choose 2) / (n choose 2))
		where:
			n_i is count of allele i in sample
			n is the sum of n_i (allele1 plus allele2 in this case, as we assume bi-allelic sites
			
	inputs:
		allele_count_list is a list of the counts of the different alleles present at a site
	assumptions:
		snps
		sympy installed and imported for binomial coefficient calculation 
	"""
	n_tot = 0
	for i in allele_count_list:
		n_tot += i	
	if n_tot == 0:
		return 0
	elif n_tot == 1:
		return 0
	else:
		pi_denom = float(sympy.binomial(n_tot,2))
		pi_numer = 0.0
		for i in allele_count_list:
			pi_numer += float(sympy.binomial(i,2))
		return (1.0 - (pi_numer / pi_denom))
		
def count_diffs(allele_list):
	""" Takes a list or string of alleles to count the number of differences among chromosomes.
	Example: For an input site from 3 diploid individuals with genotypes (A/T), (T/T), and (C/A),
	appropriate inputs would be either ["A","T","T","T","C","A"] or "ATTTCA".
	
	Returns:
		The number of pairwise differences as an integer
	"""
	diffs = 0
	for index,i in enumerate(allele_list):
		for j in allele_list[index + 1:]:
			if i != j:
				diffs += 1
	return diffs
	
def rand_sample(data_input, n_sites):
	"""Outputs a numpy array of resampled (with replacement) values for bootstrap function
	Inputs:
		data_input is the list of values to resample
		n_vals is the number of values to be resampled
	Requirements:
		numpy installed and imported
		
	Note: IF THE NUMBER OF SITES TO BE RESAMPLED IS LONGER THAN THE DATA_INPUT, ALL
		ADDITIONAL SITES ARE ASSUMED TO HAVE A VALUE OF 0
	"""
	dim = len(data_input)
	# create random indices to grab the random sample
	indices = np.random.random_integers(0, n_sites, n_sites)
	# return all random values with indices less than the length of the input table
	#		all other sites are assumed to be zero
	return np.take(data_input, indices[ indices < dim ])
	

def bootstrap_pi_distribution(data_input, n_sites, replicates):
	""" Returns a bootstrapped distribution (with replacement) as a list.
	data_input is the list of values
	n_sites is the total number of callable site (will be longer than data_input if
		the input VCF only contained polymorphic sites)
	replicates is the number of bootstrap replicates
	"""
	resamples = []
	n_sites = float(n_sites)
	for i in range(replicates):
		resamples.append(np.sum(rand_sample(data_input,n_sites)) / n_sites)
	return resamples
	


###########################################################################
## Parse VCF file, calculate pi per site, and count pairwise differences ##
###########################################################################


print "Beginning diversity calculations"
counter = 0
for record in vcf_reader:
	try:
		if float(record.INFO["QD"]) >= args.QD and float(record.INFO["FS"]) <= args.FS and float(record.QUAL) >= args.QUAL and float(record.INFO["MQ"]) >= args.MAPQ:
			if record.CHROM in x_linked:
				for pop in populations:
					allele_list=[]
					male_count = 0
					total_count = 0
					for indv in pop[0]:
						if indv in males:
							male_count += 1
						total_count += 1
						call = record.genotype(indv)
						if call['GT'] != None:
							if indv in males:
								try:
									if float(call["DP"]) >= args.min_cov:
										# Randomly selects an allele to include for males
										allele_list.append(call.gt_bases[int(np.random.choice([0,2],1))])
								except TypeError:
									print "Error reading DP - site excluded"
							else:
								try:
									if float(call["DP"]) >= args.min_cov:
										# call.gt_bases returns in the format "A/T", so this grabs the A and
										#		the T, while skipping the / (or |)
										allele_list.append(call.gt_bases[0])
										allele_list.append(call.gt_bases[2])
								except TypeError:
									print "Error reading DP - site excluded"
					# Process allele list and calculate pi and number of differences
					if len(allele_list) == male_count + (2 * (total_count - male_count)):
						pop[3][0].append(pi_overall(count_diffs(allele_list), len(allele_list), 1.0))
						allele_count = collections.Counter(allele_list)
						pop[3][1].append(pi_site([allele_count[x] for x in allele_count]))
			else:
				for pop in populations:
					allele_list = []
					total_count = 0
					for indv in pop[0]:
						call = record.genotype(indv)
						total_count += 1
						if call['GT'] != None:
							try:
								if float(call["DP"]) >= args.min_cov:
									allele_list.append(call.gt_bases[0])
									allele_list.append(call.gt_bases[2])
							except TypeError:
								print "Error reading DP - site excluded"
					if len(allele_list) == total_count * 2:
						pop[2][0].append(pi_overall(count_diffs(allele_list), len(allele_list), 1.0))
						allele_count = collections.Counter(allele_list)
						pop[2][1].append(pi_site([allele_count[x] for x in allele_count]))
	except KeyError:
		print "KeyError - site removed for poorly formated INFO"
		print record
	counter += 1
	if counter % 10000 == 0:
		print "%d records complete..." % (counter)

print "VCF traversal complete"
				

# Close the vcf file
vcf_file.close()

				

###########################################################
## Calculate mean diversity (and bootstrap, if selected) ##
###########################################################

print "Calculating diversity statistics..."
print ""
if args.bootstrap > 0:
	print "Preparing to calculate %d bootstrap replicates with replacement" % (args.bootstrap)
else:
	print "No bootstrap resampling"

#### Charlesworth Pi
print ""
print "Now computing mean nucleotide diversity using equation from Charlesworth and Charlesworth (2010)"
print ""

for pop in populations:
	print ""
	print "%s: " % pop[1]
	auto_seq = 0
	x_seq = 0
	for i in pop[-1]:
		if i in x_linked:
			x_seq += pop[-1][i]
		else:
			auto_seq += pop[-1][i]
	print "Total callable sequence:"
	print "    Autosomes: %d" % auto_seq
	print "    X-linked: %d" % x_seq
	
	print "Mean autosomal diversity: %f" % (float(np.sum(pop[2][0])) / float(auto_seq))
	if args.bootstrap > 0:
		auto_dist = bootstrap_pi_distribution(pop[2][0], auto_seq, args.bootstrap)
		print "     Bootstrap confidence intervals:"
		print "          2.5 = %f" % np.percentile(auto_dist, 2.5)
		print "          Median = %f" % np.percentile(auto_dist,50)
		print "          97.5 = %f" % np.percentile(auto_dist,97.5)
	print "Mean X chromosome diversity: %f" % (float(np.sum(pop[3][0])) / float(x_seq))
	if args.bootstrap > 0:
		x_dist = bootstrap_pi_distribution(pop[3][0], x_seq, args.bootstrap)
		print "     Bootstrap confidence intervals:"
		print "          2.5 = %f" % np.percentile(x_dist, 2.5)
		print "          Median = %f" % np.percentile(x_dist,50)
		print "          97.5 = %f" % np.percentile(x_dist,97.5)
	if args.bootstrap > 0:
		x_a_dist = np.asarray(x_dist) / np.asarray(auto_dist)
		print "Mean X/A diversity ratio: %f" % np.mean(x_a_dist)
		print "     Bootstrap confidence intervals:"
		print "          2.5 = %f" % np.percentile(x_a_dist, 2.5)
		print "          Median = %f" % np.percentile(x_a_dist,50)
		print "          97.5 = %f" % np.percentile(x_a_dist,97.5)

	
#### Hohenlohe Pi
print ""
print "Now computing mean nucleotide diversity using equation from Hohenlohe et al. (2010)"
print ""

for pop in populations:
	print ""
	print "%s: " % pop[1]
	auto_seq = 0
	x_seq = 0
	for i in pop[-1]:
		if i in x_linked:
			x_seq += pop[-1][i]
		else:
			auto_seq += pop[-1][i]
	print "Total callable sequence:"
	print "    Autosomes: %d" % auto_seq
	print "    X-linked: %d" % x_seq
	
	print "Mean autosomal diversity: %f" % (float(np.sum(pop[2][1])) / float(auto_seq))
	if args.bootstrap > 0:
		auto_dist = bootstrap_pi_distribution(pop[2][1], auto_seq, args.bootstrap)
		print "     Bootstrap confidence intervals:"
		print "          2.5 = %f" % np.percentile(auto_dist, 2.5)
		print "          Median = %f" % np.percentile(auto_dist,50)
		print "          97.5 = %f" % np.percentile(auto_dist,97.5)
	print "Mean X chromosome diversity: %f" % (float(np.sum(pop[3][1])) / float(x_seq))
	if args.bootstrap > 0:
		x_dist = bootstrap_pi_distribution(pop[3][1], x_seq, args.bootstrap)
		print "     Bootstrap confidence intervals:"
		print "          2.5 = %f" % np.percentile(x_dist, 2.5)
		print "          Median = %f" % np.percentile(x_dist,50)
		print "          97.5 = %f" % np.percentile(x_dist,97.5)
	if args.bootstrap > 0:
		x_a_dist = np.asarray(x_dist) / np.asarray(auto_dist)
		print "Mean X/A diversity ratio: %f" % np.mean(x_a_dist)
		print "     Bootstrap confidence intervals:"
		print "          2.5 = %f" % np.percentile(x_a_dist, 2.5)
		print "          Median = %f" % np.percentile(x_a_dist,50)
		print "          97.5 = %f" % np.percentile(x_a_dist,97.5)











	

