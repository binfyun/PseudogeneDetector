
from __future__ import print_function
import time
from datetime import datetime
import pysam
import sys
import os
import math
import numpy as np
import subprocess
#from termcolor import colored, cprint
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
#import click

year = sys.argv[1]
sampleList = "/dmp/hot/zehira/All_Results/MASTER_SAMPLE_LIST.txt"
ref = "/dmp/data/pubdata/hg-fasta/production/Homo_sapiens_assembly19.fasta"
impact_gene = "/dmp/data/mskdata/interval-lists/production/genelist.with_aa.interval_list"
in_file = open("%s_list_of_pseudoGenesBams_new.txt"%year, "r")
out_file = open("%s_candidate_pseudo.txt"%year, "w")

print("Run""\t""MRN""\t""Gene""\t""N_fc""\t""T_fc""\t""Tumor""\t""Region""\t""SplitReads""\t""ProperReads""\t""AF", file=out_file)


def estimateInsertSizeDistribution(prefix, bamfile, pos, chr, gene):#, alignments=50000):
	"estimate insert size from first alignments in bam file. returns mean and stddev of insert sizes."
	#assert isPaired(bamfile, start, end, tcr), 'can only estimate insert size from paired bam files'
	#mappedFile = pysam.AlignmentFile(bamfile,"rb")
	#end = int(end)
	#start = int(start)
	#tcr = str(tcr)
	pseudo = []
	e_i_junc = {}
	pos1 = pos

	for i, j in zip(pos1, range(len(pos1))):
		i = int(i)
		if (j % 2) == 0: 
			region = '%s:%s-%s'%(chr, i-1,i)
		else:
			region = '%s:%s-%s'%(chr, i,i+1)
		pseudo.append(region)
		samfile = bamfile.fetch(region=region)
		# only get positive to avoid double counting
		#inserts = np.array([read.template_length for read in samfile if read.is_proper_pair and read.template_length > 0])
		#inserts = np.array([read.template_length for read in samfile if read.template_length < 1000])
		split_inserts = np.array([])
		inserts = np.array([])
		all_reads = np.array([])
		for read in samfile:
			if not read.is_duplicate and read.template_length != 0:
				all_reads = np.append(all_reads, abs(read.template_length))
				if abs(read.template_length) >= 1000:
					split_inserts = np.append(split_inserts, read.template_length)
				else:
					inserts = np.append(inserts, read.template_length)
		#split_inserts, inserts = np.array([read.template_length for read in samfile if read.template_length >= 1000 else])
		#print colored("Split reads:, %s, %s", "red")%(len(split_inserts), abs(np.mean(abs(split_inserts))))
		#print colored("Proper reads:, %s, %s", "green")%(len(inserts), abs(np.mean(abs(inserts))))
		#print>>out_file, prefix, (region, j), ("Split reads: %s %s")%(len(split_inserts), abs(np.mean(abs(split_inserts)))), ("Proper reads: %s %s")%(len(inserts), abs(np.mean(abs(inserts))))
		if float((len(inserts)+len(split_inserts))) == 0:
			print(prefix, (region, j), "\t", (len(split_inserts)), "\t", (len(inserts)), "\t", 0, file=out_file)
		else:
			print(prefix, (region, j), "\t", (len(split_inserts)), "\t", (len(inserts)), "\t", round(float(len(split_inserts))/float((len(inserts)+len(split_inserts))),3), file=out_file)

############# Get some stats and histogram => plot script needs reworks #######################
		# plt.hist(all_reads, alpha=0.5, histtype='stepfilled', bins="auto")
		# #plt.xlim([0, 8000])
		# plt.savefig("%s_region_readDist.png", dpi=300)%gene
		# #plt.tight_layout()
		# plt.close('all')
		# #ax = sns.distplot(all_reads, fit=norm, kde=False) # - gaussian distribution
		# ax = sns.distplot(all_reads) # - kde
		# #plt.xlim([0, 8000])
		# fig = ax.get_figure()
		# fig.savefig("%s_kde.png", dpi=300)%gene
		# plt.close('all')

def prepareINPUTS(smplist, infile):
	cancer = open(sampleList, "r")
	sample_cancer = {}
	for sample in cancer.readlines():
		samples = sample.strip().split("\t")
		run = samples[2]
		tumor = samples[10]
		sample_cancer[run]=tumor
	cancer.close()

	gene_inter = {}
	with open("/dmp/data/mskdata/interval-lists/production/genelist.with_aa.interval_list", "r") as intervals:
		for lines in intervals.readlines():
			line = lines.strip().split("\t")
			Gene = line[4].split(":")[0]
			Chrom = line[0]
			upos = line[1]
			dpos = line[2]
			gene_inter.setdefault(Gene, []).append([upos,dpos])

	debug = open("debg.txt", "w")
	for items in in_file.readlines():
		item = items.strip().split(" ")
		bam = item[0]
		gene = item[1].strip()
		chrom = item[2]
		n_fc = item[3]
		t_fc = item[4]
		run_n = bam.split("/")[6]
		sample = item[0].split("/")[8].split("_")[0].split("-")[0]
		print(sample, file=debug)
		bamfile = pysam.AlignmentFile(bam, "rb")
		if sample in sample_cancer:
			print(sample, sample_cancer[sample], file=debug)
			if gene == "MLL3":
				gene = "KMT2C"
			prefix = run_n + "\t" + sample + "\t" + gene + "\t" + n_fc + "\t" + t_fc + "\t" + sample_cancer[sample] + "\t"
			if gene in gene_inter:
				l = gene_inter[gene]
				pos = [item for sublist in l for item in sublist]
				print(gene, pos, file=debug)
				estimateInsertSizeDistribution(prefix, bamfile, pos, chrom, gene)
		else:
			continue
	bamfile.close()
	debug.close()


def main():
	prepareINPUTS(sampleList, in_file)
	in_file.close()
	out_file.close()


if __name__ == '__main__':
    main()