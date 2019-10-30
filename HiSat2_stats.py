#!/usr/bin/python
# coding: utf8

# Imports:
import re, os, sys, time, glob
import datetime as datetime

##################################################################
# Script to summarise HiSat2 reports
# Data will be saved
# By Jens Theine (jens.theine at uni-bielefeld .de), 30.10.2019
##################################################################


mapper_out_path = "/path/to/mapping/hisat2_reports/"
mapper_out = "/path/to/output_file"

mapper_file_names = []
proc = 0

mapper_reads = []
unpaired = []
aligned0 = []
aligned1 = []
alignedover1 = []
overall_alignment_rate = []
tmp = ""

mapper_file_names = glob.glob(mapper_out_path + "*-out_log.txt") + glob.glob(mapper_out_path + "*/*-out_log.txt")

for name in mapper_file_names:
	with open(name, "r") as mfn:
		for line in mfn:
			tmp = ""
			if "reads; of these:" in line:
				reads = line
				reads = int(reads.replace(" reads; of these:\n",""))
				mapper_reads.append(reads)
			if "were unpaired; of these:" in line:
				start = "'  "
				end = " ("
				tmp = line[line.find(start)+len(start):line.rfind(end)]
				unpaired.append(tmp)
			if "aligned 0 times" in line:
				start = "    "
				end = " ("
				tmp = line[line.find(start)+len(start):line.rfind(end)]
				aligned0.append(tmp)
			if "aligned exactly 1 time" in line:
				start = "    "
				end = " ("
				tmp = line[line.find(start)+len(start):line.rfind(end)]
				aligned1.append(tmp)
			if "aligned >1 times" in line:
				start = "    "
				end = " ("
				tmp = line[line.find(start)+len(start):line.rfind(end)]
				alignedover1.append(tmp)
			proc = proc + 1
			if "overall alignment rate" in line:
				tmp = line[0:5]
				overall_alignment_rate.append(tmp)

i = 0
o = 0

for i in mapper_reads:
	o = int(o) + int(i)
all_reads = o
res_mapper_reads = str(o / len(mapper_reads))
i = 0
o = 0


for i in unpaired:
	o = float(o) + float(i)
res_unpaired = str(o / len(unpaired))
i = 0
o = 0


for i in aligned0:
	o = float(o) + float(i)
res_aligned0 = str(o / len(aligned0))	
i = 0
o = 0


for i in aligned1:
	o = float(o) + float(i)
res_aligned1 = str(o / len(aligned1))
i = 0
o = 0


for i in alignedover1:
	o = float(o) + float(i)
res_alignedover1 = str(o / len(alignedover1))
i = 0
o = 0


for i in overall_alignment_rate:
	o = float(o) + float(i)
res_overall_alignment_rate = str(o / len(overall_alignment_rate))


with open(mapper_out, "w") as out:
	out.write("mean reads" + "\t" + "all reads" + "\t" + "unpaired reads" + "\t" + "aligned 0 times" + "\t" + "aligned exactly 1 time" + "\t" + "aligned >1" + "\t" + "overall alignment rate" + "\n")
	out.write(str(res_mapper_reads) + "\t" + str(all_reads) + "\t" + str(res_unpaired) + "\t" + str(res_aligned0) + "\t" + str(res_aligned1) + "\t" + str(res_alignedover1) + "\t" + str(res_overall_alignment_rate) + "\n")

print "Number of detected FeatureCount summaries (lines processed / files processed): (" + str(proc) + "/" + str( len( mapper_file_names ) ) + ")"
