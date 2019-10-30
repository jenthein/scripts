#!/usr/bin/python
# coding: utf8

# Imports:
import re, os, sys, time, glob
import datetime as datetime

##################################################################
# Script to summarise Featurecounts reports
# Data will be saved
# By Jens Theine (jens.theine at uni-bielefeld .de), 30.10.2019
##################################################################


counter_out_path = "/path/to/mapping/featurecounts_reports/"
fc_out = "/path/to/output_file"

fc_assigned_list = []
fc_unassigned_unmapped_list = []
fc_unassigned_multimapping_list = []
fc_unassigned_nofeatures_list = []
fc_unassigned_ambiguity_list = []
fc_file_names = []
proc = 0

fc_file_names = glob.glob(counter_out_path + "*cat.sorted.bam.gz.countTable.txt.summary") + glob.glob(counter_out_path + "*/*cat.sorted.bam.gz.countTable.txt.summary")

for name in fc_file_names:
	with open(name, "r") as fc:
		for line in fc:
			if line.split("\t")[0] == "Assigned":
				assigned = line.split("\t")[1]
				assigned = assigned.replace("\n","")
				fc_assigned_list.append(assigned)
			if line.split("\t")[0] == "Unassigned_Unmapped":
				unassigned_unmapped = line.split("\t")[1]
				unassigned_unmapped = unassigned_unmapped.replace("\n","")
				fc_unassigned_unmapped_list.append(unassigned_unmapped)
			if line.split("\t")[0] == "Unassigned_MultiMapping":
				unassigned_multimapping = line.split("\t")[1]
				unassigned_multimapping = unassigned_multimapping.replace("\n","")
				fc_unassigned_multimapping_list.append(unassigned_multimapping)
			if line.split("\t")[0] == "Unassigned_NoFeatures":
				unassigned_nofeatures = line.split("\t")[1]
				unassigned_nofeatures = unassigned_nofeatures.replace("\n","")
				fc_unassigned_nofeatures_list.append(unassigned_nofeatures)
			if line.split("\t")[0] == "Unassigned_Ambiguity":
				unassigned_ambiguity = line.split("\t")[1]
				unassigned_ambiguity = unassigned_ambiguity.replace("\n","")
				fc_unassigned_ambiguity_list.append(unassigned_ambiguity)
			proc = proc + 1

i = 0
o = 0

for i in fc_assigned_list:
	o = int(o) + int(i)
	res_fc_assigned = str(o / len(fc_assigned_list))
i = 0
o = 0

for i in fc_unassigned_unmapped_list:
	o = int(o) + int(i)
	res_fc_unassigned_unmapped = str(int(o) / len(fc_unassigned_unmapped_list))
i = 0
o = 0

for i in fc_unassigned_multimapping_list:
	o = int(o) + int(i)
	res_fc_unassigned_multimapping = str(int(o) / len(fc_unassigned_multimapping_list))
i = 0
o = 0

for i in fc_unassigned_nofeatures_list:
	o = int(o) + int(i)
	res_fc_unassigned_nofeatures = str(int(o) / len(fc_unassigned_nofeatures_list))
i = 0
o = 0

for i in fc_unassigned_ambiguity_list:
	o = int(o) + int(i)
	res_fc_unassigned_ambiguity = str(int(o) / len(fc_unassigned_ambiguity_list))


with open(fc_out, "w") as out:
	out.write("assigned" + "\t" + "unassigned_unmapped" + "\t" + "unassigned_multimapping" + "\t" + "unassigned_nofeatures" + "\t" + "unassigned_ambiguity" + "\n")
	out.write(res_fc_assigned + "\t" + res_fc_unassigned_unmapped + "\t" + res_fc_unassigned_multimapping + "\t" + res_fc_unassigned_nofeatures + "\t" + res_fc_unassigned_ambiguity)

print "Number of detected FeatureCount summaries (lines processed / files processed): (" + str(proc) + "/" + str( len( fc_file_names ) ) + ")"
