#!/usr/bin/python
# coding: utf8

# Imports:
import csv
import numpy as np
from scipy.stats import pearsonr

##################################################################
# Script to - pearson - correlate expression profiles of two genotypes along two years
# Data will be saved (csv)
# from 4 sourcefiles (2 years, 2 genotypes; avg. triplicates)
# By Jens Theine (jens.theine at uni-bielefeld .de), 30.10.2019
##################################################################

gois = []

############################# CONFIG #############################

# IDs + TPMs:
FB17_TPMs = "/path/to/counttable_FB17.csv"
SB17_TPMs = "/path/to/counttable_SB17.csv"
FB14_TPMs = "/path/to/counttable_FB14.csv"
SB14_TPMs = "/path/to/counttable_SB14.csv"

# Output file:
outfile = "/path/to/output/file.csv"

# Gene(s) of interest: (not obligatory)
#goi = ['VIT_207s0001g00300','VIT_208s0001g00300']


##################################################################

# Get all geneids of all genes from first counttable (if gois not stated):
if len(gois) == 0:
	print "Holen aller geneids..."
	with open( FB17_TPMs, "r" ) as temp_in:
		for line in temp_in:
			if "VIT" in line:
				geneid = line.split(";")[0]
				gois.append(geneid[0:18])

# Prepare csv
fields=['geneid','SB pearson-correlation (14 vs. 17)','SB pvalue (14 vs. 17)','FB pearson-correlation (14 vs. 17)','FB pvalue (14 vs. 17)']
with open('/vol/gf-grape/project_NoViSys/members/jenthein/_RV-SBFB/__Timeshift_Test_Phasenangleichung+Corr/out/corr_FBSB1417.csv', 'a') as f:
	writer = csv.writer(f)
	writer.writerow(fields)

counter = 0

### Iterate over all geneids::
for i in gois:
	goi = str(i)
	counter = counter + 1
	#print goi
	
	SB14 = []
	FB14 = []
	SB17 = []
	FB17 = []
	
	# Get first TPMs from counttable:
	with open( SB14_TPMs, "r" ) as SB14_TPMs_tmp:
		for line in SB14_TPMs_tmp:
			geneid = line.split(";")[0]
			if goi == geneid:
				line = line.replace(",",".")
				SB14.append(float(line.split(";")[1]))
				SB14.append(float(line.split(";")[2]))
				#SB14.append(float(line.split(";")[3]))
				SB14.append(float(line.split(";")[4]))
				SB14.append(float(line.split(";")[5]))
				SB14.append(float(line.split(";")[6]))
				SB14.append(float(line.split(";")[7]))
				#SB14.append(float(line.split(";")[8])) 
				SB14.append(float(line.split(";")[9]))
				SB14.append(float(line.split(";")[10]))
				#SB14.append(float(line.split(";")[11]))
				#SB14.append(float(line.split(";")[12]))
			
	# Get second TPMs from counttable:
	with open( FB14_TPMs, "r" ) as FB14_TPMs_tmp:
		for line in FB14_TPMs_tmp:
			geneid = line.split(";")[0]
			if goi == geneid:
				line = line.replace(",",".")
				FB14.append(float(line.split(";")[1]))
				FB14.append(float(line.split(";")[2]))
				#FB14.append(float(line.split(";")[3]))
				FB14.append(float(line.split(";")[4]))
				FB14.append(float(line.split(";")[5]))
				FB14.append(float(line.split(";")[6]))
				FB14.append(float(line.split(";")[7]))
				#FB14.append(float(line.split(";")[8]))
				FB14.append(float(line.split(";")[9]))
				FB14.append(float(line.split(";")[10]))
				#FB14.append(float(line.split(";")[11]))
				#FB14.append(float(line.split(";")[12]))
	
	#  Get third TPMs from counttable:
	with open( SB17_TPMs, "r" ) as SB17_TPMs_tmp:
		for line in SB17_TPMs_tmp:
			geneid = line.split(";")[0]
			if goi == geneid:
				line = line.replace(",",".")
				SB17.append(float(line.split(";")[1]))
				#SB17.append(float(line.split(";")[2]))
				#SB17.append(float(line.split(";")[3]))
				SB17.append(float(line.split(";")[4]))
				#SB17.append(float(line.split(";")[5]))
				SB17.append(float(line.split(";")[6]))
				SB17.append(float(line.split(";")[7]))
				SB17.append(float(line.split(";")[8]))
				SB17.append(float(line.split(";")[9]))
				#SB17.append(float(line.split(";")[10]))
				SB17.append(float(line.split(";")[11]))
				#SB17.append(float(line.split(";")[12]))
				SB17.append(float(line.split(";")[13]))

	#  Get fourths TPMs from counttable:
	with open( FB17_TPMs, "r" ) as FB17_TPMs_tmp:
		for line in FB17_TPMs_tmp:
			geneid = line.split(";")[0]
			if goi == geneid:
				line = line.replace(",",".")
				FB17.append(float(line.split(";")[1]))
				#FB17.append(float(line.split(";")[2]))
				#FB17.append(float(line.split(";")[3]))
				FB17.append(float(line.split(";")[4]))
				#FB17.append(float(line.split(";")[5]))
				FB17.append(float(line.split(";")[6]))
				FB17.append(float(line.split(";")[7]))
				FB17.append(float(line.split(";")[8]))
				FB17.append(float(line.split(";")[9]))
				#FB17.append(float(line.split(";")[10]))
				FB17.append(float(line.split(";")[11]))
				#FB17.append(float(line.split(";")[12]))
				FB17.append(float(line.split(";")[13]))

	if len(SB14) != 0:

		# Calculate Pearsons r and p-value for each geneid between years:
		SBcorr, SBpvalue = pearsonr(SB14, SB17)
		FBcorr, FBpvalue = pearsonr(FB14, FB17)
		print str(counter) + "/" + str(len(gois))
	
		# Save results in csv:
		fields=[goi,SBcorr,SBpvalue,FBcorr,FBpvalue]
		with open(outfile, 'a') as f:
			writer = csv.writer(f)
			writer.writerow(fields)
			print fields
	
print "..done"
	
