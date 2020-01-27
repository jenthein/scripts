######################################################################################
# Script to estimate the temporal x-shift between two similar expression profiles
# Data will be saved (plot + csv)
# from a counttable/datamatrix (2 genotypes (G1/G2); triplicates will be averaged)
# By Jens Theine (jens.theine_at_uni-bielefeld.de), 27.01.2020
######################################################################################

# Imports:
import csv
import time
import xlrd
import scipy.spatial.distance
from scipy.stats import pearsonr
from scipy.spatial.distance import directed_hausdorff
import numpy as np
import matplotlib.pyplot as plt

############################# CONFIG #############################

# Working directory (an "out" directory is expected in this folder):
wd = "/vol/gf-grape/project_NoViSys/members/jenthein/_RV-G1G2/__Timeshift_Test_Phasenangleichung+Corr/"

# Workbook (IDs + counts):
wb = xlrd.open_workbook(wd + "data_matrix.xlsx")

# Timepoints (the same for both genotypes):
tps_x = [0,7,13,21,27,35,41,49,56,63,70,77]

# Gene(s) of interest: (expects a list f.ex.: ['gene1','gene2']; Wildcard is allowed: f.ex.: ['gene'])
goi = ['VIT'] #default: "VIT"

### General settings:
comp_xend = -28			# Threshold of comparison (-28 days), phenotype differs at a maximum of 28 days..
steps = 0.0001			# x-shifts steps (0.0001 days)

# Perform additional tests (just a print output):
hptp_check = 0                  # Basic tests (highest and lowest points) - TRUE / FALSE
add_tests_check = 0             # Additional tests (Pearsonr, Hausdorff, cdist functions) - TRUE / FALSE

# Debugging controls on (just a print output):
debug = 0                       # TRUE / FALSE

# Plot (each) transcript:
plot = 0			# turn on/off all plots
plot_G2 = 1			# Plot G2
plot_G1 = 1			# Plot G1
plot_shifted = 1                # plot shifted plot

# Save the plot?:
plot_save = 1                   # Should the plot be saved? (default = 1)

# Save the csv?:
csv_save = 1                    # Should the csv file be saved? (default = 1)

########################### CONFIG END ############################

sh = wb.sheet_by_index(0)

# Creation of a .csv file:

if csv_save == 1:
        fields=['goi','x-shift']
        with open(wd + "out/set.csv", "w") as f:
                writer = csv.writer(f)
                writer.writerow(fields)


### Main script:

for i in goi:
	
	goi = str(i)
	
	for row_number in xrange(sh.nrows):
	
		geneid_xls = str(sh.cell(row_number,0))
		geneid_xls = geneid_xls[7:-1]
		
		if goi in geneid_xls:	# Take just defined genes (and skip of header line, etc.)
			
			geneid = geneid_xls
			
			if debug == 1:
                            print "\n\n\nTranscript: ", geneid
				
			### File handling, calculation of average replicate values and saving of results for each gene:
		
			G1_y_means_tps = []
			G2_y_means_tps = []
			
			col_numbers = sh.ncols -6
		
			for col_number in range(1,col_numbers):
				if col_number == 1 or ((col_number % 6) == 0 and col_number >= 6):
					
					if col_number != 1:
						col_number = col_number +1
					
					tripl_G1 = []
					tripl_G2 = []
		
					val1 = str(sh.cell(row_number,col_number))[7:]
					tripl_G1.append(float(val1))
					val2 = str(sh.cell(row_number,col_number+1))[7:]
					tripl_G1.append(float(val2))
					val3 = str(sh.cell(row_number,col_number+2))[7:]
					tripl_G1.append(float(val3))
					
					val4 = str(sh.cell(row_number,col_number+3))[7:]
					tripl_G2.append(float(val4))
					val5 = str(sh.cell(row_number,col_number+4))[7:]
					tripl_G2.append(float(val5))
					val6 = str(sh.cell(row_number,col_number+5))[7:]
					tripl_G2.append(float(val6))
					
					# Control of each identified values from the datamatrix:
					if debug == 1:
						print val1, val2, val3, val4, val5, val6
					
					if len(tripl_G2) == 3:
						G2_y_means_tps.append(((tripl_G2[0]+tripl_G2[1]+tripl_G2[2])/3))
					if len(tripl_G1) == 3:
						G1_y_means_tps.append(((tripl_G1[0]+tripl_G1[1]+tripl_G1[2])/3))
		
			# Conversion of data into np.arrays:
			G2_data = np.column_stack((G2_y_means_tps, tps_x))
			G1_data = np.column_stack((G1_y_means_tps, tps_x))
		
			# # Control of test input data (means and timepoints):
			if debug == 1:
				print "\n"	
				print "Gene-ID:\n",geneid
				print "\n"	
				print "Spaetburgunder_y_means(" + str(len(G1_y_means_tps)) + "):\n",G1_data
				print "\n"
				print "Fruehburgunder_y_means(" + str(len(G2_y_means_tps)) + "):\n",G2_data
				print "\n"
				print "Timepoints_x(" + str(len(tps_x)) + "):\n", tps_x
				print "\n"
			
			
			### Manipulation of x values for G1 and initial cdist test:
			
			test_dist = sum(sum(scipy.spatial.distance.cdist(G1_data, G2_data, 'euclidean')))
			
			# Initial Shift calculation:
			
			tps_x_temp = []
			
			for i in tps_x:
				tps_x_temp.append(i - steps) # x-shift is performed
			G1_data_temp = np.column_stack((G1_y_means_tps, tps_x_temp))
			temp_test_dist_shifted = sum(sum(scipy.spatial.distance.cdist(G1_data_temp, G2_data, 'euclidean')))
			
			# Approximation to the optimal x-shift coefficient. In a perfect situation the test should be 0 (cdist)! 

			test_dist_shifted_res = []
			coeff_res = []
			j = 0
			coeff = steps
			
			while temp_test_dist_shifted <= test_dist:		# Iterating over the x axis (while the test gets better)
                                
				if coeff < comp_xend:                           # ..break if it reaches the threshold = 28)
					break
				
				coeff -= steps                                  # x steps
				
				if debug == 1:                                  # Control of the iterating  shift
					print coeff
				
				j += 1
				tps_x_temp = []
				for i in tps_x:
					tps_x_temp.append(i + coeff) # x-shift is performed
				G1_data_temp = np.column_stack((G1_y_means_tps, tps_x_temp))
				temp_test_dist_shifted = sum(sum(scipy.spatial.distance.cdist(G1_data_temp, G2_data, 'euclidean')))
				test_dist_shifted_res.append(temp_test_dist_shifted)
				coeff_res.append(coeff)
			
				
			### Plotting:
			# Plot (G2):	
			if plot_G2 == 1:			
				plt.scatter(tps_x,G2_y_means_tps, c="b", label='G2')
				plt.plot(tps_x,G2_y_means_tps, '-b')
			# Plot (G1):
			if plot_G1 == 1:				
				plt.scatter(tps_x,G1_y_means_tps, c="r", label='G1')
				plt.plot(tps_x,G1_y_means_tps, '-r')
			# Plot (altered state / x-shift):			
			if plot_shifted == 1:
				temp_label = 'G1(x-shift=' + str(coeff) + ')'
				plt.scatter(tps_x_temp,G1_y_means_tps, c="y", label=temp_label)
				plt.plot(tps_x_temp,G1_y_means_tps, '-y')
			plt.legend(loc='upper left')
			
			if plot_save == 1:
                            plt.savefig(wd + "out/" + geneid + '.png') # Save the Plot!
			
			if plot == 1:
				plt.show()
			plt.close()
			

			if coeff == comp_xend:	# If the test isnt getting better
				coeff = "NA"
				
			# Write x-shift in csv:				
			if csv_save == 1:
                                fields=[geneid,coeff]
                                with open(wd + "/out/set.csv", "a") as f:
                                        writer = csv.writer(f)
                                        writer.writerow(fields)
			
			if add_tests_check == 1:
			
                            # Pearsonr:
                            corr, value = pearsonr(G1_y_means_tps, G2_y_means_tps)
                            print "pearsonr = " + str(corr)
                            
                            # Hausdorff-Test:
                            print "\nHausdorff-Test(original): ",directed_hausdorff(G1_data, G2_data)[0]
                            print "\nHausdorff-Test(manipulated): ",directed_hausdorff(G1_data_temp, G2_data)[0]
                            print "\n"
                            
                            # Further tests (functions of cdist):
                            tests = ["euclidean", "braycurtis", "canberra", "chebyshev", "cityblock", "correlation", "cosine", "dice", "hamming", "jaccard", "kulsinski", "mahalanobis", "matching", "minkowski", "rogerstanimoto", "russellrao", "seuclidean", "sokalmichener", "sokalsneath", "sqeuclidean"]
                            for test in tests:
                                print "\n" + test + "(original): ",sum(sum(scipy.spatial.distance.cdist(G1_data, G2_data, test)))
                                print "\n" + test + "(manipulated): ",sum(sum(scipy.spatial.distance.cdist(G1_data_temp, G2_data, test)))
                                print "\n"

			if hptp_check == 1:
			
				# Peak comparison:
				G1_TP = G1_y_means_tps.index(min(G1_y_means_tps))
				G2_TP = G2_y_means_tps.index(min(G2_y_means_tps))
				diff_TP = abs(tps_x[G1_TP] - tps_x[G2_TP])
			
				G1_HP = G1_y_means_tps.index(max(G1_y_means_tps))
				G2_HP = G2_y_means_tps.index(max(G2_y_means_tps))
				diff_HP = abs(tps_x[G1_HP] - tps_x[G2_HP])
			
				# Delta calculation (amount of the number):
				diff = abs(diff_HP - diff_TP)
		
				## Controls:
				if debug == 1:
					#Values for TPs und HPs:
					print "\nG1_TP", G1_TP
					print "G2_TP", G2_TP
					print "G1_HP", G1_HP
					print "G2_HP", G2_HP
					# Deltas:
					print "\n"
					print "Diff(HP-TP):",diff
					print "diff_HP",diff_HP
					print "diff_TP",diff_TP
