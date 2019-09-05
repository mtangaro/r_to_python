#################################### REDI OUT TABLE ########################################################
#Region	Position	Reference	Strand	Coverage-q30	MeanQ	BaseCount[A,C,G,T]\	   	   #
#AllSubs	Frequency	gCoverage-q30	gMeanQ	gBaseCount[A,C,G,T]	gAllSubs	gFrequency #
############################################################################################################

import os, sys, argparse
from scipy import stats
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-c", action = 'store', dest = 'min_coverage', 
		type = int, default=10,  help='Coverage-q30')
parser.add_argument("-cpval", action = 'store', dest = 'pvalue_correction',
                type = int, default = 0, help = '1 --> Bonferroni correction / 2 --> Benjamini hochberg')
parser.add_argument("-input_file", action = 'store', dest = 'samples_informations_file',
		type = str, default= 'empty', help = 'Comma separated file  e.g: SRR3306830,Control \
		SRR3306829,Healthy...etc')
parser.add_argument("-f", action = 'store', dest = 'min_edit_frequency',
		type = float, default=0.1, help='Editing Frequency')
parser.add_argument("-mts", action = 'store', dest = 'min_sample_testing',
		type = float, default=50.0, help="min percentage of each sample category")
parser.add_argument("-sig", action = 'store', dest = 'only_significant',
		type = str, default = 'no', help = 'Return only significant editing events')
parser.add_argument("-linear", action = 'store_true', help = 'Enable linear model')

args = parser.parse_args()
min_coverage = args.min_coverage
min_edit_frequency = args.min_edit_frequency
min_sample_testing = args.min_sample_testing
only_significants = args.only_significant
pvalue_correction = args.pvalue_correction
samples_informations_file = args.samples_informations_file
enable_linear_model = args.linear

if args.samples_informations_file == 'empty':
	parser.error('sample_informations_file is MISSING!' + '\n' + \
	'Please type "python reditools_output_parser.py -h" for more details on usage of this script')

def Set_Chr_Nr(Chr):
    """ Sort by chromosome """
    if Chr: 
	New = Chr.lstrip('chr').split('_')[0]
        if New == 'X': New = 23
        elif New == 'Y': New = 24
        elif New == 'M': New = 25
        else: New = int(New)
    else:
        New = 0
    return New

def Sample_percentage(row):
	"""Percentage of samples from each type"""
	percentage = (len(filter(lambda x: x!= '-', row))/float(len(row)))*100
	return round(percentage)

def Sample_count(row):
	"""Number of samples from each type"""
	count = len(filter(lambda x: x!= '-', row))
	return count 

def get_bh(pvalue,siglevel):
	"""B-H correction """
	pvalue.sort()
	x=1
	y=0
	p=0
	for i in pvalue:
		nf=i[0]*len(pvalue)
		fdr=nf/x
		if fdr<=siglevel:
			i[1].append('True')
			p=i[0]
			y+=1
		else: i[1].append('False')
		x+=1
	return pvalue,y,p
	
def get_b(pvalue,siglevel): 
	"""Bonferroni correction"""
	pvalue.sort()
	y=0
	pp=1.0
	for i in pvalue:
		p=i[0]*len(pvalue)
		if p<=siglevel:
			i[1].append('True')
			y+=1
			if p<pp: pp=p
		else: i[1].append('False')
	return pvalue,y,pp

def only_sig(row):
	"""Returns only significant events"""
	if(row[-1] != '-' and row[-1] != 0.0 and row[-1] <= 0.05):
		print '\t'.join(map(str,row))

# If required the linear model is run through call_differential_editing_sites.py script
if enable_linear_model:

    cmd = 'python ./call_differential_editing_sites.py'
    os.system(cmd)

else:

    sample_informations = {}
    with open(samples_informations_file, 'r') as f:
        for line in f:
            if line.startswith('Sample'): #TODO REVERT BACK TO SRR!!!
                line = map(str.strip, line.split('\t'))
                sample_informations.setdefault(line[0], line[1])
    
    print sample_informations #TODO REMOVE!!!

    #cwd = filter(os.path.isdir, os.listdir(os.getcwd()))
    #all_available_sites = []
    #sample_edited_sites = {}
    #for directory in cwd:
    #    print directory
    #    path = list(os.walk(directory + '/editing/'))
    #    print path[0][0]
    #    print path[0][2][0]
    #    #table = path[1][0] + '/' + path[1][-1][-1] 
    #    table = str(path[0][0]) + str(path[0][2][0])
    #    print table
    #    with open(table,'r') as a:
    #        for line in a:
    #            if line.startswith('chr'):
    #                s = map(str.strip, line.split("\t"))
    #                site, freq, coverage = s[0] + "_" + s[1], s[8], s[4]
    #                if site not in all_available_sites: all_available_sites.append(site)
    #                if (int(coverage) >= min_coverage) and (float(freq) >= min_edit_frequency):
    #                	sample_edited_sites.setdefault((directory, site), []).append(freq)

    import pandas as pd
    import math

    all_available_sites = []
    sample_edited_sites = {}
    editing_file= '../data/editing_sites.txt'
    editing_table = pd.read_csv(editing_file,sep='\t')

    # iterating over rows using iterrows() function  
    for index, row in editing_table.iterrows():
        site = row['chromosome'] + "_" + str(row['position'])
        for col in editing_table.columns:
            # skip first 3 columns
            if col == 'chromosome': continue
            if col == 'position': continue
            if col == 'type_editing': continue

            # remove nan values if there
            if not isinstance(row[col], basestring):
                if math.isnan(float(row[col])): continue

            values =  str(row[col]).split("^")
            coverage =  int(values[2])
            freq = float(values[1])/float(values[2])
 
            if site not in all_available_sites: all_available_sites.append(site)
            if (int(coverage) >= min_coverage) and (float(freq) >= min_edit_frequency):
                sample_edited_sites.setdefault((col, site), []).append(freq)

    print sample_edited_sites #TODO REMOVE!!!

    table_columns = map(lambda x: x + '_' + sample_informations[x], sorted(sample_informations.keys()))
    
    disease = [i for i in table_columns if i.upper().find('AUTISM') != -1] #TODO REVERT BACK TO DIS!!!
    controls = [i for i in table_columns if i.upper().find('NORMAL') != -1] #TODO REVERT BACK TO CTRL!!!
 
    header = ['SITES'] + controls + disease + ['[num_controls/num_disease]'] + ['delta_diff'] + \
    ['pvalue (Mannwhitney)'] 
    
    if pvalue_correction == 1:
    	header += ['pvalue Bonferroni corrected']
    if pvalue_correction == 2:
    	header += ['pvalue BH corrected']
    
    print '\t'.join(header)
    
    for chrom in sorted(all_available_sites, key = lambda x: Set_Chr_Nr(x)):
    	row = [chrom]
    	for col in header[1:header.index('[num_controls/num_disease]')]:
    		row.append(sample_edited_sites.get((col.split('_')[0],chrom), ['-'])[0])
    	ctrls = zip(*(zip(controls,row[1:])))[1]
    	dss = zip(*(zip(disease,row[len(ctrls)+1:])))[1] 
    	row.append(str([Sample_count(ctrls), Sample_count(dss)]))
    	if Sample_percentage(ctrls) >= min_sample_testing and Sample_percentage(dss) >= min_sample_testing:
    		ctrls_mean = sum(map(float, filter(lambda x : x!= '-', ctrls)))/len(filter(lambda x : x!= '-', ctrls))
    		dss_mean = sum(map(float, filter(lambda x : x!= '-', dss)))/len(filter(lambda x: x!= '-', dss))
    		delta_diff =  abs(ctrls_mean - dss_mean)
    		pvalue=stats.mannwhitneyu(ctrls,dss, alternative='two-sided')
    		row.append(round(delta_diff, 3))
    		row.append(str(round(pvalue[1], 3)))
    		correction_argmnt = [(pvalue[1], map(float, filter(lambda x : x!= '-', ctrls+dss)))]
    		if pvalue_correction == 1:
    			row.append(round(get_b(correction_argmnt, 0.05)[-1], 6))
    		elif pvalue_correction == 2:
    			row.append(round(get_bh(correction_argmnt, 0.05)[-1], 6))
    	else:
    		if pvalue_correction == 0:
    			row += ['-', '-']
    		else:
    			row += ['-', '-', '-']
    	if pvalue_correction != 0 and only_significants == 'yes':
    		only_sig(row)
    	else:
    		print '\t'.join(map(str,row))
    
    
    
