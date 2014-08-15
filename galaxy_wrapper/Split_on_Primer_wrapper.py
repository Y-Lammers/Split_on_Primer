#!/usr/bin/env python

# Small python wrapper that passes the xml commands to the HTS-barcode-checker script
# This requires the HTS-barcode-checker script to be present in the $PATH variable

# import the used modules
import sys, os
from subprocess import Popen, PIPE

# call the Split_on_Primer script with the parameters and commands from
# the Split_on_Primer.xml galaxy file. Catch the stdout containing the filenames.
# open the stat output file, which will contain the # reads per marker
path = Popen(['Split_on_Primer.py'] + sys.argv[1:-4], stdout=PIPE, stderr=PIPE)
filelist = path.communicate()[0].strip().split('\n')
stat_file = open(sys.argv[-3], 'w')

# parse through the output filenames and rename them so they can be picked up
# by the galaxy workflow.
for file in filelist:
	base = os.path.splitext(os.path.basename(file))[0]
	stat_file.write('{0}\t{1}\n'.format(base, sum(1 for line in open(file))/4))
	dest = '{0}_{1}_{2}_{3}_{4}'.format('primary', sys.argv[-2], base, 'visible', 'fastq')
	dest = os.path.join(sys.argv[-1], dest)
	os.rename(file, dest)
stat_file.close()
