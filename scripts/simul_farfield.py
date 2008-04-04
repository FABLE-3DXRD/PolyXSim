#!/usr/bin/env python

# Modules to import 
import sys
from Simul_farfield import check_input
from Simul_farfield import find_refl
from Simul_farfield import generate_grains
from Simul_farfield import make_image
import logging
logging.basicConfig(level=logging.INFO,format='\n%(levelname)s: %(message)s')

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--input", action="store",
                  dest="filename", type="string",
                  help="Name of the file containing the input parameters")
parser.add_option("-d", "--debug", action="store_true",
                  dest="debug",default =False,
                  help="Run in debug mode")
options , args = parser.parse_args()

if options.filename == None:
    parser.print_help()
    print "\nNo input file supplied [-i filename]\n"
    sys.exit()



# Is the input file available?

# Read and check input
print '\n'
myinput = check_input.parse_input(input_file=options.filename)  # Make instance of parse_input class
try:
    myinput.read()                                # read input file
except:
    sys.exit()

myinput.check()                               # check validity of input
if myinput.missing == True:                   # if problem exit
    logging.info('MISSING ITEMS')
    sys.exit()
myinput.initialize()                            # if ok initialize

#print myinput.param
generate_grains.generate_grains(myinput.param)
generate_grains.save_grains(myinput.param)
generate_grains.save_ubi(myinput.param)

# Determine the reflection parameters for grains
graindata = find_refl.find_refl(myinput.param)
graindata.frameinfo = myinput.frameinfo
graindata.run()
graindata.save()
graindata.write_gve()

if myinput.param['make_image'] != 0:
	image = make_image.make_image(graindata)
	image.make_image()
