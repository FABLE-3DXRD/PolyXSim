#!/usr/bin/env python
#This module has been made by Gaelle for the GUI. This is a copy of the script.
# Modules to import 
import sys,os
from polyxsim import check_input
from polyxsim import file_io
from polyxsim import find_refl
from polyxsim import generate_grains
from polyxsim import make_image
from polyxsim import make_imagestack
from polyxsim import reflections
from polyxsim import help_input
#import logging
#logging.basicConfig(level=logging.INFO,format='\n%(levelname)s: %(message)s')

from optparse import OptionParser

def get_options(parser):

    parser = OptionParser()
    parser.add_option("-i", "--input", action="store",
                      dest="filename", type="string",
                      help="Name of the file containing the input parameters")
    parser.add_option("-d", "--debug", action="store_true",
                      dest="debug",default =False,
                      help="Run in debug mode")
    parser.add_option("-p", "--print", action="store_true",
                      dest="print_input",default =False,
                      help="Show input parameters and syntax")
    parser.add_option("-k","--killfile", action="store",
                      dest="killfile", default=None, type="string",
                      help="Name of file to create halt PolyXSim")
    
    options , args = parser.parse_args()
    return options

def run(options):
    # Check for print_input
    try:
        options.print_input
    except:
        options.print_input = None;
    if options.print_input:
        print help_input.show_input()
        sys.exit()

    # Check if filename is specified
    try:
        options.filename
    except:
        options.filename = None;
    if options.filename == None:
        print "\nNo input file supplied [-i filename]\n"
        #Gaelle comment : ? sys.exit() and add raise error instead
        sys.exit()
    #print 'options = ',options
    #print '\n'

    # Check killfile does not exist
    try:
        options.killfile
    except:
        options.killfile = None;
    if options.killfile is not None and os.path.exists(options.killfile):
        print "The purpose of the killfile option is to create that file"
        print "only when you want PolyXsim to stop"
        print "If the file already exists when you start PolyXsim, it is"
        print "stopped immediately"
        raise ValueError("Your killfile "+options.killfile+" already exists")
    
    # Is the input file available?
    
    # Read and check input
    
    # Make instance of parse_input class
    print('Reading input\n')
    
    myinput = check_input.parse_input(input_file=options.filename)
    
    myinput.read()                                # read input file
    
    print('Checking input\n')
    myinput.check()                               # check validity of input
    check_input.interrupt(options.killfile)
    
    
#     if myinput.missing == True:                   # if problem exit
#         print('MISSING ITEMS')
#         sys.exit()
    
    if len(myinput.errors) > 0:
        myinput.show_errors()
        sys.exit()

    print('Initialize parameters etc\n')
    myinput.initialize()                            # if ok initialize
    check_input.interrupt(options.killfile)
    
    # Generate reflections
    hkl = []
    
    for phase in myinput.param['phase_list']:
        if  ('structure_phase_%i' %phase) in myinput.param:
            xtal_structure = reflections.open_structure(myinput.param,phase)
            #print 'UNIT CELL', myinput.param['unit_cell_phase_%i' %phase]
            print('Generating miller indices')
            hkl_tmp = reflections.gen_miller(myinput.param,phase)
            if myinput.param['structure_factors'] != 0:
                print('Structure factor calculation')
                hkl.append(reflections.calc_intensity(hkl_tmp,
                                                      xtal_structure,
                                                      options.killfile))
            else:
                hkl.append(reflections.add_intensity(hkl,myinput.param))
                print('No structure factor calculation')
        else:
            hkl_tmp = reflections.gen_miller(myinput.param,phase)
            hkl.append(reflections.add_intensity(hkl_tmp,myinput.param))
    
        check_input.interrupt(options.killfile)
    #    if options.killfile is not None and os.path.exists(options.killfile):
    #        raise KeyboardInterrupt()
    
    generate_grains.generate_grains(myinput.param)
    check_input.interrupt(options.killfile)
    print('Write grains file')
    file_io.write_grains(myinput.param)
    check_input.interrupt(options.killfile)
    print('Write res file')
    file_io.write_res(myinput.param)
    check_input.interrupt(options.killfile)
    
    if '.hkl' in myinput.param['output']:
        print('Write hkl file')
        file_io.write_hkl(myinput.param,hkl)
    if '.fcf' in myinput.param['output']:
        print('Write fcf file')
        file_io.write_fcf(myinput.param,hkl)
    if '.ubi' in myinput.param['output']:
        print('Write UBI file')
        file_io.write_ubi(myinput.param)
    if '.par' in myinput.param['output']:
        print('Write detector.par file')
        file_io.write_par(myinput.param)
    check_input.interrupt(options.killfile)
    
    
    # Determine the reflection parameters for grains
    graindata = find_refl.find_refl(myinput.param,hkl,options.killfile)
    graindata.frameinfo = myinput.frameinfo
    print('Determine reflections positions')
    graindata.run()
    if '.ref' in myinput.param['output']:
        print('Write reflection file')
        graindata.save()
    if '.gve' in myinput.param['output']:
        print('Write g-vector file')
        graindata.write_gve()
    if '.ini' in myinput.param['output']:
        print('Write GrainSpotter ini file - Remember it is just a template')
        graindata.write_ini()
    if '.flt' in myinput.param['output']:
        print('Write filtered peaks file')
        graindata.write_flt()

    if myinput.param['make_image'] == 1:
        if  myinput.param['peakshape'][0] == 2:
            image = make_imagestack.make_image(graindata,options.killfile)
            image.setup_odf()
            image.make_image_array()
            image.make_image()
            image.correct_image()
        else:
            image = make_image.make_image(graindata,options.killfile)
            image.make_image()


