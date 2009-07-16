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
import logging
logging.basicConfig(level=logging.INFO,format='\n%(levelname)s: %(message)s')

from optparse import OptionParser

def get_options(parser):

    parser = OptionParser()
    parser.add_option("-i", "--input", action="store",
                      dest="filename", type="string",
                      help="Name of the file containing the input parameters")
    parser.add_option("-d", "--debug", action="store_true",
                      dest="debug",default =False,
                      help="Run in debug mode")
    parser.add_option("-k","--killfile", action="store",
                      dest="killfile", default=None, type="string",
                      help="Name of file to create halt PolyXSim")
    
    options , args = parser.parse_args()
    return options

def run(options):
    if options.filename == None:
        parser.print_help()
        print "\nNo input file supplied [-i filename]\n"
        #Gaelle comment : ? sys.exit() and add raise error instead
        raise ValueError("No input file supplied [-i filename]")
    #print 'options = ',options
    #print '\n'
    if options.killfile is not None and os.path.exists(options.killfile):
        print "The purpose of the killfile option is to create that file"
        print "only when you want PolyXsim to stop"
        print "If the file already exists when you start PolyXsim, it is"
        print "stopped immediately"
        raise ValueError("Your killfile "+options.killfile+" already exists")
    
    
    # Is the input file available?
    
    # Read and check input
    
    # Make instance of parse_input class
    logging.info('Reading input\n')
    
    myinput = check_input.parse_input(input_file=options.filename)
    
    myinput.read()                                # read input file
    #except:
    #    logging.error('Cannot read input - exit')
    #sys.exit()
    
    logging.info('Checking input\n')
    myinput.check()                               # check validity of input
    check_input.interrupt(options.killfile)
    
    
    if myinput.missing == True:                   # if problem exit
        logging.info('MISSING ITEMS')
        sys.exit()
    
    logging.info('Initialize parameters etc\n')
    myinput.initialize()                            # if ok initialize
    check_input.interrupt(options.killfile)
    
    # Generate reflections
    hkl = []
    
    for phase in myinput.param['phase_list']:
        if  ('structure_phase_%i' %phase) in myinput.param:
            xtal_structure = reflections.open_structure(myinput.param,phase)
            #print 'UNIT CELL', myinput.param['unit_cell_phase_%i' %phase]
            logging.info('Generating miller indices')
            hkl_tmp = reflections.gen_miller(myinput.param,phase)
            if myinput.param['structure_factors'] != 0:
                logging.info('Structure factor calculation')
                hkl.append(reflections.calc_intensity(hkl_tmp,
                                                      xtal_structure,
                                                      options.killfile))
            else:
                hkl.append(reflections.add_intensity(hkl,myinput.param))
                logging.info('No structure factor calculation')
        else:
            hkl_tmp = reflections.gen_miller(myinput.param,phase)
            hkl.append(reflections.add_intensity(hkl_tmp,myinput.param))
    
        check_input.interrupt(options.killfile)
    #    if options.killfile is not None and os.path.exists(options.killfile):
    #        raise KeyboardInterrupt()
    
    generate_grains.generate_grains(myinput.param)
    check_input.interrupt(options.killfile)
    logging.info('Write grains file')
    file_io.write_grains(myinput.param)
    check_input.interrupt(options.killfile)
    logging.info('Write res file')
    file_io.write_res(myinput.param)
    check_input.interrupt(options.killfile)
    
    if '.hkl' in myinput.param['output']:
        logging.info('Write hkl file')
        file_io.write_hkl(myinput.param,hkl)
    if '.fcf' in myinput.param['output']:
        logging.info('Write fcf file')
        file_io.write_fcf(myinput.param,hkl)
    if '.ubi' in myinput.param['output']:
        logging.info('Write UBI file')
        file_io.write_ubi(myinput.param)
    if '.par' in myinput.param['output']:
        logging.info('Write detector.par file')
        file_io.write_par(myinput.param)
    check_input.interrupt(options.killfile)
    
    
    
    # Determine the reflection parameters for grains
    graindata = find_refl.find_refl(myinput.param,hkl,options.killfile)
    graindata.frameinfo = myinput.frameinfo
    logging.info('Determine reflections positions')
    graindata.run()
    logging.info('Save reflections')
    graindata.save()
    if '.gve' in myinput.param['output']:
        logging.info('Write g-vector file')
        graindata.write_gve()
    if '.ini' in myinput.param['output']:
        logging.info('Write GrainSpotter ini file - Remember it is just a template')
        graindata.write_ini()
    if '.flt' in myinput.param['output']:
        logging.info('Write filtered peaks file')
        graindata.write_flt()
    #Gaelle to Henning ? image should be indented with my editor. Is it correct?
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


