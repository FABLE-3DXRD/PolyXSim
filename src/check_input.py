#!/usr/bin/env python

#
# Checking input  
#

from string import split
import sys, os 
from fabio import deconstruct_filename, jump_filename
import variables

import numpy as N
import logging

logging.basicConfig(level=logging.DEBUG,format='%(levelname)s %(message)s')

class parse_input:
    def __init__(self,input_file = None):
        self.filename = input_file
        self.entries = {}
        self.grainno = 0
        self.entries['pos_grains'] = []
        self.no_pos = 0 # keeping track of no of grain position
        # Experimental setup
        self.needed_items = {
                    'wavelength' : 'Missing input: wavelenght [wavelength in angstrom]',
                    'distance'   : 'Missing input: distance [sample-detector distance in mm)]',
                    'dety_center': 'Missing input: dety_center [beamcenter, y in pixel coordinatees]',
                    'detz_center': 'Missing input: detz_center [beamcenter, z in pixel coordinatees]',
                    'y_size'     : 'Missing input: y_size [Pixel size y in mm]',
                    'z_size'     : 'Missing input: z_size [Pixel size z in mm]',
                    'dety_size'  : 'Missing input: dety_size [detector y size in pixels]',
                    'detz_size'  : 'Missing input: detz_size [detector z size in pixels]',
                    'omega_start'       : 'Missing input: omega_start [Omega start in degrees]',
                    'omega_end'       : 'Missing input: omega_end [Omega end in degrees]',
                    'omega_step'      : 'Missing input: omega_step [Omega step size in degrees]',
                    #'sgno'       : 'Missing input: sgno [space group number has to be given when sysconditions are not]',
                    'unit_cell'  : 'Missing input: unit_cell [unit cell parameters: a,b,c,alpha,beta, gamma]',
                    'no_grains'  : 'Missing input: no_grains [number of grains]',
                    #'U_grains'   : 'Missing input: U_grains [crystal orientations]',
                    'direc'      : 'Missing input: direc [directory to save output]',
                    'theta_min'   : 'Missing input: theta_min [Minimum theta angle for generation of reflections in degrees]',
                    'theta_max'   : 'Missing input: theta_max [Maximum theta angle for generation of reflections in degrees]',
                    'imagefile' :  'Missing input: imagefile [the first image file to be processed]'
                    }
        self.optional_items = {
            'sgno': 1,
            'pos_grains' : [[0, 0, 0]],
            'tilt_x'     : 0,
            'tilt_y'     : 0,
            'tilt_z'     : 0,
            'beampol_factor' : 1,
            'beampol_angle' : 0.0,
            'spatial' : None,
            'flood' : None,
            'dark' : None,
            'darkoffset' : None
            }

        
    def read(self):     
        try:
            f = open(self.filename,'r')
        except:
            print 'No file named %s' %self.filename 
            return False
        
        self.input = f.readlines()
        f.close()

        for lines in self.input:
            if lines[0] != '#' and lines[0] != '%':
                if lines.find('=') > -1:
                    line = split(lines,'=')
                    if line[1].find('#') > -1:
                        line[1] = split(line[1],'#')[0]
                    if line[1].find('%') > -1:
                        line[1] = split(line[1],'%')[0]

                    key = line[0].strip()
                    val = line[1].strip()
                    if key == 'U_grains':
                        key = 'U_grains_%i' %self.grainno
                        self.grainno = self.grainno + 1
                    if key == 'pos_grains':
                        spos = self.entries[key]
                        spos.append(eval(val))
                        self.entries[key] = spos
                    else:
                        self.entries[key] = eval(val)
                

                
    def check(self):
        self.missing = False

        for item in self.needed_items:
            if item not in self.entries:
                print self.needed_items[item]
                self.missing = True
        assert self.grainno == self.entries['no_grains']
                
#        for i in range(self.entries['no_grains'])
#        
    def initialize(self):
        fileinfo = deconstruct_filename(self.entries['imagefile'])
        self.entries['filetype'] = fileinfo.format
        self.entries['stem'] = fileinfo.stem
        if 'start_frame' not in self.entries:
            
            self.entries['start_frame'] = fileinfo.num
        self.entries['filetype'] = fileinfo.format
        for item in self.optional_items:
            if (item not in self.entries):
                self.entries[item] = self.optional_items[item]
            if (self.entries[item] == []):
                self.entries[item] = self.optional_items[item]*self.entries['no_grains']

        # Generate FILENAME of frames
        omega_step = self.entries['omega_step']
        omega_start  = self.entries['omega_start']
        omega_end  = self.entries['omega_end']
        omega_sign = self.entries['omega_sign']
        start_frame = self.entries['start_frame']
        omegalist = omega_sign*N.arange(omega_start,omega_end+omega_step,omega_step)
        nframes = int((omega_end-omega_start)/omega_step)
        omegalist.sort()
        i=0
        logging.info("Generating frame data...")

        #Initialize frameinfo container
        self.frameinfo = [] 
        
        if omega_sign > 0:
            filerange = N.arange(start_frame,start_frame+nframes)
        else:
            filerange = N.arange((start_frame-1)+nframes,(start_frame-1),omega_sign)
            # reverse omega_start/omega_end
            self.entries['omega_end'] = omega_start*omega_sign 
            self.entries['omega_start'] = omega_end*omega_sign

        for no in filerange:
            self.frameinfo.append(variables.frameinfo_cont(no))
            self.frameinfo[i].name = jump_filename(self.entries['imagefile'],no)
            self.frameinfo[i].omega = omegalist[i];
            self.frameinfo[i].nrefl = 0 # Initialize number of reflections on frame
            self.frameinfo[i].refs = [] # Initialize number of reflections on frame
            i += 1
        logging.debug("Printing frameinfo...")

        # Does output directory exist?
        if not os.path.exists(self.entries['direc']):
            os.mkdir(self.entries['direc'])
            
# if exist('theta_max') == 0
#         disp('Warning - missing input: theta_max [Maximum theta angle for generation of reflections in degrees]')
#         % Find maximum theta for generation of all possible reflections on
#         % the detector from the detector specs
#         dety_center_mm = dety_center * y_size;
#         detz_center_mm = detz_center * y_size;
#         dety_size_mm = dety_size * y_size;
#         detz_size_mm = detz_size * y_size;
#         c2c(1) = sqrt((dety_center_mm-dety_size_mm)^2 + (detz_center_mm-det_zsize_mm)^2);
#         c2c(2) = sqrt((dety_center_mm-dety_size_mm)^2 + (detz_center_mm-0)^2);
#         c2c(3) = sqrt((dety_center_mm-0)^2 + (detz_center_mm-det_zsize_mm)^2);
#         c2c(4) = sqrt((dety_center_mm-0)^2 + (detz_center_mm-0)^2);
#         c2c_max = max(c2c);
#         theta_max = atan(c2c_max/distance)/2 * 180/pi;
#         disp(['NOTICE: To make full detector coverage sets theta_max = ',num2str(theta_max)])
#         clear dety_center_mm detz_center_mm dety_size_mm detz_size_mm c2c c2c_max
# end



# % Generate FILENAME of frame
# for no=1:nframes
#     file(no).name = sprintf('%s/%s%0.4d.tif',direc,fileprefix,no);
# end


if __name__=='__main__':

    #import check_input
    try:
        filename = sys.argv[1] 
    except:
        print 'Usage: check_input.py  <input.inp>'
        sys.exit()

    myinput = parse_input(input_file = filename)
    myinput.read()
    print myinput.entries
    myinput.check() 
    if myinput.missing == True:
        print 'MISSING ITEMS'
    myinput.evaluate()
    print myinput.entries
