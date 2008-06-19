#!/usr/bin/env python

#
# Checking input  
#

from string import split
import sys, os 
import variables
from xfab import tools

import numpy as n
import logging

logging.basicConfig(level=logging.DEBUG,format='%(levelname)s %(message)s')

class parse_input:
    def __init__(self,input_file = None):
        self.filename = input_file
        self.param = {}
        self.grainno = 0
#        self.no_pos = 0 # keeping track of no of grain position
        # Experimental setup
        self.needed_items = {
                    'wavelength' : 'Missing input: wavelength [wavelength in angstrom]',
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
                    'no_grains'  : 'Missing input: no_grains [number of grains]',
                    'direc'      : 'Missing input: direc [directory to save output]',
                                        }
        self.optional_items = {
            'sgno': None,
            'sgname': None,
            'tilt_x'     : 0,
            'tilt_y'     : 0,
            'tilt_z'     : 0,
            'wedge': 0.0,
            'beampol_factor' : 1,
            'beampol_angle' : 0.0,
            'start_frame': 0,
            'omega_sign': 1,
            'noise' : 0,
            'psf': 0,
            'make_image': 1,
            'bg': 0,
            'peakshape': [0,0],
            'spatial' : None,
            'flood' : None,
            'dark' : None,
            'darkoffset' : None,
            'gen_U'   : 0,
            'gen_pos' : [0,0],
            'gen_eps' : [0,0,0,0,0],   
            'gen_size': [0,0,0,0],
            'sample_xyz': None,
            'sample_cyl': None,
            'direc': '.',
            'prefix': 'test',
            'odf_type' : 1,
            'odf_scale' : 0.02,
            'mosaicity' : 0.2,
            'theta_min' : 0.0,
            'theta_max' : None,
            'unit_cell' : None,
            'structure_file': None,
            'odf_file': None, 
            'output': None
            }

        
    def read(self):     
        try:
            f = open(self.filename,'r')
        except IOError:
            logging.error('No file named %s' %self.filename)
            raise IOError
        
        self.input = f.readlines()
        f.close()

        for lines in self.input:
            if lines.find('#') != 0:
                if lines.find('#') > 0:
                    lines = split(lines,'#')[0]
                line = split(lines)
                if len(line) != 0:
                    key = line[0]
                    val = line[1:]
                    valtmp = '['
                    if len(val) > 1:
                        for i in val:
                            valtmp = valtmp + i +','
							
                        val = valtmp + ']'
                    else:
                        val = val[0]
                    self.param[key] = eval(val)
					
           

                
    def check(self):
        self.missing = False

# check that all needed items are present
        for item in self.needed_items:
            if item not in self.param:
                print self.needed_items[item]
                self.missing = True
		
# set all non-read items to defaults
        for item in self.optional_items:
            if (item not in self.param):
                self.param[item] = self.optional_items[item]
            if (self.param[item] == []):
                self.param[item] = self.optional_items[item]*self.param['no_grains']

# assert that the correct number of arguments are given
        for key in self.param:
            val = self.param[key] 
            if val != None and key != 'output':
                if key == 'peakshape':
                    assert len(val) <= 3, 'Wrong number of arguments for %s' %key
                elif key == 'sample_cyl' or key == 'gen_pos':
                    assert len(val) == 2 , 'Wrong number of arguments for %s' %key
                elif key == 'sample_xyz' or 'pos_grains' in key:
                    assert len(val) == 3, 'Wrong number of arguments for %s' %key
                elif key == 'gen_size':
                    assert len(val) == 4, 'Wrong number of arguments for %s' %key
                elif key == 'gen_eps':
                    assert len(val) == 5, 'Wrong number of arguments for %s' %key
                elif key == 'unit_cell' or 'eps_grains' in key:
                    assert len(val) == 6, 'Wrong number of arguments for %s' %key
                elif 'U_grains' in key:
                    assert len(val) == 9, 'Wrong number of arguments for %s' %key
                    # reshape U-matrices
                    self.param[key] = n.array(self.param[key])
                    self.param[key].shape = (3,3)
                else:
                    assert type(val) != list, 'Wrong number of arguments for %s' %key

# read U, pos, eps and size for all grains		
        grain_list_U = []
        grain_list_pos = []
        grain_list_eps = []
        grain_list_size = []
        no_grains = self.param['no_grains']

        for item in self.param:
            if '_grains_' in item:
                if 'U' in item:
                    grain_list_U.append(eval(split(item,'_grains_')[1]))
                elif 'pos' in item:
                    grain_list_pos.append(eval(split(item,'_grains_')[1]))
                elif 'eps' in item:
                    grain_list_eps.append(eval(split(item,'_grains_')[1]))
                elif 'size' in item:
                    grain_list_size.append(eval(split(item,'_grains_')[1]))
					
# assert that input U, pos, eps size are correct in format
# (same number of grains and same specifiers or else not input) 
        grain_list_U.sort()
        grain_list_pos.sort()
        grain_list_eps.sort()
        grain_list_size.sort()
        if len(grain_list_U) != 0 and self.param['gen_U'] == 0:
            assert len(grain_list_U) == no_grains, \
                'Input number of grains does not agree with number\n' +\
                ' of U_grains, check for multiple names'
            self.param['grain_list'] = grain_list_U
            if len(grain_list_pos) != 0 and self.param['gen_pos'][0] == 0:
                assert grain_list_U == grain_list_pos, \
                    'Specified grain number for U_grains and pos_grains disagree'
                if len(grain_list_eps) != 0 and self.param['gen_eps'][0] == 0:
                    assert grain_list_U == grain_list_eps, \
                        'Specified grain number for U_grains and eps_grains disagree'
                if len(grain_list_size) != 0 and self.param['gen_size'][0] == 0:
                    assert grain_list_U == grain_list_size, \
                        'Specified grain number for U_grains and size_grains disagree'
        else:
            if len(grain_list_pos) != 0 and self.param['gen_pos'][0] == 0:
                assert len(grain_list_pos) == no_grains, \
                    'Input number of grains does not agree with number\n'+\
                    ' of pos_grains, check for multiple names'
                self.param['grain_list'] = grain_list_pos
                if len(grain_list_eps) != 0 and self.param['gen_eps'][0] == 0:
                    assert grain_list_pos == grain_list_eps, \
                        'Specified grain number for pos_grains and eps_grains disagree'
                    if len(grain_list_size) != 0 and self.param['gen_size'][0] == 0:
                        assert grain_list_pos == grain_list_size, \
                            'Specified grain number for pos_grains and size_grains disagree'
            elif len(grain_list_eps) != 0 and self.param['gen_eps'][0] == 0:
                assert len(grain_list_eps) == no_grains, \
                    'Input number of grains does not agree with number'+\
                    ' of eps_grains, check for multiple names'
                self.param['grain_list'] = grain_list_eps
                if len(grain_list_size) != 0 and self.param['gen_size'][0] == 0:
                    assert grain_list_eps == grain_list_size, \
                        'Specified grain number for eps_grains and size_grains disagree'
            elif len(grain_list_size) != 0 and self.param['gen_size'][0] == 0:
                assert len(grain_list_size) == no_grains, \
                    'Input number of grains does not agree with number\n'+\
                    ' of size_grains, check for multiple names'
                self.param['grain_list'] = grain_list_size
            else:
                self.param['grain_list'] = range(no_grains)

# assert that all information needed to generate grains is present	
        assert len(grain_list_U) != 0 or self.param['gen_U'] != 0,\
            'Information on U generations missing'
        assert len(grain_list_pos) != 0 or self.param['gen_pos'][0] != 0,\
            'Information on position generation missing'
        assert len(grain_list_eps) != 0 or self.param['gen_eps'][0] != 0,\
            'Information on strain tensor generation missing'
        assert len(grain_list_size) != 0 or self.param['gen_size'][0] != 0,\
            'Information on grain size generation missing'
			
#assert that not both sample_xyz and sample_cyl are given
        if self.param['sample_xyz'] != None:
            assert self.param['sample_cyl'] == None, 'Both sample_xyz and sample_cyl are given'
			
# assert that mean grain size != 0 and if mean > 0 then min < mean < max, assure that min non-negative
        if self.param['gen_size'][0] != 0:
            assert self.param['gen_size'][1] != 0, 'Invalid gen_size command, mean size 0'
            if self.param['gen_size'][1] > 0:
                assert self.param['gen_size'][2] < self.param['gen_size'][1], \
                    'grain_min larger than grain_size'
                assert self.param['gen_size'][3] > self.param['gen_size'][1], \
                    'grain_max smaller than grain_size'
                if self.param['gen_size'][2] < 0:
                    self.param['gen_size'][2] = 0
                    
#check that the given grain_size and no_grains are consistent with sample_vol, adjust max to sample size
        if self.param['sample_xyz'] != None:
            assert self.param['sample_xyz'][0] > 0 and self.param['sample_xyz'][1] > 0 and self.param['sample_xyz'][2] > 0,\
                'Invalid sample_xyz <= 0'
            self.param['sample_vol'] = self.param['sample_xyz'][0]*\
                                       self.param['sample_xyz'][1]*\
                                       self.param['sample_xyz'][2]
            diam_limit = (6*self.param['sample_vol']/\
                         (n.exp(.5)*n.pi*self.param['no_grains']))**(1/3.)
            assert abs(self.param['gen_size'][1]) < diam_limit, \
                'The sample volume is too small to contain the '+\
                'specified number of grains with the given grain size'
            if len(grain_list_size) > 0:
                for i in range(len(grain_list_size)):
                    assert self.param['size_grains_%s' %(self.param['grain_list'][i])] < diam_limit, \
                'The sample volume is too small to contain the '+\
                'specified number of grains with the given grain size'
            self.param['gen_size'][3] = min(self.param['sample_xyz'][0],
                                            self.param['sample_xyz'][1],
                                            self.param['sample_xyz'][2])
        elif self.param['sample_cyl'] != None:
            assert self.param['sample_cyl'][0] > 0 and self.param['sample_cyl'][1] > 0,\
                'Invalid sample_cyl <= 0'
            self.param['sample_vol'] = n.pi*self.param['sample_cyl'][0]*\
                                       self.param['sample_cyl'][0]*self.param['sample_cyl'][1]
            diam_limit = (6*self.param['sample_vol']/\
                         (n.exp(.5)*n.pi*self.param['no_grains']))**(1/3.)
            assert abs(self.param['gen_size'][1]) < diam_limit, \
                'The sample volume is too small to contain the\n '+\
                'specified number of grains with the given grain size'
            if len(grain_list_size) > 0:
                for i in range(len(grain_list_size)):
                    assert self.param['size_grains_%s' %(self.param['grain_list'][i])] < diam_limit, \
                'The sample volume is too small to contain the '+\
                'specified number of grains with the given grain size'
            self.param['gen_size'][3] = min(2*self.param['sample_cyl'][0],
                                            self.param['sample_cyl'][1])
        else:					
            self.param['sample_vol'] = None

#check that a file name with the odf file is input is odf_type chosen to be 2.
        if self.param['peakshape'][0] == 2:
            if self.param['odf_type'] == 2:
                assert self.param['odf_file'] != None, 'No filename given for ODF'

#If no structure file is given - unit_cell should be               
        if self.param['structure_file'] == None:
            print 'NO structure file'
            assert self.param['unit_cell'] != None, \
                'Missing input: structure_file or unit_cell' 
            assert self.param['sgno'] != None or self.param['sgname'] != None , \
                'Missing input: no space group information, please input either sgno or sgname' 
            from xfab import sg
            if self.param['sgno'] == None:
                self.param['sgno'] = sg.sg(sgname = self.param['sgname']).no
            else:
                self.param['sgname'] = sg.sg(sgno = self.param['sgno']).name
                
			
    def initialize(self):
        # Frame generation
        if self.param['make_image'] != 0:
            if self.param['output'] == None:
                self.param['output'] = '.edf'
            if ('.edf' not in self.param['output'] and '.tif' not in self.param['output']):
                self.param['output'].append('.edf')
     
        # Does output directory exist?
        if not os.path.exists(self.param['direc']):
            os.mkdir(self.param['direc'])

	    # Generate FILENAME of frames
        omega_step = self.param['omega_step']
        omega_start  = self.param['omega_start']
        omega_end  = self.param['omega_end']
        omega_sign = self.param['omega_sign']
        start_frame = self.param['start_frame']
        omegalist = omega_sign*n.arange(omega_start,omega_end+omega_step,omega_step)
        nframes = int((omega_end-omega_start)/omega_step)
        omegalist.sort()
#        i=0
#        logging.info("Generating frame data...")
        #Initialize frameinfo container
        self.frameinfo = [] 
        
        if omega_sign > 0:
            filerange = n.arange(start_frame,start_frame+nframes)
        else:
            filerange = n.arange((start_frame-1)+nframes,(start_frame-1),omega_sign)
            # reverse omega_start/omega_end
            self.param['omega_end'] = omega_start*omega_sign 
            self.param['omega_start'] = omega_end*omega_sign

        for no in filerange:
            self.frameinfo.append(variables.frameinfo_cont(no))
            self.frameinfo[no].name = '%s/%s_frame%0.4d' \
                %(self.param['direc'],self.param['prefix'],no)
            self.frameinfo[no].omega = omegalist[no];
            self.frameinfo[no].nrefl = 0 # Initialize number of reflections on frame
            self.frameinfo[no].refs = [] # Initialize number of reflections on frame
#            i += 1
#        logging.debug("Printing frameinfo...")
            
        if self.param['theta_max'] == None:
            # Find maximum theta for generation of all possible reflections on
            # the detector from the detector specs
            dety_center_mm = self.param['dety_center'] * self.param['y_size']
            detz_center_mm = self.param['detz_center'] * self.param['z_size']
            dety_size_mm = self.param['dety_size'] * self.param['y_size']
            detz_size_mm = self.param['detz_size'] * self.param['z_size']
            c2c = n.zeros((4))
            c2c[0] = (dety_center_mm-dety_size_mm)**2 + (detz_center_mm-detz_size_mm)**2
            c2c[1] = (dety_center_mm-dety_size_mm)**2 + (detz_center_mm-0)**2
            c2c[2] = (dety_center_mm-0)**2 + (detz_center_mm-detz_size_mm)**2
            c2c[3] = (dety_center_mm-0)**2 + (detz_center_mm-0)**2
            c2c_max = n.max(n.sqrt(c2c))
            theta_max = n.arctan(c2c_max/self.param['distance'])/2.0 * 180./n.pi
            logging.info('To make full detector coverage sets theta_max: %f' %theta_max)
            self.param['theta_max'] = theta_max
			

if __name__=='__main__':

    #import check_input
    try:
        filename = sys.argv[1] 
    except:
        print 'Usage: check_input.py  <input.inp>'
        sys.exit()

    myinput = parse_input(input_file = filename)
    myinput.read()
    print myinput.param
    myinput.check() 
    if myinput.missing == True:
        print 'MISSING ITEMS'
    myinput.evaluate()
    print myinput.param
