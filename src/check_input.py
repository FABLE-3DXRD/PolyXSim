#!/usr/bin/env python

#
# Checking input  
#

from string import split
import sys, os 
import variables
from xfab import tools

import numpy as N
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
                    'unit_cell'  : 'Missing input: unit_cell [unit cell parameters: a,b,c,alpha,beta, gamma]',
                    'no_grains'  : 'Missing input: no_grains [number of grains]',
                    'direc'      : 'Missing input: direc [directory to save output]',
                    'theta_min'   : 'Missing input: theta_min [Minimum theta angle for generation of reflections in degrees]',
                    'theta_max'   : 'Missing input: theta_max [Maximum theta angle for generation of reflections in degrees]',
                                        }
        self.optional_items = {
            'sgno': 1,
            'tilt_x'     : 0,
            'tilt_y'     : 0,
            'tilt_z'     : 0,
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
            'gen_U'   : None,
            'gen_pos' : None,
            'gen_eps' : None,   
            'sample_xyz': None,
            'sample_cyl': None,
            'grain_size': None,
            'grain_min_max': None,
            'direc': '.',
            'prefix': 'test',
            'format' : '.edf'
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

# assert that the correct number of arguments are given
                    if key == 'peakshape':
                        assert len(val) == 2 or len(val) == 3, 'Wrong number of arguments for %s' %key
                    elif key == 'sample_cyl'or key == 'grain_min_max':
                        assert len(val) == 2 , 'Wrong number of arguments for %s' %key
                    elif key == 'sample_xyz' or 'pos_grains' in key:
                        assert len(val) == 3, 'Wrong number of arguments for %s' %key
                    elif key == 'gen_eps':
                        assert len(val) == 4, 'Wrong number of arguments for %s' %key
                    elif key == 'unit_cell' or 'eps_grains' in key:
                        assert len(val) == 6, 'Wrong number of arguments for %s' %key
                    elif 'U_grains' in key:
                        assert len(val) == 9, 'Wrong number of arguments for %s' %key
                    else:
                        assert len(val) == 1, 'Wrong number of arguments for %s' %key

# evaluate and store 
                    valtmp = '['
                    if len(val) > 1:
                        for i in val:
                            valtmp = valtmp + i +','
							
                        val = valtmp + ']'
                    else:
                        val = val[0]

                    self.param[key] = eval(val)
					
                    if 'U_grains' in key:
                        self.param[key] = N.array(self.param[key])
                        self.param[key].shape = (3,3)
           

                
    def check(self):
		self.missing = False

		for item in self.needed_items:
			if item not in self.param:
				print self.needed_items[item]
				self.missing = True
		
		grain_list_U = []
		grain_list_pos = []
		grain_list_eps = []
		grain_list_size = []
		no_grains = self.param['no_grains']

# read U, pos, eps and size for all grains		
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
					
# assert that input U, pos, eps size are correct in format (same number of grains and same specifiers or else not input) 
		grain_list_U.sort()
		grain_list_pos.sort()
		grain_list_eps.sort()
		grain_list_size.sort()
		if len(grain_list_U) != 0 and 'gen_U' not in self.param:
			assert len(grain_list_U) == no_grains, 'Input number of grains does not agree with number of U_grains, check for multiple names'
			self.param['grain_list'] = grain_list_U
			if len(grain_list_pos) != 0 and 'gen_pos' not in self.param:
				assert grain_list_U == grain_list_pos, 'Specified grain number for U_grains and pos_grains disagree'
			if len(grain_list_eps) != 0 and 'gen_eps' not in self.param:
				assert grain_list_U == grain_list_eps, 'Specified grain number for U_grains and eps_grains disagree'
			if len(grain_list_size) != 0 and 'grain_size' not in self.param:
				assert grain_list_U == grain_list_size, 'Specified grain number for U_grains and size_grains disagree'
		else:
			if len(grain_list_pos) != 0 and 'gen_pos' not in self.param:
				assert len(grain_list_pos) == no_grains, 'Input number of grains does not agree with number of pos_grains, check for multiple names'
				self.param['grain_list'] = grain_list_pos
				if len(grain_list_eps) != 0 and 'gen_eps' not in self.param:
					assert grain_list_pos == grain_list_eps, 'Specified grain number for pos_grains and eps_grains disagree'
				if len(grain_list_size) != 0 and 'grain_size' not in self.param:
					assert grain_list_pos == grain_list_size, 'Specified grain number for pos_grains and size_grains disagree'
			elif len(grain_list_eps) != 0 and 'gen_eps' not in self.param:
				assert len(grain_list_eps) == no_grains, 'Input number of grains does not agree with number of eps_grains, check for multiple names'
				self.param['grain_list'] = grain_list_eps
				if len(grain_list_size) != 0 and 'grain_size' not in self.param:
					assert grain_list_eps == grain_list_size, 'Specified grain number for eps_grains and size_grains disagree'
			elif len(grain_list_size) != 0 and 'grain_size' not in self.param:
				assert len(grain_list_size) == no_grains, 'Input number of grains does not agree with number of size_grains, check for multiple names'
				self.param['grain_list'] = grain_list_size
			else:
				self.param['grain_list'] = range(no_grains)

# give default values for generation if no info is read				
		if len(grain_list_U) == 0 and 'gen_U' not in self.param:
			self.param['gen_U'] = 0
		if len(grain_list_pos) == 0 and 'gen_pos' not in self.param:
			self.param['gen_pos'] = 0
		if len(grain_list_eps) == 0 and 'gen_eps' not in self.param:
			self.param['gen_eps'] = [0,0,0,0]
		if len(grain_list_size) == 0 and 'grain_size' not in self.param:
			self.param['grain_size'] = -0.075
		if 'grain_size' in self.param and 'grain_min_max' not in self.param:
			self.param['grain_min_max'] = [0,10*abs(self.param['grain_size'])]			
			
#assert that not both sample_xyz and sample_cyl are given
		if 'sample_xyz' in self.param:
			assert 'sample_cyl' not in self.param, 'Both sample_xyz and sample_cyl are given'
			
# assert that mean grain size != 0 and if mean > 0 then min < mean < max, assure that min non-negative
		if 'grain_size' in self.param:
			assert self.param['grain_size'] != 0, 'grain_size 0 is an invalid command'
			if self.param['grain_size'] > 0:
				assert self.param['grain_min_max'][0] < self.param['grain_size'], 'grain_min larger than grain_size'
				assert self.param['grain_min_max'][1] > self.param['grain_size'], 'grain_max smaller than grain_size'
				if self.param['grain_min_max'][0] < 0:
					self.param['grain_min_max'][0] = 0
		
#check that the given grain_size and no_grains are consistent with sample_vol, adjust max to sample size
			if 'sample_xyz' in self.param:
				self.param['sample_vol'] = self.param['sample_xyz'][0]*self.param['sample_xyz'][1]*self.param['sample_xyz'][2]
				diam_limit = (6*self.param['sample_vol']/(N.exp(.5)*N.pi*self.param['no_grains']))**(1/3.)
				assert abs(self.param['grain_size']) < diam_limit, 'The sample volume is too small to contain the specified number of grains with the given grain size'
				self.param['grain_min_max'][1] = min(self.param['sample_xyz'][0],self.param['sample_xyz'][1],self.param['sample_xyz'][2])
			elif 'sample_cyl' in self.param:
				self.param['sample_vol'] = N.pi*self.param['sample_cyl'][0]*self.param['sample_cyl'][0]*self.param['sample_cyl'][1]
				diam_limit = (6*self.param['sample_vol']/(N.exp(.5)*N.pi*self.param['no_grains']))**(1/3.)
				assert abs(self.param['grain_size']) < diam_limit, 'The sample volume is too small to contain the specified number of grains with the given grain size'
				self.param['grain_min_max'][1] = min(2*self.param['sample_cyl'][0],self.param['sample_cyl'][1])
			else:					
				self.param['sample_vol'] = None

			
    def initialize(self): 
		# set all non-read items to defaults
        for item in self.optional_items:
            if (item not in self.param):
                self.param[item] = self.optional_items[item]
            if (self.param[item] == []):
                self.param[item] = self.optional_items[item]*self.param['no_grains']
 
        # Does output directory exist?
        if not os.path.exists(self.param['direc']):
            os.mkdir(self.param['direc'])

			# Generate FILENAME of frames
        omega_step = self.param['omega_step']
        omega_start  = self.param['omega_start']
        omega_end  = self.param['omega_end']
        omega_sign = self.param['omega_sign']
        start_frame = self.param['start_frame']
        omegalist = omega_sign*N.arange(omega_start,omega_end+omega_step,omega_step)
        nframes = int((omega_end-omega_start)/omega_step)
        omegalist.sort()
#        i=0
#        logging.info("Generating frame data...")
        #Initialize frameinfo container
        self.frameinfo = [] 
        
        if omega_sign > 0:
            filerange = N.arange(start_frame,start_frame+nframes)
        else:
            filerange = N.arange((start_frame-1)+nframes,(start_frame-1),omega_sign)
            # reverse omega_start/omega_end
            self.param['omega_end'] = omega_start*omega_sign 
            self.param['omega_start'] = omega_end*omega_sign

        for no in filerange:
            self.frameinfo.append(variables.frameinfo_cont(no))
            self.frameinfo[no].name = '%s/%s_frame%0.4d%s' %(self.param['direc'],self.param['prefix'],no,self.param['format'])
            self.frameinfo[no].omega = omegalist[no];
            self.frameinfo[no].nrefl = 0 # Initialize number of reflections on frame
            self.frameinfo[no].refs = [] # Initialize number of reflections on frame
#            i += 1
#        logging.debug("Printing frameinfo...")
            
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
    print myinput.param
    myinput.check() 
    if myinput.missing == True:
        print 'MISSING ITEMS'
    myinput.evaluate()
    print myinput.param
