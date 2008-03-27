#!/usr/bin/env python

#
# Checking input  
#

from string import split
import sys, os 
from fabio import deconstruct_filename, jump_filename
import variables
from xfab import tools

import numpy as N
import logging

logging.basicConfig(level=logging.DEBUG,format='%(levelname)s %(message)s')

class parse_input:
    def __init__(self,input_file = None):
        self.filename = input_file
        self.entries = {}
        self.grainno = 0
#        self.entries['pos_grains'] = []
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
#            'pos_grains' : [[0, 0, 0]],
            'tilt_x'     : 0,
            'tilt_y'     : 0,
            'tilt_z'     : 0,
            'beampol_factor' : 1,
            'beampol_angle' : 0.0,
            'spatial' : None,
            'flood' : None,
            'dark' : None,
            'darkoffset' : None,
#			'gen_U'   : 0,
#			'gen_pos' : 0,
#			'gen_eps' : [0,0],	
#			'sample_xyz': [0,0,0],
#			'sample_cyl': [0,0]
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
                    self.entries[key] = eval(val)
                

                
    def check(self):
		self.missing = False

		for item in self.needed_items:
			if item not in self.entries:
				print self.needed_items[item]
				self.missing = True
		
		grain_list_U = []
		grain_list_pos = []
		grain_list_eps = []
		no_grains = self.entries['no_grains']

# read U, pos and eps for all grains		
		for item in self.entries:
			if '_grains_' in item:
				if 'U' in item:
					grain_list_U.append(eval(split(item,'_grains_')[1]))
				elif 'pos' in item:
					grain_list_pos.append(eval(split(item,'_grains_')[1]))
				elif 'eps' in item:
					grain_list_eps.append(eval(split(item,'_grains_')[1]))
					
# assert that input U, pos and eps are correct in format (same number of grains and same specifiers or else not input) 
		grain_list_U.sort()
		grain_list_pos.sort()
		grain_list_eps.sort()
		if len(grain_list_U) != 0 and 'gen_U' not in self.entries:
			assert len(grain_list_U) == no_grains, 'Input number of grains does not agree with number of U_grains, check for multiple names'
			self.entries['grain_list'] = grain_list_U
			if len(grain_list_pos) != 0 and 'gen_pos' not in self.entries:
				assert grain_list_U == grain_list_pos, 'Specified grain number for U_grains and pos_grains disagree'
			if len(grain_list_eps) != 0 and 'gen_eps' not in self.entries:
				assert grain_list_U == grain_list_eps, 'Specified grain number for U_grains and eps_grains disagree'
		else:
			if len(grain_list_pos) != 0 and 'gen_pos' not in self.entries:
				assert len(grain_list_pos) == no_grains, 'Input number of grains does not agree with number of pos_grains, check for multiple names'
				self.entries['grain_list'] = grain_list_pos
				if len(grain_list_eps) != 0 and 'gen_eps' not in self.entries:
					assert grain_list_pos == grain_list_eps, 'Specified grain number for pos_grains and eps_grains disagree'
			elif len(grain_list_eps) != 0 and 'gen_eps' not in self.entries:
				assert len(grain_list_eps) == no_grains, 'Input number of grains does not agree with number of eps_grains, check for multiple names'
				self.entries['grain_list'] = grain_list_eps
			else:
				self.entries['grain_list'] = range(no_grains)
							
# Generate U if gen_U exists or the values are not input
		if len(grain_list_U) == 0 or 'gen_U' in self.entries: 
			print 'Grain orientations will be randomly generated\n'
			for i in self.entries['grain_list']:
				phi1 = N.random.rand()*2*N.pi
				phi2 = N.random.rand()*2*N.pi
				PHI = N.random.rand()*N.pi
				self.entries['U_grains_%s' %(i)] = tools.euler2U(phi1,PHI,phi2)
				
# Generate pos if gen_pos exists
# If gen_pos != 0 use the sample shape and size if this is input, else default to (0,0,0)
# Finally, if the positions are not input or gen_pos == 0, default to (0,0,0)	
		if 'gen_pos' in self.entries:
			if self.entries['gen_pos'] != 0:
				if 'sample_xyz' in self.entries:
					assert 'sample_cyl' not in self.entries, 'sample_xyz and sample_cyl are input simultaneously'
					sample_x = self.entries['sample_xyz'][0]
					sample_y = self.entries['sample_xyz'][1]
					sample_z = self.entries['sample_xyz'][2]
					print 'Grain positions will be randomly generated within a box of ', sample_x, sample_y, sample_z, 'mm\n'
					for i in self.entries['grain_list']:
						x = (N.random.rand()-0.5)*sample_x
						y = (N.random.rand()-0.5)*sample_y
						z = (N.random.rand()-0.5)*sample_z
						self.entries['pos_grains_%s' %(i)] = [x, y, z]

				elif 'sample_cyl'in self.entries:
					sample_r = self.entries['sample_cyl'][0]
					sample_z = self.entries['sample_cyl'][1]
					print 'Grain positions will be randomly generated within a cylinder of radius ', sample_r, ' and length ', sample_z, 'mm\n'
					for i in self.entries['grain_list']:
						r = N.random.rand()*sample_r
						w = N.random.rand()*2*N.pi
						z = (N.random.rand()-0.5)*sample_z
						self.entries['pos_grains_%s' %(i)] = [r*N.cos(w), r*N.sin(w), z]

				else:
					print 'Grain positions will be set to (0,0,0)\n'
					for i in self.entries['grain_list']:
						self.entries['pos_grains_%s' %(i)] = [0, 0, 0]

			else:
				print 'Grain positions will be set to (0,0,0)\n'
				for i in self.entries['grain_list']:
					self.entries['pos_grains_%s' %(i)] = [0, 0, 0]
					
		elif len(grain_list_pos) == 0:
			print 'Grain positions will be set to (0,0,0)\n'
			for i in self.entries['grain_list']:
				self.entries['pos_grains_%s' %(i)] = [0, 0, 0]

# Generate strain tensors if gen_eps exists
# Use a normal distribution with the specified mean and spread for diagonal and off-diagonal elements
# Finally, if the strain tensors are not input default to (0,0,0,0,0,0)	
		if 'gen_eps' in self.entries: 
			mean_diag = self.entries['gen_eps'][0]
			spread_diag = self.entries['gen_eps'][1]
			mean_offdiag = self.entries['gen_eps'][2]
			spread_offdiag = self.entries['gen_eps'][3]
			print 'Grain strain tensors will be randomly generated using a normal distribution'
			print 'For diagonal elements, mean: ', mean_diag, ' and spread: ', spread_diag
			print 'For off-diagonal elements, mean: ', mean_offdiag, ' and spread: ', spread_offdiag
			
			for i in self.entries['grain_list']:
				eps = tools.geneps(mean_diag,spread_diag,mean_offdiag,spread_offdiag)
				self.entries['eps_grains_%s' %(i)] = eps
				
		elif len(grain_list_eps) == 0: 
			print 'Grain strain tensors will be set to (0,0,0,0,0,0) - no strain'
			for i in self.entries['grain_list']:
				self.entries['eps_grains_%s' %(i)] = [0, 0, 0, 0, 0, 0]
	
	
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
