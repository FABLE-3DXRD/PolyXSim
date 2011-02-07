#!/usr/bin/env python

#
# Checking input  
#

from string import split
from copy import copy
import sys, os 
import variables
from xfab import tools,sg

import numpy as n
#import logging
#logging.basicConfig(level=logging.DEBUG,format='%(levelname)s %(message)s')

def interrupt(killfile):
    if killfile is not None and os.path.exists(killfile):
        raise KeyboardInterrupt


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
                    'omega_start': 'Missing input: omega_start [Omega start in degrees]',
                    'omega_end'  : 'Missing input: omega_end [Omega end in degrees]',
                    'omega_step' : 'Missing input: omega_step [Omega step size in degrees]',
                    'no_grains'  : 'Missing input: no_grains [number of grains]',
                    'direc'      : 'Missing input: direc [directory to save output]',
                                        }
        self.optional_items = {
            'sgno'  : None,
            'sgname': None,
            'tilt_x': 0,
            'tilt_y': 0,
            'tilt_z': 0,
            'wedge': 0.0,
            'beam_width': None,
            'beampol_apply' : 1,
            'beampol_factor': 1,
            'beampol_direct' : 0.0,
            'lorentz_apply' : 1,
            'start_frame': 0,
            'omega_sign': 1,
            'noise' : 0,
            'psf': 0,
            'make_image': 1,
            'o11': 1,
            'o12': 0,
            'o21': 0,
            'o22': 1,
            'bg' : 0,
            'peakshape': [0,0],
            'spatial' : None,
            'flood' : None,
            'dark' : None,
            'darkoffset' : None,
            'gen_phase': [0],
            'gen_U' : 0,
            'gen_pos' : [0,0],
            'gen_eps' : [0,0,0,0,0],   
            'gen_size' : [0,0,0,0],
            'sample_xyz': None,
            'sample_cyl': None,
            'direc' : '.',
            'stem': 'test',
            'odf_type' : 1,
            'odf_scale' : 0.02,
            'odf_cut' : None,
            'odf_sub_sample': 1,
            'mosaicity' : 0.2,
            'theta_min' : 0.0,
            'theta_max' : None,
            'unit_cell' : None,
            'structure_file': None,
            'odf_file': None, 
            'output': None,
            'structure_factors': 1,
            'no_phases': 1,
            'intensity_const': 0
            }
        self.output_types = ['.edf', 
                             '.tif', 
                             '.ref',
                             '.par', 
                             '.flt',
                             '.gve',
                             '.ini',
                             '.ubi']
        
    def read(self):     
        try:
            f = open(self.filename,'r')
        except IOError:
            raise IOError, 'No file named %s' %self.filename
        
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
                    # This ensures that file names with space can be handled 
                    if key == 'direc' or key == 'stem' or key == 'spatial' or 'structure' in key:
                        valtmp = ''
                        valend = ''
                        sepa = ' '
                    else:   
                        valtmp = '['
                        valend = ']'
                        sepa = ','
                    
                    if len(val) > 1:
                        for i in val:
                            valtmp = valtmp + i + sepa
                        # remove last separator
                        val = valtmp[:-len(sepa)] + valend
                    else:
                        val = val[0]
                        
                    # Problems using the Windows path separator as they can 
                    # cause part of a string to be interpreted as a special charactor
                    # like \n - newline.
                    if key == 'direc' or key == 'stem' or key == 'spatial' or ('structure' in key and 'factors' not in key):
                        # Hack to remove a final backslash in the directory path
                        if str(val).endswith('\\\''):
                            val = val[:-2]+"'"
                        # Unfortunately the escape character can be changed like this
                        val = val.replace('\\x','/x')
                        # before taking care of the rest
                        val = eval(val)
                        val = val.replace('\t','\\t')
                        val = val.replace('\n','\\n')
                        val = val.replace('\c','\\c')
                        val = val.replace('\b','\\b')
                        val = val.replace('\f','\\f')
                        val = val.replace('\r','\\r')
                        val = val.replace('\v','\\v')
                        # Added value to key
                        self.param[key] = val                        
                    else:
                        self.param[key] = eval(val)
           
    def check(self):
        self.missing = False
        self.errors = {}

# check that all needed items are present
        for item in self.needed_items:
            if item not in self.param:
                #print self.needed_items[item]
                #self.missing = True
                self.errors[item] =  self.needed_items[item]
		
# set all non-read items to defaults
        for item in self.optional_items:
            if (item not in self.param):
                self.param[item] = self.optional_items[item]
            if (self.param[item] == []):
                self.param[item] = self.optional_items[item]*self.param['no_grains']

# assert that the correct number of arguments are given
        for key in self.param:
            val = self.param[key] 
            if val != None and key != 'output' and key != 'grain_list':
                if key == 'peakshape':
                    if type(val) == type(1):
                        #assert val == 0, 'Wrong number of arguments for %s' %key
                        if val != 0: 
                            self.errors[key] = 'Wrong number of arguments' 
                        else:
                            val = [0, 0]
                            self.param[key] = val
                    else:
                        #assert len(val) <= 3, 'Wrong number of arguments for %s' %key
                        if len(val) > 3: 
                            self.errors[key] = 'Wrong number of arguments' 
                elif key == 'sample_cyl' or key == 'gen_pos':
                    #assert len(val) == 2 , 'Wrong number of arguments for %s' %key
                    if len(val) != 2: 
                        self.errors[key] = 'Wrong number of arguments' 
                elif key == 'sample_xyz' or 'pos_grains' in key:
                    #assert len(val) == 3, 'Wrong number of arguments for %s' %key
                    if len(val) != 3: 
                        self.errors[key] = 'Wrong number of arguments' 
                elif 'gen_size' in key:
                    #assert len(val) == 4, 'Wrong number of arguments for %s' %key
                    if len(val) != 4: 
                        self.errors[key] = 'Wrong number of arguments' 
                elif 'gen_eps' in key:
                    if type(val) == type(1):
                        #assert val == 0, 'Wrong number of arguments for %s' %key
                        if val != 0: 
                            self.errors[key] = 'Wrong number of arguments' 
                    else:
                        #assert len(val) == 5, 'Wrong number of arguments for %s' %key
                        if len(val) != 5: 
                            self.errors[key] = 'Wrong number of arguments' 
                elif key == 'gen_phase':
                    try: 
                        dummy = len(val)
                    except:
                        val = [val]
                        self.param[key] = val
                    #assert len(val) > 0, 'Wrong number of arguments for %s' %key
                    if len(val) < 1: 
                        self.errors[key] = 'Wrong number of arguments' 
                elif key == 'unit_cell' or 'eps_grains' in key:
                    #assert len(val) == 6, 'Wrong number of arguments for %s' %key
                    if len(val) != 6: 
                        self.errors[key] = 'Wrong number of arguments' 
                elif 'U_grains' in key:
                    if len(val) != 3:
                        #assert len(val) == 9, 'Wrong number of arguments for %s' %key
                        if len(val) != 9: 
                            self.errors[key] = 'Wrong number of arguments' 
                        else:
                            self.param[key] = n.array(self.param[key])
                            self.param[key].shape = (3,3)
                    else:
                        #assert val.shape == (3,3), 'Wrong number of arguments for %s' %key
                        if  val.shape != (3,3): 
                            self.errors[key] = 'Wrong number of arguments' 
                    # reshape U-matrices
            elif key == 'output':
                for item in val:
                    if item not in self.output_types:
                        self.errors[key] = 'Output type given %s is not an option' %item 

# Check no of phases
        no_phases = self.param['no_phases']

        phase_list_structure = []
        phase_list_unit_cell = []
        phase_list_sgno = []
        phase_list_sgname = []
        phase_list_gen_size = []
        phase_list_gen_eps = []
        phase_list = []

        for item in self.param:
            if '_phase_' in item:
                if 'structure' in item:
                    phase_list_structure.append(eval(split(item,'_phase_')[1]))
                elif 'unit_cell' in item:
                    phase_list_unit_cell.append(eval(split(item,'_phase_')[1]))
                elif 'sgno' in item:
                    phase_list_sgno.append(eval(split(item,'_phase_')[1]))
                elif 'sgname' in item:
                    phase_list_sgname.append(eval(split(item,'_phase_')[1]))
                elif 'gen_size' in item:
                    phase_list_gen_size.append(eval(split(item,'_phase_')[1]))
                elif 'gen_eps' in item:
                    phase_list_gen_eps.append(eval(split(item,'_phase_')[1]))
                    

        phase_list_structure.sort()
        phase_list_unit_cell.sort()
        phase_list_sgno.sort()
        phase_list_sgname.sort()
        phase_list_gen_size.sort()
        phase_list_gen_eps.sort()
                       
        if len(phase_list_structure) != 0:
#             assert len(phase_list_structure) == no_phases, \
#                 'Input number of structural phases does not agree with number\n' +\
#                 ' of structure_phase, check for multiple names or missing files.'
#             phase_list = phase_list_structure
            if len(phase_list_structure) != no_phases:
                self.errors['phase_list_structure'] = \
                    'Input number of structural phases does not agree with number\n' +\
                    ' of structure_phase, check for multiple names or missing files.'
            else:
                phase_list = phase_list_structure
        elif  len(phase_list_unit_cell) != 0:
#             assert  len(phase_list_unit_cell) == no_phases, \
#                 'Input number of structural phases does not agree with number\n' +\
#                 ' of unit_cell, check for multiple names or missing linies.'
#             phase_list = phase_list_unit_cell
            if len(phase_list_unit_cell) != no_phases:
                self.errors['phase_list_unit_cell'] = \
                    'Input number of structural phases does not agree with number\n' +\
                    ' of unit_cell, check for multiple names or missing linies.'
            else:
                phase_list = phase_list_unit_cell
            if len(phase_list_sgno) == 0:
#                 assert len(phase_list_sgname) == no_phases, \
#                     'Input number of structural phases does not agree with number\n' +\
#                     'of space group information given (sgno or sgname),\n' +\
#                     'check for multiple names or missing linies.'
                if len(phase_list_sgname) != no_phases:
                    self.errors['phase_list_sgname'] = \
                        'Input number of structural phases does not agree with number\n' +\
                        'of space group information given (sgno or sgname),\n' +\
                        'check for multiple names or missing linies.'
#                 assert phase_list_sgname == phase_list_unit_cell, \
#                     'The phase numbers given to unit_cell does not match those in sgname'
                if phase_list_sgname != phase_list_unit_cell:
                    self.errors['phase_list_sgname_2'] = \
                        'The phase numbers given to unit_cell does not match those in sgname'

                # add sgno for phase to parameter list
                for phase in phase_list_sgname:
                    self.param['sgno_phase_%i' %phase] = sg.sg(sgname = self.param['sgname_phase_%i' %phase]).no

            elif len(phase_list_sgname) == 0:
#                 assert len(phase_list_sgno) == no_phases, \
#                     'Input number of structural phases does not agree with number\n' +\
#                     'of space group information given (sgno or sgname),\n' +\
#                     'check for multiple names or missing linies.'
                if len(phase_list_sgno) != no_phases:
                    self.errors['phase_list_sgno'] = \
                        'Input number of structural phases does not agree with number\n' +\
                        'of space group information given (sgno or sgname),\n' +\
                        'check for multiple names or missing linies.'
#                 assert phase_list_sgno == phase_list_unit_cell, \
#                     'The phase numbers given to unit_cell does not match those in sgno.'
                if phase_list_sgno != phase_list_unit_cell:
                    self.errors['phase_list_sgno_2'] = \
                        'The phase numbers given to unit_cell does not match those in sgno.'

                # add sgname for phase to parameter list
                for phase in phase_list_sgno:
                    self.param['sgname_phase_%i' %phase] = sg.sg(sgno = self.param['sgno_phase_%i' %phase]).name
            else:
                # both sg numbers and names in input check if they point at the same space group
                for phase in phase_list_sgno:
#                     assert self.param['sgname_phase_%i' %phase] == sg.sg(sgno = self.param['sgno_phase_%i' %phase]).name, \
#                         '\nSpace group is specified both as space group name and number - \n' + \
#                         'and they do not correspond to the same space group. \n' + \
#                         'Please sort this out in the input file.' 
                    if self.param['sgname_phase_%i' %phase] != sg.sg(sgno = self.param['sgno_phase_%i' %phase]).name:
                        self.errors['sgname_phase_list_sgno_2'] = \
                            '\nSpace group is specified both as space group name and number - \n' + \
                            'and they do not correspond to the same space group. \n' + \
                            'Please sort this out in the input file for phase %i.' %phase 
                
                
        if len(phase_list_gen_size) != 0:
#             assert len(phase_list_gen_size) == no_phases, \
#                 'Input number of structural phases does not agree with number\n' +\
#                 'of gen_size_phase_ keywords given'
            if len(phase_list_gen_size) != no_phases:
                self.errors['phase_list_gen_size_1'] = \
                    'Input number of structural phases does not agree with number\n' +\
                    'of gen_size_phase_ keywords given'
#            assert phase_list_gen_size == phase_list, \
#                    'The phase numbers given to gen_size_phase does not match those\n' +\
#                    'in crystallographic part - structure_phase_X or unit_cell_phase_X.'
            
            if phase_list_gen_size != phase_list:
                self.errors['phase_list_gen_size_2'] = \
                    'The phase numbers given to gen_size_phase does not match those\n' +\
                    'in crystallographic part - structure_phase_X or unit_cell_phase_X.'
            self.param['gen_size'][0] = 1
        else:
            if len(phase_list) > 0:
                for phase in phase_list:
                    self.param['gen_size_phase_%i' %phase] = copy(self.param['gen_size'])
            else:
                phase = 0 
                self.param['gen_size_phase_%i' %phase] = copy(self.param['gen_size'])
                
        if len(phase_list_gen_eps) != 0:
#             assert len(phase_list_gen_eps) == no_phases, \
#                 'Input number of structural phases does not agree with number\n' +\
#                 'of gen_size_phase_ keywords given'
            if len(phase_list_gen_eps) != no_phases:
                self.errors['phase_list_gen_eps_1'] = \
                    'Input number of structural phases does not agree with number\n' +\
                    'of gen_size_phase_ keywords given'
#             assert phase_list_gen_eps == phase_list, \
#                     'The phase numbers given to gen_size_phase does not match those\n' +\
#                     'in crystallographic part - structure_phase_X or unit_cell_phase_X.'
            if phase_list_gen_eps != phase_list:
                self.errors['phase_list_gen_eps_2'] = \
                    'The phase numbers given to gen_size_phase does not match those\n' +\
                    'in crystallographic part - structure_phase_X or unit_cell_phase_X.'
            self.param['gen_eps'][0] = 1
        else:
            if len(phase_list) > 0:
                for phase in phase_list:
                    self.param['gen_eps_phase_%i' %phase] = copy(self.param['gen_eps'])
            else:
                phase = 0 
                # If strain is not provided and no generation of strains have be asked for
                # "set" all eps to zero. 
                if self.param['gen_eps'][0] == 0:
                    self.param['gen_eps'] = [1, 0.0, 0.0, 0.0, 0.0]
                self.param['gen_eps_phase_%i' %phase] = copy(self.param['gen_eps'])
        if self.param['gen_phase'][0] != 0:
#             assert len(self.param['gen_phase'][1:]) == no_phases*2, 'Missing info for  -  gen_phase'
            if len(self.param['gen_phase'][1:]) != no_phases*2:
                self.errors['gen_phase'] = 'Missing info for  -  gen_phase'

        # Make sure both sgname and sgno exist
        if len(phase_list_sgno) == 0:
            for phase in phase_list_sgno:
                    self.param['sgno_phase_%i' %phase] = sg.sg(sgname = self.param['sgname_phase_%i' %phase]).no
        if len(phase_list_sgname) == 0:
            for phase in phase_list_sgname:
                    self.param['sgname_phase_%i' %phase] = sg.sg(sgno = self.param['sgno_phase_%i' %phase]).name


# Init no of grains belonging to phase X if not generated
        if self.param['gen_phase'][0] != 1:
            if len(phase_list) == 0:
                self.param['no_grains_phase_0'] = self.param['no_grains']
            else:
                for phase in phase_list:
                    self.param['no_grains_phase_%i' %phase] = 0
        else:
            for i in range(self.param['no_phases']):
                phase = self.param['gen_phase'][i*2+1]
                no_grains_phase = int(self.param['gen_phase'][i*2+2])
                self.param['no_grains_phase_%i' %phase] = no_grains_phase
            
# read U, pos, eps and size for all grains		
        grain_list_U = []
        grain_list_pos = []
        grain_list_eps = []
        grain_list_size = []
        grain_list_phase = []
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
                elif 'phase' in item[:5]:
                    grain_list_phase.append(eval(split(item,'_grains_')[1]))
                    self.param['no_grains_phase_%i' %self.param[item]] += 1

#assert that the number of grains in all match 

        sum_of_grains = 0
        #print 'no of phases', self.param['no_phases']
        if self.param['no_phases'] > 1:
            for phase in phase_list:
                sum_of_grains += self.param['no_grains_phase_%i' %phase]
#             assert sum_of_grains == no_grains, \
#                 'Input number of grains (%i) does not agree ' %no_grains +\
#                 'with number of phase_grains_ keywords (%i)' %sum_of_grains
            if sum_of_grains != no_grains:
                self.errors['no_grains_2'] = \
                    'Input number of grains (%i) does not agree ' %no_grains +\
                    'with number of phase_grains_ keywords (%i)' %sum_of_grains
        else:
            self.param['no_grains_phase_0'] = no_grains

# assert that input U, pos, eps size are correct in format
# (same number of grains and same specifiers or else not input) 
        grain_list_U.sort()
        grain_list_pos.sort()
        grain_list_eps.sort()
        grain_list_size.sort()
        grain_list_phase.sort()
        if len(grain_list_U) != 0 and self.param['gen_U'] == 0:
#             assert len(grain_list_U) == no_grains, \
#                 'Input number of grains does not agree with number\n' +\
#                 ' of U_grains, check for multiple names'
            if len(grain_list_U) != no_grains:
                self.errors['grain_list_U'] = \
                    'Input number of grains does not agree with number\n' +\
                    ' of U_grains, check for multiple names'
            self.param['grain_list'] = grain_list_U
            if len(grain_list_pos) != 0 and self.param['gen_pos'][0] == 0:
#                 assert grain_list_U == grain_list_pos, \
#                     'Specified grain numbers for U_grains and pos_grains disagree'
                if grain_list_U != grain_list_pos:
                    self.errors['grain_list_U_pos'] = \
                        'Specified grain numbers for U_grains and pos_grains disagree'
                if len(grain_list_eps) != 0 and self.param['gen_eps'][0] == 0:
#                     assert grain_list_U == grain_list_eps, \
#                         'Specified grain numbers for U_grains and eps_grains disagree'
                    if grain_list_U != grain_list_eps:
                        self.errors['grain_list_U_eps'] = \
                            'Specified grain numbers for U_grains and eps_grains disagree'
                if len(grain_list_size) != 0 and self.param['gen_size'][0] == 0:
#                     assert grain_list_U == grain_list_size, \
#                         'Specified grain numbers for U_grains and size_grains disagree'
                    if grain_list_U != grain_list_size:
                        self.errors['grain_list_U_size'] = \
                            'Specified grain numbers for U_grains and size_grains disagree'
                if len(grain_list_phase) != 0 and self.param['gen_phase'][0] == 0:
#                     assert grain_list_U == grain_list_phase, \
#                         'Specified grain numbers for U_grains and phase_grains disagree'
                    if grain_list_U != grain_list_phase:
                        self.errors['grain_list_U_phase'] = \
                            'Specified grain numbers for U_grains and phase_grains disagree'
        else:
            if len(grain_list_pos) != 0 and self.param['gen_pos'][0] == 0:
#                 assert len(grain_list_pos) == no_grains, \
#                     'Input number of grains does not agree with number\n'+\
#                     ' of pos_grains, check for multiple names'
                if len(grain_list_pos) != no_grains:
                    self.errors['grain_list_pos_nograins'] = \
                        'Input number of grains does not agree with number\n'+\
                        ' of pos_grains, check for multiple names'
                self.param['grain_list'] = grain_list_pos
                if len(grain_list_eps) != 0 and self.param['gen_eps'][0] == 0:
#                     assert grain_list_pos == grain_list_eps, \
#                         'Specified grain number for pos_grains and eps_grains disagree'
                    if grain_list_pos != grain_list_eps:
                        self.errors['grain_list_pos_eps'] = \
                            'Specified grain number for pos_grains and eps_grains disagree'
                if len(grain_list_size) != 0 and self.param['gen_size'][0] == 0:
#                     assert grain_list_pos == grain_list_size, \
#                         'Specified grain number for pos_grains and size_grains disagree'
                    if grain_list_pos != grain_list_size:
                        self.errors['grain_list_pos_size'] = \
                            'Specified grain number for pos_grains and size_grains disagree'
                if len(grain_list_phase) != 0 and self.param['gen_phase'][0] == 0:
#                     assert grain_list_pos == grain_list_phase, \
#                         'Specified grain numbers for pos_grains and phase_grains disagree'
                    if grain_list_pos != grain_list_phase:
                        self.errors['grain_list_pos_phase'] = \
                            'Specified grain number for pos_grains and phase_grains disagree'
            elif len(grain_list_eps) != 0 and self.param['gen_eps'][0] == 0:
#                 assert len(grain_list_eps) == no_grains, \
#                     'Input number of grains does not agree with number'+\
#                     ' of eps_grains, check for multiple names'
                if len(grain_list_eps) != no_grains:
                    self.errors['grain_list_eps_nograins'] = \
                        'Input number of grains does not agree with number'+\
                        ' of eps_grains, check for multiple names'
                self.param['grain_list'] = grain_list_eps
                if len(grain_list_size) != 0 and self.param['gen_size'][0] == 0:
#                     assert grain_list_eps == grain_list_size, \
#                         'Specified grain number for eps_grains and size_grains disagree'
                    if grain_list_eps != grain_list_size:
                        self.errors['grain_list_eps_size'] = \
                            'Specified grain number for eps_grains and size_grains disagree'
                if len(grain_list_phase) != 0 and self.param['gen_phase'][0] == 0:
#                     assert grain_list_eps == grain_list_phase, \
#                         'Specified grain numbers for eps_grains and phase_grains disagree'
                    if grain_list_eps != grain_list_phase:
                        self.errors['grain_list_eps_phase'] = \
                            'Specified grain number for eps_grains and phase_grains disagree'
            elif len(grain_list_size) != 0 and self.param['gen_size'][0] == 0:
#                 assert len(grain_list_size) == no_grains, \
#                     'Input number of grains does not agree with number\n'+\
#                     ' of size_grains, check for multiple names'
                if len(grain_list_size) != no_grains:
                    self.errors['grain_list_size_nograins'] = \
                        'Input number of grains does not agree with number\n'+\
                        ' of size_grains, check for multiple names'
                self.param['grain_list'] = grain_list_size
                if len(grain_list_phase) != 0 and self.param['gen_phase'][0] == 0:
#                     assert grain_list_size == grain_list_phase, \
#                         'Specified grain numbers for size_grains and phase_grains disagree'
                    if grain_list_size != grain_list_phase:
                        self.errors['grain_list_size_phase'] = \
                            'Specified grain numbers for size_grains and' +\
                            'phase_grains disagree'
            elif len(grain_list_phase) != 0 and self.param['gen_phase'][0] == 0:
#                 assert len(grain_list_phase) == no_grains, \
#                     'Input number of grains does not agree with number\n'+\
#                     ' of phase_grains, check for multiple names'
                if len(grain_list_phase) != no_grains:
                    self.errors['grain_list_phase_nograins'] = \
                        'Input number of grains does not agree with number\n'+\
                        ' of phase_grains, check for multiple names'
                self.param['grain_list'] = grain_list_phase
            else:
                self.param['grain_list'] = range(no_grains)

# assert that all information needed to generate grains is present	
#         assert len(grain_list_U) != 0 or self.param['gen_U'] != 0,\
#             'Information on U generations missing'
#         assert len(grain_list_pos) != 0 or self.param['gen_pos'][0] != 0,\
#             'Information on position generation missing'
#         assert len(grain_list_eps) != 0 or self.param['gen_eps'][0] != 0,\
#             'Information on strain tensor generation missing'
#         assert len(grain_list_size) != 0 or self.param['gen_size'][0] != 0,\
#             'Information on grain size generation missing'
        if len(grain_list_U) == 0 and self.param['gen_U'] == 0:
            self.errors['grain_list_gen_U'] = \
                'Information on U generations missing'
        if len(grain_list_pos) == 0 and self.param['gen_pos'][0] == 0:
            self.errors['grain_list_gen_pos'] = \
                'Information on position generation missing'
# Not used anymore as zero strains are given if not strains are provide 
# and no gen_eps's are given. Might be used to write warnings later.
#        if len(grain_list_eps) == 0 and self.param['gen_eps'][0] == 0:
#            self.errors['grain_list_gen_eps'] = \
#                'Information on strain tensor generation missing'
        if len(grain_list_size) == 0 and self.param['gen_size'][0] == 0:
            self.errors['grain_list_gen_size'] = \
                'Information on grain size generation missing'
			
#If no structure file is given - unit_cel should be               
        if len(phase_list) == 0:
            # This is a monophase simulation probably using the "old" keywords
            if self.param['structure_file'] == None:
                #print 'NO structure file'
#                 assert self.param['unit_cell'] != None, \
#                     'Missing input: structure_file or unit_cell' 
                if self.param['unit_cell'] == None:
                    self.errors['unit_cell'] = \
                        'Missing input: either structure_file or unit_cell has to be specified' 
                
                # rename keyword
                self.param['unit_cell_phase_0'] = self.param['unit_cell']
                # and delete old one
                del self.param['unit_cell']
#                 assert self.param['sgno'] != None or self.param['sgname'] != None , \
#                     'Missing input: no space group information, please input either sgno or sgname' 
                if self.param['sgno'] == None and self.param['sgname'] == None:
                    self.errors['sgno'] = \
                        'Missing input: no space group information, '+\
                        'please input either sgno or sgname' 
                if self.param['sgno'] == None:
                    self.param['sgno_phase_0'] = sg.sg(sgname = self.param['sgname']).no
                    # rename keyword
                    self.param['sgname_phase_0'] = self.param['sgname']
                    # and delete old one
                    del self.param['sgname']
                else:
                    self.param['sgname_phase_0'] = sg.sg(sgno = self.param['sgno']).name
                    # rename keyword
                    self.param['sgno_phase_0'] = self.param['sgno']
                    # and delete old one
                    del self.param['sgno']
            else:
                # rename keyword
                self.param['structure_phase_0'] = self.param['structure_file']
                # and delete old one
                del self.param['structure_file']
            phase_list = [0]
        self.param['phase_list'] = phase_list

        # make old inp files work
        if len(grain_list_phase) == 0 and self.param['no_phases'] == 1:
            self.param['grain_list_phase_%i' %self.param['phase_list'][0]] = self.param['grain_list']



#assert that not both sample_xyz and sample_cyl are given
#         if self.param['sample_xyz'] != None:
#             assert self.param['sample_cyl'] == None, 'Both sample_xyz and sample_cyl are given'
        if self.param['sample_xyz'] != None and self.param['sample_cyl'] != None:
            self.errors['sample_dim'] = \
                'Both sample_xyz and sample_cyl are given'
			
# assert that mean grain size != 0 and if mean > 0 then min < mean < max, assure that min non-negative
        if self.param['gen_size'][0] != 0:
            for phase in  phase_list:
                phase_key = 'gen_size_phase_%i' %phase
#                 assert  self.param[phase_key][1] != 0, 'Invalid gen_size command, mean size 0'
                if self.param[phase_key][1] == 0:
                    self.errors['phase_key_1'] = 'Invalid gen_size command, mean size 0'
                if self.param[phase_key][1] > 0:
#                     assert self.param[phase_key][2] < self.param[phase_key][1], \
#                         'grain_min larger than grain_size for phase %i' %phase
#                     assert self.param[phase_key][3] > self.param[phase_key][1], \
#                         'grain_max smaller than grain_size for phase %i' %phase
                    if self.param[phase_key][2] > self.param[phase_key][1]:
                        self.errors['phase_key_2'] = \
                            'gen_size (phase %i): the minimum grain size is made larger than the mean grain size - it should be smaller' %phase
                    if self.param[phase_key][3] < self.param[phase_key][1]:
                        self.errors['phase_key_3'] = \
                            'gen_size (phase %i): the maximum grain size is made smaller than the mean grain size - it should be larger' %phase
                    if self.param[phase_key][2] < 0:
                        self.param[phase_key][2] = 0
                    
#check that the given grain_size and no_grains are consistent with sample_vol, adjust max to sample size
        if self.param['sample_xyz'] != None:
#             assert self.param['sample_xyz'][0] > 0 and self.param['sample_xyz'][1] > 0 and self.param['sample_xyz'][2] > 0,\
#                 'Invalid sample_xyz <= 0'
            if self.param['sample_xyz'][0] < 0 or \
                    self.param['sample_xyz'][1] < 0 or \
                    self.param['sample_xyz'][2] < 0:
                self.errors['sample_xyz'] = 'Invalid sample_xyz all values should be positive'
                self.param['sample_vol'] = None
            else:
                self.param['sample_vol'] = self.param['sample_xyz'][0]*\
                    self.param['sample_xyz'][1]*\
                    self.param['sample_xyz'][2]
                sample_min_dim = min(self.param['sample_xyz'])
        elif self.param['sample_cyl'] != None:
#             assert self.param['sample_cyl'][0] > 0 and self.param['sample_cyl'][1] > 0,\
#                 'Invalid sample_cyl <= 0'
            if self.param['sample_cyl'][0] <= 0 or \
                    self.param['sample_cyl'][1] <= 0:
                self.errors['sample_cyl'] = \
                    'Invalid sample_cyl <= 0'
                self.param['sample_vol'] = None
            else:
                self.param['sample_vol'] = n.pi*self.param['sample_cyl'][0]*\
                    self.param['sample_cyl'][0]*self.param['sample_cyl'][1]/4.
            sample_min_dim = min(self.param['sample_cyl'])
        else:					
            self.param['sample_vol'] = n.inf

        if self.param['sample_vol'] != None and self.param['gen_size'][0] != 0:
            diam_limit = (6*self.param['sample_vol']/\
                         (n.exp(.5)*n.pi*self.param['no_grains']))**(1/3.)
            mean_diam = 0
            vol = []
            for phase in phase_list:

                weight = self.param['no_grains_phase_%i' %phase]/self.param['no_grains']
                vol.append(abs(self.param['gen_size_phase_%i' %phase][1])**3 *\
                               n.pi/6.* self.param['no_grains_phase_%i' %phase])
                mean_diam += abs(self.param['gen_size_phase_%i' %phase][1])*weight
            for i in range(self.param['no_phases']):
                self.param['vol_frac_phase_%i' %phase_list[i]] = vol[i]/n.sum(vol)
            

#             assert mean_diam <= diam_limit, \
#                 'The sample volume is too small to contain the '+\
#                 'specified number of grains with the given grain size'
            if mean_diam > diam_limit:
                self.errors['mean_diam'] = \
                    'The sample volume is too small to contain the '+\
                    'specified number of grains with the given grain size'
            if len(grain_list_size) > 0:
                for i in range(len(grain_list_size)):
#                     assert self.param['size_grains_%s' %(self.param['grain_list'][i])] < diam_limit, \
#                         'The sample volume is too small to contain the '+\
#                         'specified number of grains with the given grain size'
                    if self.param['size_grains_%s' %(self.param['grain_list'][i])] >= diam_limit:
                        self.errors['size_grains_%s' %(self.param['grain_list'][i])] = \
                            'The sample diameter is too small to contain the size of the grain ' +\
                            'by size_grains_%s' %(self.param['grain_list'][i])


#check that a file name with the odf file is input is odf_type chosen to be 2.
        if self.param['peakshape'][0] == 2:
            if self.param['odf_type'] == 2:
                assert self.param['odf_file'] != None, 'No filename given for ODF'

    def show_errors(self):
        if len(self.errors) > 0:
            print 'List of errors and/or inconsistencies found in input: '
            print '----------------------------------------------------- '
            no = 0
            for i in self.errors:
                no += 1
                print 'Error %3i : ' %no, self.errors[i]
            print '----------------------------------------------------- \n'
            


    def initialize(self):
        # Frame generation
        if self.param['make_image'] != 0:
            if self.param['output'] == None:
                self.param['output'] = '.edf'
            if ('.edf' not in self.param['output']) and \
                    ('.tif' not in self.param['output']) and \
                    ('.tif16bit' not in self.param['output']):
                self.param['output'].append('.edf')
     
        # Does output directory exist?
        if not os.path.exists(self.param['direc']):
            os.makedirs(self.param['direc'])
        
	    # Generate FILENAME of frames
        omega_step = self.param['omega_step']
        omega_start  = self.param['omega_start']
        omega_end  = self.param['omega_end']

        modulus = n.abs(omega_end-omega_start)%omega_step
        if  modulus > 1e-9:
            if omega_step-modulus > 1e-9:
                raise ValueError, 'The omega range does not match an integer number of omega steps' 

        # print omega_start,omega_end,omega_step, (n.abs(omega_end-omega_start)+1e-19)%omega_step
        omega_sign = self.param['omega_sign']
        start_frame = self.param['start_frame']
        omegalist = omega_sign*n.arange(omega_start,omega_end+omega_step+1e-19,omega_step)
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
        #print filerange,len(filerange)
        #print omegalist,len(omegalist)
        
        i = 0
        for no in filerange:
            self.frameinfo.append(variables.frameinfo_cont(no))
            self.frameinfo[i].name = '%s/%s%0.4d' \
                %(self.param['direc'],self.param['stem'],filerange[no])
            self.frameinfo[i].omega = omegalist[no];
            self.frameinfo[i].nrefl = 0 # Initialize number of reflections on frame
            self.frameinfo[i].refs = [] # Initialize number of reflections on frame
            i += 1
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
#            print('To make full detector coverage sets theta_max: %f' %theta_max)
            self.param['theta_max'] = theta_max
			

    def init_values(self):
        #Gaelle duplicate code from check_input for the GUI while loading the file to avoid a check input and all the errors launched if file is wrong
        # set all non-read items to defaults
        for item in self.optional_items:
            if (item not in self.param):
                self.param[item] = self.optional_items[item]
            if (self.param[item] == []):
                self.param[item] = self.optional_items[item]*self.param['no_grains']
        #Gaelle : create a new parameter to get all grains id (contains a list of keys with U_grains_Y
        self.param['UgrainsId']=[]  
        self.param['posgrainsId']=[]      
        self.param['epsgrainsId']=[]
        self.param['sizegrainsId']=[]
        # assert that the correct number of arguments are given
        for key in self.param:
            val = self.param[key] 
            if val != None and key != 'output' and key != 'grain_list':
                if key == 'peakshape':
                    if type(val) == type(1):
                        assert val == 0, 'Wrong number of arguments for %s' %key
                        val = [0, 0]
                        self.param[key] = val
                    else:
                        assert len(val) <= 3, 'Wrong number of arguments for %s' %key
                elif key == 'sample_cyl' or key == 'gen_pos':
                    assert len(val) == 2 , 'Wrong number of arguments for %s' %key
                elif key == 'sample_xyz' or 'pos_grains' in key:
                    if 'pos_grains' in key:
                      self.param['posgrainsId'].append(key)
                    assert len(val) == 3, 'Wrong number of arguments for %s' %key
                elif 'gen_size' in key:
                    assert len(val) == 4, 'Wrong number of arguments for %s' %key
                elif 'gen_eps' in key:
                    if type(val) == type(1):
                        assert val == 0, 'Wrong number of arguments for %s' %key
                    else:
                        assert len(val) == 5, 'Wrong number of arguments for %s' %key
                elif key == 'gen_phase':
                    try: 
                        dummy = len(val)
                    except:
                        val = [val]
                        self.param[key] = val
                    assert len(val) > 0, 'Wrong number of arguments for %s' %key
                elif key == 'unit_cell' or 'eps_grains' in key:
                    if 'eps_grains' in key:
                        self.param['epsgrainsId'].append(key)
                    assert len(val) == 6, 'Wrong number of arguments for %s' %key
                elif 'U_grains' in key:
                    #Gaelle : add this key to my new parameter so that it' s easier to get in the GUI.
                    self.param['UgrainsId'].append(key) 
                    #end add Gaelle
                    if len(val) != 3:
                        assert len(val) == 9, 'Wrong number of arguments for %s' %key
                    else:
                        assert val.shape == (3,3), 'Wrong number of arguments for %s' %key
                    # reshape U-matrices
                    self.param[key] = n.array(self.param[key])
                    self.param[key].shape = (3,3)
                elif 'size_grains' in key:
                    #added by gaelle
                    self.param['sizegrainsId'].append(key)
#                else:
#                    assert type(val) != list, 'Wrong number of arguments for %s' %key


        # init no of phases
        no_phases = self.param['no_phases']

        phase_list_structure = []
        phase_list_unit_cell = []
        phase_list_sgno = []
        phase_list_sgname = []
        phase_list_gen_size = []
        phase_list_gen_eps = []
        phase_list = []
        
        for item in self.param:
            if '_phase_' in item:
                if 'structure' in item:
                    phase_list_structure.append(eval(split(item,'_phase_')[1]))
                elif 'unit_cell' in item:
                    phase_list_unit_cell.append(eval(split(item,'_phase_')[1]))
                elif 'sgno' in item:
                    phase_list_sgno.append(eval(split(item,'_phase_')[1]))
                elif 'sgname' in item:
                    phase_list_sgname.append(eval(split(item,'_phase_')[1]))
                elif 'gen_size' in item:
                    phase_list_gen_size.append(eval(split(item,'_phase_')[1]))
                elif 'gen_eps' in item:
                    phase_list_gen_eps.append(eval(split(item,'_phase_')[1]))
                    

        phase_list_structure.sort()
        phase_list_unit_cell.sort()
        phase_list_sgno.sort()
        phase_list_sgname.sort()
        phase_list_gen_size.sort()
        phase_list_gen_eps.sort()
    
        # Init no of grains belonging to phase X if not generated
        if self.param['gen_phase'][0] != 1:
            if len(phase_list) == 0:
                self.param['no_grains_phase_0'] = self.param['no_grains']
            else:
                for phase in phase_list:
                    self.param['no_grains_phase_%i' %phase] = 0
        else:
            for i in range(self.param['no_phases']):
                phase = self.param['gen_phase'][i*2+1]
                no_grains_phase = int(self.param['gen_phase'][i*2+2])
                self.param['no_grains_phase_%i' %phase] = no_grains_phase
            
        #If no structure file is given - unit_cell should be               
        if len(phase_list) == 0:
            # This is a monophase simulation probably using the "old" keywords
            if self.param['structure_file'] == None:
                #print 'NO structure file'
                assert self.param['unit_cell'] != None, \
                    'Missing input: structure_file or unit_cell' 
                # rename keyword
                self.param['unit_cell_phase_0'] = self.param['unit_cell']
                # and delete old one
                del self.param['unit_cell']
                assert self.param['sgno'] != None or self.param['sgname'] != None , \
                    'Missing input: no space group information, please input either sgno or sgname' 
                if self.param['sgno'] == None:
                    self.param['sgno_phase_0'] = sg.sg(sgname = self.param['sgname']).no
                    # rename keyword
                    self.param['sgname_phase_0'] = self.param['sgname']
                    # and delete old one
                    del self.param['sgname']
                else:
                    self.param['sgname_phase_0'] = sg.sg(sgno = self.param['sgno']).name
                    # rename keyword
                    self.param['sgno_phase_0'] = self.param['sgno']
                    # and delete old one
                    del self.param['sgno']
            else:
                # rename keyword
                self.param['structure_phase_0'] = self.param['structure_file']
                # and delete old one
                del self.param['structure_file']
            phase_list = [0]
        self.param['phase_list'] = phase_list

       #make old input file work
        #if len(grain_list_phase) == 0 and self.param['no_phases'] == 1:
           # self.param['grain_list_phase_%i' %self.param['phase_list'][0]] = self.param['grain_list']
            

if __name__=='__main__':

    #import check_input
    try:
        filename = sys.argv[1] 
    except:
        print 'Usage: check_input.py  <input.inp>'
        sys.exit()

    myinput = parse_input(input_file = filename)
    myinput.read()
    myinput.check() 
    if myinput.missing == True:
        print 'MISSING ITEMS'
    myinput.evaluate()
    print myinput.param
