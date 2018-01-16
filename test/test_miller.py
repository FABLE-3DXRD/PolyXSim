from __future__ import absolute_import
import unittest
from polyxsim import check_input
from polyxsim import reflections
import numpy as n

class test_gen_miller(unittest.TestCase):
    def test_(self):  ## test method names begin 'test*'
        param = {}
        param['theta_min'] = 0 
        param['theta_max'] = 5
        param['wavelength'] = 0.2647
        param['unit_cell_phase_0'] = [8.5312,4.8321,10.125,90.00,92.031,90.00]
        param['sgno'] = 4
        param['sgname_phase_0'] = 'P21'
        param['cell_choice_phase_0'] = 'standard' 
        hkl = reflections.gen_miller(param,0)
        self.assertEqual(len(hkl),498)

    def test_open_structure(self):
        param = {}
        param['structure_phase_0'] = 'oPPA.cif'
#        param['structure_datablock'] = 'oPPA'
        structure = reflections.open_structure(param,0)
        self.assertEqual(param['sgno_phase_0'],4)
        self.assertEqual([8.5312,4.8321,10.125,90.00,92.031,90.00],
                          param['unit_cell_phase_0'])
        

    def test_calc_intensity(self):
        myinput = check_input.parse_input(input_file='simul.inp')
        myinput.read()
        myinput.check()
        myinput.initialize()
        myinput.param['structure_phase_0'] = 'oPPA.cif'
        myinput.param['structure_datablock'] = 'oPPA'

        xtal_structure = reflections.open_structure(myinput.param,0)
        hkl = reflections.gen_miller(myinput.param,0)
        hkl = reflections.calc_intensity(hkl,xtal_structure)

    def test_add_intensity(self): 
        param = {}
        param['theta_min'] = 0 
        param['theta_max'] = 3
        param['wavelength'] = 0.2647
        param['unit_cell_phase_0'] = [8.5312,4.8321,10.125,90.00,92.031,90.00]
        param['cell_choice_phase_0']='standard'
        param['sgno_phase_0'] = 4
        param['sgname_phase_0'] = 'p21'
        hkl = reflections.gen_miller(param,0)
        hkl2 = reflections.add_intensity(hkl,param)
        self.assertEqual(n.sum(hkl2[:,3]),len(hkl2)*2**15)
        param['structure_int'] = 2**14
        hkl2 = reflections.add_intensity(hkl,param)
        self.assertEqual(n.sum(hkl2[:,3]),len(hkl2)*2**14)


if __name__ == '__main__':
    unittest.main()
