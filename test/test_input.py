from __future__ import absolute_import
import unittest
import numpy as n
from polyxsim import check_input

class test_input(unittest.TestCase):
    def test_reading(self):  ## test method names begin 'test*'
        myinput = check_input.parse_input(input_file='simul.inp')
        myinput.read()
    def test_checking(self):
        myinput = check_input.parse_input(input_file='simul.inp')
        myinput.read()
        myinput.check()
    def test_initialize(self):
        myinput = check_input.parse_input(input_file='simul.inp')
        myinput.read()
        myinput.check()
        myinput.initialize()
    def test_single_phase(self):
        myinput = check_input.parse_input(input_file='simul_single_phase.inp')
        myinput.read()
        gen_eps_before_check = myinput.param['gen_eps'].copy()
        myinput.check()
        myinput.initialize()
        gen_eps_before_after_check = myinput.param['gen_eps'].copy()
        for gen_eps1,gen_eps2 in zip(gen_eps_before_check, gen_eps_before_after_check):
            self.assertEqual(gen_eps1, gen_eps2)


if __name__ == '__main__':
    unittest.main()
