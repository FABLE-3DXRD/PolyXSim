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

if __name__ == '__main__':
    unittest.main()
