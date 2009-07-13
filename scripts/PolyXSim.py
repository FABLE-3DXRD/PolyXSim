#!/usr/bin/env python

# Modules to import 
from polyxsim import gopolyxsim

import logging
logging.basicConfig(level=logging.INFO,format='\n%(levelname)s: %(message)s')

if __name__=="__main__":
    options = None
    try:
        from optparse import OptionParser
        parser = OptionParser()
        options  = gopolyxsim.get_options(parser)
        print options
        gopolyxsim.run(options)
    except:
        if options != None:
            parser.print_help()
        raise 