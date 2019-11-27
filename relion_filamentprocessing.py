import numpy as np
import os
import argparse

import parse_star
import filament_functions

parser = argparse.ArgumentParser()

parser.add_argument('--input', '--i', required = True, type= str, help = 'Input the path to a RELION starfile' )
parser.add_argument('--reset_tilt', action = 'store_true', help = 'Option to fit the tilt values for each filament')
parser.add_argument('--make_superparticles', type = int, default = 7, metavar = '*Specify windowing range*', help = 'Option to generate superparticles')
parser.add_argument('--plot_changes', '--p', action = 'store_true', help = 'Option to save pdf plots showing the old and updated angles for each filament')

args=parser.parse_args()

do_plots = args.plot_changes

if args.reset_tilt:
    filament_functions.reset_tilt(args.input, do_plots)

if args.make_superparticles:
    filament_functions.make_superparticles(args.input, args.make_superparticles)
