import numpy as np
import os
import argparse

import unifyparticles
import superparticles
import plotparticles

parser = argparse.ArgumentParser()

parser.add_argument('--input', '--i', required = True, type= str, help = 'Input the path to a RELION starfile' )
parser.add_argument('--reset_tilt', action = 'store_true', help = 'Option to fit the tilt values for each filament')
parser.add_argument('--make_superparticles', type = int, default = 7, metavar = '*Specify windowing range*', help = 'Option to generate superparticles')
parser.add_argument('--plot_changes_only', action = 'store_true', help = 'Option to save pdf plots showing the old and updated angles for each filament')
parser.add_argument('--plot_pdf', '--p', action = 'store_true', help = 'Plot the particle data from filaments into a pdf file')
args=parser.parse_args()

do_plots = args.plot_changes

if args.reset_tilt:
    unifyparticles.reset_tilt(args.input, do_plots)

if args.make_superparticles:
    superparticles.make_superparticles(args.input, args.make_superparticles)

if args.plot_pdf:
    plotparticles.plot_pdf(args.input)
