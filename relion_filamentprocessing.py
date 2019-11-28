import numpy as np
import os
import argparse

from filtools import unifyparticles, superparticles, plotparticles

parser = argparse.ArgumentParser()

parser.add_argument('--input', '--i', required = True, type= str, help = 'Input the path to a RELION starfile' )

parser.add_argument('--unify_tilt', '--tilt', action = 'store_true', help = 'Unify the tilt values for filaments')
parser.add_argument('--unify_rot', '--rot', type = float, help = 'Unify the rot values for filaments')
parser.add_argument('--unify_psi', '--psi', action = 'store_true', help = 'Unify the psi values for filaments')

parser.add_argument('--reset_tilt', action = 'store_true', help = 'Option to fit the tilt values for each filament')

parser.add_argument('--make_superparticles', type = int, metavar = '*Specify windowing range*', help = 'Option to generate superparticles')

parser.add_argument('--plot_changes', action = 'store_true', help = 'Option to save pdf plots showing the old and updated angles for each filament')
parser.add_argument('--plot_pdf', '--p', action = 'store_true', help = 'Plot the particle data from filaments into a pdf file')
parser.add_argument('--make_fibretwist_movie', action = 'store_true')

args=parser.parse_args()

do_plots = args.plot_changes

if args.reset_tilt:
    unifyparticles.reset_tilt(args.input, do_plots)

if args.unify_tilt:
    unifyparticles.unify_tilt(args.input, do_plots)

if args.unify_rot:
    unifyparticles.unify_rot(args.input, args.unify_rot, plot_changes = do_plots)

if args.unify_psi:
    unifyparticles.unify_psi(args.input, do_plots)

if args.make_superparticles:
    superparticles.make_superparticles(args.input, args.make_superparticles)

if args.plot_pdf:
    plotparticles.plot_filament_pdf(args.input)

if args.make_fibretwist_movie:
    plotparticles.makeTurningFilamentVideo(args.input, 60)
