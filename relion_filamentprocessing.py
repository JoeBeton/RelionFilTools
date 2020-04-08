import numpy as np
import os
import argparse

from filtools import unifyparticles, plotparticles, get_helixinimodel2d_angles

parser = argparse.ArgumentParser()

parser.add_argument('--input', '--i', nargs = '+', help = 'Input star file(s)')

parser.add_argument('--unify_tilt', '--tilt', action = 'store_true', help = 'Unify the tilt values for filaments')
parser.add_argument('--unify_rot', '--rot', type = float, help = 'Unify the rot values for filaments')
parser.add_argument('--unify_psi', '--psi', action = 'store_true', help = 'Unify the psi values for filaments')

parser.add_argument('--select_angles', '--sel', nargs = 3, metavar = '[Starfile Header] [Lower limit] [Upper limit]', help = 'Function to select particles with alignment angles that fall within the specified range')

parser.add_argument('--reset_tilt', action = 'store_true', help = 'Option to fit the tilt values for each filament')
parser.add_argument('--remove_shortfils', type = int, help = 'Option to remove the particles from filaments shorter than the stated value')

parser.add_argument('--make_superparticles', type = int, metavar = '*Specify windowing range*', help = 'Option to generate superparticles')

parser.add_argument('--plot_changes', action = 'store_true', help = 'Option to save pdf plots showing the old and updated angles for each filament')
parser.add_argument('--plot_pdf', '--p', action = 'store_true', help = 'Plot the particle data from filaments into a pdf file')
parser.add_argument('--plot_fillenhist', action = 'store_true', help = 'Plot a histogram of the filament lengths')
parser.add_argument('--compare_starfiles', action = 'store_true', help = 'Plot a histogram of the filament lengths')

parser.add_argument('--get_helixinimodel2d_angles', nargs = 1, help = '[angpix] Specialised function for me')

args=parser.parse_args()


do_plots = args.plot_changes

if args.reset_tilt:
    for starfile in args.input:
        unifyparticles.reset_tilt(starfile, do_plots)

if args.remove_shortfils:
    for starfile in args.input:
        unifyparticles.removeShortFilsFromStarfile(starfile, args.remove_shortfils)

if args.select_angles:
    for starfile in args.input:
        unifyparticles.selectParticlesbyAlignmentAngleRange(starfile, args.select_angles[0], int(args.select_angles[1]), int(args.select_angles[2]))

if args.unify_tilt:
    for starfile in args.input:
        unifyparticles.unify_tilt(starfile, do_plots)

if args.unify_rot:
    for starfile in args.input:
        unifyparticles.unify_rot(starfile, args.unify_rot, plot_changes = do_plots)

if args.unify_psi:
    for starfile in args.input:
        unifyparticles.unify_psi(starfile, do_plots)


if args.make_superparticles:
    from filtools import superparticles
    for starfile in args.input:
        superparticles.make_superparticles(starfile, args.make_superparticles)


if args.plot_pdf:
    for starfile in args.input:
        plotparticles.plot_filament_pdf(starfile)
if args.plot_fillenhist:
    for starfile in args.input:
        plotparticles.plotFilamentLengthHistogram(starfile)
if args.compare_starfiles:
    if len(args.input) != 2:
        quit('Please provide two starfiles for this function')
    else:
        plotparticles.compareFilamentNumbers(args.input[0], args.input[1])

if args.get_helixinimodel2d_angles:
    if len(args.input) != 2:
        quit('Please provide two starfiles for this function')
    else:
        get_helixinimodel2d_angles.getAngles(args.input[0], args.input[1], args.get_helixinimodel2d_angles[0])
