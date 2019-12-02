import numpy as np
from math import isnan

from filtools import parse_star
from utils.plotfit import *

def unify_tilt(starfile_path, plot_changes = False, save_changes = True):
    '''
    Super basic function that resets all tilt angles such that they equal the
    median tilt angle
    '''

    filament_data = parse_star.readFilamentsFromStarFile(starfile_path)

    for filament_no in range(filament_data.number_of_filaments):

        original_tilt_angles = filament_data.getNumpyFilamentColumn(filament_no, 'rlnAngleTilt')
        tilt_median = np.full(len(original_tilt_angles), np.median(original_tilt_angles), dtype = 'float16')

        filament_data.addFilamentDataColumn(filament_no, tilt_median, 'rlnAngleTilt')

    filament_data.writeFilamentsToStarFile()


def reset_tilt(starfile_path, plot_changes = False, save_changes = True):
    '''
    Function to reset the tilt values for all particles to 90 degrees
    '''

    filament_data = parse_star.readFilamentsFromStarFile(starfile_path)

    for filament_no in range(filament_data.number_of_filaments):

        tilt_original = filament_data.getNumpyFilamentColumn(filament_no, 'rlnAngleTilt')

        reset_tilt = np.full(len(tilt_original), 90, dtype = 'float16')
        filament_data.addFilamentDataColumn(filament_no, reset_tilt, 'rlnAngleTilt')

    print('The tilt angles for all filaments have been fitted')
    filament_data.writeFilamentsToStarFile()

    if plot_changes:
        plot_changes()

def unify_rot(starfile_path, twist, rise = 4.75, apix = 1.05, plot_changes = False):
    '''
    Function that fits the rot angles for all filaments - VERY slow at the moment
    '''
    initial_gradient = twist/(rise/apix)

    filament_data = parse_star.readFilamentsFromStarFile(starfile_path)
    print('Filaments with less than 5 particles will be removed')
    filament_data.removeShortFilaments(5)
    print('There are %s filaments to process' % filament_data.number_of_filaments)

    m_list = np.linspace((initial_gradient - (initial_gradient/2)), (initial_gradient+(initial_gradient/2)), 40)

    #np_lin_reg = np.vectorize(linearRegression, otypes = [float])

    for filament_no in range(filament_data.number_of_filaments):

        rot_angles = filament_data.getNumpyFilamentColumn(filament_no, 'rlnAngleRot')
        rln_score = filament_data.getNumpyFilamentColumn(filament_no, 'rlnMaxValueProbDistribution')
        tracklength = filament_data.getNumpyFilamentColumn(filament_no, 'rlnHelicalTrackLengthAngst')
        psi_angles = filament_data.getNumpyFilamentColumn(filament_no, 'rlnAnglePsi')
        xsh = filament_data.getNumpyFilamentColumn(filament_no, 'rlnOriginXAngst')
        ysh = filament_data.getNumpyFilamentColumn(filament_no, 'rlnOriginYAngst')

        #This makes everything take ages but is necessary
        number_of_particles = len(ysh)
        search_limit = number_of_particles*(tracklength[1]*0.1*2)

        #find the best linear plot with a gradient of 4 by changing the y-intercept
        search_limit = number_of_particles*(tracklength[1]*0.1*2)
        y_intercepts = np.arange((-180 - search_limit),(180 + search_limit),0.1)
        #Exhaustive search for gradient 4 and -4:
        pos_m4_search = exhaustiveLinearSearch(initial_gradient, y_intercepts, tracklength, rot_angles, rln_score, 4)
        neg_m4_search = exhaustiveLinearSearch(-initial_gradient, y_intercepts, tracklength, rot_angles, rln_score, 4)

        top_scored_plots = np.concatenate((pos_m4_search[-20:], neg_m4_search[-20:]))
        top_y_intercepts = np.hsplit(top_scored_plots,2)[0]

        best_plot = optimiseLinearGradient(m_list,top_y_intercepts, tracklength, rot_angles, rln_score, 4)

        best_gradient = best_plot[0]
        y_intercept = best_plot[1]

        line_of_best_fit = linRegress(tracklength,y_intercept,best_gradient)

        adjusted_rot_angles = adjustAngletoLOBF(line_of_best_fit, rot_angles, 4)

        filament_data.addFilamentDataColumn(filament_no, adjusted_rot_angles, 'rlnAngleRot')

        if filament_no % 100 == 0:
            print('Tube %s completed' % (filament_no))

    filament_data.writeFilamentsToStarFile()

def unifyXY(starfile_path, plot_changes = False, rot_outliers = []):
    '''
    Fairly basic function to fit the x and y shifts using a polynomial

    Need to implement selective changes i.e. only changing particles that are
    outliers
    '''

    filament_data = parse_star.readFilamentsFromStarFile(starfile_path)

    for filament_no in range(filament_data.number_of_filaments):

        xsh = filament_data.getNumpyFilamentColumn(filament_no, 'rlnOriginXAngst')
        ysh = filament_data.getNumpyFilamentColumn(filament_no, 'rlnOriginYAngst')
        tracklength = filament_data.getNumpyFilamentColumn(filament_no, 'rlnHelicalTrackLengthAngst')

        if len(rot_outliers) > 1:
            #Do the selective changes
            pass
        else:
            uniX, uniY = fitPolynomial(Xsh,Ysh,xax)

        #Checks if the "flatten and cluster" function failed - in which case the shifts are all 'nan'
        if math.isnan(uniX[1]):
            uniX = Xsh
        if math.isnan(uniY[1]):
            uniY = Ysh

        filament_data.addFilamentDataColumn(filament_no, uniX, 'rlnOriginXAngst')
        filament_data.addFilamentDataColumn(filament_no, uniY, 'rlnOriginYAngst')

    filament_data.writeFilamentsToStarFile()


def unify_psi(starfile_path, plot_changes = False):

    filament_data = parse_star.readFilamentsFromStarFile(starfile_path)

    for filament_no in range(filament_data.number_of_filaments):

        psi_original = filament_data.getNumpyFilamentColumn(filament_no, 'rlnAnglePsi')
        relion_score = filament_data.getNumpyFilamentColumn(filament_no, 'rlnMaxValueProbDistribution')
        tracklength = filament_data.getNumpyFilamentColumn(filament_no, 'rlnHelicalTrackLengthAngst')

        psi_low_bound = int(sorted(psi_original)[0])
        psi_high_bound = int(sorted(psi_original)[-1])
        psi_low_bound_search = np.arange(psi_low_bound - 8, psi_low_bound + 8, 0.5)
        psi_high_bound_search = np.arange(psi_high_bound - 8, psi_high_bound + 8, 0.5)
        y_intercepts = np.concatenate((psi_low_bound_search, psi_high_bound_search))

        sorted_plots = exhaustiveLinearSearch(0, y_intercepts, tracklength, psi_original, relion_score, 5)

        best_plot = sorted_plots[-1]
        #line_of_best_fit = linRegress(tracklength, np.full(len(tracklength),best_plot[0]),np.full(len(tracklength),0))
        line_of_best_fit = linRegress(tracklength, best_plot[0],0)
        unified_psi_angles = adjustAngletoLOBF(line_of_best_fit, psi_original, 5)

        filament_data.addFilamentDataColumn(filament_no, unified_psi_angles, 'rlnAnglePsi')

    filament_data.writeFilamentsToStarFile()
