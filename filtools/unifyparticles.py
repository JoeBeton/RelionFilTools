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
    #print('Filaments with less than 5 particles will be removed')
    #filament_data.removeShortFilaments(5)
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

def removeShortFilsFromStarfile(starfile_path, minimum_length):

    print('Removing short filaments from %s' % starfile_path)

    filament_data = parse_star.readFilamentsFromStarFile(starfile_path)

    '''Remove filaments with less than specified number of particles '''
    print('There are %i filaments in the input star file' % filament_data.number_of_filaments)

    for fil_no in range(filament_data.number_of_filaments):
        #This is just to ensure that the starfile name is correctly updated and data structure is maintained
        filament_data.addFilamentDataColumn(fil_no, filament_data.getNumpyFilamentColumn(fil_no,'rlnAngleRot'), 'noShortFilaments')

        if len(filament_data.getNumpyFilamentColumn(fil_no, 'rlnAngleRot')) < minimum_length:
            del filament_data.filaments[fil_no]

    #remake the filaments dictionary with sequential keys
    temp_filaments = {}
    filament_data.number_of_filaments = 0
    for num, key in enumerate(sorted(filament_data.filaments.keys())):
        temp_filaments[num] = filament_data.filaments[key]
        filament_data.number_of_filaments += 1

    filament_data.filaments = temp_filaments

    print('There are %i filaments in the saved star file' % filament_data.number_of_filaments)

    filament_data.writeFilamentsToStarFile()

def removeShortFilsFromObject(filament_object, minimum_length):

    '''Remove filaments with less than specified number of particles '''
    print('There are %i filaments in the input star file' % filament_object.number_of_filaments)

    for fil_no in range(filament_object.number_of_filaments):
        #This is just to ensure that the starfile name is correctly updated and data structure is maintained
        filament_object.addFilamentDataColumn(fil_no, filament_object.getNumpyFilamentColumn(fil_no,'rlnAngleRot'), 'noShortFilaments')

        if len(filament_object.getNumpyFilamentColumn(fil_no, 'rlnAngleRot')) < minimum_length:
            del filament_object.filaments[fil_no]

    #remake the filaments dictionary with sequential keys
    temp_filaments = {}
    filament_object.number_of_filaments = 0
    for num, key in enumerate(sorted(filament_object.filaments.keys())):
        temp_filaments[num] = filament_object.filaments[key]
        filament_object.number_of_filaments += 1

    filament_object.filaments = temp_filaments

    if verbose:
        print('There are %i filaments in the saved star file' % filament_object.number_of_filaments)


def selectParticlesbyAlignmentAngleRange(starfile_path, rln_header_identifier, lower_limit, upper_limit):

    particle_data = parse_star.readBlockDataFromStarfile(starfile_path)

    particle_data.selectAngularRange(rln_header_identifier, lower_limit, upper_limit)
    particle_data.writeBlockDatatoStar()

def orderFilaments(starfile):

    fil_data = parse_star.readFilamentsFromStarFile(starfile)
    fil_data.writeFilamentsToStarFile()

def removeDuplicates(starfile):

    fil_data = parse_star.readFilamentsFromStarFile(starfile)
    starting_particles = fil_data.number_of_particles

    for fil_no in sorted(fil_data.filaments.keys()):

        hel_track_lengths = fil_data.getStringListFilamentColumn(fil_no, 'rlnHelicalTrackLengthAngst')
        number_of_duplicates = len(hel_track_lengths) - len(set(hel_track_lengths))

        if number_of_duplicates > 0:
            for duplicate in range(number_of_duplicates):
                for i, track_length in enumerate(hel_track_lengths):
                    try:
                        duplicate_position = hel_track_lengths.index(track_length, i+1)
                        hel_track_lengths.pop(i)
                        fil_data.removeParticleData(fil_no, i)
                        break
                    except ValueError:
                        continue

    print('%i duplicate particles were removed from the starfile' % (starting_particles - fil_data.number_of_particles))
    fil_data.writeFilamentsToStarFile()
