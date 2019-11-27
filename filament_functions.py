import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import parse_star

def fit_tilt(starfile_path, plot_changes = False, save_changes = True):

    filament_data = parse_star.readFilamentsFromStarFile(starfile_path)

    for filament_number in range(filament_data.number_of_filaments):

        original_tilt_angles = filament_data.getNumpyFilamentColumn(filament_number, 'rlnAngleTilt')

        tilt_median = np.full(len(original_tilt_angles), np.median(original_tilt_angles), dtype = 'float16')

        filament_data.addFilamentDataColumn(tilt_median, 'rlnAngletilt')

        star_savefilename = args.fit_tilt[:-5] + '_fittilt'

    write_star(star_savefilename,MetaDataLabels,[ptcl for MT in MTs for ptcl in MT])


def reset_tilt(starfile_path, plot_changes = False, save_changes = True):

    loaded_data = parse_star.readFilamentsFromStarFile(starfile_path)

    for filament_no in range(loaded_data.number_of_filaments):

        tilt_original = loaded_data.getNumpyFilamentColumn(filament_no, 'rlnAngleTilt')

        reset_tilt = np.full(len(tilt_original), 90, dtype = 'float16')
        loaded_data.addFilamentDataColumn(filament_no, reset_tilt, 'rlnAngleTilt')

    print('The tilt angles for all filaments have been fitted')
    loaded_data.writeFilamentsToStarFile()

    if plot_changes:
        plot_changes()

def make_superparticles(starfile_path, window_size):
    print('This function will make superparticles from the particles in this directory: ' + starfile_path+', and will use a windowing size of ' + str(window_size)+' particles')

    particles = parse_star.readBlockDataFromStarfile(starfile_path)

    particles_ordered_on_rot = sorted(particles.particle_data_block, key = lambda x:x[particles.headers['rlnAngleRot']])

def plot_changes():
    pass
