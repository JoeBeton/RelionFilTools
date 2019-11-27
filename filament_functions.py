import numpy as np

import parse_star

def reset_tilt(starfile_path, plot_changes = False):

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
