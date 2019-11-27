import numpy as np

import parse_star

def make_superparticles(starfile_path, window_size):
    print('This function will make superparticles from the particles in this directory: ' + starfile_path+', and will use a windowing size of ' + str(window_size)+' particles')

    particles = parse_star.readBlockDataFromStarfile(starfile_path)

    particles_ordered_on_rot = sorted(particles.particle_data_block, key = lambda x:x[particles.headers['rlnAngleRot']])
