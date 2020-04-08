import numpy as np

from filtools import parse_star

def getAngles(stack_create_starfile, helix_inimodel2d_starfile, crossover_distance, pixel_size):

    stack_create_data = parse_star.readBlockDataFromStarfile(stack_create_starfile)
    helix_ini2d_output = parse_star.readBlockDataFromStarfile(helix_inimodel2d_starfile)

    #Add data columns for all the info needed for recconstructon:
    stack_create_data.addEmptyDataColumn('rlnAnglePsi')
    stack_create_data.addEmptyDataColumn('rlnAngleRot')
    stack_create_data.addEmptyDataColumn('rlnAngleTilt')
    stack_create_data.addEmptyDataColumn('rlnOriginXAngst')
    stack_create_data.addEmptyDataColumn('rlnOriginYAngst')

    superparticle_width = getParticleSpecificData(hel2dmod_particle_no, 'rlnImageSize')

    #iterate through the images in helix_ini2d starfile and calculate alignement angles
    for hel2dmod_particle_no in range(helix_ini2d_output.number_of_particles):

        image_name = getParticleSpecificData(hel2dmod_particle_no, 'rlnImageName')
        psi = getParticleSpecificData(hel2dmod_particle_no, 'rlnAnglePsi')
        y_origin = getParticleSpecificData(hel2dmod_particle_no, 'rlnAnglePsi')
        #Have to use these ridiculous placeholder names for now as I cba to properly edit RELION code
        position_in_superparticle = getParticleSpecificData(hel2dmod_particle_no, 'rlnCurrentIteration')

        relative_rot_angle = position_in_superparticle/superparticle_width * 180

        particle_no_in_orig_star = stack_create_data.getParticlePositionsBasedOnMetaData('rlnImageName', image_name)[0]

        stack_create_data.updateParticleDataNewHeader(particle_no_in_orig_star, 'rlnAnglePsi', psi)
        stack_create_data.updateParticleDataNewHeader(particle_no_in_orig_star, 'rlnAngleRot', relative_rot_angle)
        stack_create_data.updateParticleDataNewHeader(particle_no_in_orig_star, 'rlnAngleTilt', 90)
        stack_create_data.updateParticleDataNewHeader(particle_no_in_orig_star, 'rlnOriginXAngst', 0)
        stack_create_data.updateParticleDataNewHeader(particle_no_in_orig_star, 'rlnOriginYAngst', y_origin)

    stack_create_data.writeBlockDatatoStar()
