from filtools import(
            unifyparticles,
            plotparticles,
            get_helixinimodel2d_angles,
            parse_star
            )

import numpy as np


test_starfile1 = 'tests/test_data/test_star1.star'
test_starfile2 = 'tests/test_data/test_star2.star'

class TestFunctions:

    def test_updateAlignments(self):

        unifyparticles.updateAlignments(test_starfile1, test_starfile2)
        saved_star_path = 'tests/test_data/test_star1_updatedAlignments.star'
        saved_star = parse_star.readBlockDataFromStarfile(saved_star_path)
        updating_angles_star = parse_star.readBlockDataFromStarfile(test_starfile2)

        rot_saved_file = saved_star.getNumpyDataColumn('rlnAngleRot')
        rot_orig_file = updating_angles_star.getNumpyDataColumn('rlnAngleRot')

        assert np.isclose(rot_saved_file, rot_orig_file).all()

    def test_updateCTF(self):

        unifyparticles.updateCTF(test_starfile1, test_starfile2)
        saved_star_path = 'tests/test_data/test_star1_updatedCTF.star'
        saved_star = parse_star.readBlockDataFromStarfile(saved_star_path)
        updating_angles_star = parse_star.readBlockDataFromStarfile(test_starfile2)

        rot_saved_file = saved_star.getNumpyDataColumn('rlnDefocusU')
        rot_orig_file = updating_angles_star.getNumpyDataColumn('rlnDefocusU')

        assert np.isclose(rot_saved_file, rot_orig_file).all()
