import numpy as np

from filtools import parse_star

def unify_tilt(starfile_path, plot_changes = False, save_changes = True):

    filament_data = parse_star.readFilamentsFromStarFile(starfile_path)

    for filament_no in range(filament_data.number_of_filaments):

        original_tilt_angles = filament_data.getNumpyFilamentColumn(filament_no, 'rlnAngleTilt')
        tilt_median = np.full(len(original_tilt_angles), np.median(original_tilt_angles), dtype = 'float16')

        filament_data.addFilamentDataColumn(filament_no, tilt_median, 'rlnAngletilt')

    filament_data.writeFilamentsToStarFile()


def reset_tilt(starfile_path, plot_changes = False, save_changes = True):

    filament_data = parse_star.readFilamentsFromStarFile(starfile_path)

    for filament_no in range(filament_data.number_of_filaments):

        tilt_original = filament_data.getNumpyFilamentColumn(filament_no, 'rlnAngleTilt')

        reset_tilt = np.full(len(tilt_original), 90, dtype = 'float16')
        filament_data.addFilamentDataColumn(filament_no, reset_tilt, 'rlnAngleTilt')

    print('The tilt angles for all filaments have been fitted')
    filament_data.writeFilamentsToStarFile()

    if plot_changes:
        plot_changes()

#This functions works but is very slow and too much is hardcoded in
def unify_rot(starfile_path, plot_changes = False, save_changes = True):

    filament_data = parse_star.readFilamentsFromStarFile(starfile_path)

    print('There are %s filaments to process' % filament_data.number_of_filaments)

    np_lin_reg = np.vectorize(linearRegression, otypes = [float])

    for filament_no in range(filament_data.number_of_filaments):

        #My horrible attempt to optimise a plot for the Rot angles
        mgphID=re.search(mgre,MT[0][mgphIDX]).group(0)[:-4]
        ptclNum=len(MT)

        rot_angles = filament_data.getNumpyFilamentColumn(filament_no, 'rlnAngleRot')
        rln_score = filament_data.getNumpyFilamentColumn(filament_no, 'rlnMaxValueProbDistribution')
        tracklength = filament_data.getNumpyFilamentColumn(filament_no, 'rlnHelicalTrackLengthAngst')
        psi_angles = filament_data.getNumpyFilamentColumn(filament_no, 'rlnAnglePsi')
        xsh = filament_data.getNumpyFilamentColumn(filament_no, 'rlnOriginXAngst')
        ysh = filament_data.getNumpyFilamentColumn(filament_no, 'rlnOriginYAngst')

        #find the best linear plot with a gradient of 4 by changing the y-intercept
        search_limit = ptclNum*(vector[1]*0.1*2)
        y_intercepts = np.arange((-180-search_limit),(180+search_limit),0.1)

        #Exhaustive search for gradient 4 and -4:
        pos_m4_search = exhaustiveLinearSearch(0.1, y_intercepts, vector, rot_angles, rln_score, 4)
        neg_m4_search = exhaustiveLinearSearch(-0.1, y_intercepts, vector, rot_angles, rln_score, 4)

        top_scored_plots = np.concatenate((pos_m4_search[-20:], neg_m4_search[-20:]))
        top_y_intercepts = np.hsplit(top_scored_plots,2)[0]

        best_plot = optimiseLinearGradient(m_list,top_y_intercepts, vector, rot_angles, rln_score, 4)

        best_gradient = best_plot[0]
        y_intercept = best_plot[1]

        line_of_best_fit = np_lin_reg(vector,y_intercept,best_gradient)

        adjusted_rot_angles = adjustAngletoLOBF(line_of_best_fit, rot_angles)

        ##### Using the updated rot angles to adjust the x and y shifts
        twist, rise, pixel_size = np.full(len(rot_angles),(-0.5)), np.full(len(rot_angles), (4.7)), np.full(len(rot_angles), 1.05)

        delta_rot = np.subtract(np.absolute(rot_angles), np.absolute(adjusted_rot_angles))

        delta_shift_x = np.divide(delta_rot, (np.multiply(twist, rise, np.cos(np.radians(psi_angles)))))
        delta_shift_y = np.divide(delta_rot, (np.multiply(twist, rise, np.cos(np.radians(psi_angles)))))

        scaled_delta_shift_x = np.divide(delta_shift_x, pixel_size)
        scaled_delta_shift_y = np.divide(delta_shift_y, pixel_size)

        fixed_x = np.add(xsh, scaled_delta_shift_x)
        fixed_y = np.add(ysh, scaled_delta_shift_y)

        for num, i in enumerate(sort_array):
            MTs[mtIDX][num][rotIDX]=adjusted_rot_angles[i]
            MTs[mtIDX][num][xshIDX]=fixed_x[i]
            MTs[mtIDX][num][yshIDX]=fixed_y[i]
            #MTs[mtIDX][i][xshIDX]=uniX[i]
            #MTs[mtIDX][i][yshIDX]=uniY[i]
            ptclID=re.search(pcre,MTs[mtIDX][i][imgeIDX]).group(0)[:-1]

        if filament_no % 100 == 0:
            print('Tube %s completed' % (filament_no))

    filament_data.writeFilamentsToStarFile()
