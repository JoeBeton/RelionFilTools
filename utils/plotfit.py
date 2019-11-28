import numpy as np


def linRegress(x,y,m):
    return np.add(y, np.multiply(m,x))

def removeBigGaps(difference_calculation, search_range):
    return [x for x in gaps if x < search_range]

def correctAngles(fitted_angle): # think this is redundant now as I don't bother making everything positive
    if fitted_angle > 180:
        fitted_angle = fitted_angle - 360
    elif fitted_angle < -180:
        fitted_angle = fitted_angle + 360

    return fitted_angle

def exhaustiveLinearSearch(m, y_intercepts, x_values, angles, rln_score, search_range):

    #np_remove_gaps = np.vectorize(removeBigGaps, otypes = [float])

    for y in y_intercepts:
        plot_line = linRegress(x_values,y,m)
        distance_between_points = np.absolute(np.subtract(angles,plot_line))
        remove_big_gaps = removeBigGaps(distance_between_points, search_range)
        invert_score = np.reciprocal(remove_big_gaps)
        scored_errors = np.multiply(invert_score, np.sqrt(rln_score))
        sum_error = np.nansum(scored_errors) * len(scored_errors[~np.isnan(scored_errors)])

        if sum_error != 0:
            try:
                score_y = np.concatenate((score_y, [[y,sum_error]]), axis = 0)
            except NameError:
                score_y = np.array([[y,sum_error]])

    sorted_score = score_y[np.argsort(score_y[:,1])]

    return sorted_score

def optimiseLinearGradient(m_list, y_intercepts, vector, rot_angles, rln_score, search_range):

    for m in m_list:
        for y in y_intercepts:
            plot_line = np_lin_reg(vector,y,m)
            error = np.absolute(np.subtract(rot_angles,plot_line))
            remove_big_gaps = np_remove_gaps(error, search_range)
            invert_score = np.reciprocal(remove_big_gaps)
            scored_errors = np.multiply(remove_big_gaps, np.sqrt(rln_score))
            sum_error = np.nansum(scored_errors) * len(scored_errors[~np.isnan(scored_errors)])

            if sum_error != 0:

                try:
                    score_y = np.concatenate((score_y, [[m,y,sum_error]]), axis = 0)
                except NameError:
                    score_y = np.array([[m,y,sum_error]])

    sorted_score = score_y[np.argsort(score_y[:,2])]
    best_m_score = sorted_score[-1]
    return best_m_score

def fitPolynomial(x_shifts, y_shifts, p_nums, plot = False):

    #try to polynomial fit the x and y shifts - if polyfit fails then Alex's function is run instead
    with warnings.catch_warnings():
        try:
            coefs = np.polyfit(filt_p_nums,filt_x_shifts,2)
            fitt_object = np.poly1d(coefs)
            uniX = [fitt_object(i) for i in xax]
        except Warning:
            print('Polynomial fitting failed')
            uniX=flatten_and_cluster_shifts(Xsh,xax)
        try:
            coefs = np.polyfit(filt_p_nums,filt_y_shifts,2)
            fitt_object = np.poly1d(coefs)
            uniY = [fitt_object(i) for i in xax]
        except Warning:
            print('Polynomial fitting failed')
            uniY=flatten_and_cluster_shifts(Ysh,xax)
    if plot == True:
        if mtID % 50 == 1:

            plt.plot(xax,Ysh, 'o')
            plt.plot(np.arange(xax[0],xax[-1],0.01),[fitt_object(i) for i in np.arange(xax[0],xax[-1],0.01)])
            plt.show()

    return uniX, uniY
