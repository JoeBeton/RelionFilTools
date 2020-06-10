import numpy as np
import mrcfile as mrc
import math

class makeDecorationCylinder(object):

    def __init__(self, angpix, boxsize, rise, cyl_outer_diameter, cyl_inner_diameter, layer_thickness, output_filename):
        self.angpix = angpix
        self.boxsize = boxsize
        self.rise = rise
        self.cyl_outer_diameter = cyl_outer_diameter
        self.cyl_inner_diameter = cyl_inner_diameter
        self.layer_thickness = layer_thickness
        self.output_filename = output_filename

        self.cylinder_vol = np.ones((self.boxsize,self.boxsize,self.boxsize), dtype = 'float16')

        if self.boxsize % 2 != 0:
            quit('Please use an even number for box size')

        self.prepareVariables()

        self.layer_thickness += self.cos_edge*2

        self.outer_dnaJ_ring = self.make2Dcircle((self.outer_radius - self.layer_thickness), self.outer_radius)
        self.inner_dnaJ_ring = self.make2Dcircle(self.inner_radius, (self.inner_radius + self.layer_thickness))

        self.makeHelicalCylinder()

        self.saveAsMRC()

    def prepareVariables(self):

        self.z_step = int(self.rise/self.angpix)

        self.cos_edge = 4 #pixels

        self.outer_radius = int((self.cyl_outer_diameter/self.angpix) / 2)
        self.inner_radius = int((self.cyl_inner_diameter/self.angpix) / 2)

        self.each_layer_thickness = int(self.outer_radius - self.inner_radius)/2

        self.center_point = (self.boxsize/2, self.boxsize/2)

    def euclidDistance(self, x, y):
        return math.sqrt(((self.center_point[0] - x)**2) + (self.center_point[1] - y)**2)

    def make2Dcircle(self, inner_radius, outer_radius):

        #make a 2D hollow circle which is used to build the 3D map
        x_edge_lower = (self.boxsize/2) - (outer_radius + self.cos_edge)
        x_edge_upper = (self.boxsize/2) + (outer_radius + self.cos_edge)

        y_edge_lower = (self.boxsize/2) - (outer_radius + self.cos_edge)
        y_edge_upper = (self.boxsize/2) + (outer_radius + self.cos_edge)

        circle_2d = np.ones((self.boxsize, self.boxsize), dtype = 'float16')

        for x in range(self.boxsize):
            #Skips over obvious points
            if x < x_edge_lower or x > x_edge_upper:
                continue
            for y in range(self.boxsize):
                #Skips over obvious points
                if y < y_edge_lower or y > y_edge_upper:
                    continue

                #is the point within the outer diameter?
                radial_distance = self.euclidDistance(x, y)
                if radial_distance <= outer_radius and radial_distance >= inner_radius:
                    circle_2d[x,y] = 0
                    continue
                #is the point within the inner diameter?
                elif radial_distance < inner_radius - self.cos_edge:
                    continue
                #is the point in the outer cosine region?
                elif radial_distance > outer_radius and radial_distance < outer_radius + self.cos_edge:
                    circle_2d[x,y] = math.cos(math.pi * ((((outer_radius + self.cos_edge) - radial_distance)/self.cos_edge)/2))
                    continue
                elif radial_distance < inner_radius and radial_distance > inner_radius - self.cos_edge:
                    circle_2d[x,y] = math.cos(math.pi * ((((inner_radius-self.cos_edge) - radial_distance)/self.cos_edge)/2))
                    continue

        return circle_2d

    def makeHelicalCylinder(self):

        no_of_repeats = int((self.boxsize - (self.boxsize % self.z_step)) / self.z_step)

        #apply the inner layer of decoration
        for repeat in range(no_of_repeats + 1):
            offset = int((self.z_step) * repeat + self.cos_edge)

            for z in range(offset, offset + int((self.z_step/2)) ):
                cos_edge_end = (offset + int((self.z_step/2))) - self.cos_edge

                if z - offset < self.cos_edge:
                    cos_blur_factor = math.cos(math.pi * ((((self.cos_edge - (offset - z)) / self.cos_edge)/2)))
                    inverted_map = np.multiply(-1, np.subtract(self.inner_dnaJ_ring, 1))

                    self.cylinder_vol[z,:,:] = np.add(self.cylinder_vol[z,:,:], np.multiply(inverted_map, cos_blur_factor))
                elif z  > cos_edge_end:
                    try:
                        cos_blur_factor = math.cos(math.pi * ((self.cos_edge/(self.cos_edge + (z - cos_edge_end)))))
                        inverted_map = np.multiply(-1, np.subtract(self.inner_dnaJ_ring, 1))

                        self.cylinder_vol[z,:,:] = np.add(self.cylinder_vol[z,:,:], np.multiply(inverted_map, cos_blur_factor))
                    except IndexError: #happens at the end of the image
                        break
                else:
                    self.cylinder_vol[z,:,:] = self.inner_dnaJ_ring

        #apply the outer layer of decoration
        outer_offset = self.z_step/2
        for repeat in range(no_of_repeats):

            offset = int(((self.z_step) * repeat) + outer_offset + self.cos_edge)

            for z in range(offset, offset + int((self.z_step/2)) ):
                cos_edge_end = (offset + int((self.z_step/2))) - self.cos_edge

                if z - offset < self.cos_edge:
                    cos_blur_factor = math.cos(math.pi * ((((self.cos_edge - (offset - z)) / self.cos_edge)/2)))
                    inverted_map = np.multiply(-1, np.subtract(self.outer_dnaJ_ring, 1))

                    self.cylinder_vol[z,:,:] = np.add(self.cylinder_vol[z,:,:], np.multiply(inverted_map, cos_blur_factor))

                elif z  > cos_edge_end:
                    try:
                        cos_blur_factor = math.cos(math.pi * ((self.cos_edge/(self.cos_edge + (z - cos_edge_end)))))
                        inverted_map = np.multiply(-1, np.subtract(self.outer_dnaJ_ring, 1))

                        self.cylinder_vol[z,:,:] = np.add(self.cylinder_vol[z,:,:], np.multiply(inverted_map, cos_blur_factor))
                    except IndexError: #happens at the end of the image
                        break
                else:
                    self.cylinder_vol[z,:,:] = self.outer_dnaJ_ring#self.cylinder_vol[z,:,:] + self.outer_dnaJ_ring

    def saveAsMRC(self):

        with mrc.new(self.output_filename, overwrite = True) as mrc_file:
            mrc_file.set_data(self.cylinder_vol)
            mrc_file.voxel_size = self.angpix
        print('Saved DNAJ decoration mask as %s' % (self.output_filename))


class makeUnevenMask(object):

    def __init__(self, angpix, boxsize, cyl_outer_diameter, cyl_inner_diameter, output_filename):
        self.angpix = angpix
        self.boxsize = boxsize
        self.cyl_outer_diameter = int(cyl_outer_diameter/self.angpix)
        self.cyl_inner_diameter = int(cyl_inner_diameter/self.angpix)
        self.output_filename

        print(self.cyl_inner_diameter, self.cyl_outer_diameter)

        self.cos_edge = 5
        self.cylinder_vol = np.zeros((self.boxsize,self.boxsize,self.boxsize), dtype = 'float16')
        self.center_point = (self.boxsize/2, self.boxsize/2)

        #self.inner_circle = self.make2Dcircle(0.25, self.cyl_inner_diameter/2)
        self.outer_circle = self.make2Dcircle(1, self.cyl_outer_diameter/2)

        self.makeMask()
        self.saveAsMRC()

    def euclidDistance(self, x, y):
        return math.sqrt(((self.center_point[0] - x)**2) + (self.center_point[1] - y)**2)

    def makeMask(self):

        for z in range(self.boxsize):
            self.cylinder_vol[z,:,:] = self.outer_circle#np.subtract(self.outer_circle, self.inner_circle)


    def make2Dcircle(self, density_value, outer_radius):

        #make a 2D hollow circle which is used to build the 3D map
        x_edge_lower = (self.boxsize/2) - (outer_radius + self.cos_edge)
        x_edge_upper = (self.boxsize/2) + (outer_radius + self.cos_edge)

        y_edge_lower = (self.boxsize/2) - (outer_radius + self.cos_edge)
        y_edge_upper = (self.boxsize/2) + (outer_radius + self.cos_edge)

        circle_2d = np.zeros((self.boxsize, self.boxsize), dtype = 'float16')

        for x in range(self.boxsize):
            #Skips over obvious points
            if x < x_edge_lower or x > x_edge_upper:
                continue
            for y in range(self.boxsize):
                #Skips over obvious points
                if y < y_edge_lower or y > y_edge_upper:
                    continue

                #is the point within the outer diameter?
                radial_distance = self.euclidDistance(x, y)
                if radial_distance <= outer_radius:
                    circle_2d[x,y] = density_value
                    continue
                #is the point in the outer cosine region?
                elif radial_distance > outer_radius and radial_distance < outer_radius + self.cos_edge:
                    circle_2d[x,y] = math.cos(math.pi * (((radial_distance - outer_radius) /self.cos_edge)/2)) * density_value
                    print(circle_2d[x,y])
                    continue
        return circle_2d

    def saveAsMRC(self):

        with mrc.new(self.output_filename, overwrite = True) as mrc_file:
            mrc_file.set_data(self.cylinder_vol)
            mrc_file.voxel_size = self.angpix
        print('Saved fibre downweight mask as %s' % (self.output_filename))

class makeSemiCircleMask(object):

    def __init__(self, angpix, boxsize, cyl_outer_diameter, angle_degrees, twist, output_filename):
        self.angpix = angpix
        self.boxsize = boxsize
        self.cyl_outer_diameter = int(cyl_outer_diameter/self.angpix)
        self.angle = angle_degrees
        self.pixel_twist = twist/(4.75/self.angpix)
        self.output_filename = output_filename

        self.half_boxsize = int(self.boxsize/2)
        self.cos_edge = 5

        if self.angle >= 180:
            self.angle -= 180

        self.ang_low_lim = self.angle - 90
        if self.ang_low_lim < -180:
            self.ang_low_lim = self.ang_high_lim
            self.ang_high_lim += 180

        self.ang_high_lim = self.angle + 90
        if self.ang_high_lim > 180:
            self.ang_high_lim = self.ang_low_lim
            self.ang_low_lim -= 180

        print(self.ang_high_lim, self.ang_low_lim)

        self.polar_coordinates_2d = np.zeros((self.boxsize, self.boxsize))
        #self.polar_coordinates_2d = np.meshgrid(range(-self.half_boxsize, self.half_boxsize), range(-self.half_boxsize, self.half_boxsize))
        self.semicircle3d_vol = np.zeros((self.boxsize, self.boxsize, self.boxsize), dtype = 'float16')
        self.angleLookUpSet = set()

        #self.make2DSemiCircle()
        self.make3DSemiCircle()
        self.saveAsMRC()


    def make2DSemiCircle(self):

        semicircle2d = np.zeros((self.boxsize, self.boxsize))

        for i, x in enumerate(range(-self.half_boxsize, self.half_boxsize)):
            for j, y in enumerate(range(-self.half_boxsize, self.half_boxsize)):

                rho = math.sqrt(x**2 + y**2) * self.angpix

                if rho > (self.cyl_outer_diameter + self.cos_edge):
                    continue

                phi = math.degrees(math.atan2(y, x))

                if self.ang_high_lim > 180:
                    if phi < self.ang_high_lim - 360 and rho <= self.cyl_outer_diameter:
                        semicircle2d[i,j] = 1
                        continue
                elif self.ang_low_lim < -180:
                    if phi > self.ang_low_lim + 360 and rho <= self.cyl_outer_diameter:
                        semicircle2d[i,j] = 1
                        continue

                if self.ang_low_lim > phi or self.ang_high_lim < phi:
                    continue
                elif rho <= self.cyl_outer_diameter:
                    semicircle2d[i,j] = 1
                    continue

                #elif rho > self.cyl_outer_diameter and rho <  (self.cyl_outer_diameter + self.cos_edge):
                #    semicircle2d[i,j] = math.cos(math.pi * ((rho - self.cyl_outer_diameter)/self.cos_edge)/2)

        return semicircle2d

    def make3DSemiCircle(self):

        for z in range(self.boxsize):
            self.semicircle3d_vol[z,:,:] = self.make2DSemiCircle()
            self.ang_low_lim += self.pixel_twist
            self.ang_high_lim += self.pixel_twist

    def saveAsMRC(self):

        with mrc.new(self.output_filename, overwrite = True) as mrc_file:
            mrc_file.set_data(self.semicircle3d_vol)
            mrc_file.voxel_size = self.angpix
        print('Saved twisting semi-circle mask as %s. IMPORTANT: THIS MASK HAS NO COSINE EDGE - REMEMVER TO MAKE ON IN RELION' % (self.output_filename))


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--angpix', nargs = '?', const = 1.7875, type = float, help = 'The pixel size of the output mask')
    parser.add_argument('--box_size', '--box' nargs = '?', type = int, const = 256, help = 'The box size of the output mask')
    parser.add_argument('--out', '--o', nargs = '?', const = 'temp.mrc', help = 'Filename for output map, default is temp.mrc'  )

    parser.add_argument('--dnaJ_mask', '--J', action = 'store_true', help = 'Generate a oscillating density mask to add J seperation to map')
    parser.add_argument('--reduce_fibre', '--R', action = 'store_true', help = 'Generate a mask which downweights the fibre density')
    parser.add_argument('--semi_circle', '--S', action = 'store_true', help = 'Generate a semi-circle shaped mask - IMPORTANT: THIS MASK HAS NO SOFT EDGE')

    parser.add_argument('--cyl_outer_diameter', '--outer' nargs = 1, type = int, required = True, help = 'The outer diameter of respective masks in angstroms')
    parser.add_argument('--fibre_diameter', '--fwidth' nargs = 1, type = int, required = False, help = 'The inner diameter of the mask i.e. diameter of fibre in angstroms')
    parser.add_argument('--helical_rise', '--rise' nargs = '?', type = float, const = 40, help = 'The helical rise of DnaJ seperation in angstroms')
    parser.add_argument('--helical_twist', '--twist' nargs = 1, type = float, required = False, help = 'The helical twist for one 4.75A step in degrees')
    parser.add_argument('--layer_thickness', '--layer', nargs = '?', const = 16, type = int, help = 'Thickness of DNAJ layer (angstroms) - default 16 A'  )
    parse.add_argument('--start_angle', '--ang', nargs = '?', const = 60, type = int, help = 'The starting angle for plotting a twisting semi circle mask')

    args = parser.parse_args()

    if args.dnaJ_mask:
        #makeDecorationCylinder(1.7875, 256, 40, 310, 115, 16)

        if not args.helical_rise:
            quit('Please provide the helical rise using (probably ~ 40) using --rise')
        elif not args.fibre_diameter:
            quit('Please provide the diameter of the fibre (normally ~140 A) using --fibre_diameter')

        makeDecorationCylinder(args.angpix, args.boxsize, args.helical_rise, args.cyl_outer_diameter, args.cyl_inner_diameter, args.layer_thickness, args.out)

    elif args.reduce_fibre:

        if not args.fibre_diameter:
            quit('Please provide the diameter of the fibre - normally ~140 A')

        makeUnevenMask(args.angpix, args.boxsize, args.cyl_outer_diameter, args.cyl_inner_diameter, args.out)

    elif args.semi_circle:

        if not args.helical_twist:
            quit('Please provide the helical twist for one 4.75 A step using --twist e.g.: --twist 0.82')

        makeSemiCircleMask(args.angpix, args.boxsize, args.cyl_outer_diameter, args.start_angle, args.helical_twist, args.out)
