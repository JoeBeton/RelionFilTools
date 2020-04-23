import numpy as np
import mrcfile as mrc
import math

class makeDecorationCylinder(object):

    def __init__(self, angpix, boxsize, rise, cyl_outer_diameter, cyl_inner_diameter, layer_thickness):
        self.angpix = angpix
        self.boxsize = boxsize
        self.rise = rise
        self.cyl_outer_diameter = cyl_outer_diameter
        self.cyl_inner_diameter = cyl_inner_diameter
        self.layer_thickness = layer_thickness

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

        with mrc.new('temp.mrc', overwrite = True) as mrc_file:
            mrc_file.set_data(self.cylinder_vol)
            mrc_file.voxel_size = self.angpix


if __name__ == '__main__':

    makeDecorationCylinder(1.7875, 256, 40, 340, 130, 25)
