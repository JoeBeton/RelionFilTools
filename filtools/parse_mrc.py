import numpy as np
import struct

class readMRCfileNumpy(object):

    def __init__(self, filename, mrc_data = [], apix = 0, box_size = 0, nr_of_open_images = 0,
    stack_position = 0, binary_filename = ''):
        self.filename = filename
        self.mrc_data = mrc_data
        self.apix = apix
        self.box_size = box_size
        self.nr_of_open_images = nr_of_open_images
        self.stack_position = stack_position
        self.binary_filename = binary_filename

        self.readMRCheader()
        self.loadMRCintoNumpyArray()

    def readMRCheader(self):

        self.binary_filename = self.filename[7:]
        self.stack_position = int(self.filename[0:6])

        print(self.binary_filename, self.stack_position)

        with open(self.binary_filename, 'rb') as open_mrc:

            struct_format_string = '<'+(10*'l')+(6*'f')+(3*'l')+(3*'f')+(27*'l')+(3*'f')+(4*'c')+'lfl'
            mrc_header = list(struct.unpack(struct_format_string, open_mrc.read(224)))

            self.box_size = mrc_header[7]
            self.apix = mrc_header[7]/mrc_header[10]


    def loadMRCintoNumpyArray(self):

        read_offset = 1024 + (((self.box_size**2)*16)*(self.stack_position-1))

        with open(self.binary_filename, 'rb') as open_mrc:
            open_mrc.read(read_offset)
            self.mrc_data = np.fromfile(open_mrc, dtype = 'float32', count = (self.box_size**2))

        print(len(self.mrc_data))
