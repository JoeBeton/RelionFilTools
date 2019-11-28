import itertools
import numpy as np


class readFilamentsFromStarFile(object):

    '''Makes an object which arranges the data from a star file such that the
    information for each filament are grouped together in individual entries to
    a dictionary which can be accessed by their filament number

    Also includes functions to access specific data and edit and save the loaded
    particle data'''

    def __init__(self, filename, number_of_filaments = 0, headers = {}, optics_info = [],
    filaments = {}, new_data_headers = {}, star_comments = [], number_updated_columns = 0):
        self.filename = filename
        self.number_of_filaments = number_of_filaments
        self.optics_info = optics_info
        self.headers = headers
        self.filaments = filaments
        self.new_data_headers = new_data_headers
        self.star_comments = star_comments
        self.number_updated_columns = number_updated_columns

        self.loadFilamentsFromStar()

    def loadFilamentsFromStar(self):

        '''Reads in the a starfile and sorts the particles from each filament into
        a single dictionary entry into the "filaments" dictionary'''

        with open(self.filename, 'r') as starfile:
            lines = [i.strip() for i in starfile if len(i.strip()) !=0]

        optics = False
        star_data = []
        full_data_dict = {}

        for line in lines:
            if line[0] == '#':
                self.star_comments.append(line)
            elif line == 'data_optics':
                optics = True
                self.optics_info.append(line)
                pass
            elif line == 'data_particles':
                optics = False
                pass
            elif optics == True:
                self.optics_info.append(line)
                pass
            elif line == 'loop_':
                pass
            elif line[0] == '_':
                header_name = line.split()[0][1:]
                header_position = int(line.split()[1][1:])
                self.headers[header_name] = header_position -1
                pass
            else:
                #Construct a dictionary with micrograph names as keys and all particle info from each micrograph as data
                try:
                    full_data_dict[line.split()[self.headers['rlnMicrographName']]].append(line.split())
                except KeyError:
                    full_data_dict[line.split()[self.headers['rlnMicrographName']]] = [line.split()]

        tube_id_column_number = self.headers['rlnHelicalTubeID']

        #Construct the "filaments" dictionary with the filament "number" as a key and all the particles from that filament as entries
        for micrograph_key in sorted(full_data_dict.keys()):
            micrograph_data = full_data_dict[micrograph_key]
            filament_positions = {}
            #Make dictionaries which contain the positions for each filament in the original array/starfile
            for num, particle in enumerate(micrograph_data):
                try:
                    filament_positions[particle[tube_id_column_number]].append(num)
                except KeyError:
                    filament_positions[particle[tube_id_column_number]] = [num]

            #Use the dictionaries with the array positions to find the particles
            #and enter them into the "filaments" dictionary
            for filament_key in sorted(filament_positions.keys()):
                particle_position_list = filament_positions[filament_key]

                #Sort the particles based on helical track tracklength
                #Note: the zip() function rearranges the list into a more numpy like
                #structure: i.e. columns can be easily indexed using array[i] for i'th column
                sorted_particles = sorted([micrograph_data[i] for i in particle_position_list], key = lambda x: float(x[self.headers['rlnHelicalTrackLengthAngst']]))
                #particles = [micrograph_data[i] for i in particle_position_list]

                sorted_structured_particles = list(zip(*sorted_particles))

                self.filaments[self.number_of_filaments] = sorted_structured_particles
                self.number_of_filaments += 1


    def getAllFilamentData(self, filament_number):

        '''Returns all the particle data for a specific filament'''

        return self.filaments[filament_number]

    def getNumpyFilamentColumn(self, filament_number, header_name):

        '''Returns the specified data column as a 1D numpy array from a
        single filament '''

        try:
            return np.array(self.filaments[filament_number][self.headers[header_name]], dtype = 'float16')
        except ValueError:
            return np.array(self.filaments[filament_number][self.headers[header_name]])

    def getStringListFilamentColumn(self, filament_number, header_name):

        '''Returns a list which contains the column data as a list of strings'''

        return [str(i) for i in self.filaments[filament_number][self.headers[header_name]]]

    def addFilamentDataColumn(self, filament_number, new_data_column, name_of_altered_data_column):

        '''Appends a new column of data on a single filament particle stack and
        updates the updated data headers '''

        new_column_number = len(sorted(self.headers.keys())) + len(sorted(self.new_data_headers.keys())) - 1
        filament_data = self.filaments[filament_number]

        filament_data.append(new_data_column)

        self.new_data_headers[name_of_altered_data_column] = new_column_number
        self.filaments[filament_number] = filament_data

    def removeShortFilaments(self, minimum_filament_length, verbose = False):

        '''Remove filaments with less than specified number of particles '''

        if verbose:
            print('There are %i filaments in the input star file' % self.number_of_filaments)

        for fil_no in range(self.number_of_filaments):
            if len(self.getNumpyFilamentColumn(fil_no, 'rlnAngleRot')) < minimum_filament_length:
                del self.filaments[fil_no]

        #remake the filaments dictionary with sequential keys
        temp_filaments = {}
        self.number_of_filaments = 0
        for num, key in enumerate(sorted(self.filaments.keys())):
            temp_filaments[num] = self.filaments[key]
            self.number_of_filaments += 1

        self.filaments = temp_filaments

        if verbose:
            print('There are %i filaments in the saved star file' % self.number_of_filaments)

    def writeFilamentsToStarFile(self, save_updated_data = True):

        '''Writes the data from all the filaments to a starfile, updating columns
        of edited data as specified

        Need to do a fair bit of faffing around for this due to the way the
        data is loaded and handled'''

        save_file_name = self.filename[:-5] + '_updated'

        #Make an ordered list of the original headers
        for key in self.headers.keys():
            try:
                ordered_header_list.append([self.headers[key], key])
            except NameError:
                ordered_header_list = [[self.headers[key],key]]
        ordered_header_list.sort(key = lambda x:x[0])

        #Updates the column positions in headers to write the edited data columns rather than original data
        #Also updates savefilename to include all the edited columns
        if save_updated_data and len(self.new_data_headers.keys()) > 0:
            for key in self.new_data_headers.keys():
                self.headers[key] = self.new_data_headers[key]
                save_file_name = save_file_name + key

        with open(save_file_name + '.star', 'w') as write_star:

            if len(self.optics_info) > 0:
                write_star.write('\n' + self.star_comments[0] + '\n\n')
                for i in self.optics_info:
                    write_star.write(str(i + '\n'))
                write_star.write(str('\n'))

            write_star.write(str('\n ' + self.star_comments[0] + ' \n\n_data_particles\n\n_loop\n'))

            #Write out the header info using the ordered_header_list
            for number, header in enumerate(ordered_header_list):
                write_star.write('_%s #%i\n' % (header[1], number + 1))

            #Write out the data
            for filament_number in range(self.number_of_filaments):
                all_filament_data = []

                #Can't just shunt the data directly into text file due to abstract method of writing new data
                for header in ordered_header_list:
                    all_filament_data.append(self.getStringListFilamentColumn(filament_number, header[1]))

                write_data = list(zip(*all_filament_data))
                for i in write_data:
                    [write_star.write('%s\t' % j) for j in i]
                    write_star.write('\n')
        print('New starfile saved as ' + save_file_name + '.star')


class readBlockDataFromStarfile(object):

    '''Reads in a starfile as a block of data (rather than seperating out individual
    filaments) which is helpful for functions like making superparticles

    This would also be used for standard single particle projects'''

    def __init__(self, filename, headers = {}, optics_info = [],
    particle_data_block = [], new_data_headers = {}, number_of_particles = 0):
        self.filename = filename
        self.headers = headers
        self.optics_info = optics_info
        self.particle_data_block = particle_data_block
        self.new_data_headers = new_data_headers
        self.number_of_particles = number_of_particles

        self.loadBlockDataFromStar()

    def loadBlockDataFromStar(self):

        ''' Loads the particle data from a starfile into a single big list.

        Similar code to loadFilamentsFromStar but doesn't seperate particles into
        individual filaments '''

        with open(self.filename, 'r') as starfile:
            lines = [i.strip() for i in starfile if len(i.strip()) !=0 and i[0] !='#']

        optics = False
        star_data = []
        full_data_dict = {}

        for line in lines:
            if line == 'data_optics':
                optics = True
                self.optics_info.append(line)
                pass
            elif line == 'data_particles':
                optics = False
                pass
            elif optics == True:
                self.optics_info.append(line)
                pass
            elif line == 'loop_':
                pass
            elif line[0] == '_':
                header_name = line.split()[0][1:]
                header_position = int(line.split()[1][1:])
                self.headers[header_name] = header_position -1
                pass
            else:
                #Load a simple multidimensional list for the particles
                try:
                    temp_data_block.append(line.split())
                except NameError:
                    temp_data_block = [line.split()]
                self.number_of_particles += 1

        #Makes a multidimensional array that can be easily indexed for sepecific columns
        self.particle_data_block = list(zip(*temp_data_block))

    def getSpecificDataColumn(self, header_name):

        '''Retrieves a specific column of data as a numpy 1D array'''

        try:
            return np.array(self.particle_data_block[self.headers[header_name]], dtype = 'float16')
        except ValueError:
            return np.array(self.particle_data_block[self.headers[header_name]])

    def addColumntoBlockData(self, new_data_column, header_name):

        '''Adds a new column to the data and updates the header info'''

        new_column_number = len(sorted(self.headers.keys())) + len(sorted(self.new_data_headers.keys())) - 1
        self.new_data_headers[header_name] = new_column_number

        self.particle_data_block.append(new_data_column)

    def writeBlockDatatoStar(self, new_starfile_name, updated_data = True):

        save_file_name = self.filename[:-5] + '_updated'

        #Make an ordered list of the original headers
        for key in self.headers.keys():
            try:
                ordered_header_list.append([self.headers[key], key])
            except NameError:
                ordered_header_list = [[self.headers[key],key]]
        ordered_header_list.sort(key = lambda x:x[0])

        #Updates the column positions in headers to write the edited data columns rather than original data
        #Also updates savefilename to include all the edited columns
        if save_updated_data and len(self.new_data_headers.keys()) > 0:
            for key in self.new_data_headers.keys():
                self.headers[key] = self.new_data_headers[key]
                save_file_name = save_file_name + key

        with open(save_file_name + '.star', 'w') as write_star:

            if len(self.optics_info) > 0:
                write_star.write('\n' + self.star_comments[0] + '\n\n')
                for i in self.optics_info:
                    write_star.write(str(i + '\n'))
                write_star.write(str('\n'))

            write_star.write(str('\n ' + self.star_comments[0] + ' \n\n_data_particles\n\n_loop\n'))

            #Write out the header info using the ordered_header_list
            for number, header in enumerate(ordered_header_list):
                write_star.write('_%s #%i\n' % (header[1], number + 1))

            #Write out the data
            #Can't just shunt the data directly into text file due to abstract method of writing new data
            for header in ordered_header_list:
                all_filament_data.append(self.getStringListFilamentColumn(filament_number, header[1]))

            write_data = list(zip(*all_filament_data))
            for i in write_data:
                [write_star.write('%s\t' % j) for j in i]
                write_star.write('\n')


if __name__ == '__main__':
    #run tests
    quit('I still need to write tests')
