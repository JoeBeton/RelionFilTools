import numpy as np
import os


class readFilamentsFromStarFile(object):

    '''Makes an object which arranges the data from a star file such that the
    information for each filament are grouped together in individual entries to
    a dictionary which can be accessed by their filament number

    Also includes functions to access specific data and edit and save the loaded
    particle data'''

    def __init__(self, filename):
        self.filename = filename
        self.number_of_filaments = 0
        self.optics_info = []
        self.headers = {}
        self.filaments = {}
        self.filament_no_of_particles = {}
        self.new_data_headers = {}
        self.star_comments = []
        self.number_updated_columns = 0
        self.number_of_particles = 0
        self.fil_no_in_micrograph = {}
        self.rln_fil_no_in_micrograph = {}

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

        #Construct the "filaments" dictionary with the filament "number" as a key
        #and all the particles from that filament as entries
        for micrograph_key in sorted(full_data_dict.keys()):
            micrograph_data = full_data_dict[micrograph_key]
            filament_positions = {}

            #Make dictionaries which contain the positions for each filament in the original array/starfile
            for num, particle in enumerate(micrograph_data):

                try:
                    filament_positions[int(particle[tube_id_column_number])].append(num)
                    self.number_of_particles += 1
                except KeyError:
                    filament_positions[int(particle[tube_id_column_number])] = [num]
                    self.number_of_particles += 1

            #Use the dictionaries with the array positions to find the particles
            #and enter them into the "filaments" dictionary
            for filament_key in sorted(filament_positions.keys()):

                #Keeps track of the actual number of filaments in an image
                try:
                    self.fil_no_in_micrograph[micrograph_key] += 1
                except KeyError:
                    self.fil_no_in_micrograph[micrograph_key] = 1

                particle_position_list = filament_positions[filament_key]

                #Sort the particles based on helical track tracklength
                sorted_particles = sorted([micrograph_data[i] for i in particle_position_list], key = lambda x: float(x[self.headers['rlnHelicalTrackLengthAngst']]))

                #The zip() function rearranges the list into a more numpy like
                #structure: i.e. columns can be easily indexed using array[i] for i'th column
                sorted_structured_particles = list(zip(*sorted_particles))

                #Keeps track of the RELION tube numbering - important if new filaments need to be added to the object
                try:
                    self.rln_fil_no_in_micrograph[micrograph_key]
                except KeyError: #catches when a new dictionary entry needs to be made
                    self.rln_fil_no_in_micrograph[micrograph_key] = int(sorted_structured_particles[self.headers['rlnHelicalTubeID']][0])
                else:
                    #complex looking code that just checks if the new tube number is the biggest encountered so far
                    if int(sorted_structured_particles[self.headers['rlnHelicalTubeID']][0]) > self.rln_fil_no_in_micrograph[micrograph_key]:
                        self.rln_fil_no_in_micrograph[micrograph_key] = int(sorted_structured_particles[self.headers['rlnHelicalTubeID']][0])

                self.filaments[self.number_of_filaments] = sorted_structured_particles
                self.filament_no_of_particles[self.number_of_filaments] = len(particle_position_list)
                self.number_of_filaments += 1


    def getAllFilamentData(self, filament_number):

        '''Returns all the particle data for a specific filament'''

        return self.filaments[filament_number]

    def getNumpyFilamentColumn(self, filament_number, header_name):

        '''Returns the specified data column as a 1D numpy array from a
        single filament '''

        try:
            return np.array(self.filaments[filament_number][self.headers[header_name]], dtype = 'float32')
        except ValueError: #happens when data column contains strings e.g. micrograph name
            return np.array(self.filaments[filament_number][self.headers[header_name]])

    def getStringListFilamentColumn(self, filament_number, header_name):

        '''Returns the column data as a list of strings'''

        return [str(i) for i in self.filaments[filament_number][self.headers[header_name]]]

    def getHelicalTrackLengthList(self, filament_number):

        '''Returns the helical track lengths for one filament as a list of strings

        Need a specific function for removing duplicates as some particles have
        track length 0.000000 and some -0.000000 in RELION files'''

        #return [str(i) if i != '-0.000000' else str('0.000000') for i in self.filaments[filament_number][self.headers['rlnHelicalTrackLengthAngst']]]
        return [float(i) for i in self.filaments[filament_number][self.headers['rlnHelicalTrackLengthAngst']]]

    def getTupleFilamentDataColumnSpecificParticles(self, header, fil_no, particle_array):

        return [self.filaments[filament_number][self.headers[header]][i] for i in particle_array]

    def addFilamentDataColumn(self, filament_number, new_data_column, name_of_altered_data_column):

        '''Appends a new column of data on a single filament particle stack and
        updates the updated data headers '''

        new_column_number = len(sorted(self.headers.keys())) + len(sorted(self.new_data_headers.keys())) - 1

        filament_data = self.filaments[filament_number]

        filament_data.append(new_data_column)

        self.new_data_headers[name_of_altered_data_column] = new_column_number
        self.filaments[filament_number] = filament_data

    def addNewFilamentFromOtherStar(self, other_star_object, old_fil_no):

        ''' Adds a filament from another star file as a new filament - crucially
        not inserting the particles with the filaments of

        Uses the parse_star object of the other star file

        This function updates the filament number for the merged filaments to avoid
        annoying problems with RELION errors: its designed to join subtracted particles
        from the same fibre'''

        #mic_name = other_star_object.filaments[old_fil_no][other_star_object.headers['rlnMicrographName']][0]
        mic_name = self.getRlnFilamentNumberandMicrograph(old_fil_no)[0]

        #Work out where to start with numbering the new filaments
        try:
            new_fil_no = self.rln_fil_no_in_micrograph[mic_name] + 1
        except KeyError: #possible error where particles from one image are only in new starfile
            self.rln_fil_no_in_micrograph[mic_name] = 0
            new_fil_no = 1

        #print(self.fil_no_in_micrograph[mic_name])
        #print(new_fil_no)

        #make sure the header information for the new filaments is in the correct format
        temp_particle_block = other_star_object.getAllFilamentData(old_fil_no)

        for header in self.headers.keys():
            #Give the filaments from the new data an original tube number - should help avoid RELION errors/bugs
            if header == 'rlnHelicalTubeID':
                temp_particle_block[self.headers[header]] = [new_fil_no for _ in range(other_star_object.filament_no_of_particles[old_fil_no])]
                continue
            #reorganise the new data so it matches the "old" star file
            try:
                other_star_object.filaments[other_star_object.headers[header]]
            except KeyError:
                raise KeyError('The header option %s is present in the %s starfile but not the %s starfile' % (header, self.filename, other_star_object.filename))
            else:
                temp_particle_block[self.headers[header]] = other_star_object.getStringListFilamentColumn(old_fil_no, header)

        #print(temp_particle_block[self.headers['rlnHelicalTubeID']])

        #Update the various filament numbers and information
        self.filaments[self.number_of_filaments] = temp_particle_block
        self.filament_no_of_particles[self.number_of_filaments] = len(other_star_object.getStringListFilamentColumn(old_fil_no, 'rlnOriginXAngst'))
        self.fil_no_in_micrograph[mic_name] += 1
        self.number_of_filaments += 1
        self.rln_fil_no_in_micrograph[mic_name] += 1
        self.number_of_particles += other_star_object.filament_no_of_particles[old_fil_no]

    def fixExpansionOneFilament(self, fil_no, reference_star, expansion_factor):

        ref_fil_particles = reference_star.getAllFilamentData(fil_no)
        expanded_fil_particles = self.getAllFilamentData(fil_no)

        ref_fil_rot = reference_star.getNumpyFilamentColumn(fil_no, 'rlnAngleRot')
        expanded_fil_rot = self.getNumpyFilamentColumn(fil_no, 'rlnAngleRot')

        #identify the change in rot angle for expanded particles
        first_image_name = reference_star.getStringListFilamentColumn(fil_no, 'rlnImageName')[0]
        expanded_image_name_list = self.getStringListFilamentColumn(fil_no, 'rlnImageName')
        mic_name = self.getStringListFilamentColumn(fil_no, 'rlnImageName')[0]

        #Not sure what this is supposed to achieve
        for position, img_name in enumerate(expanded_image_name_list):
            if img_name == first_image_name:
                try:
                    particle_position_list.append(position)
                except NameError:
                    particle_position_list = [position]
        particle_position_list.sort()

        #get the minimum rot_difference between expanded particles - assumes star is not in order
        #This should probably go into its own function somewhere
        p_count = 0
        for i, img_name in expanded_image_name_list:
            min_rot_angle_diff = 1e9
            if first_image_name == img_name:
                p_count =+ 1
                try:
                    rot_angle_diff = abs(rot - expanded_fil_rot[i])
                except NameError:
                    rot = expanded_fil_rot[i]
                    continue
                else:
                    if rot_angle_diff < min_rot_angle_diff:
                        min_rot_angle_diff = rot_angle_diff
            if pcount == expansion_factor:
                break
        rot_difference = min_rot_angle_diff

        new_filaments[0] = ref_fil_particles

        for expand_factor in range(1, expansion_factor):
            for particle_rot in ref_fil_rot:
                for i, expanded_particle_rot in enumerate(expanded_fil_rot):
                    if abs(expanded_particle_rot - particle_rot)/(expand_factor-1) == rot_difference:
                        try:
                            new_filament_positions[expand_factor].append(i)
                        except KeyError:
                            new_filaments[expand_factor] = [i]

        #make a new dictionary new_filaments which contains all the particle data for each expanded set
        for expanded_fil_keys in sorted(new_filament.keys()):
            new_tube_no = self.fil_no_in_micrograph[mic_name] + expanded_fil_keys
            temp_particle_block = reference_star.getAllFilamentData(fil_no)

            for header in self.headers.keys():
                #Give the filaments from the new data an original tube number - should help avoid RELION errors/bugs
                if header == 'rlnHelicalTubeID':
                    temp_particle_block[self.headers[header]] = tuple([new_fil_no for _ in range(other_star_object.filament_no_of_particles[old_fil_no])])
                    continue

                #reorganise the new data so it matches the "old" star file
                temp_particle_block[self.headers[header]] = self.getTupleFilamentDataColumnSpecificParticles(header, fil_no, new_filaments[expanded_fil_keys])

            new_filaments[expanded_fil_keys] = temp_particle_block


    def removeFilamentDuplicateParticles(self, fil_no):

        hel_track_lengths = self.getHelicalTrackLengthList(fil_no)
        number_of_duplicates = self.filament_no_of_particles[fil_no] - len(set(hel_track_lengths))
        #number_of_duplicates = len(hel_track_lengths) - len(set(hel_track_lengths))

        if number_of_duplicates > 0:
            for duplicate in range(number_of_duplicates):
                for i, track_length in enumerate(hel_track_lengths):
                    try:
                        duplicate_position = hel_track_lengths.index(track_length, i+1)
                        hel_track_lengths.pop(i)
                        self.removeParticleData(fil_no, i)
                        break
                    except ValueError:
                        continue
        else:
            return 1

    def removeFilament(self, fil_no):

        '''Remove all the particles from one filament and update all the relevant
        variables to account for this '''

        self.number_of_particles =- self.filament_no_of_particles[fil_no]

        self.filaments.pop(fil_no)
        self.filament_no_of_particles.pop(fil_no)

        #rename each of the dictionary entries to account for deleted filament
        for new_fil_no in range(fil_no, self.number_of_filaments):
            self.filaments[new_fil_no] = self.filaments[new_fil_no + 1][:]
            self.filament_no_of_particles[new_fil_no] = int(self.filament_no_of_particles)

        self.filaments.pop(self.number_of_filaments)
        #self.filament_no_of_particles.pop(self.number_of_filaments)
        self.number_of_filaments =- 1

    def removeParticleData(self, fil_no, particle_no):
        particles_from_filament = list(zip(*self.getAllFilamentData(fil_no)))
        particles_from_filament.pop(particle_no)
        self.number_of_particles -= 1
        self.filament_no_of_particles[fil_no] -= 1
        self.filaments[fil_no] = list(zip(*particles_from_filament))

    def getNumberofParticlesinFilament(self, filament_no):
        '''Defunct function as now the filament_no_of_particles dict keeps track of this

        will delete soon'''
        return self.filament_no_of_particles[fil_no]

    def getRlnFilamentNumberandMicrograph(self, filament_no):

        micrograph_name = self.getStringListFilamentColumn(filament_no, 'rlnMicrographName')[0]
        rln_tube_number = int(self.getStringListFilamentColumn(filament_no, 'rlnHelicalTubeID')[0])

        return (micrograph_name, rln_tube_number)

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

            write_star.write(str('\n ' + self.star_comments[0] + ' \n\ndata_particles\n\nloop_\n'))

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

    def __init__(self, filename):
        self.filename = filename
        self.headers = {}
        self.optics_info = []
        self.particle_data_block = []
        self.particles = []
        self.new_data_headers = {}
        self.number_of_particles = 0
        self.star_comments = []

        self.loadBlockDataFromStar()

    def loadBlockDataFromStar(self):

        ''' Loads the particle data from a starfile into a single big list.

        Similar code to loadFilamentsFromStar but doesn't seperate particles into
        individual filaments '''

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
            elif line == 'data_':
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
        self.particles = tuple(temp_data_block)
        self.particle_data_block = list(zip(*temp_data_block))

    def getNumpyDataColumn(self, header_name):

        '''Retrieves a specific column of data as a numpy 1D array'''

        try:
            return np.array(self.particle_data_block[self.headers[header_name]], dtype = 'float32')
        except ValueError:
            return np.array(self.particle_data_block[self.headers[header_name]])

    def getStringDataColumn(self, header_name):

        '''Retreves a specific column of data as a list of strings - useful for
        saving star files '''

        return [str(x) for x in self.particle_data_block[self.headers[header_name]]]

    def addColumntoBlockData(self, new_data_column, header_name):

        '''Adds a new column to the data and updates the header info

        Raises an error if the new data column is not the correct type or shape'''

        if type(new_data_column) is not list or not np.ndarray:
            raise TypeError('Only lists or numpy arrays can be used as data columns')

        if len(new_data_column) != self.number_of_particles or len(new_data_column[0]) != 1:
            raise ValueError('The new data column is an incorrect length or shape')

        new_column_number = len(sorted(self.headers.keys())) + len(sorted(self.new_data_headers.keys())) - 1
        self.new_data_headers[header_name] = new_column_number

        self.particle_data_block.append(new_data_column)


    def getOneParticleData(self, particle_no):

        '''Return all the data for one particle - use a list of strings to ensure
        a predictable result rather than a mixed list'''

        return list(self.particles[particle_no])

    def getParticleSpecificDataString(self, particle_no, header_name):
        return self.particle_data_block[self.headers[header_name]][particle_no]

    def getParticleSpecificDataFloat(self, particle_no, header_name):
        return float(self.particle_data_block[self.headers[header_name]][particle_no])

    def getParticleSpecificNewData(self, particle_no, header_name):
        self.particle_data_block[self.new_data_headers[header_name]][particle_no] = new_data

    def addEmptyDataColumn(self, new_header_name):

        ''' Function to add an empty new data column to the particle stack
        which can be edited particle by particle '''

        new_column_number = len(sorted(self.headers.keys())) + len(sorted(self.new_data_headers.keys()))

        ### Check new header name doesn't already exist
        try:
            temp = self.new_data_headers[new_header_name]
            raise ValueError('Nnew header already exists')
        except KeyError:
            self.new_data_headers[new_header_name] = new_column_number

        list_of_zeros = [0] * self.number_of_particles

        self.particle_data_block.append(list_of_zeros)

    def updateParticleDataNewHeader(self, particle_no, header_name, new_data):
        self.particle_data_block[self.new_data_headers[header_name]][particle_no] = new_data

    def updateParticleData(self, particle_no, header_name, new_data):
        self.particle_data_block[self.headers[header_name]][particle_no] = new_data

    def getParticlePositionsBasedOnMetaData(self, header_name, metadata_value):

        ''' Returns a list of the indexes for particles which match the metatdata
        value '''

        return [i for i, n in enumerate(self.particle_data_block[self.headers[header_name]]) if n == metadata_value]


    def updateColumnsWithNewData(self):

        '''Updates the particle_data_block with new data as specified in the
        new_data_headers dictionary '''

        for key in self.new_data_headers.keys():
            self.particle_data_block[self.headers[key]] = self.particle_data_block[self.new_data_headers[key]]

    def selectAngularRange(self, header_name, lower_limit, upper_limit):

        '''Selects the particles within the stated range for a given alignment angle '''

        print('There are ' + str(self.number_of_particles) + ' particles to process')

        anglelist = self.getNumpyDataColumn(header_name)

        for particle_no in range(self.number_of_particles):

            if lower_limit < anglelist[particle_no] and upper_limit > anglelist[particle_no]:
                try:
                    particles_within_angular_range.append(self.getOneParticleData(particle_no))
                except NameError:
                    particles_within_angular_range = [self.getOneParticleData(particle_no)]
            elif (lower_limit - 180) < anglelist[particle_no] and (upper_limit-180) > anglelist[particle_no]:
                try:
                    particles_within_angular_range.append(self.getOneParticleData(particle_no))
                except NameError:
                    particles_within_angular_range = [self.getOneParticleData(particle_no)]

            #if particle_no % 1000 == 0:
            #    print('Particle number ' + str(particle_no) + ' has been processed')

        self.paticles = particles_within_angular_range
        self.particle_data_block = list(zip(*particles_within_angular_range))

        self.new_data_headers[header_name + 'Range' + str(lower_limit) + 'to' +str(upper_limit)] = 0

    def writeBlockDatatoStar(self, save_updated_data = True, save_new_data = False):

        save_file_name = self.filename[:-5] + '_updated'

        if save_new_data:
            for key in self.new_data_headers.keys():
                self.headers[key] = self.new_data_headers[key]

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

            write_star.write(str('\n ' + self.star_comments[0] + ' \n\n_data_particles\n\nloop_\n'))

            #Write out the header info using the ordered_header_list
            for number, header in enumerate(ordered_header_list):
                write_star.write('_%s #%i\n' % (header[1], number + 1))

            #Write out the data
            #Can't just shunt the data directly into text file due to abstract method of writing new data
            for header in ordered_header_list:
                try:
                    all_bulk_data.append(self.getStringDataColumn(header[1]))
                except NameError:
                    all_bulk_data = [self.getStringDataColumn(header[1])]

            write_data = list(zip(*all_bulk_data))

            for i in write_data:
                [write_star.write('%s\t' % j) for j in i]
                write_star.write('\n')

            print('New starfile saved as ' + save_file_name + '.star')

class parseAllModelsinDirectory(object):

    def __init__(self, directory):
        self.directory = directory

        if self.directory[-1] != '/':
            self.directory += '/'

        self.model_star_list = []
        self.number_of_classes = 0
        self.number_of_iterations = 0
        self.classes = {}
        self.headers = {}

        self._getAllModelStarFilenames()
        self._loadAllModelStars()

    def _getAllModelStarFilenames(self):

        dir_list = os.listdir(self.directory)

        for filename in dir_list:
            if filename[-10:] == 'model.star':
                self.model_star_list.append(filename)

    def _loadAllModelStars(self):

        #put the star files in the right order i.e. it001 -> it025
        self.model_star_list = sorted(self.model_star_list, key = lambda x:int(x[-14:-11]))

        for filename in self.model_star_list:
            self._readOneModelStar(filename)
            self.number_of_iterations += 1

        #format the data for easy access
        for key in self.classes.keys():
            self.classes[key] = list(zip(*self.classes[key]))

    def _readOneModelStar(self, filename):

        self.number_of_classes = 0

        with open(self.directory + filename, 'r') as model_star:
            line = model_star.readline()
            while line != 'data_model_classes\n':
                line = model_star.readline()

            while line != 'data_model_class_1\n':
                line = model_star.readline()
                try:
                    temp_data.append(line)
                except NameError:
                    temp_data = [line]

        #Sort out the temp_data into a labelled thing
        for line in temp_data:
            if line[0] == '_':
                header_name = line.split()[0][1:]
                header_position = int(line.split()[1][1:])
                self.headers[header_name] = header_position -1
            elif line.startswith('#'):
                continue
            elif line.strip() =='loop_':
                continue
            elif line.strip() == 'data_model_class_1':
                continue
            elif len(line.strip()) == 0:
                continue
            else:
                try:
                    data.append(line)
                except NameError:
                    data = [line]

        for num, data_line in enumerate(data):
            try:
                self.classes[num].append(data_line.split())
            except KeyError:
                self.classes[num] = [data_line.split()]
            self.number_of_classes += 1

    def getClassDistForOneClass(self, class_number):
        return [float(i) for i in self.classes[class_number][self.headers['rlnClassDistribution']]]


if __name__ == '__main__':
    #run tests
    quit('I still need to write tests')
