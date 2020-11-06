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

    def reloadFilamentObject(self):

        '''This function reassembles the self object which can be necessary after
        making big changes to the object - e.g. merging two big star files '''

        pass

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

    def getFilamentColumnSpecificDecimalPlaces(self, filament_number, header_name, decimal_places):

        return [round(float(i), decimal_places) for i in self.filaments[filament_number][self.headers[header_name]]]

    def getStringListFilamentColumn(self, filament_number, header_name):

        '''Returns the column data as a list of strings'''

        return [str(i) for i in self.filaments[filament_number][self.headers[header_name]]]

    def getHelicalTrackLengthList(self, filament_number):

        '''Returns the helical track lengths for one filament as a list of strings

        Need a specific function for removing duplicates as some particles have
        track length 0.000000 and some -0.000000 in RELION files'''

        #return [str(i) if i != '-0.000000' else str('0.000000') for i in self.filaments[filament_number][self.headers['rlnHelicalTrackLengthAngst']]]
        return [float(i) for i in self.filaments[filament_number][self.headers['rlnHelicalTrackLengthAngst']]]

    def getFilamentDataColumnSpecificParticles(self, header, fil_no, particle_positions):

        return [self.filaments[fil_no][self.headers[header]][i] for i in particle_positions]

    def addFilamentDataColumn(self, filament_number, new_data_column, name_of_altered_data_column):

        '''Appends a new column of data on a single filament particle stack and
        updates the updated data headers '''

        new_column_number = len(sorted(self.headers.keys())) + len(sorted(self.new_data_headers.keys())) - 1

        filament_data = self.filaments[filament_number]

        filament_data.append(new_data_column)

        self.new_data_headers[name_of_altered_data_column] = new_column_number
        self.filaments[filament_number] = filament_data

    def addNewFilament(self, new_filament_data, mic_name):

        #some test needs to be inserted to check its all in the correct format

        self.filaments[self.number_of_filaments] = new_filament_data
        self.filament_no_of_particles[self.number_of_filaments] = len(new_filament_data[0])
        self.number_of_particles += len(new_filament_data[0])
        try:
            self.fil_no_in_micrograph[mic_name] += 1
        except KeyError:
            self.fil_no_in_micrograph[mic_name] = 1
        try:
            self.rln_fil_no_in_micrograph[mic_name] += 1
        except KeyError:
            self.rln_fil_no_in_micrograph[mic_name] = 1
        self.number_of_filaments += 1

    def addNewFilamentFromOtherStar(self, other_star_object, other_fil_no):

        ''' Adds a filament from another star file as a new filament - crucially
        not inserting the particles with the filaments of

        Uses the parse_star object of the other star file

        This function updates the filament number for the merged filaments to avoid
        annoying problems with RELION errors: its designed to join subtracted particles
        from the same fibre'''

        mic_name = other_star_object.getRlnFilamentNumberandMicrograph(other_fil_no)[0]

        #Work out where to start with numbering the new filaments
        try:
            new_fil_no = self.rln_fil_no_in_micrograph[mic_name] + 1
        except KeyError: #possible error where particles from one image are only in new starfile
            self.rln_fil_no_in_micrograph[mic_name] = 0
            self.fil_no_in_micrograph[mic_name] = 0
            new_fil_no = 1

        #make sure the header information for the new filaments is in the correct format
        temp_particle_block = other_star_object.getAllFilamentData(other_fil_no)

        for header in self.headers.keys():
            #Give the filaments from the new data an original tube number - should help avoid RELION errors/bugs
            if header == 'rlnHelicalTubeID':
                temp_particle_block[self.headers[header]] = [new_fil_no for _ in range(other_star_object.filament_no_of_particles[other_fil_no])]
                continue
            #reorganise the new data so it matches the "old" star file
            try:
                other_star_object.headers[header]
            except KeyError:
                raise KeyError('The header option %s is present in the %s starfile but not the %s starfile' % (header, self.filename, other_star_object.filename))
            else:
                temp_particle_block[self.headers[header]] = other_star_object.getStringListFilamentColumn(other_fil_no, header)

        #Update the various filament numbers and information
        self.filaments[self.number_of_filaments] = temp_particle_block
        self.filament_no_of_particles[self.number_of_filaments] = other_star_object.filament_no_of_particles[other_fil_no]
        self.fil_no_in_micrograph[mic_name] += 1
        self.number_of_filaments += 1
        self.rln_fil_no_in_micrograph[mic_name] += 1
        self.number_of_particles += other_star_object.filament_no_of_particles[other_fil_no]

    def fixExpansionOneFilament(self, fil_no, reference_star_object, expansion_factor):

        ref_fil_particles = reference_star_object.getAllFilamentData(fil_no)
        expanded_fil_particles = self.getAllFilamentData(fil_no)

        ref_fil_rot = reference_star_object.getFilamentColumnSpecificDecimalPlaces(fil_no, 'rlnAngleRot', 2)
        expanded_fil_rot = self.getFilamentColumnSpecificDecimalPlaces(fil_no, 'rlnAngleRot', 2)

        #identify the change in rot angle for expanded particles
        #first_image_name = reference_star_object.getStringListFilamentColumn(fil_no, 'rlnImageName')[0]
        reference_image_names = reference_star_object.getStringListFilamentColumn(fil_no, 'rlnImageName')
        first_image_name = reference_image_names[0]
        expanded_image_name_list = self.getStringListFilamentColumn(fil_no, 'rlnImageName')
        mic_name = self.getStringListFilamentColumn(fil_no, 'rlnMicrographName')[0]

        #build a dictionary that points to where each set of expanded particles are
        expanded_particle_position_dict = {}
        for image_name in reference_image_names:
            for i, expanded_image_name in enumerate(expanded_image_name_list):
                if image_name == expanded_image_name:
                    try:
                        expanded_particle_position_dict[image_name].append(i)
                    except KeyError:
                        expanded_particle_position_dict[image_name] = [i]

        #find rot angle change between expanded particles
        min_rot_angle_diff = 1e9
        for particle_position in expanded_particle_position_dict[first_image_name]:
            try:
                rot_difference = abs(prev_rot - expanded_fil_rot[particle_position])
            except NameError:
                prev_rot = expanded_fil_rot[particle_position]
            else:
                if rot_difference < min_rot_angle_diff:
                    min_rot_angle_diff = round(rot_difference, 2)

        #Following code assembled dictionaries that points to where each expanded particle is in the original filament
        new_filament_positions = {} #dictionary containing the positions of each set of expanded filaments
        for i, image_name in enumerate(reference_image_names):
            for p_position in expanded_particle_position_dict[image_name]:

                rot_change = round(ref_fil_rot[i] - expanded_fil_rot[p_position], 2)
                expanded_particle_set = int(rot_change/min_rot_angle_diff)

                try:
                    new_filament_positions[expanded_particle_set].append(p_position)
                except KeyError:
                    new_filament_positions[expanded_particle_set] = [p_position]

        #make a new dictionary now containing all the particle data for each expanded set
        new_filaments = {}
        for i, expansion_coeff in enumerate(sorted(new_filament_positions.keys())):
            #if expansion_coeff != 0:
            new_tube_ID = self.rln_fil_no_in_micrograph[mic_name] + i
            #print(new_tube_ID)
            temp_particle_block = reference_star_object.getAllFilamentData(fil_no)

            for header in self.headers.keys():
                #Give the filaments from the new data an original tube number - should help avoid RELION errors/bugs
                if header == 'rlnHelicalTubeID':
                    temp_particle_block[self.headers[header]] = [new_tube_ID for _ in range(reference_star_object.filament_no_of_particles[fil_no])]
                    continue
                temp_particle_block[self.headers[header]] = self.getFilamentDataColumnSpecificParticles(header, fil_no, new_filament_positions[expansion_coeff])

            new_filaments[expansion_coeff] = temp_particle_block[:]

        #actually add these new filaments to self.filaments etc
        for expanded_fil_key in sorted(new_filaments.keys()):
            self.addNewFilament(new_filaments[expanded_fil_key], mic_name)

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
        variables to account for this

        This function is slow if running on a large number of filaments - in which
        case use the removeMultipleFilaments function'''

        mic_name = self.getStringListFilamentColumn(fil_no, 'rlnMicrographName')[0]

        self.number_of_particles =- self.filament_no_of_particles[fil_no]
        self.filaments.pop(fil_no)
        self.filament_no_of_particles.pop(fil_no)
        #self.number_of_filaments -= 1

        #rename each of the dictionary entries to account for deleted filament
        for new_fil_no in range(fil_no, self.number_of_filaments):
            self.filaments[new_fil_no] = self.filaments[new_fil_no + 1][:]
            self.filament_no_of_particles[new_fil_no] = int(self.filament_no_of_particles[new_fil_no + 1])

        self.fil_no_in_micrograph[mic_name] =- 1
        self.rln_fil_no_in_micrograph[mic_name] =- 1

        self.filaments.pop(self.number_of_filaments)
        self.filament_no_of_particles.pop(self.number_of_filaments)
        self.number_of_filaments =- 1

    def removeMultipleFilaments(self, filament_numbers):

        for fil_no in filament_numbers:
            mic_name = self.getStringListFilamentColumn(fil_no, 'rlnMicrographName')[0]
            self.number_of_particles =- self.filament_no_of_particles[fil_no]
            self.filaments.pop(fil_no)
            self.filament_no_of_particles.pop(fil_no)
            self.number_of_filaments -= 1
            self.fil_no_in_micrograph[mic_name] =- 1
            self.rln_fil_no_in_micrograph[mic_name] =- 1

        #fixes the numbering system for consistency
        for i, fil_no in enumerate(sorted(self.filaments.keys())):
            if i != fil_no:
                self.filaments[i] = self.filaments[fil_no]
                self.filament_no_of_particles[i] = self.filaments[fil_no]

        for i, fil_no in enumerate(sorted(self.filaments.keys())):
            if fil_no < self.number_of_filaments:
                continue
            else:
                self.filaments.pop(fil_no)
                self.filament_no_of_particles.pop(fil_no)

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

    def writeFilamentsToStarFile(self, save_updated_data = True, suffix = None):

        '''Writes the data from all the filaments to a starfile, updating columns
        of edited data as specified

        Need to do a fair bit of faffing around for this due to the way the
        data is loaded and handled'''

        if not suffix:
            save_file_name = self.filename[:-5] + '_updated'
        else:
            save_file_name = self.filename[:-5] + suffix

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
            for filament_number in sorted(self.filaments.keys()):
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

    def __init__(self, filename, index_particles=False):
        self.filename = filename
        self.index_particles = index_particles
        self.headers = {}
        self.optics_info = []
        self.particle_data_block = []
        self.particles = []
        self.new_data_headers = {}
        self.number_of_particles = 0
        self.star_comments = []
        self.particle_index = {}

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

                if self.index_particles:
                    #This will bug out for expanded particles
                    p_index = {line.split()[self.headers['rlnImageName']]:self.number_of_particles}
                    try:
                        self.particle_index[line.split()[self.headers['rlnMicrographName']]].update(p_index)
                    except KeyError:
                        self.particle_index[line.split()[self.headers['rlnMicrographName']]] = p_index
                self.number_of_particles += 1

        #Makes a multidimensional array that can be easily indexed for sepecific columns
        self.particle_data_block = list(zip(*temp_data_block))

        # self.particles is tuple as this is for reference only - not for updating data
        self.particles = tuple(temp_data_block)

        # hacky bs to set each column in particle_data_block to be a list to allow assignment
        for n, col in enumerate(self.particle_data_block):
            self.particle_data_block[n] = list(col)


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

        if len(new_data_column) != self.number_of_particles:
            raise ValueError('The new data column is an incorrect length or shape')

        new_column_number = len(sorted(self.headers.keys())) + len(sorted(self.new_data_headers.keys()))
        self.new_data_headers[header_name] = new_column_number

        self.particle_data_block.append(new_data_column)


    def getOneParticleData(self, particle_no):

        '''Return all the data for one particle - use a list of strings to ensure
        a predictable result rather than a mixed list'''

        return list(self.particles[particle_no])

    def getParticleValueString(self, particle_no, header_name):
        return self.particle_data_block[self.headers[header_name]][particle_no]

    def getParticleSpecificDataFloat(self, particle_no, header_name):
        return float(self.particle_data_block[self.headers[header_name]][particle_no])

    def getParticleSpecificNewData(self, particle_no, header_name):
        self.particle_data_block[self.new_data_headers[header_name]][particle_no] = new_data

    def getParticleMicrograph(self, particle_no):
        return self.particles[particle_no][self.headers['rlnMicrographName']]

    def getParticleImageName(self, particle_no):
        return self.particles[particle_no][self.headers['rlnImageName']]

    def getParticleNumber(self, mic_name, img_name):
        return self.particle_index[mic_name][img_name]

    def setParticleValue(self, particle_no, header_name, value):
        #print(self.particle_data_block[self.headers[header_name]])
        self.particle_data_block[self.headers[header_name]][particle_no] = value

    def getParticleValue(self, particle_no, header_name):
        return self.particles[particle_no][self.headers[header_name]]

    def addEmptyDataColumn(self, new_header_name):

        ''' Function to add an empty new data column to the particle stack
        which can be edited particle by particle '''

        new_column_number = len(sorted(self.headers.keys())) + len(sorted(self.new_data_headers.keys()))

        ### Check new header name doesn't already exist
        try:
            temp = self.new_data_headers[new_header_name]
            raise ValueError('New header already exists')
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

    def writeBlockDatatoStar(self, save_updated_data=True, save_new_data =False, suffix=None):

        if not suffix:
            save_file_name = self.filename[:-5] + '_updated'
        else:
            save_file_name = self.filename[:-5] + suffix

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
                if not suffix:
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
