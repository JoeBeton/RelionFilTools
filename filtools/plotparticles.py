import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns


from filtools import parse_star, parse_mrc

def plot_changes():
    pass

def plot_filament_pdf(starfile_path):

    filament_data = parse_star.readFilamentsFromStarFile(starfile_path)

    step = int(filament_data.number_of_filaments/100)

    if step < 1:
        step = 1

    #Code to plot the graphs as a pdf file
    with PdfPages(starfile_path[:-5] + '_alignmentplot.pdf') as pdf:
        for i in range(0, filament_data.number_of_filaments, step):
            phi=filament_data.getNumpyFilamentColumn(i, 'rlnAngleRot')
            the=filament_data.getNumpyFilamentColumn(i, 'rlnAngleTilt')
            psi=filament_data.getNumpyFilamentColumn(i, 'rlnAnglePsi')
            xsh=filament_data.getNumpyFilamentColumn(i, 'rlnOriginXAngst')
            ysh=filament_data.getNumpyFilamentColumn(i, 'rlnOriginYAngst')
            tracklength = filament_data.getNumpyFilamentColumn(i, 'rlnHelicalTrackLengthAngst')

            f,(ax1,ax2,ax3,ax4)=plt.subplots(1,4,figsize=(12,6))
            plt.setp([ax1,ax2,ax3,ax4],xticks=[i for i in range(1,len(phi)+1,2)],xlabel='Particle Number')

            ax1.plot(tracklength, phi, 'o')
            ax1.set_ylim([-181,181])
            ax1.set_yticks([i for i in range(-180,180+1,40)])
            ax1.set_ylabel('Phi Angle')
            ax1.set_title('Phi')

            ax2.plot(tracklength, the, 'o')
            ax2.set_yticks([i for i in range(-180,180+1,40)])
            ax2.set_ylim([-181,181])
            ax2.set_ylabel('Theta Angle')
            ax2.set_title('Theta')

            ax3.plot(tracklength, psi, 'o')
            ax3.set_yticks([i for i in range(-180,180+1,40)])
            ax3.set_ylim([-181,181])
            ax3.set_ylabel('Psi Angle')
            ax3.set_title('Psi')

            ax4.plot(tracklength, xsh, 'o',label='Xshift')
            ax4.plot(tracklength, ysh, 'o',label='Yshift')
            ax4.set_title('X/Y Shifts')
            ax4.set_ylabel('Shift (pixels)')
            ax4.legend(loc='upper right')

            pdf.savefig()
            plt.close()

def plotFilamentLengthHistogram(starfile_path):

    filament_data = parse_star.readFilamentsFromStarFile(starfile_path)

    filament_length_array = []
    longest_filament = 0
    shortest_fil = 1e6

    for key in sorted(filament_data.filaments.keys()):
        fil_length = len(filament_data.getNumpyFilamentColumn(key, 'rlnAnglePsi'))
        filament_length_array.append(fil_length)

        if fil_length > longest_filament:
            longest_filament = fil_length
        if fil_length < shortest_fil:
            shortest_fil = fil_length

    bins = longest_filament - shortest_fil + 1

    with PdfPages(starfile_path[:-5] + '_filLengthHist.pdf') as pdf:
        plt.hist(filament_length_array, bins, histtype= 'bar')
        plt.xlabel('Number of particles per filament')
        plt.ylabel('Occurence')
        pdf.savefig()
        plt.close()

    print('Saved a histogram plot showing the filament lengths as: %s_filLengthHist.pdf' % starfile_path[:-5])

def compareFilamentNumbers(starfile1_path, starfile2_path):

    #Initialising seaborn settings
    sns.set()
    sns.set_style("white", {'font.family': ['sans-serif']})
    sns.set_context("poster", font_scale=0.5)
    sns.color_palette(palette=None)

    filament_data1 = parse_star.readFilamentsFromStarFile(starfile1_path)
    filament_data2 = parse_star.readFilamentsFromStarFile(starfile2_path)

    dataset1_uniquefils = set()
    dataset1_fil_lengths = dict()
    dataset2_uniquefils = set()
    dataset2_fil_lengths = dict()

    for key in sorted(filament_data1.filaments.keys()):
        uniquefil_identifier = filament_data1.getRlnFilamentNumberandMicrograph(key)
        dataset1_uniquefils.add(uniquefil_identifier)
        dataset1_fil_lengths[uniquefil_identifier] = filament_data1.getNumberofParticlesinFilament(key)
        #print(len(dataset1_uniquefils), len(dataset2_uniquefils))

    for key in sorted(filament_data2.filaments.keys()):
        uniquefil_identifier = filament_data2.getRlnFilamentNumberandMicrograph(key)
        dataset2_uniquefils.add(uniquefil_identifier)
        dataset2_fil_lengths[uniquefil_identifier] = filament_data2.getNumberofParticlesinFilament(key)

    ratio_array = []

    for unique_fils_ds1 in dataset1_uniquefils:
        if unique_fils_ds1 in dataset2_uniquefils:
            ds1_nofil = dataset1_fil_lengths[unique_fils_ds1]
            ds2_nofil = dataset2_fil_lengths[unique_fils_ds1]

            ratio = ds1_nofil / (ds1_nofil + ds2_nofil)

            ratio_array.append(ratio)
        else:
            ratio_array.append(1)

    for unique_fils_ds2 in dataset2_uniquefils:
        if unique_fils_ds2 in dataset1_uniquefils:
            ds1_nofil = dataset1_fil_lengths[unique_fils_ds2]
            ds2_nofil = dataset2_fil_lengths[unique_fils_ds2]

            ratio = ds1_nofil / (ds1_nofil + ds2_nofil)

            ratio_array.append(ratio)
        else:
            ratio_array.append(0)

    with PdfPages(starfile1_path[:-5] + '_starFileComparison.pdf') as pdf:
        plt.hist(ratio_array, 40, histtype= 'bar')
        plt.xlabel('Ratio between starfile 1 and starfile 2')
        plt.ylabel('Occurence')
        pdf.savefig()
        plt.close()

    print('Saved a histogram plot showing the filament lengths as: %s_starFileComparison.pdf' % starfile1_path[:-5])
