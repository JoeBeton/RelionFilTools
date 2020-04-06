import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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
