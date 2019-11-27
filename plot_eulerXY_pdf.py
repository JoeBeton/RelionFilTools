from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sys import argv
import parse_star

script,starfile=argv

filament_data = parse_star.readFilamentsFromStarFile(starfile)

#Code to plot the graphs as a pdf file
with PdfPages(starfile[:-5] + '_alignmentplot.pdf') as pdf:
	for i in range(filament_data.number_of_filaments):

		#If statement only plots one in every 50 "tubes" to prevent making an insanely large pdf file
		if i % 50 == 1:

			phi=filament_data.getSpecificFilamentData(i, 'rlnAngleRot')
			the=filament_data.getSpecificFilamentData(i, 'rlnAngleTilt')
			psi=filament_data.getSpecificFilamentData(i, 'rlnAnglePsi')
			xsh=filament_data.getSpecificFilamentData(i, 'rlnOriginXAngst')
			ysh=filament_data.getSpecificFilamentData(i, 'rlnOriginYAngst')
			tracklength = filament_data.getSpecificFilamentData(i, 'rlnHelicalTrackLengthAngst')

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
