import seaborn as sns
from filtools import parse_star
import os
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def plotClassDistributions(job_directory):

    #Initialising seaborn settings
    sns.set()
    sns.set_style("white", {'font.family': ['sans-serif']})
    sns.set_context("poster", font_scale=0.5)
    sns.color_palette(palette=None)
    sns.set(rc={'figure.figsize':(23.4,16.54)})

    job_classes = parse_star.parseAllModelsinDirectory(job_directory)

    if job_directory[-1] != '/':
        job_directory+='/'

    save_file_path = job_directory + 'class_distribution.pdf'

    with PdfPages(save_file_path) as pdf:
        for i in range(job_classes.number_of_classes):
            class_label = 'Class number ' + str(i+1)
            sns.lineplot(range(job_classes.number_of_iterations), job_classes.getClassDistForOneClass(i), label = class_label)
        plt.xlabel('Iteration Number')
        plt.ylabel('Class Distribution')
        pdf.savefig()
        plt.close()

    print('Saved a plot of the class distribution for job %s as %s' % (job_directory, save_file_path))


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '--i', nargs = '+', help = 'Class2D/3D job directory (e.g. Class2D/job010/) or directories')
    args = parser.parse_args()

    for folder in args.input:
        if not os.path.isdir(folder):
            print('%s is not a directory - please provide the path to a directory' % (folder))
            continue

        plotClassDistributions(folder)
