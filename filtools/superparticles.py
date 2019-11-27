import numpy as np
import EMAN2 as eman2

from filtools import parse_star, parse_mrc

def make_superparticles(starfile_path, window_size):

    print('This function will make superparticles from this starfile: ' + starfile_path+', and will use a windowing size of ' + str(window_size)+' particles')
    print('The make superparticles function should be run from your RELION directory to ensure folders are made in sensible places')

    particles = parse_star.readBlockDataFromStarfile(starfile_path)

    #sort the particles based on their rot angle
    sorted_particles = sorted(list(zip(*particles.particle_data_block)), key = lambda x: float(x[particles.headers['rlnAngleRot']]))
    particles.particle_data_block = list(zip(*sorted_particles))

    #Need some code to make the folders etc

    for particle_num in range(particles.number_of_particles):

        particle_names = particles.particle_data_block[particles.headers['rlnImageName']][particle_num:particle_num+window_size]

        #Extremely confusing code that just tries to load the newest particle into memory
        try:
            binary_particles.addImage(particle_names[particle_num+window_size])
        except NameError:
            for name in particle_names:
                try:
                    binary_particles.addImage(name)
                except NameError:
                    binary_particles = parse_mrc.readMRCfileNumpy(name)

        binary_particles.sumLoadedImages()
        binary_particles.saveSummedImage()
        binary_particles.closeImage(0)


        try:
            superparticle_imageNames.append(particle.headers['rlnImageName'][:-5]+'_SP.star')
        except NameError:
            superparticle_imageNames = [particle.headers['rlnImageName'][:-5]+'_SP.star']

    #particles.addColumntoBlockData(superparticle_imageNames, 'rlnImageName')

        '''
    mgphs=[mgph for mgph in group(ptcls,mgphIDX)]

    for mhIDX,mgph in enumerate(mgphs):

        mgphID=re.search(mgre,mgph[0][mgphIDX]).group(0)[:-4]
        print(mgphID)
        mgphSTK=EMData.read_images('%s%s.mrcs' % (args.extractpath,mgphID) )

        for i in range(len(mgph)):
            mgph[i].append(mgphSTK[i])
            ptclID=re.search(pcre,mgph[i][imgIDX]).group(0)[:-1]
            mgphs[mhIDX][i][imgIDX]='%s@%s/Particles/%s_SPs.mrcs' % (ptclID,args.outpath,mgphID)

        superparticles=[]

        for MTidx,MT in enumerate(group(mgph,tubeIDX)):

            for ptclIDX,ptcl in enumerate(MT):
                t=Transform()
                psi,xsh,ysh=ptcl[psiIDX],ptcl[xshIDX],ptcl[yshIDX]
                t.set_params({'type':'2d','alpha':psi,'tx':xsh,'ty':ysh})
                ptcl[-1].transform(t)

                ptclID=re.search(pcre,ptcl[imgIDX]).group(0)[:-1]

            for i in range(len(MT)):
                lw,hi=i-2,i+2

                if lw<0:
                    lw=0
                if hi>len(MT):
                    hi=len(MT)

                avg=Averagers.get('mean')

                for j in range(lw,hi):
                    img=EMData(MT[j][-1])
                    avg.add_image(img)

                avg_img=avg.finish()
                t=Transform()
                t.set_params({'type':'2d','alpha':-MT[i][psiIDX],'tx':-MT[i][xshIDX]/args.apix,'ty':-MT[i][yshIDX]/args.apix})
                avg_img.transform(t)
                avg_img.append_image('./%s/Particles/%s_SPs.mrcs' % (args.outpath,mgphID) )

    ptcls=[ptcl for mgph in mgphs for ptcl in mgph]
    write_star('%s/superparticles.star' % args.outpath,MetaDataLabels,ptcls)
    '''
