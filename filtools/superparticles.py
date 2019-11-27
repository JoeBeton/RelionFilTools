import numpy as np

from filtools import parse_star

def make_superparticles(starfile_path, window_size):

    quit()

    #from EMAN2 import *

    print('This function will make superparticles from this starfile: ' + starfile_path+', and will use a windowing size of ' + str(window_size)+' particles')
    print('The make superparticles function should be run from your RELION directory to ensure folders are made in sensible places')

    particles = parse_star.readBlockDataFromStarfile(starfile_path)

    particles_ordered_on_rot = sorted(particles.particle_data_block, key = lambda x:x[particles.headers['rlnAngleRot']])

    try:
        mkdir('%s' % args.outpath)
        mkdir('%s/Particles' % args.outpath)
    except OSError:
        rmtree('%s' % args.outpath)
        mkdir('%s' % args.outpath)
        mkdir('%s/Particles' % args.outpath)

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

def doNothing():
    pass
