#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
from ImageD11 import columnfile
from numpy import array,concatenate,floor,arange,zeros,abs,sqrt,pi,sin,int8
import sys

class determine_overlap:
    def __init__(self,direc='.',stem=None,dy=11,dz=11,dw=2.):
        self.direc = direc
        self.stem = stem
        self.setno = 0
        self.dy = dy
        self.dz = dz
        self.dw = dw
        
    def readgff(self):
        filename = '%s/%s.gff' \
            %(self.direc,self.stem)
        gff_file = columnfile.columnfile(filename)
        self.grains = gff_file.grain_id 

    def readref(self):
        omega = array(())
        tth = array(())
        dety  = array(())
        detz  = array(())
        for grainno in self.grains:
            filename = '%s/%s_gr%0.4d_set%0.4d.ref' \
                %(self.direc,self.stem,grainno,self.setno)
            file = columnfile.columnfile(filename)
            tth = concatenate((tth,file.tth))
            omega = concatenate((omega,file.omega))
            dety = concatenate((dety,file.dety))
            detz = concatenate((detz,file.detz))
        order = tth.argsort()
        self.omega = omega[order]
        self.tth = tth[order]
        self.dety = dety[order]
        self.detz = detz[order]

    def find_overlap(self):
        import time
        dtth = .5  # Don't compare position of refs further apart than dtth 
        dpix = (self.dy+self.dz)/2.0
        nrefl = self.tth.shape[0]
        
        nover=zeros((nrefl))
        print('Ready to compare all %i reflections' %nrefl)
        print('')
        t1 = time.time()
        for i in range(1,nrefl):
#            if i%1000 == 0:
#                print 'Comparing reflection %i' %i
            t2 = time.time()-t1
            print('\rCompared %4.1f pct in %5.1f sec' %(100.*i/nrefl,t2), end=' ')                        
            sys.stdout.flush()


            j=i-1
            while j > -1 and self.tth[i]-self.tth[j] < dtth :
                if abs(self.omega[i]-self.omega[j]) < self.dw:
                    peak_distance = sqrt((self.dety[i]-self.dety[j])**2+\
                        (self.detz[i]-self.detz[j])**2)
                    if peak_distance < dpix:
#                    if abs(self.dety[i]-self.dety[j]) < self.dy:
#                      if abs(self.detz[i]-self.detz[j]) < self.dz:
                        nover[i] = 1
                        nover[j] = 1
                j = j - 1
        print('Number of overlaps %i out of %i refl.' %(nover.sum(),nrefl))

        self.nover = nover

    def overlap_new(self):
        """
        New routine for identifying possible spot overlaps. This routine is much faster 
        than the "overlaps", but does not provide information about which specific spots
        are overlapped, but merely flags a spot as overlapped or not.

        This routine sets up an array over all images and adds the value 
        of 1 at the position of each spot
        A filter (elliptic, user defines the size in pixels and image length)
        is made with center at the peak position.
        if the sum of filter * image is larger than 1 it meens the the center 
        of another (or more) spots is found within the area considered to give overlaps.

        """

        from scipy import sparse
        
        # build one big array of reflection info of all grains
#         A = self.grain[0].refs
#         for grainno in range(1,self.param['no_grains']):
#             A = concatenate((A,self.grain[grainno].refs))
        
        #initialize
        nrefl = self.tth.shape[0]
#         nrefl = A.shape[0]
        self.nover=zeros((nrefl))
        overlaps = dict([(i,[]) for i in range(nrefl)])
        print('Ready to compare all %i reflections' %nrefl)

        # Make filter 
#         y_max  = self.param['overlap_specs'][1]*2+1
#         z_max  = self.param['overlap_specs'][2]*2+1
#         om_max = self.param['overlap_specs'][0]*2+1
        y_max = self.dy*2+1
        z_max = self.dz*2+1
        om_max = self.dw

        y_cen = (y_max-1)/2
        z_cen = (z_max-1)/2
        om_cen = 1 #(om_max-1)/2
        print('om_cen',om_cen)

        filter = zeros((y_max,z_max))

        for i in range(y_max):
            for j in range(z_max):
                if (i-y_cen)**2/float(y_cen)**2 + (j-z_cen)**2/float(z_cen)**2 <= 1:
                    filter[i,j] = 1

        # Loop over rotation ranges
        
        #images = range(self.param['frame_range'][0],
         #              self.param['frame_range'][1]+1)

#         # make A array for only including reflections within this rotation range 
#         Atest = A[A[:,22] >= self.param['frame_range'][0]]
#         Atest = Atest[Atest[:,22] <= self.param['frame_range'][1]]
#         nrefl = Atest.shape[0]

        # make stack of empty images as a dictionary of sparse matrices
        omega_sign = 1
        omega_start = 0.
        omega_end = 90.
        omega_step = 0.5
        self.dety_size = 3072
        self.detz_size = 3072


        omegalist = omega_sign*arange(omega_start,omega_end+omega_step+1e-19,omega_step)
        stacksize = len(omegalist)
        frames= zeros((stacksize,int(self.dety_size),int(self.detz_size)),int8)

        # add one in center position of reflections into the image stack. 
        for i in range(nrefl):
#             frameindex =  images.index(int(Atest[i,22]))
            frameindex = omega_sign*\
                floor((self.omega[i]-omega_start)/omega_step)

            frames[frameindex,int(self.dety[i]), int(self.detz[i])] += 1
         
        # determine if reflections are overlapped
        for i in range(nrefl):
            filter_use = filter.copy()
            frameindex =  omega_sign*\
                floor((self.omega[i]-omega_start)/omega_step)
#            print frameindex
            yc = int(self.dety[i])
            zc = int(self.detz[i])
            om_1 = frameindex-om_cen
            om_2 = frameindex+om_cen+1
            y_1 = yc - y_cen
            y_2 = yc + y_cen + 1
            z_1 = zc - z_cen
            z_2 = zc + z_cen + 1

            # Check if filterbox is extending outside detector images 
            # if so change filterbox accordingly
            if om_1 < 0 : 
                om_1 = 0
            if om_2 > stacksize+1 : 
                om_2 = stacksize+1
            if y_1 < 0: 
                filter_use = filter_use[abs(y_1):,:]
                y_1 = 0
            if y_2 > self.dety_size+1: 
                filter_use = filter_use[:self.dety_size+1-y_2 ,:]
                self.dety_size+1
            if z_1 < 0: 
                filter_use = filter_use[:,abs(z_1):]
                z_1 = 0
            if z_2 > self.detz_size + 1:
                filter_use = filter_use[:,:self.detz_size+1-z_2]
                z_2 = self.detz_size+1

#            print om_1, om_2
            box = frames[om_1:om_2, y_1:y_2, z_1:z_2]

            # Calculate no of overlaps with this reflection
            no_over = (box*filter_use).sum()-1

#            print no_over
            # if any put overlap flag in reflection list
            if no_over > 0:
                self.nover[i] += 1
            


    def shells(self,nshells=10,ntype='VOL',wavelength=0.71073):
        stl = sin(pi/180.*self.tth/2)/wavelength
        stlmax = stl[-1]
        vol_shell = (4./3.*pi*stlmax**3)/nshells
        stlmax_shell = zeros(nshells+1)
        n_in_shell = len(stl)/nshells
        
        for i in range(1,nshells+1):

            #equal volume shells
            if ntype == 'VOL':
                V=i*vol_shell
                stlmax_shell[i]= (3*V/(4*pi))**(1/3.)
            #equal stl shells
            if ntype == 'RES':
                stlmax_shell[i] = stlmax/nshells*i
            # equal no refl shells
            if ntype == 'REF':
                stlmax_shell[i] = stl[n_in_shell*i]
            
        refl_shell = zeros(nshells)
        nover_shell = zeros(nshells)
        for i in range(1,nshells+1):
            lt = (stl <= stlmax_shell[i])
            gt = (stl > stlmax_shell[i-1])
            ref_range = lt*gt
            nref_range = ref_range.sum()
            nover_range = self.nover[ref_range]
            nover_shell[i-1] = nover_range.sum()
            refl_shell[i-1] = nref_range
        self.refl_shell = refl_shell
        self.overlap_shell = nover_shell
        self.stlmax_shell = stlmax_shell[1:]
        
    def output(self):
        dmin_shell = 1./(2*self.stlmax_shell)
        overlap_frac = self.overlap_shell*100./self.refl_shell
        print('Table of overlapping reflections as a function of resolution')
        print('')
        print('%10s%10s%10s%10s%15s' %('stl_max',
                                       'd_min',
                                       'reflec',
                                       'overlap',
                                       'overlap ratio'))
        print('%10s%10s%10s%10s%15s' %('-------',
                                       '-----',
                                       '------',
                                       '-------',
                                       '-------------'))
        for i in range(len(self.refl_shell)):
            print('%10.4f%10.4f%10i%10i%15.2f' %(self.stlmax_shell[i],
                                                 dmin_shell[i],
                                                 self.refl_shell[i],
                                                 self.overlap_shell[i],
                                                 overlap_frac[i]))
        print('%10s%10s%10s%10s%15s' %('-------',
                                       '-----',
                                       '------',
                                       '-------',
                                       '-------------'))
        print('%10s%10s%10i%10i%15.2f' %('  All    ',
                                         '',
                                         self.refl_shell.sum(),
                                         self.overlap_shell.sum(),
                                         self.overlap_shell.sum()/self.refl_shell.sum()*100))
        print('')

    def plot(self):
        import pylab as p
        overlap_frac = self.overlap_shell*1./self.refl_shell
        p.plot(self.stlmax_shell,overlap_frac)
        p.xlabel('$\sin(\ttheta)/\lambda$')
        p.ylabel('Fraction of overlap')
        p.show()


if __name__=='__main__':
    
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-n", "--stem", action="store",
                      dest="stem", type="string",
                      help="Stem of the harvest files")
    parser.add_option("-D", "--dir", action="store",
                      dest="direc", type="string",default='.',
                      help="directory of the harvest files")
    parser.add_option("-y", "--ddety", action="store",
                      dest="dy", type="int",default=11,
                      help="distance limit in dety (pixels)")
    parser.add_option("-z", "--ddetz", action="store",
                      dest="dz", type="int",default=11,
                      help="distance limit in detz (pixels)")
    parser.add_option("-w", "--domega", action="store",
                      dest="dw", type="float",default=2.,
                      help="distance limit in omega (degrees)")
    parser.add_option("-l", "--lambda", action="store",
                      dest="wavelength", type="float",default=2.,
                      help="distance limit in omega (degrees)")
    parser.add_option("-s", "--nshells", action="store",
                      dest="nshells", type="int",default=9,
                      help="Show overlap fractions divided into resolution shells")
    parser.add_option("-p", "--plot", action="store_true",
                      dest="doplot",default=False,
                      help="Plot the overlap fractions")

    parser.add_option("-t", "--type", action="store",
                      dest="type", type="string",default='VOL',
                      help="Resolution shells divided to have equal volume (VOL), resolution (RES), or no. of reflections (REF)")
    options , args = parser.parse_args()
    if options.stem == None:
        parser.print_help()
        print("No stem (part of file name before the file number) provided [-n stem]\n")
        sys.exit()
        
    print(options)

    overlap = determine_overlap(options.direc,
                                options.stem,
                                options.dy,
                                options.dz,
                                options.dw)
    overlap.readgff()
    overlap.readref()
    overlap.find_overlap()
#    overlap.overlap_new()
    overlap.shells(nshells=options.nshells,ntype=options.type,wavelength=options.wavelength)
    overlap.output()
    if options.doplot:
        overlap.plot()
