#!/usr/bin/env python

from ImageD11 import columnfile
from numpy import array,concatenate,zeros,abs,sqrt,pi,sin
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
        print 'Ready to compare all %i reflections' %nrefl
        print ''
        t1 = time.clock()
        for i in range(1,nrefl):
#            if i%1000 == 0:
#                print 'Comparing reflection %i' %i
            t2 = time.clock()-t1
            print '\rCompared %4.1f pct in %5.1f sec' %(100.*i/nrefl,t2),                        
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
        print 'Number of overlaps %i out of %i refl.' %(nover.sum(),nrefl)

        self.nover = nover

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
        print'Table of overlapping reflections as a function of resolution'
        print ''
        print '%10s%10s%10s%10s%15s' %('stl_max',
                                       'd_min',
                                       'reflec',
                                       'overlap',
                                       'overlap ratio')
        print '%10s%10s%10s%10s%15s' %('-------',
                                       '-----',
                                       '------',
                                       '-------',
                                       '-------------')
        for i in range(len(self.refl_shell)):
            print '%10.4f%10.4f%10i%10i%15.2f' %(self.stlmax_shell[i],
                                                 dmin_shell[i],
                                                 self.refl_shell[i],
                                                 self.overlap_shell[i],
                                                 overlap_frac[i])
        print '%10s%10s%10s%10s%15s' %('-------',
                                       '-----',
                                       '------',
                                       '-------',
                                       '-------------')
        print '%10s%10s%10i%10i%15.2f' %('  All    ',
                                         '',
                                         self.refl_shell.sum(),
                                         self.overlap_shell.sum(),
                                         self.overlap_shell.sum()/self.refl_shell.sum()*100)
        print ''

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
    parser.add_option("-t", "--type", action="store",
                      dest="type", type="string",default='VOL',
                      help="Resolution shells divided to have equal volume (VOL), resolution (RES), or no. of reflections (REF)")
    options , args = parser.parse_args()
    if options.stem == None:
        parser.print_help()
        print "No stem (part of file name before the file number) provided [-n stem]\n"
        sys.exit()
        
    print options

    overlap = determine_overlap(options.direc,
                                options.stem,
                                options.dy,
                                options.dz,
                                options.dw)
    overlap.readgff()
    overlap.readref()
    overlap.find_overlap()
    overlap.shells(nshells=options.nshells,ntype=options.type,wavelength=options.wavelength)
    overlap.output()
    overlap.plot()
