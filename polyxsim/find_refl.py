import numpy as n
from xfab import tools
from xfab import sg
from xfab import detector
from xfab.structure import int_intensity
import variables,check_input,file_io
import sys
from ImageD11 import blobcorrector
import logging
logging.basicConfig(level=logging.DEBUG,format='%(levelname)s %(message)s')

A_id = variables.refarray().A_id


class find_refl:
    def __init__(self,param,hkl,options):
        self.killfile = options.killfile
        self.param = param
        self.hkl = hkl
        self.grain = []
    
        # Simple transforms of input and set constants
        self.K = -2*n.pi/self.param['wavelength']
        self.S = n.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
        
        # Detector tilt correction matrix
        self.R = tools.detect_tilt(self.param['tilt_x'],
                                   self.param['tilt_y'],
                                   self.param['tilt_z'])

        # Spatial distortion
        if self.param['spatial'] != None:
            self.spatial = blobcorrector.correctorclass(self.param['spatial'])

        # %No of images
        self.nframes = (self.param['omega_end']-self.param['omega_start'])/self.param['omega_step']
        
        # Generate Miller indices for reflections within a certain resolution
        logging.info('Generating reflections')


        print 'Finished generating reflections\n'
    
    def run(self):
        spot_id = 0
        # Generate orientations of the grains and loop over all grains
        print 'no of grains ',self.param['no_grains']
        for grainno in range(self.param['no_grains']):
            A = []
            U = self.param['U_grains_%s' %(self.param['grain_list'][grainno])]
            if len(self.param['phase_list']) == 1:
                phase = self.param['phase_list'][0]
            else:
                phase = self.param['phase_grains_%s' %(self.param['grain_list'][grainno])]
            unit_cell = self.param['unit_cell_phase_%i' %phase]

            self.grain.append(variables.grain_cont(U))
            gr_pos = n.array(self.param['pos_grains_%s' %(self.param['grain_list'][grainno])])
            gr_eps = n.array(self.param['eps_grains_%s' %(self.param['grain_list'][grainno])])
            # Calculate the B-matrix based on the strain tensor for each grain
            B = tools.epsilon2B(gr_eps,unit_cell) 
            # add B matrix to grain container
            self.grain[grainno].B = B
            V = tools.CellVolume(unit_cell)
            grain_vol = n.pi/6 * self.param['size_grains_%s' %self.param['grain_list'][grainno]]**3 

#            print 'GRAIN NO: ',self.param['grain_list'][grainno]
#            print 'GRAIN POSITION of grain ',self.param['grain_list'][grainno],': ',gr_pos
#            print 'STRAIN TENSOR COMPONENTS (e11 e12 e13 e22 e23 e33) of grain ',self.param['grain_list'][grainno],':\n',gr_eps
#            print 'U of grain ',self.param['grain_list'][grainno],':\n',U
            nrefl = 0
  
            # Calculate these values:
            # totalnr, grainno, refno, hkl, omega, 2theta, eta, dety, detz
            # For all reflections in Ahkl that fulfill omega_start < omega < omega_end.
            # All angles in Grain are in degrees
            for hkl in self.hkl[self.param['phase_list'].index(phase)]:
                check_input.interrupt(self.killfile)
                Gc = n.dot(B,hkl[0:3])
                Gw =   n.dot(self.S,n.dot(U,Gc))
                tth = tools.tth2(Gw,self.param['wavelength'])
                costth = n.cos(tth)
                (Omega, Eta) = tools.find_omega_wedge(Gw,
                                                      tth,
                                                      self.param['wedge'])
                if len(Omega) > 0:
                    for solution in range(len(Omega)):
                        omega = Omega[solution]
                        eta = Eta[solution]
                        if  (self.param['omega_start']*n.pi/180) < omega and\
                                omega < (self.param['omega_end']*n.pi/180):
                            # form Omega rotation matrix
                            Om = tools.OMEGA(omega)
                            Gt = n.dot(Om,Gw)
  
                            # Calc crystal position at present omega
                            [tx,ty]= n.dot(Om[:2,:2],gr_pos[:2])
                            tz = gr_pos[2]
                            
                            # Calc detector coordinate for peak 
                            (dety, detz) = detector.det_coor(Gt, 
                                                             costth,
                                                             self.param['wavelength'],
                                                             self.param['distance'],
                                                             self.param['y_size'],
                                                             self.param['z_size'],
                                                             self.param['dety_center'],
                                                             self.param['detz_center'],
                                                             self.R,
                                                             tx,ty,tz)

                            #If shoebox extends outside detector exclude it
                            if (-0.5 > dety) or\
                               (dety > self.param['dety_size']-0.5) or\
                               (-0.5 > detz) or\
                               (detz > self.param['detz_size']-0.5):
                                continue


                            if self.param['spatial'] != None :
                                # To match the coordinate system of the spline file
                                (x,y) = detector.detyz2xy([dety,detz],
                                                          self.param['o11'],
                                                          self.param['o12'],
                                                          self.param['o21'],
                                                          self.param['o22'],
                                                          self.param['dety_size'],
                                                          self.param['detz_size'])
                                # Do the spatial distortion
                                (xd,yd) = self.spatial.distort(x,y)
                                # transform coordinates back to dety,detz
                                (detyd,detzd) = detector.xy2detyz([xd,yd],
                                                          self.param['o11'],
                                                          self.param['o12'],
                                                          self.param['o21'],
                                                          self.param['o22'],
                                                          self.param['dety_size'],
                                                          self.param['detz_size'])
                            else:
                                detyd = dety
                                detzd = detz

                            if self.param['beampol_apply'] == 1:
                                #Polarization factor (Kahn, J. Appl. Cryst. (1982) 15, 330-337.)
                                rho = n.pi/2.0 + eta + self.param['beampol_direct']*n.pi/180.0 
                                P = 0.5 * (1 + costth*costth +\
                                         self.param['beampol_factor']*n.cos(2*rho)*n.sin(tth)**2)
                            else:
                                P = 1.0

                            #Lorentz factor
                            if self.param['lorentz_apply'] == 1:
                                if eta != 0:
                                    L=1/(n.sin(tth)*abs(n.sin(eta)))
                                else:
                                    L=n.inf;
                            else:
                                L = 1.0
 
                            overlaps = 0 # set the number overlaps to zero
                            
                            if self.param['intensity_const'] != 1: 
                                intensity = int_intensity(hkl[3],
                                                          L,
                                                          P,
                                                          self.param['beamflux'],
                                                          self.param['wavelength'],
                                                          V,
                                                          grain_vol)
                            else:
                                intensity = hkl[3]

                            A.append([self.param['grain_list'][grainno],
                                      nrefl,spot_id,
                                      hkl[0],hkl[1],hkl[2],
                                      tth,omega,eta,
                                      dety,detz,
                                      detyd,detzd,
                                      Gw[0],Gw[1],Gw[2],
                                      L,P,hkl[3],intensity,overlaps])
                            nrefl = nrefl+1
                            spot_id = spot_id+1

#           print 'Length of Grain', len(self.grain[0].refl)
            A = n.array(A)
            if len(A) > 0:
                # sort rows according to omega
                A = A[n.argsort(A,0)[:,A_id['omega']],:]
                
                # Renumber the reflections  
                A[:,A_id['ref_id']] = n.arange(nrefl)

                # Renumber the spot_id
                A[:,A_id['spot_id']] = n.arange(n.min(A[:,A_id['spot_id']]),
                                            n.max(A[:,A_id['spot_id']])+1)
            else:
                A = n.zeros((0,len(A_id)))
            # save reflection info in grain container
            self.grain[grainno].refs = A 
            print '\rDone %3i grain(s) of %3i' %(grainno+1,self.param['no_grains']),
            sys.stdout.flush()

        print '\n'


    def overlap(self):
        
        dtth = 1*n.pi/180.  # Don't compare position of refs further apart than dtth 

        # build one big array of reflection info of all grains
        A = self.grain[0].refs
        for grainno in range(1,self.param['no_grains']):
            A = n.concatenate((A,self.grain[grainno].refs))
        logging.debug('Finished concatenating ref arrays')
        A = A[n.argsort(A,0)[:,A_id['tth']],:] # sort rows according to tth
        logging.debug('Sorted full ref array after twotheta')
        nrefl = A.shape[0]
        
        nover=n.zeros((nrefl))
        logging.debug('Ready to compare all %i reflections',nrefl)
        overlaps = dict([(i,[]) for i in range(nrefl)])
        for i in range(1,nrefl):
            if i%1000 == 0:
                logging.debug('Comparing reflection %i', i)
            j=i-1
            while j > -1 and A[i,A_id['tth']]-A[j,A_id['tth']] < dtth :
                if abs(A[i,A_id['omega']]-A[j,A_id['omega']]) \
                        < n.pi/180.0*self.param['omega_step']*self.param['sbox_omega']:
                    peak_distance = n.sqrt((A[i,A_id['detyd']]-A[j,A_id['detyd']])**2+\
                        (A[i,A_id['detzd']]-A[j,A_id['detzd']])**2)
                    if peak_distance < (self.param['sbox_y']+self.param['sbox_z'])/2.0:
                            overlaps[A[i,A_id['spot_id']]].append([A[j,A_id['grain_id']],
                                                                   A[j,A_id['ref_id']]])
                            overlaps[A[j,A_id['spot_id']]].append([A[i,A_id['grain_id']],
                                                                   A[i,A_id['ref_id']]])
                            self.grain[int(A[i,A_id['grain_id']])].refs[A[i,A_id['ref_id']],
                                                                        A_id['overlaps']] += 1
                            self.grain[int(A[j,A_id['grain_id']])].refs[A[j,A_id['ref_id']],
                                                                        A_id['overlaps']] += 1
                j = j - 1
        print 'Number of overlaps %i out of %i refl.' %(n.sum(nover),nrefl)
        co = 0
        # How to find the info for reflection with spot_id
        #refl_with_spotid = A[(A[:,A_id['spot_id']]==spot_id),:]
        
        for i in range(nrefl):
            if len(overlaps[i]) > 0:
                co +=1
                print i, overlaps[i]
        print co

    def save(self,grainno=None):
        """
        write PolyXSim ref file 
        
        """
        file_io.write_ref(self.param,self.grain,grainno)
            
    
    def write_gve(self):
        """
        Write gvector (gve) file, for format see
        http://fable.wiki.sourceforge.net/imaged11+-+file+formats
        
        Henning Osholm Sorensen, RisoeDTU, 2008.
        python translation: Jette Oddershede, Risoe DTU, March 31 2008
        """
        file_io.write_gve(self.param,self.grain,self.hkl)


    def write_ini(self):
        """
        write input file for GrainSpotter
        """
        file_io.write_ini(self.param,self.hkl)

        
    def write_flt(self):
        """
         Write filtered peaks (flt) file, for format see
         http://fable.wiki.sourceforge.net/imaged11+-+file+formats
                 
        """
        file_io.write_flt(self.param,self.grain)
