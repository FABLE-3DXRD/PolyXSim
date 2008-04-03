import numpy as N
from xfab import tools
from xfab import sg
import variables
import sys
from ImageD11 import blobcorrector
import logging
logging.basicConfig(level=logging.DEBUG,format='%(levelname)s %(message)s')

A_id = variables.refarray().A_id


class find_refl:
    def __init__(self,param):
        self.param = param
        self.grain = []
        # determine position of reflections
    
        # Simple transforms of input and set constants
        sintlmin = N.sin(self.param['theta_min']*N.pi/180)/self.param['wavelength']
        sintlmax = N.sin(self.param['theta_max']*N.pi/180)/self.param['wavelength']
        self.K = -2*N.pi/self.param['wavelength']
        self.S = N.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
        self.A0inv = N.array(tools.FormAinv(self.param['unit_cell']))
        self.V = N.array(tools.CellVolume(self.param['unit_cell']))
    #    self.sbox_y_half = (self.param['sbox_y']-1)/2
    #    self.sbox_z_half = (self.param['sbox_z']-1)/2
        # Detector tilt correction matrix
        if self.param['tilt_x'] != 0 or self.param['tilt_y'] != 0 or self.param['tilt_z'] != 0:
             tilt_x = self.param['tilt_x']
             tilt_y = self.param['tilt_y']
             tilt_z = self.param['tilt_z']
             Rx = N.array([[            1,            0,            0],
                           [            0,  cos(tilt_x), -sin(tilt_x)],
                           [            0,  sin(tilt_x),  cos(tilt_x)]])
             Ry = N.array([[  cos(tilt_y),            0,  sin(tilt_y)],
                           [            0,            1,            0],
                           [ -sin(tilt_y),            0,  cos(tilt_y)]])
             Rz = N.array([[  cos(tilt_z), -sin(tilt_z),            0],
                           [  sin(tilt_z),  cos(tilt_z),            0],
                           [            0,            0,            1]])
             self.R = N.dot(Rx,N.dot(Ry,Rz))
        else:
            self.R = [0]
        # Spatial distortion
        if self.param['spatial'] != None:
            self.spatial = blobcorrector.correctorclass(self.param['spatial'])

        # %No of images
        self.nframes = (self.param['omega_end']-self.param['omega_start'])/self.param['omega_step']
        
        # Generate Miller indices for reflections within a certain resolution
        print 'Generating reflections'
        self.hkl  = tools.genhkl(self.param['unit_cell'],self.param['sysconditions'],sintlmin,sintlmax)
        #print self.hkl
        print 'Finished generating reflections\n'
    
    def run(self):
        spot_id = 0
        # Generate orientations of the grains and loop over all grains
        for grainno in range(self.param['no_grains']):
            A = []
            U = self.param['U_grains_%s' %(self.param['grain_list'][grainno])]
            self.grain.append(variables.grain_cont(U))
            gr_pos = N.array(self.param['pos_grains_%s' %(self.param['grain_list'][grainno])])
            gr_eps = N.array(self.param['eps_grains_%s' %(self.param['grain_list'][grainno])])
            B = tools.epsilon2B(gr_eps,self.A0inv)  # Calculate the B-matrix based on the strain tensor for each grain
#            print 'GRAIN NO: ',self.param['grain_list'][grainno]
#            print 'GRAIN POSITION of grain ',self.param['grain_list'][grainno],': ',gr_pos
#            print 'STRAIN TENSOR COMPONENTS (e11 e12 e13 e22 e23 e33) of grain ',self.param['grain_list'][grainno],':\n',gr_eps
#            print 'U of grain ',self.param['grain_list'][grainno],':\n',U
            nrefl = 0
  
            # Calculate these values:
            # totalnr, grainno, refno, hkl, omega, 2theta, eta, dety, detz
            # For all reflections in Ahkl that fulfill omega_start < omega < omega_end.
            # All angles in Grain are in degrees
            for hkl in self.hkl:
                Gtmp = N.dot(B,hkl)
                Gtmp = N.dot(U,Gtmp)
                Gw =   N.dot(self.S,Gtmp)
                #Gw = self.S*U*self.B*hkl
                #print G
                Glen = N.sqrt(N.dot(Gw,Gw))
                theta = N.arcsin(Glen/(2*abs(self.K)))
                costth = N.cos(2*theta)

                Omega = tools.find_omega(Gw,costth)
  
                if len(Omega) > 0:
                    for omega in Omega:
                        if  (self.param['omega_start']*N.pi/180) < omega and omega < (self.param['omega_end']*N.pi/180):
                            Om = N.array([[N.cos(omega), -N.sin(omega), 0],
                                        [N.sin(omega),  N.cos(omega), 0],
                                        [  0       ,    0       , 1]])
                            Gt = N.dot(Om,Gw)
                            eta = N.arctan2(-Gt[1],Gt[2])
                            if eta < 0.0:  # We want eta to be [0,2pi] not [-pi,pi]
                                eta = eta +2*N.pi 
  
                            # Calc crystal position at present omega
                            tx = 0
                            ty = 0
                            tz = 0
                            [tx,ty]= N.dot(Om[:2,:2],gr_pos[:2])
                            tz = gr_pos[2]
                            
                            # NOTE THE + signs here means that the detector coord system adheres to standard
                            if len(self.R) == 3: # If detector tilt calc position like this 
                                R = self.R
                                v = N.array([costth, 
                                             self.param['wavelength']/(2*N.pi)*Gt[1],
                                             self.param['wavelength']/(2*N.pi)*Gt[2]])
                                t = R[0,0]*self.param['distance']/N.sum(R[:,0]*v)
                                Ltv = N.array([tx-self.param['distance'], ty, tz])+ t*v
                                dety = N.sum(R[:,1]*Ltv)/self.param['y_size'] + self.param['dety_center']
                                detz = N.sum(R[:,2]*Ltv)/self.param['z_size'] + self.param['detz_center']
                            else:
                                konst = self.param['wavelength']*(self.param['distance'] - tx)/(2*N.pi*costth)
                                dety = self.param['dety_center'] + (ty + Gt[1]*konst)/self.param['y_size']
                                detz = self.param['detz_center'] + (tz + Gt[2]*konst)/self.param['z_size']

                            if self.param['spatial'] != None :
                                # To match the coordinate system of the spline file
                                # SPLINE(i,j): i = detz; j = (dety_size-1)-dety
                                x = detz 
                                y = self.param['dety_size']-1-dety
                                
                                (xd,yd) = self.spatial.distort(x,y)
                                detyd = self.param['dety_size']-1-yd
                                detzd = x
                            else:
                                detyd = dety
                                detzd = detz
                             #If shoebox extends outside detector exclude it
 #                           if ( self.param['sbox_y'] > detyd) or (detyd > self.param['dety_size']-self.param['sbox_y']) or (self.param['sbox_z'] > detzd) or (detzd > self.param['detz_size']-self.param['sbox_z']):
 #                                continue
                            
 #                           frame_center = N.floor((omega*180/N.pi-self.param['omega_start'])/self.param['omega_step'])
 #                           delta_sbox_omega =  int((self.param['sbox_omega']-1)/2)
 #                           frame_limits = [frame_center - delta_sbox_omega, frame_center + delta_sbox_omega]
                            #Polarization factor (Kahn et al, J. Appl. Cryst. (1982) 15, 330-337.)
                            rho = N.pi/2.0 + eta + self.param['beampol_direct']*N.pi/180.0 
                            P = 0.5 * (1 + costth*costth +\
                                        self.param['beampol_factor']*N.cos(2*rho)*N.sin(2*theta)*N.sin(2*theta))
                            #Lorentz factor
                            if eta != 0:
                                L=1/(N.sin(2*theta)*abs(N.sin(eta)))
                            else:
                                L=N.inf;
 
                            overlaps = 0 # set the number overlaps to zero
                            #logging.debug("frame_center: %i, omega: %f" %(frame_center,omega*180/N.pi))
                            #logging.debug("frame_limits: %i, %i" %(frame_limits[0],frame_limits[1]))
                            A.append([grainno,nrefl,spot_id,
                                      hkl[0],hkl[1],hkl[2],
                                      2*theta,omega,eta,
                                      dety,detz,
                                      detyd,detzd,
                                      Gw[0],Gw[1],Gw[2],
                                      L,P])
 #                                     frame_limits[0],frame_limits[1],
  #                                    overlaps])
                            nrefl = nrefl+1
                            spot_id = spot_id+1

#           print 'Length of Grain', len(self.grain[0].refl)
            A = N.array(A)
            A = A[N.argsort(A,0)[:,A_id['omega']],:] # sort rows according to omega
            A[:,A_id['ref_id']] = N.arange(nrefl)     # Renumber the reflections  
            A[:,A_id['spot_id']] = N.arange(N.min(A[:,A_id['spot_id']]),N.max(A[:,A_id['spot_id']])+1) # Renumber the spot_id
            self.grain[grainno].refs = A
            

    def overlap(self):
        
        dtth = 1*N.pi/180.  # Don't compare position of refs further apart than dtth 

        # build one big array of reflection info of all grains
        A = self.grain[0].refs
        for grainno in range(1,self.param['no_grains']):
            A = N.concatenate((A,self.grain[grainno].refs))
        logging.debug('Finished concatenating ref arrays')
        A = A[N.argsort(A,0)[:,A_id['tth']],:] # sort rows according to tth
        logging.debug('Sorted full ref array after twotheta')
        nrefl = A.shape[0]
        
        nover=N.zeros((nrefl))
        logging.debug('Ready to compare all %i reflections',nrefl)
        overlaps = dict([(i,[]) for i in range(nrefl)])
        for i in range(1,nrefl):
            if i%1000 == 0:
                logging.debug('Comparing reflection %i', i)
            j=i-1
            while j > -1 and A[i,A_id['tth']]-A[j,A_id['tth']] < dtth :
                if abs(A[i,A_id['omega']]-A[j,A_id['omega']]) < N.pi/180.0*self.param['omega_step']*self.param['sbox_omega']:
#                    if abs(A[i,A_id['detyd']]-A[j,A_id['detyd']]) < self.param['sbox_y']:
#                        if abs(A[i,A_id['detzd']]-A[j,A_id['detzd']]) < self.param['sbox_z']:
                    peak_distance = N.sqrt((A[i,A_id['detyd']]-A[j,A_id['detyd']])**2+\
                        (A[i,A_id['detzd']]-A[j,A_id['detzd']])**2)
                    if peak_distance < (self.param['sbox_y']+self.param['sbox_z'])/2.0:
                            overlaps[A[i,A_id['spot_id']]].append([A[j,A_id['grain_id']],A[j,A_id['ref_id']]])
                            overlaps[A[j,A_id['spot_id']]].append([A[i,A_id['grain_id']],A[i,A_id['ref_id']]])
                            self.grain[int(A[i,A_id['grain_id']])].refs[A[i,A_id['ref_id']],A_id['overlaps']] += 1
                            self.grain[int(A[j,A_id['grain_id']])].refs[A[j,A_id['ref_id']],A_id['overlaps']] += 1
                j = j - 1
        print 'Number of overlaps %i out of %i refl.' %(N.sum(nover),nrefl)
        co = 0
        # How to find the info for reflection with spot_id
        #refl_with_spotid = A[(A[:,A_id['spot_id']]==spot_id),:]
        
        for i in range(nrefl):
            if len(overlaps[i]) > 0:
                co +=1
                print i, overlaps[i]
        print co

    def save(self,grainno=None):
        if grainno == None:
            savegrains = range(len(self.grain))
        else:
            savegrains = grainno
        for grainno in savegrains:
            A = self.grain[grainno].refs
            setno = 0
            filename = '%s/%s_gr%0.4d_set%0.4d.ref' \
                %(self.param['direc'],self.param['prefix'],grainno,setno)
            f = open(filename,'w')
            format = "%d "*6 + "%f "*12 + "%d "*3 + "\n"
            ( nrefl, ncol ) = A.shape
#            print nrefl, ncol
            out = "#"
            A_col = dict([[v,k] for k,v in A_id.items()])
            for col in A_col:
                out = out + ' %s' %A_col[col]
            out = out +"\n"

            f.write(out)
            for i in range(nrefl):
                out = format %(A[i,A_id['grain_id']],
                               A[i,A_id['ref_id']],
                               A[i,A_id['spot_id']],   
                               A[i,A_id['h']],
                               A[i,A_id['k']],
                               A[i,A_id['l']],
                               A[i,A_id['tth']]*180/N.pi,
                               A[i,A_id['omega']]*180/N.pi,
                               A[i,A_id['eta']]*180/N.pi,
                               A[i,A_id['dety']],
                               A[i,A_id['detz']],
                               A[i,A_id['detyd']],
                               A[i,A_id['detzd']],
                               A[i,A_id['gv1']],
                               A[i,A_id['gv2']],
                               A[i,A_id['gv3']],
                               A[i,A_id['L']],
                               A[i,A_id['P']],
                               A[i,0],
                               A[i,0],
                               A[i,0]
#                               A[i,A_id['frame_start']],
#                               A[i,A_id['frame_end']],
#                               A[i,A_id['overlaps']]
                           )
                f.write(out)
        
            f.close()   
            
    
    def write_gve(self):
#  Write gvector (gve) file, for format see
# http://fable.wiki.sourceforge.net/imaged11+-+file+formats
# 
#Henning Osholm 2008,
# python translation: Jette Oddershede, Risoe DTU, March 31 2008
#
        filename = '%s/%s.gve' %(self.param['direc'],self.param['prefix'])
        f = open(filename,'w')
        lattice = sg.sg(sgno=self.param['sgno']).name[0]
        format = "%f "*6 + "%s "*1 +"\n"
        out = format %(self.param['unit_cell'][0],self.param['unit_cell'][1],
		               self.param['unit_cell'][2],self.param['unit_cell'][3],
					   self.param['unit_cell'][4],self.param['unit_cell'][5], lattice)
        f.write(out)
        out = "# wavelength = %s\n" %(self.param['wavelength'])
        f.write(out)
        out = "# wedge = 0\n"
        f.write(out)
        out = "# ds h k l\n" 
        f.write(out)
		
        A = self.grain[0].refs
        A = A[N.argsort(-A,0)[:,A_id['tth']],:] # sort rows according to tth, descending
        format = "%f "*1 + "%d "*3 +"\n"
        for i in range(A.shape[0]):
            out = format %(self.param['wavelength']/(2*N.sin(.5*A[i,A_id['tth']])),
                           A[i,A_id['h']],
                           A[i,A_id['k']],
                           A[i,A_id['l']]
                            )
            f.write(out)

        out = "# xr yr zr dety detz ds eta omega\n" 
        f.write(out)
        format = "%f "*8 + "\n"
        for grainno in range(1,self.param['no_grains']):
            A = N.concatenate((A,self.grain[grainno].refs))
			
        nrefl = A.shape[0]
        for i in range(nrefl):
            out = format %(A[i,A_id['gv1']],
			               A[i,A_id['gv2']],
			               A[i,A_id['gv3']],
			               A[i,A_id['dety']],
			               A[i,A_id['detz']],
						   self.param['wavelength']/(2*N.sin(.5*A[i,A_id['tth']])),
			               A[i,A_id['eta']]*180/N.pi,
			               A[i,A_id['omega']]*180/N.pi
                           )
            f.write(out)
		
        f.close()   
            
    
