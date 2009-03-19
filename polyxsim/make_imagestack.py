import numpy as n

from xfab import tools
from xfab import detector
from fabio import edfimage,tifimage
from scipy import ndimage
import variables,check_input
import generate_grains
import time 
import sys

A_id = variables.refarray().A_id

class make_image:
    def __init__(self,graindata,killfile):
        self.graindata = graindata
        self.killfile = killfile

        # wedge NB! wedge is in degrees
        self.wy = self.graindata.param['wedge']*n.pi/180.
        self.wx = 0.

    def setup_odf(self):
		
            odf_scale = self.graindata.param['odf_scale'] 
            if self.graindata.param['odf_type'] == 1:
                odf_spread = self.graindata.param['mosaicity']/4
                odf_spread_grid = odf_spread/odf_scale
                sigma = odf_spread_grid*n.ones(3)
                r1_max = n.ceil(3*odf_spread_grid)
                r1_range = r1_max*2 + 1
                r2_range = r1_max*2 + 1
                r3_range = r1_max*2 + 1
                mapsize = r1_range*n.ones(3)
                odf_center = r1_max*n.ones(3)
                print 'size of ODF map', mapsize
                self.odf = generate_grains.gen_odf(sigma,odf_center,mapsize)
                #from pylab import *
                #imshow(self.odf[:,:,odf_center[2]])
                #show()
            elif self.graindata.param['odf_type'] == 3:
                odf_spread = self.graindata.param['mosaicity']/4
                odf_spread_grid = odf_spread/odf_scale
                r1_max = n.ceil(3*odf_spread_grid)
                r2_max = n.ceil(3*odf_spread_grid)
                r3_max = n.ceil(3*odf_spread_grid)
                r1_range = r1_max*2 + 1
                r2_range = r2_max*2 + 1
                r3_range = r3_max*2 + 1
                print 'size of ODF map', r1_range*n.ones(3)
                odf_center = r1_max*n.ones(3)
                self.odf= n.zeros((r1_range,r2_range,r3_range))
                # Makes spheric ODF for debug purpuses
                for i in range(self.odf.shape[0]):
                    for j in range(self.odf.shape[1]):
                        for k in range(self.odf.shape[2]):
                            r = [i-(r1_max), j-(r2_max), k-(r3_max)]
                            if n.linalg.norm(r) > r1_max:
                                 self.odf[i,j,k] = 0
                            else:
                                 self.odf[i,j,k] = 1
                #from pylab import *
                #imshow(self.odf[:,:,r3_max],interpolation=None)
                #show()
            elif self.graindata.param['odf_type'] == 2:
                file = self.graindata.param['odf_file']
                print 'Read ODF from file_ %s' %file
                file = open(file,'r')
                (r1_range, r2_range, r3_range) = file.readline()[9:].split()
                r1_range = int(r1_range)
                r2_range = int(r2_range)
                r3_range = int(r3_range)
                odf_scale = float(file.readline()[10:])
                oneD_odf = n.fromstring(file.readline(),sep=' ')
                elements = r1_range*r2_range*r3_range
                self.odf = oneD_odf[:elements].reshape(r1_range,r2_range,r3_range)
                if self.graindata.param['odf_sub_sample'] > 1:
                    sub =self.graindata.param['odf_sub_sample']
                    print 'subscale =',sub
                    r1_range_sub = r1_range * self.graindata.param['odf_sub_sample']
                    r2_range_sub = r2_range * self.graindata.param['odf_sub_sample']
                    r3_range_sub = r3_range * self.graindata.param['odf_sub_sample']
                    odf_fine = n.zeros((r1_range_sub,r2_range_sub,r3_range_sub))
                    for i in range(r1_range):
                        for j in range(r2_range):
                            for k in range(r3_range):
                                odf_fine[i*sub:(i+1)*sub,
                                         j*sub:(j+1)*sub,
                                         k*sub:(k+1)*sub] = self.odf[i,j,k]
                    self.odf = odf_fine.copy()/(sub*sub*sub)
                    r1_range = r1_range_sub 
                    r2_range = r2_range_sub 
                    r3_range = r3_range_sub
                    odf_scale = odf_scale/sub
                    print 'odf_scale', odf_scale

                #[r1_range, r2_range, r3_range] = self.odf.shape
                odf_center = [(r1_range)/2, r2_range/2, r3_range/2]
                print odf_center
                #self.odf[:,:,:] = 0.05
                print self.odf.shape
                #from pylab import *
                #imshow(self.odf[:,:,odf_center[2]])
                #show()
            self.Uodf = n.zeros(r1_range*r2_range*r3_range*9).\
                reshape(r1_range,r2_range,r3_range,3,3)
            if self.graindata.param['odf_cut'] != None:
                self.odf_cut = self.odf.max()*self.graindata.param['odf_cut']
            else:
                self.odf_cut = 0.0
            for i in range(self.odf.shape[0]):
                for j in range(self.odf.shape[1]):
                    for k in range(self.odf.shape[2]):
                        r = odf_scale*n.pi/180.*\
                            n.array([i-odf_center[0],
                                     j-odf_center[1],
                                     k-odf_center[2]])
                        self.Uodf[i,j,k,:,:] = tools.rod_to_u(r)
	    
            if self.graindata.param['odf_type'] !=  2:
                file = open(self.graindata.param['stem']+'.odf','w')
                file.write('ODF size: %i %i %i\n' %(r1_range,r2_range,r3_range))
                file.write('ODF scale: %f\n' %(odf_scale))
                for i in range(int(r1_range)):
                    self.odf[i,:,:].tofile(file,sep=' ',format='%f')
                    file.write(' ')
                file.close()
	    
            return self.Uodf


    def make_image_array(self):
        from scipy import sparse
        #make stack of empty images as a dictionary of sparse matrices
        print 'Build sparse image stack'
        stacksize = len(self.graindata.frameinfo)
        self.frames = {}
        for i in range(stacksize):
            self.frames[i]=sparse.lil_matrix((int(self.graindata.param['dety_size']),
                                              int(self.graindata.param['detz_size'])))


    def make_image(self,grainno=None,refl = None):
        from scipy import ndimage
        if grainno == None:
            do_grains = range(self.graindata.param['no_grains'])
        else:
            do_grains = [grainno]

        # loop over grains
        for grainno in do_grains:
            gr_pos = n.array(self.graindata.param['pos_grains_%s' \
                    %(self.graindata.param['grain_list'][grainno])])
            B = self.graindata.grain[grainno].B
            SU = n.dot(self.graindata.S,self.graindata.grain[grainno].U)
            if refl == None:
                do_refs = range(len(self.graindata.grain[grainno].refs))
            else:
                do_refs = [refl]
            # loop over reflections for each grain

            for nref in do_refs:
			    # exploit that the reflection list is sorted according to omega
                print '\rDoing reflection %i of %i for grain %i of %i' %(nref+1,
                                                                         len(self.graindata.grain[grainno].refs),
                                                                         grainno+1,self.graindata.param['no_grains']),
                sys.stdout.flush()
                #print 'Doing reflection: %i' %nref
                if self.graindata.param['odf_type'] == 3:
                    intensity = 1
                else:
                    intensity = self.graindata.grain[grainno].refs[nref,A_id['Int']]
        
                hkl = n.array([self.graindata.grain[grainno].refs[nref,A_id['h']],
                               self.graindata.grain[grainno].refs[nref,A_id['k']],
                               self.graindata.grain[grainno].refs[nref,A_id['l']]])
                Gc  = n.dot(B,hkl)
                for i in range(self.odf.shape[0]):
                    for j in range(self.odf.shape[1]):
                        for k in range(self.odf.shape[2]):
                            check_input.interrupt(self.killfile)
                            if self.odf[i,j,k] > self.odf_cut:
                                Gtmp = n.dot(self.Uodf[i,j,k],Gc)
                                Gw =  n.dot(SU,Gtmp)
                                Glen = n.sqrt(n.dot(Gw,Gw))
                                tth = 2*n.arcsin(Glen/(2*abs(self.graindata.K)))
                                costth = n.cos(tth)
                                Qw = Gw*self.graindata.param['wavelength']/(4.*n.pi)
                                (Omega, eta) = tools.find_omega_general(Qw,
                                                                        tth,
                                                                        self.wx,
                                                                        self.wy)
                                try:
                                    minpos = n.argmin(n.abs(Omega-self.graindata.grain[grainno].refs[nref,A_id['omega']]))
                                except:
                                    print Omega
                                if len(Omega) == 0:
                                    continue
                                omega = Omega[minpos]
                                # if omega not in rotation range continue to next step
                                if (self.graindata.param['omega_start']*n.pi/180) > omega or\
                                    omega > (self.graindata.param['omega_end']*n.pi/180):
                                    continue
                                Om = tools.form_omega_mat_general(omega,self.wx,self.wy) 
                                Gt = n.dot(Om,Gw)
						
                                # Calc crystal position at present omega
                                [tx,ty,tz]= n.dot(Om,gr_pos) 

                                (dety, detz) = detector.det_coor(Gt, 
                                                                 costth,
                                                                 self.graindata.param['wavelength'],
                                                                 self.graindata.param['distance'],
                                                                 self.graindata.param['y_size'],
                                                                 self.graindata.param['z_size'],
                                                                 self.graindata.param['dety_center'],
                                                                 self.graindata.param['detz_center'],
                                                                 self.graindata.R,
                                                                 tx,ty,tz)


                                if self.graindata.param['spatial'] != None :
                                    # To match the coordinate system of the spline file
                                    # SPLINE(i,j): i = detz; j = (dety_size-1)-dety
                                    # Well at least if the spline file is for frelon2k
                                    (x,y) = detector.detyz_to_xy([dety,detz],
                                                              self.graindata.param['o11'],
                                                              self.graindata.param['o12'],
                                                              self.graindata.param['o21'],
                                                              self.graindata.param['o22'],
                                                              self.graindata.param['dety_size'],
                                                              self.graindata.param['detz_size'])
                                    # Do the spatial distortion
                                    (xd,yd) = self.spatial.distort(x,y)

                                    # transform coordinates back to dety,detz
                                    (dety,detz) = detector.xy_to_detyz([xd,yd],
                                                                    self.graindata.param['o11'],
                                                                    self.graindata.param['o12'],
                                                                    self.graindata.param['o21'],
                                                                    self.graindata.param['o22'],
                                                                    self.graindata.param['dety_size'],
                                                                    self.graindata.param['detz_size'])
								   
                                if dety > -0.5 and dety <= self.graindata.param['dety_size']-0.5 and\
                                    detz > -0.5 and detz <= self.graindata.param['detz_size']-0.5:
                                    dety = int(round(dety))
                                    detz = int(round(detz))
                                    frame_no = int(n.floor((omega*180/n.pi-self.graindata.param['omega_start'])/\
                                                self.graindata.param['omega_step']))
                                    self.frames[frame_no][dety,detz] =  self.frames[frame_no][dety,detz]+ intensity*self.odf[i,j,k]

    def correct_image(self):
        no_frames = len(self.graindata.frameinfo)
        print '\nGenerating ', no_frames, 'frames'
        for frame_no in self.frames:
            t1 = time.clock()

            frame = self.frames[frame_no].toarray()
            if self.graindata.param['bg'] > 0:
                frame = frame + self.graindata.param['bg']*n.ones((self.graindata.param['dety_size'],
								       self.graindata.param['detz_size']))
            # add noise
            if self.graindata.param['noise'] != 0:
                frame = n.random.poisson(frame)
            # apply psf
            if self.graindata.param['psf'] != 0:
                frame = ndimage.gaussian_filter(frame,self.graindata.param['psf']*0.5)
	    # limit values above 16 bit to be 16bit
	    frame = n.clip(frame,0,2**16-1)
	    # convert to integers
	    frame = n.uint16(frame)
	    #flip detector orientation according to input: o11, o12, o21, o22
	    frame = detector.trans_orientation(frame,
	      				 self.graindata.param['o11'],
	      				 self.graindata.param['o12'],
	      				 self.graindata.param['o21'],
	      				 self.graindata.param['o22'],
	      				 'inverse')
	    # Output frames 
	    if '.edf' in self.graindata.param['output']:
	            self.write_edf(frame_no,frame)
	    if '.tif' in self.graindata.param['output']:
	            self.write_tif(frame_no,frame)
	    print '\rDone frame %i took %8f s' %(frame_no+1,time.clock()-t1),
	    sys.stdout.flush()
				
    def write_edf(self,framenumber,frame):
		e=edfimage.edfimage()
		e.data=frame
		e.dim2,e.dim1=frame.shape
		e.header = {}
		e.header['Dim_1']=e.dim1
		e.header['Dim_2']=e.dim2
		e.header['col_end']=e.dim1-1
		e.header['row_end']=e.dim2-1
		e.header['DataType']='UnsignedShort'
		e.header['Image']=1
		e.header['ByteOrder']='Low'
		e.header['time']=time.asctime()
		e.header['Omega']= self.graindata.frameinfo[framenumber].omega +\
		    self.graindata.param['omega_step']/2.0
		e.header['OmegaStep']=self.graindata.param['omega_step']
		e.header['grainfile']='%s/%s_%0.4dgrains.txt' \
			%(self.graindata.param['direc'],self.graindata.param['stem'],self.graindata.param['no_grains'])
		e.write('%s%s' %(self.graindata.frameinfo[framenumber].name,'.edf'))
				
    def write_tif(self,framenumber,frame):
		e=tifimage.tifimage()
		e.data=frame
		e.write('%s%s' %(self.graindata.frameinfo[framenumber].name,'.tif'))


