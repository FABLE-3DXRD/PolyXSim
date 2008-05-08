import numpy as n
import sys
from xfab import tools
import variables
from fabio import edfimage,tifimage
import time
from scipy import ndimage
from scipy.stats import norm

A_id = variables.refarray().A_id

class make_image:
	def __init__(self,graindata):
		self.graindata = graindata

	def make_image(self):
		"""
		makeimage script produces edf diffraction images using the reflection information
		
                Henning Osholm Sorensen, June 23, 2006.
		python translation Jette Oddershede, March 31, 2008
		"""
		if self.graindata.param['make_image'] != 0:
	
			peakshape = self.graindata.param['peakshape']
	
			if peakshape[0] == 0: # spike peak, (2peak_add+1)x(2peak_add+1) pixels
				peak_add = int(round(peakshape[1]))
				frame_add = 0
				peakwsig = 0
			elif peakshape[0] == 1: # 3d Gaussian peak
				peak_add = 0
				frame_add = int(round(peakshape[1]))
				peakwsig = peakshape[2]
			
			framedimy = self.graindata.param['dety_size']+2*frame_add
			framedimz = self.graindata.param['detz_size']+2*frame_add
			
			totalrefl = 0
			no_frames = len(self.graindata.frameinfo)
			print 'Generating ', no_frames, 'frames'
	
			for i in range(no_frames):
				t1 = time.clock()
				nrefl = 0
				frame = n.zeros((framedimy,framedimz))
				omega = self.graindata.frameinfo[i].omega
				omega_step = self.graindata.param['omega_step']
				for j in range(self.graindata.param['no_grains']):  # loop over grains
					for k in range(len(self.graindata.grain[j].refs)): # loop over reflections for each grain
						# exploit that the reflection list is sorted according to omega
						if self.graindata.grain[j].refs[k,A_id['omega']]*180/n.pi > omega+omega_step+2*peakwsig:
							break
						elif self.graindata.grain[j].refs[k,A_id['omega']]*180/n.pi < omega-2*peakwsig:
							continue
						dety = int(round(self.graindata.grain[j].refs[k,A_id['dety']]))
						detz = int(round(self.graindata.grain[j].refs[k,A_id['detz']]))
						yrange = range(dety+frame_add-peak_add,dety+frame_add+peak_add+1)
						zrange = range(detz+frame_add-peak_add,detz+frame_add+peak_add+1)
						intensity = int(round(self.graindata.grain[j].refs[k,A_id['Int']]))
						nrefl = nrefl + 1
						totalrefl = totalrefl + 1
						# Gaussian along omega
						if peakshape[0] == 1:
							fraction = norm.cdf((omega-self.graindata.grain[j].refs[k,A_id['omega']]*180/n.pi+omega_step)/(0.5*peakwsig))\
							          -norm.cdf((omega-self.graindata.grain[j].refs[k,A_id['omega']]*180/n.pi)/(0.5*peakwsig))
						else:
							fraction = 1

						# Generate spikes, possibly more than 1x1 pixel
						for y in yrange:
							for z in zrange:
								if y > 0 and y < framedimy and z > 0 and z < framedimz:
									frame[y,z] = frame[y,z] + fraction*intensity/(len(yrange)*len(zrange))
				# 2D Gaussian on detector					
				if peakshape[0] == 1:
					frame = ndimage.gaussian_filter(frame,peakshape[1]*0.5)
					frame = frame[frame_add:framedimy-frame_add,frame_add:framedimz-frame_add]
				# add background
				if self.graindata.param['bg'] > 0:
					frame = frame + self.graindata.param['bg']*n.ones\
						((self.graindata.param['dety_size'],self.graindata.param['detz_size']))
				# add noise
				if self.graindata.param['noise'] != 0:
					frame = n.random.poisson(frame)
				# apply psf
				if self.graindata.param['psf'] != 0:
					frame = ndimage.gaussian_filter(frame,self.graindata.param['psf']*0.5)
				# convert to integers and flip to same orientation as experimental frames
				frame = n.transpose(n.flipud(n.int16(frame)))
				# Output frames 
				if self.graindata.param['format'] == '.edf':
					self.write_edf(i,frame)
				elif self.graindata.param['format'] == '.tif':
					self.write_tif(i,frame)
				print '\rDone frame %i took %8f s' %(i+1,time.clock()-t1),
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
			%(self.graindata.param['direc'],self.graindata.param['prefix'],self.graindata.param['no_grains'])
		e.write(self.graindata.frameinfo[framenumber].name)
				
	def write_tif(self,framenumber,frame):
		e=tifimage.tifimage()
		e.data=frame
		e.write(self.graindata.frameinfo[framenumber].name)
