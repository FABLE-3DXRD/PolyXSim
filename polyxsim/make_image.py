import numpy as n
import sys
from xfab import tools,detector

import variables,check_input
from fabio import edfimage,tifimage
import time
from scipy import ndimage
from scipy.stats import norm

from PIL import Image

A_id = variables.refarray().A_id

class make_image:
	def __init__(self,graindata,options):
		self.graindata = graindata
		self.killfile = options.killfile

	def make_image(self):
		"""
		makeimage script produces edf diffraction images using the reflection information
		
                Henning Osholm Sorensen, June 23, 2006.
		python translation Jette Oddershede, March 31, 2008
		"""
		if self.graindata.param['make_image'] != 0:
	
			peakshape = self.graindata.param['peakshape']
	
			if peakshape[0] == 0: # spike peak, 2x2 pixels
				peak_add = 1
				frame_add = 1
				peakwsig = 0
			elif peakshape[0] == 1: # 3d Gaussian peak
				peak_add = max(1,int(round(peakshape[1])))
				frame_add = max(1,int(round(peakshape[1])))
				peakwsig = peakshape[2]
			
			framedimy = self.graindata.param['dety_size']+2*frame_add
			framedimz = self.graindata.param['detz_size']+2*frame_add
			
			totalrefl = 0
			no_frames = len(self.graindata.frameinfo)
			print 'Generating ', no_frames, 'frames'
	
			for i in range(no_frames):
				check_input.interrupt(self.killfile)
				t1 = time.clock()
				nrefl = 0
				frame = n.zeros((framedimy,framedimz))
				omega = self.graindata.frameinfo[i].omega
				omega_step = self.graindata.param['omega_step']
				# loop over grains
				for j in range(self.graindata.param['no_grains']):
					# loop over reflections for each grain
					for k in range(len(self.graindata.grain[j].refs)):
						# exploit that the reflection list is sorted according to omega
						if self.graindata.grain[j].refs[k,A_id['omega']]*180/n.pi > \
							    omega+omega_step+2*peakwsig:
							break
						elif self.graindata.grain[j].refs[k,A_id['omega']]*180/n.pi < \
							    omega-2*peakwsig:
							continue
						dety = self.graindata.grain[j].refs[k,A_id['dety']]
						detz = self.graindata.grain[j].refs[k,A_id['detz']]
						ndety = int(round(dety))
						ndetz = int(round(detz))
						yrange = range(ndety+frame_add-peak_add,ndety+frame_add+peak_add+1)
						zrange = range(ndetz+frame_add-peak_add,ndetz+frame_add+peak_add+1)
						intensity = int(round(self.graindata.grain[j].refs[k,A_id['Int']]))
						nrefl = nrefl + 1
						totalrefl = totalrefl + 1
						# Gaussian along omega
						if peakshape[0] == 1:
							fraction = norm.cdf((omega-self.graindata.grain[j].refs[k,A_id['omega']]*180/n.pi+omega_step)/(0.5*peakwsig))\
									  -norm.cdf((omega-self.graindata.grain[j].refs[k,A_id['omega']]*180/n.pi)/(0.5*peakwsig))
						else:
							fraction = 1.                       
						# Generate spikes, 2x2 pixels
						for y in yrange:
							for z in zrange:
								if y > 0 and y < framedimy and z > 0 and z < framedimz and abs(dety+frame_add-y) < 1 and abs(detz+frame_add-z) < 1:
#									frame[y-1,z] = frame[y-1,z] + fraction*intensity*(1-abs(dety+frame_add-y))*(1-abs(detz+frame_add-z))
									frame[y,z] = frame[y,z] + fraction*intensity*(1-abs(dety+frame_add-y))*(1-abs(detz+frame_add-z))

				# 2D Gaussian on detector					
				if peakshape[0] == 1:
					frame = ndimage.gaussian_filter(frame,peakshape[1]*0.5)
				# add background
				if self.graindata.param['bg'] > 0:
					frame = frame + self.graindata.param['bg']*n.ones((framedimy,framedimz))
				# add noise
				if self.graindata.param['noise'] != 0:
					frame = n.random.poisson(frame)
				# apply psf
				if self.graindata.param['psf'] != 0:
					frame = ndimage.gaussian_filter(frame,self.graindata.param['psf']*0.5)
				# resize, convert to integers and flip to same orientation as experimental frames
				frame = frame[frame_add:framedimy-frame_add,frame_add:framedimz-frame_add]

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
					self.write_edf(i,frame)
				if '.tif' in self.graindata.param['output']:
					self.write_tif(i,frame)
				if '.tif16bit' in self.graindata.param['output']:
					self.write_tif16bit(i,frame)
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
			%(self.graindata.param['direc'],self.graindata.param['stem'],self.graindata.param['no_grains'])
		e.write('%s%s' %(self.graindata.frameinfo[framenumber].name,'.edf'))
				
	def write_tif(self,framenumber,frame):
		e=tifimage.tifimage()
		e.data=frame
		e.write('%s%s' %(self.graindata.frameinfo[framenumber].name,'.tif'))
	def write_tif16bit(self,framenumber,frame):
		size = frame.shape[:2][::-1]
		pilimage = Image.frombuffer('I',size,frame.tostring(),"raw",'I;16',0,1)
		pilimage.save('%s%s' %(self.graindata.frameinfo[framenumber].name,'.tif'))



