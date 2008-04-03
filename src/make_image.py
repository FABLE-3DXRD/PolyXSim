import numpy as n
import sys
from xfab import tools
import variables
from fabio import edfimage,tifimage
import time

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
	
			peakshape = 0
			intensity = 1000
			bgint = 10
	
			if peakshape == 0: # spike peak, 3x3 pixels 
				peakwsig = 0
				pixellimit = 0
			elif peakshape == 1: # 3d Gaussian peak
				peakwsig = 0.2
				pixellimit = 2
	
			totalrefl = 0
			no_frames = len(self.graindata.frameinfo)
			print 'Generating ', no_frames, 'frames'
	
			for i in range(no_frames):
				t1 = time.clock()
				nrefl = 0
				frame = n.zeros((self.graindata.param['dety_size'],self.graindata.param['detz_size']))
				omega = self.graindata.frameinfo[i].omega
				for j in range(self.graindata.param['no_grains']):  # loop over grains
					for k in range(len(self.graindata.grain[j].refs)): # loop over reflections for each grain
						# exploit that the reflection list is sorted according to omega
						if self.graindata.grain[j].refs[k,A_id['omega']]*180/n.pi > omega+self.graindata.param['omega_step']+5*peakwsig:
							break
						elif self.graindata.grain[j].refs[k,A_id['omega']]*180/n.pi < omega-5*peakwsig:
							continue
						dety = int(round(self.graindata.grain[j].refs[k,A_id['dety']]))+1
						detz = int(round(self.graindata.grain[j].refs[k,A_id['detz']]))+1
						nrefl = nrefl + 1
						totalrefl = totalrefl + 1
						if peakshape == 0:
							for y in range(dety-1,dety+2):
								for z in range(detz-1,detz+2):
									if y > 0 and y < self.graindata.param['dety_size'] and z > 0 and z < self.graindata.param['detz_size']:
										frame[y,z] = frame[y,z] + intensity
										
						elif peakshape == 1:
							if -5*pixellimit >= dety or\
							   dety >= self.graindata.param['dety_size']+5*pixellimit or\
							   -5*pixellimit >= detz or\
							   detz >= self.graindata.param['detz_size']+5*pixellimit:
								continue
							else:
								pass
						
				# print nrefl			
				frame = frame + bgint*n.ones\
				    ((self.graindata.param['dety_size'],self.graindata.param['detz_size']))
				if self.graindata.param['noise'] == True:
					frame = n.random.poisson(frame)

				#
				frame = n.transpose(n.flipud(frame))
				# Output frames 
				if self.graindata.param['format'] == '.edf':
					self.write_edf(i,frame)
				elif self.graindata.param['format'] == '.tif':
					self.write_tif(i,frame)
				print '\rDone frame %i took %8f' %(i+1,time.clock()-t1),
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
		e.write(self.graindata.frameinfo[framenumber].name)
				
	def write_tif(self,framenumber,frame):
		e=tifimage.tifimage()
		e.data=frame
		e.write(self.graindata.frameinfo[framenumber].name)
