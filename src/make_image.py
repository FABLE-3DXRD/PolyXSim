import numpy as n
from xfab import tools
import variables
from fabio import edfimage
import time

A_id = variables.refarray().A_id


def make_image(param,frameinfo,graindata):
# makeimage script produces edf diffraction images using the reflection information
#
#Henning Osholm Sorensen, June 23, 2006.
# python translation Jette Oddershede, March 31, 2008
	if param['make_image'] != 0:
	
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
		no_frames = len(frameinfo)
		print 'Generating ', no_frames, 'frames'
	
		for i in range(no_frames):
			nrefl = 0
			frame = n.zeros((param['dety_size'],param['detz_size']))
			omega = frameinfo[i].omega
			for j in range(param['no_grains']):  # loop over grains
				for k in range(len(graindata.grain[j].refs)): # loop over reflections for each grain
					# exploit that the reflection list is sorted according to omega
					if graindata.grain[j].refs[k,A_id['omega']]*180/n.pi > omega+param['omega_step']+5*peakwsig:
						break
					elif graindata.grain[j].refs[k,A_id['omega']]*180/n.pi < omega-5*peakwsig:
						continue
					dety = int(round(graindata.grain[j].refs[k,A_id['dety']]))+1
					detz = int(round(graindata.grain[j].refs[k,A_id['detz']]))+1
					nrefl = nrefl + 1
					totalrefl = totalrefl + 1
					if peakshape == 0:
						for y in range(dety-1,dety+2):
							for z in range(detz-1,detz+2):
								if y > 0 and y < param['dety_size'] and z > 0 and z < param['detz_size']:
									frame[y,z] = frame[y,z] + intensity
									
					elif peakshape == 1:
						if -5*pixellimit >= dety or dety >= param['dety_size']+5*pixellimit or -5*pixellimit >= detz or detz >= param['detz_size']+5*pixellimit:
							continue
						else:
							pass
						
#			print nrefl			
			frame = frame + bgint*n.ones((param['dety_size'],param['detz_size']))
			filename = frameinfo[i].name
			write(filename,frame)
				
def write(fname,frame):
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
#		e.header['Omega']=
#		e.header['OmegaStep']=
		e.write(fname)
				