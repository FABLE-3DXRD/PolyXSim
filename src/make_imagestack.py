import numpy as n
from xfab import tools
import variables
from fabio import edfimage,tifimage
import time
from polyxsim import generate_grains

#A_id = variables.refarray().A_id

class make_image:
	def __init__(self,graindata):
            self.graindata = graindata

        def setup_odf(self):

            if self.graindata.param['odf_type'] == 1:
                odf_scale = 0.02
                odf_spread = self.graindata.param['mosaicity']/4
                odf_spread_grid = odf_spread/odf_scale
                sigma = odf_spread_grid*n.ones(3)
                r1_max = n.ceil(3*odf_spread_grid)
                r1_range = r1_max*2 + 1
                r2_range = r1_max*2 + 1
                r3_range = r1_max*2 + 1
                mapsize = r1_range*n.ones(3)
                odf_center = r1_max*n.ones(3)
                odf = generate_grains.gen_odf(sigma,odf_center,mapsize)
            else:
                [r1_range, r2_range, r3_range] = odf.shape
                odf_center = [(r1_range)/2, r2_range/2, r3_range/2]


            self.Uodf = n.zeros(r1_range*r2_range*r3_range*9).\
		reshape(r1_range,r2_range,r3_range,3,3)
	    print self.Uodf
            for i in range(odf.shape[0]):
                for j in range(odf.shape[1]):
                    for k in range(odf.shape[2]):
                        r = odf_scale*n.pi/180.*\
                            n.array([i-odf_center[0],
                                     j-odf_center[1],
                                     k-odf_center[2]])

                        self.Uodf[i,j,k,:,:] = tools.rod2U(r)

	    return self.Uodf

# % Makes spheric ODF for debug purpuses
# %for i = 1:size(odf,1)
# %    for j = 1:size(odf,2)
# %        for k = 1:size(odf,3)
# %            r = [i-(r1_max+1), j-(r2_max+1), k-(r3_max+1)];
# %            if norm(r) > r1_max
# %               odf(i,j,k) = 0;
# %            else
# %               odf(i,j,k) = 1;
# %            end
# %        end
# %    end
# %end



	def make_image(self):
            from scipy import sparse

	    #make stack of empty images as a dictionary of sparse matrices
            stacksize = len(self.graindata.frameinfo)
	    frames = {}
	    for i=1:stacksize:
		    frames[i]=sparse.lil_matrix(self.graindata.param['dety_size']
						self.graindata.param['detz_size'])


	    # nr = size(A,1);
	    # A = sortrows(A,14);

            # loop over grains
	    for grainno in range(self.graindata.param['no_grains']):
		    gr_pos = n.array(self.graindata.param['pos_grains_%s' \
				     %(self.graindata.param['grain_list'][grainno])])
		    gr_eps = n.array(self.graindata.param['eps_grains_%s' \
                                     %(self.graindata.param['grain_list'][grainno])])
		    B = self.graindata.grain[grainno].B
		    SU = n.dot(self.graindata.S,self.graindata.grain[grainno].U)
		    # loop over reflections for each grain
		    for nref in range(len(self.graindata.grain[grainno].refs)):
			    # exploit that the reflection list is sorted according to omega
			    hkl = n.array([self.graindata.grain[grainno].refs[nref,A_id['h']],
					   self.graindata.grain[grainno].refs[nref,A_id['k']],
					   self.graindata.grain[grainno].refs[nref,A_id['l']]])
			    Gc  = n.dot(B,hkl)
			    for i in range(odf.shape[0]):
				    for j in range(odf.shape[1]):
					    for k in range(odf.shape[2]):
                                                Gtmp = n.dot(self.Uodf[i,j,k],Gc)
						Gw =   n.dot(SU,Gtmp)
						Glen = n.sqrt(n.dot(Gw,Gw))
						tth = 2*n.arcsin(Glen/(2*abs(self.K)))
						costth = n.cos(tth)
						Omega = tools.find_omega(Gw,tth)
						minpos = n.argmin(n.abs(Omega*180/pi-\
						omega = Omega[minpos]
						# if omega not in rotation range continue to next step
						if (self.graindata.param['omega_start']*n.pi/180) > omega or\
						       omega > (self.param['omega_end']*n.pi/180):
						     continue
						Om = tools.OMEGA(Omega[minpos])
						Gt = n.dot(Om,Gw)
						
	                                        # Calc crystal position at present omega
						[tx,ty]= n.dot(Om[:2,:2],gr_pos[:2])
						tz = gr_pos[2]

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


						if self.param['spatial'] != None :
							# To match the coordinate system of the spline file
							# SPLINE(i,j): i = detz; j = (dety_size-1)-dety
							# Well at least if the spline file is for frelon2k
							x = detz 
							y = self.param['dety_size']-1-dety
                                
							(xd,yd) = self.spatial.distort(x,y)
							dety = self.param['dety_size']-1-yd
							detz = x
								   
						if dety > -0.5 & dety <= self.param['dety_size']-0.5 &\
						   detz > -0.5 & detz <= self.param['detz_size']-0.5:
						      dety = int(round(dety))
						      detz = int(round(detz))
						      


# for reflex=1:nr
#     int=A(reflex,21); % Integrated intensity
#     calc_pos_detector
#     print 'Reflex = %i,  Grain %i, hkl: %i %i %i' %(reflex,A(reflex,2),A(reflex,4),A(reflex,5),A(reflex,6)


# % Finish up and save frame
# for i=1:stacksize
#     fullframe = full(frames{i}); %make frame non-sparse matrix
#     % Add background counts
#     fullframe = fullframe + bgint*ones(detzsize,detysize);

#     if psf ~= 0 %Apply point-spread-function
#         fullframe = conv2(fullframe,psf_map,'same');
#     end
    
#     if addnoise == 1
#         framenoise = noise(fullframe,'poisson',1);
#         fullframe = dip_array(framenoise);
#     end

#     saveframe(direc,prefix, i, fullframe);
# end
