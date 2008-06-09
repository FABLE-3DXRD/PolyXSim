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


            Uodf = n.zeros(r1_range*r2_range*r3_range*9).\
		reshape(r1_range,r2_range,r3_range,3,3)
	    print Uodf
            for i in range(odf.shape[0]):
                for j in range(odf.shape[1]):
                    for k in range(odf.shape[2]):
                        r = odf_scale*n.pi/180.*\
                            n.array([i-odf_center[0],
                                     j-odf_center[1],
                                     k-odf_center[2]])

                        Uodf[i,j,k,:,:] = tools.rod2U(r)

	    return Uodf

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
		    gr_pos = n.array(self.grainno.param['pos_grains_%s' %(self.param['grain_list'][grainno])])
		    gr_eps = n.array(self.grainno.param['eps_grains_%s' %(self.param['grain_list'][grainno])])
		    B = tools.epsilon2B(gr_eps,self.param['unit_cell'])
		    # loop over reflections for each grain
		    for k in range(len(self.graindata.grain[j].refs)):
			    # exploit that the reflection list is sorted according to omega
			    hkl = self.graindata.grain[j].refs[k,A_id['omega']]
			    Gtmp = n.dot(B,hkl)
			    Gtmp = n.dot(U,Gtmp)
			    Gw =   n.dot(self.S,Gtmp)

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
