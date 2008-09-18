import numpy as n
from xfab import tools
from xfab import sg
from xfab import detector
import variables


def generate_U(no_grains):
	# generate random U (orientations) for no_grains grains
	#
	# INPUT: no_grains
	# OUTPUT: [[U11 U12 U13],[U21 U22 U23],[U31 U32 U33]] for each grain
	#
	# Henning Osholm Sorensen, Risoe DTU,
	# Adapted by Jette Oddershede, RISOE DTU, March 27, 2008
	
	print 'Grain orientations will be randomly generated\n'
	U = n.zeros((no_grains,3,3))
	for i in range(no_grains):
		phi1 = n.random.rand()*2*n.pi
		phi2 = n.random.rand()*2*n.pi
		PHI  = n.random.rand()*n.pi
		U[i] = tools.euler2U(phi1,PHI,phi2)
		
	return U
		
def generate_pos(no_grains,gen_pos,sample_xyz=None,sample_cyl=None):
	# generate random positions for no_grains grains for a given sample shape, 
	# EITHER a box (sample_xyz != None) OR a cylinder (sample_cyl != None)
	# If no shape is given or gen_pos == 0 default to (0,0,0)
	#
	# INPUT: 	no_grains
	#              	gen_pos (==0 default to (0,0,0), else random positions)
	#		sample_xyz (optional) sample dimensions along x, y and z
	#		sample_cyl (optional) sample cylinder radius and length (along z)
	# OUTPUT: [x y z] for each grain
	#
	# Jette Oddershede, RISOE DTU, March 27, 2008
	
	pos = n.zeros((no_grains,3))
	if sample_xyz != None and gen_pos != 0:
		print 'Grain positions will be randomly generated within a box of ', sample_xyz[0], sample_xyz[1], sample_xyz[2], 'mm\n'
		for i in range(no_grains):
			for j in range(3):
				pos[i,j] = (n.random.rand()-0.5)*sample_xyz[j]

	elif sample_cyl != None and gen_pos != 0:
		print 'Grain positions will be randomly generated'
		print 'within a cylinder of radius ', sample_cyl[0], ' and length ', sample_cyl[1], 'mm\n'
		for i in range(no_grains):
			r = n.random.rand()*sample_cyl[0]
			w = n.random.rand()*2*n.pi
			z = (n.random.rand()-0.5)*sample_cyl[1]
			pos[i] = [r*n.cos(w), r*n.sin(w), z]

	else:
		print 'Grain positions will be set to (0,0,0)\n'

	return pos

def generate_eps(no_grains,gen_eps):
	# generate strain tensor components for no_grains grains using a normal distribution with 
	# mean and spread for diagonal and off-diagonal elements as specified by gen_eps
	#
	# INPUT: no_grains, gen_eps =[mean_diag,spread_diag,mean_offdiag,spread_offdiag]
	# OUTPUT: [eps11, eps12, eps13, eps22, eps23, eps33] for each grain
	#
	# Jette Oddershede, RISOE DTU, March 27 2008
	
	eps = n.zeros((no_grains,6))
	diag = [0,3,5]
	offdiag = [1,2,4]
	
	if gen_eps == [0,0,0,0]:
		print 'No strain\n'
	else:
		print 'Grain strain tensors will be randomly generated using a normal distribution'
		print 'For diagonal elements, mean: ', gen_eps[0], ' and spread: ', gen_eps[1]
		print 'For off-diagonal elements, mean: ', gen_eps[2], ' and spread: ', gen_eps[3],'\n'

		for i in range(no_grains):
			if gen_eps[1] == 0:   # non-zero spread of diagonal elements
				for j in diag:
					eps[i,j] = gen_eps[0]
				
			else:
				for j in diag:
					eps[i,j] = n.random.normal(gen_eps[0],gen_eps[1])

			if gen_eps[3] == 0:   # non-zero spread of off-diagonal elements
				for j in offdiag:
					eps[i,j] = gen_eps[2]
				
			else:
				for j in offdiag:
					eps[i,j] = n.random.normal(gen_eps[2],gen_eps[3])
			
	return eps

def grain_size(no_grains,grain_size,grain_min_max,sample_vol=None):
	# generate grain sizes (mean diameters) from a lognormal distribution
	# with a specified mean (grain_size) and upper and lower cut-offs given by grain_min_max
	# For specified sample shape and size calculate the space filling fraction
	#
	# INPUT: no_grains, gen_size (mean size), 
	#               grain_min_max (optional)[min,max]
	#               sample_xyz (optional) sample volume 
	# OUTPUT: mean diameter (in mm) for each grain
	#
	# Jette Oddershede, RISOE DTU, April 2 2008
	
	size = n.zeros((no_grains))
	grain_vol = 0
	if grain_size < 0:
		print 'All grains have the following mean diameter: ', abs(grain_size), ' mm\n'
		for i in range(no_grains):
			size[i] = abs(grain_size)
			grain_vol = grain_vol + n.pi/6*size[i]**3
			
	else:
	# NB the standard lognormal distribution
        # lognormal() = exp(normal(mean=0, spread=1) has a mean of exp(.5)
		print 'Grain sizes from a lognormal distribution with a mean diameter of: ',\
		abs(grain_size), ' mm\n'
		for i in range(no_grains):
			while size[i] <= grain_min_max[0] or size[i] > grain_min_max[1]:
				size[i] = grain_size*n.random.lognormal()/n.exp(.5)

			grain_vol = grain_vol + n.pi/6*size[i]**3
	
	if sample_vol != None:
		fraction = grain_vol/sample_vol
		print 'The generated grains cover the following fraction of the sample volume: %6f' %fraction
	return size
	

def generate_grains(param):
# Generate U if gen_U on
	if param['gen_U'] != 0: 
		U = generate_U(param['no_grains'])
		for i in range(param['no_grains']):
			param['U_grains_%s' %(param['grain_list'][i])] = U[i]
        param['gen_U'] = 0
			
# Generate pos if gen_pos on
	if param['gen_pos'][0] != 0: 
		pos = generate_pos(param['no_grains'],param['gen_pos'][1],sample_xyz=param['sample_xyz'],sample_cyl=param['sample_cyl'])
		for i in range(param['no_grains']):
			param['pos_grains_%s' %(param['grain_list'][i])] = pos[i]
        param['gen_pos'][0] = 0

	
# If pick out grain for phase
	if param['gen_phase'][0] != 0:
		#Making random list of grain numbers
		rand_grain_order = n.argsort(n.random.randn(param['no_grains']))

		no_picked_grains = 0
		for i in range(param['no_phases']):
			phase = param['gen_phase'][i*2+1]
			no_grains_phase = int(param['gen_phase'][i*2+2])
			param['no_grains_phase_%i' %phase] = no_grains_phase
			param['grain_list_phase_%i' %phase] = \
			    n.array(param['grain_list'])[rand_grain_order[no_picked_grains:no_picked_grains+no_grains_phase]]
			no_picked_grains += no_grains_phase
			for grain in param['grain_list_phase_%i' %phase]:
				param['phase_grains_%i' %grain] = phase


	param['gen_phase'][0] = 0  # Done

# Generate eps if gen_eps on
	if param['gen_eps'][0] != 0: 
		eps = generate_eps(param['no_grains'],param['gen_eps'][1:5])
		for i in range(param['no_grains']):
			param['eps_grains_%s' %(param['grain_list'][i])] = eps[i]
        param['gen_eps'][0] = 0


# Generate size if gen_size on
	for phase in param['phase_list']:
		if param['gen_size_phase_%i' %phase][0] != 0: 
			size = grain_size(param['no_grains_phase_%i' %phase],
					  param['gen_size_phase_%i' %phase][1],
					  param['gen_size_phase_%i' %phase][2:4],
					  param['sample_vol']*param['vol_frac_phase_%i' %phase])
			for i in range(param['no_grains_phase_%i' %phase]):
				param['size_grains_%s' %(param['grain_list_phase_%i' %phase][i])] = size[i]
			param['gen_size_phase_%i' %phase][0] = 0
	
	
	
def save_grains(param):
#  Save the generated grain parameters, pos, U and eps
#
# INPUT: The parameter set from the input file and the grain generator
# OUTPUT: grainno x y z phi1 PHI phi2 U11 U12 U13 U21 U22 U23 U31 U32 U33 eps11 eps12 eps13 eps22 eps23 eps33
#
# Jette Oddershede, Risoe DTU, March 31 2008
#

    filename = '%s/%s_%0.4dgrains.txt' %(param['direc'],param['stem'],param['no_grains'])
    f = open(filename,'w')
#    format = "%d "*1 + "%f "*1 + "%e"*1 + "%f"*18 + "\n"
    format = "%d "*1 + "%f "*1 + "%e "*1 + "%f "*21 + "\n"
    out = "# grainno grainsize grainvolume x y z phi1 PHI phi2 U11 U12 U13 U21 U22 U23 U31 U32 U33 eps11 eps12 eps13 eps22 eps23 eps33 \n"
    f.write(out)
    for i in range(param['no_grains']):
        euler = 180/n.pi*tools.U2euler(param['U_grains_%s' %(param['grain_list'][i])])
        out = format %(param['grain_list'][i],
                       param['size_grains_%s' %(param['grain_list'][i])],
                       n.pi/6*(param['size_grains_%s' %(param['grain_list'][i])])**3.,
                       param['pos_grains_%s' %(param['grain_list'][i])][0],
                       param['pos_grains_%s' %(param['grain_list'][i])][1],
                       param['pos_grains_%s' %(param['grain_list'][i])][2],
                       euler[0],
                       euler[1],
                       euler[2],
                       param['U_grains_%s' %(param['grain_list'][i])][0,0],
                       param['U_grains_%s' %(param['grain_list'][i])][0,1],
                       param['U_grains_%s' %(param['grain_list'][i])][0,2],
                       param['U_grains_%s' %(param['grain_list'][i])][1,0],
                       param['U_grains_%s' %(param['grain_list'][i])][1,1],
                       param['U_grains_%s' %(param['grain_list'][i])][1,2],
                       param['U_grains_%s' %(param['grain_list'][i])][2,0],
                       param['U_grains_%s' %(param['grain_list'][i])][2,1],
                       param['U_grains_%s' %(param['grain_list'][i])][2,2],
                       param['eps_grains_%s' %(param['grain_list'][i])][0],
                       param['eps_grains_%s' %(param['grain_list'][i])][1],
                       param['eps_grains_%s' %(param['grain_list'][i])][2],
                       param['eps_grains_%s' %(param['grain_list'][i])][3],
                       param['eps_grains_%s' %(param['grain_list'][i])][4],
                       param['eps_grains_%s' %(param['grain_list'][i])][5],
                           )
        f.write(out)
    f.close()   
            

def write_ubi(param):
#  Save the generated UBI's
#
# INPUT: The parameter set from the input file and the grain generator
# OUTPUT: for each grain the 3x3 UBI matrix
#
# Jette Oddershede, Risoe DTU, April 4 2008
#
    filename = '%s/%s.ubi' %(param['direc'],param['stem'])
    f = open(filename,'w')
    format = "%f "*3 + "\n"
    for i in range(param['no_grains']):
        U = param['U_grains_%s' %(param['grain_list'][i])]
        gr_eps = n.array(param['eps_grains_%s' %(param['grain_list'][i])])
	if param['no_phases'] == 1:
            phase = param['phase_list'][0]
	else:
            phase = param['phase_grains_%s' %(param['grain_list'][i])]
	# Calculate the B-matrix based on the strain tensor for each grain
        B = tools.epsilon2B(gr_eps,param['unit_cell_phase_%i' %phase])/(2*n.pi) 
        UBI = n.linalg.inv(n.dot(U,B))
        for j in range(3):
            out = format %(UBI[j,0],UBI[j,1],UBI[j,2])
            f.write(out)
			
        out = "\n"
        f.write(out)
    f.close()   



	
def write_res(param):
    """
    Save the generated grain parameters in an input file to facilitate restart of PolyXSim with same grains

    INPUT: PolyXSim input and generated grain parameters
    OUTPUT: .res file in PolyXSim input format

    Jette Oddershede, Risoe DTU, June 18 2008
    """
    filename = '%s/%s.res' %(param['direc'],param['stem'])
    f = open(filename,'w')
			
    #initialise and sort keys alphabetically
    out = "" 
    keys = param.keys()
    keys.sort()

    for item in keys:
        # rule out None entries
        if param[item] != None:
            # treat all strings, remember quotation marks
            if type(param[item]) == str:
                out += "%s '%s'\n" %(item,param[item])
            # treat all lists, special case for strings
            elif type(param[item]) == list:
                out += '%s' %item
                for i in range(len(param[item])):
                    if type(param[item][i]) == str:
                        out += " '%s'" %param[item][i]
                    else:
                        out += ' %s' %param[item][i]
                out += '\n'  
            # treat all arrays, loop over one or two dimensions
            elif type(param[item]) == n.ndarray: 
                out += '%s' %item
                dim = len(n.shape(param[item]))
                if dim == 1:
                    for i in range(len(param[item])):
                        out += ' %s' %param[item][i]
                elif dim == 2:
                    for i in range(len(param[item])):
                        for j in range(len(param[item][i])):
                            out += ' %s' %param[item][i][j]
                out += '\n'    
            # remaining entries; integers and floats
            else:
                out += "%s %s\n" %(item,param[item])

    f.write(out)
    f.close()   
	
def write_par(param):
    """
    Save the generated UBI's

    INPUT: The detector info used for the simulations
    OUTPUT: The corresponding detector.par file for ImageD11

    Jette Oddershede, Risoe DTU, June 17 2008
    """



    #Prepare detector part of p.par output

    #Calc beam center in ImageD11 coordinate system 
    (z_center, y_center) = detector.detyz2xy([param['dety_center'],param['detz_center']],
					     param['o11'],
					     param['o12'],
					     param['o21'],
					     param['o22'],
					     param['dety_size'],
					     param['detz_size'])
			

    dout = "chi 0.0\n" 
    dout = dout + "distance %f\n" %(param['distance']*1000.) 
    dout = dout + "fit_tolerance 0.5\n" 
    dout = dout + "o11 %i\n" %param['o11']
    dout = dout + "o12 %i\n" %param['o12']
    dout = dout + "o21 %i\n" %param['o21']
    dout = dout + "o22 %i\n" %param['o22']
    dout = dout + "omegasign %f\n" %param['omega_sign']
    dout = dout + "t_x 0\n" 
    dout = dout + "t_y 0\n" 
    dout = dout + "t_z 0\n" 
    dout = dout + "tilt_x %f\n" %param['tilt_x']
    dout = dout + "tilt_y %f\n" %param['tilt_y']
    dout = dout + "tilt_z %f\n" %param['tilt_z']
    dout = dout + "wavelength %f\n" %param['wavelength']
    dout = dout + "wedge 0.0\n"
    dout = dout + "y_center %f\n" %y_center
    dout = dout + "y_size %f\n" %(param['y_size']*1000.)
    dout = dout + "z_center %f\n" %z_center
    dout = dout + "z_size %f\n" %(param['z_size']*1000.)



    for phase in param['phase_list']:
	    if param['no_phases'] > 1:
                filename = '%s/%s_phase_%i.par' %(param['direc'],param['stem'],phase)
            else:
                filename = '%s/%s.par' %(param['direc'],param['stem'])
	    f = open(filename,'w')
	    
	    unit_cell = param['unit_cell_phase_%i' %phase]
	    out = "cell__a %s\n" %unit_cell[0]
	    out = out + "cell__b %s\n" %unit_cell[1]
	    out = out + "cell__c %s\n" %unit_cell[2]
	    out = out + "cell_alpha %s\n" %unit_cell[3]
	    out = out + "cell_beta %s\n" %unit_cell[4]
	    out = out + "cell_gamma %s\n" %unit_cell[5]	
	    out = out + "cell_lattice_[P,A,B,C,I,F,R] %s\n" %param['sgname_phase_%i' %phase][0]
	    out = out + dout

	    f.write(out)
	    f.close()   



def gen_odf(sigma,pos,mapsize):
    """
    Generate a 3D orietation distributio function as a 3D Gaussion
    Henning Osholm Sorensen, Risoe-DTU, June 9, 2008
    """
    odf = n.zeros(mapsize)
    
    for i in range(int(mapsize[0])):
	for j in range(int(mapsize[1])):
            for k in range(int(mapsize[2])):
                odf[i,j,k] = -0.5*(((i-pos[0])/sigma[0])**2\
			          +((j-pos[1])/sigma[1])**2\
                                  +((k-pos[2])/sigma[2])**2)

    odf = 1/(2*n.pi*n.sqrt(2*n.pi)*sigma[0]*sigma[1]*sigma[2])*n.exp(odf)
    return odf
