import numpy as n
from xfab import tools
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
	
def generate_grains(param):
# Generate U if gen_U exists
	if param['gen_U'] != None: 
		U = generate_U(param['no_grains'])
		for i in range(param['no_grains']):
			param['U_grains_%s' %(param['grain_list'][i])] = U[i]
			
# Generate pos if gen_pos exists
	if param['gen_pos'] != None: 
		pos = generate_pos(param['no_grains'],param['gen_pos'],param['sample_xyz'],param['sample_cyl'])
		for i in range(param['no_grains']):
			param['pos_grains_%s' %(param['grain_list'][i])] = pos[i]

# Generate eps if gen_eps exists
	if param['gen_eps'] != None: 
		eps = generate_eps(param['no_grains'],param['gen_eps'])
		for i in range(param['no_grains']):
			param['eps_grains_%s' %(param['grain_list'][i])] = eps[i]
			
	return param
	
	
def save_grains(param):
#  Save the generated grain parameters, pos, U and eps
#
# INPUT: The parameter set from the input file and the grain generator
# OUTPUT: grainno x y z U11 U12 U13 U21 U22 U23 U31 U32 U33 eps11 eps12 eps13 eps22 eps23 eps33
#
# Jette Oddershede, Risoe DTU, March 31 2008
#

    filename = '%s/%s_%0.4dgrains.txt' %(param['direc'],param['stem'],param['no_grains'])
    f = open(filename,'w')
    format = "%d "*1 + "%f "*18 +"\n"
    out = "# grainno x y z U11 U12 U13 U21 U22 U23 U31 U32 U33 eps11 eps12 eps13 eps22 eps23 eps33 \n"
    f.write(out)
    for i in range(param['no_grains']):
        out = format %(param['grain_list'][i],
                       param['pos_grains_%s' %(param['grain_list'][i])][0],
                       param['pos_grains_%s' %(param['grain_list'][i])][1],
                       param['pos_grains_%s' %(param['grain_list'][i])][2],
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
            


