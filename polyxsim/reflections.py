import numpy as n
from xfab import tools,structure,sg
import check_input

def gen_miller(param,phase):
	"""
	Generate set of miller indices.
	Henning Osholm Sorensen, Risoe DTU,
	"""

        sintlmin = n.sin(param['theta_min']*n.pi/180)/param['wavelength']
        sintlmax = n.sin(param['theta_max']*n.pi/180)/param['wavelength']

	hkl  = tools.genhkl(param['unit_cell_phase_%i' %phase],
				 sg.sg(sgno=param['sgno_phase_%i' %phase]).syscond,
				 sintlmin,
				 sintlmax)

	return hkl
		
def open_structure(param,phase):
	file = param['structure_phase_%i' %phase]
	if file[-3:] == 'cif':
		if ('structure_datablock_phase_%i' %phase) in param:
			datablock = param['structure_datablock_phase_%i' %phase]
		else:
			datablock = None
		struct = structure.build_atomlist()
		struct.CIFread(ciffile=file,cifblkname=datablock)
	elif file[-3:] == 'pdb':
		struct = structure.build_atomlist()
		struct.PDBread(pdbfile=file)
	else:
		raise IOError, 'Unknown structure file format'
	param['sgno_phase_%i' %phase] = sg.sg(sgname=struct.atomlist.sgname).no
	param['sgname_phase_%i' %phase] = struct.atomlist.sgname
	param['unit_cell_phase_%i' %phase] =  struct.atomlist.cell
	return struct

def calc_intensity(hkl,struct,killfile=None):
	"""
	Calculate the reflection intensities
        """
	int = n.zeros((len(hkl),1))
	for i in range(len(hkl)):
		check_input.interrupt(killfile)
		(Fr, Fi) = structure.StructureFactor(hkl[i],
						     struct.atomlist.cell,
						     struct.atomlist.sgname,
						     struct.atomlist.atom,
						     struct.atomlist.dispersion)
		int[i] = Fr**2 + Fi**2			
	hkl = n.concatenate((hkl,int),1)
	return hkl

def add_intensity(hkl,param):
	"""
	Calculate the reflection intensities
        """
	if 'structure_int' in param:
		int = param['structure_int']
	else:
		int = 2**15

	int = n.ones((len(hkl),1))*int
	hkl = n.concatenate((hkl,int),1)
	return hkl

