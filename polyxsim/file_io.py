import numpy as n
from xfab import tools,detector
import variables

A_id = variables.refarray().A_id

def write_flt(param,grain):
    """
     Write filtered peaks (flt) file, for format see
     http://fable.wiki.sourceforge.net/imaged11+-+file+formats
     
     python translation: Jette Oddershede, Risoe DTU, June 4 2008
    
    """
    filename = '%s/%s.flt' %(param['direc'],param['stem'])
    f = open(filename,'w')
    out = '#  sc  fc  omega  Number_of_pixels  avg_intensity  s_raw  f_raw  sigs  sigf  covsf  sigo  covso  covfo  sum_intensity  sum_intensity^2  IMax_int  IMax_s  IMax_f  IMax_o  Min_s  Max_s  Min_f  Max_f  Min_o  Max_o  dety  detz  onfirst  onlast  spot3d_id \n'
    f.write(out)
    	
    A = grain[0].refs
    for grainno in range(1,param['no_grains']):
        A = n.concatenate((A,grain[grainno].refs))
    A = A[n.argsort(A,0)[:,A_id['omega']],:] # sort rows according to omega
    format = "%f "*3 + "%i "*1 +"%f "*12 + "%i "*2   +"%f "*1 + "%i "*4 +"%f "*4 + "%i "*3 +"\n"

    for i in range(A.shape[0]):
        (sc, fc) = detector.detyz2xy([A[i,A_id['dety']],A[i,A_id['detz']]],
                                     param['o11'],
                                     param['o12'],
                                     param['o21'],
                                     param['o22'],
                                     param['dety_size'],
                                     param['detz_size'])
        if param['spatial'] == None:
            sr = sc
            fr = fc
        else:
            (sr, fr) = detector.detyz2xy([A[i,A_id['detyd']],A[i,A_id['detzd']]],
                                         param['o11'],
                                         param['o12'],
                                         param['o21'],
                                         param['o22'],
                                         param['dety_size'],
                                         param['detz_size'])

        out = format %(sc,
                       fc,
                       A[i,A_id['omega']]*180/n.pi,
                       25,
                       A[i,A_id['Int']]/25,
                       sr,
                       fr,
                       1,
                       1,
                       0,
                       1,
                       0,
                       0,
                       A[i,A_id['Int']],
                       (A[i,A_id['Int']]/25)**2, #best estimate of sum_i I_i^2
                       A[i,A_id['Int']]/10,
                       sc,
                       fc,
                       A[i,A_id['omega']]*180/n.pi,
                       sc-2,
                       sc+2,
                       fc-2,
                       fc+2,
                       A[i,A_id['omega']]*180/n.pi,
                       0,
                       A[i,A_id['dety']],
                       A[i,A_id['detz']],
                       0,
                       0,
                       A[i,A_id['spot_id']]
                      )
        f.write(out)

    	
    f.close()   
      
      
def write_grains(param):
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
            

def write_gve(param,grain,hkl):
    """
    Write gvector (gve) file, for format see
    http://fable.wiki.sourceforge.net/imaged11+-+file+formats
    
    Henning Osholm Sorensen, RisoeDTU, 2008.
    python translation: Jette Oddershede, Risoe DTU, March 31 2008
    changed October 1, 2008 to new .gve format including detector.par and [xl,yl,zl]
    """

    A = grain[0].refs
    for grainno in range(1,param['no_grains']):
        A = n.concatenate((A,grain[grainno].refs))
    		
    nrefl = A.shape[0]

    # from detector.par 
    (z_center, y_center) = detector.detyz2xy([param['dety_center'],param['detz_center']],
					     param['o11'],
					     param['o12'],
					     param['o21'],
					     param['o22'],
					     param['dety_size'],
					     param['detz_size'])		
    dout = "# chi 0.0\n" 
    dout = dout + "# distance %f\n" %(param['distance']*1000.) 
    dout = dout + "# fit_tolerance 0.5\n" 
    dout = dout + "# o11 %i\n" %param['o11']
    dout = dout + "# o12 %i\n" %param['o12']
    dout = dout + "# o21 %i\n" %param['o21']
    dout = dout + "# o22 %i\n" %param['o22']
    dout = dout + "# omegasign %f\n" %param['omega_sign']
    dout = dout + "# t_x 0\n" 
    dout = dout + "# t_y 0\n" 
    dout = dout + "# t_z 0\n" 
    dout = dout + "# tilt_x %f\n" %param['tilt_x']
    dout = dout + "# tilt_y %f\n" %param['tilt_y']
    dout = dout + "# tilt_z %f\n" %param['tilt_z']
    dout = dout + "# y_center %f\n" %y_center
    dout = dout + "# y_size %f\n" %(param['y_size']*1000.)
    dout = dout + "# z_center %f\n" %z_center
    dout = dout + "# z_size %f\n" %(param['z_size']*1000.)

    for phase in param['phase_list']:
        if param['no_phases'] > 1:
            filename = '%s/%s_phase_%i.gve' %(param['direc'],param['stem'],phase)
        else:
            filename = '%s/%s.gve' %(param['direc'],param['stem'])
        f = open(filename,'w')
        lattice = param['sgname_phase_%i' %phase][0]
        format = "%f "*6 + "%s "*1 +"\n"
        unit_cell = param['unit_cell_phase_%i' %phase]
        out = format %(unit_cell[0],unit_cell[1],unit_cell[2],
                       unit_cell[3],unit_cell[4],unit_cell[5],
                       lattice)
        out = out + "# wavelength = %s\n" %(param['wavelength'])
        out = out + "# wedge = %f\n" %param['wedge']
        out = out + "# axis = 0.000 0.0000 1.0000\n" 
	    # insert detector.par as comment
        out = out + "# cell__a %s\n" %unit_cell[0]
        out = out + "# cell__b %s\n" %unit_cell[1]
        out = out + "# cell__c %s\n" %unit_cell[2]
        out = out + "# cell_alpha %s\n" %unit_cell[3]
        out = out + "# cell_beta %s\n" %unit_cell[4]
        out = out + "# cell_gamma %s\n" %unit_cell[5]	
        out = out + "# cell_lattice_[P,A,B,C,I,F,R] %s\n" %param['sgname_phase_%i' %phase][0]
        out = out + dout
        # continue with gve format
        out = out +"# ds h k l\n" 
        f.write(out)
    	
        thkl = hkl[param['phase_list'].index(phase)].copy()
        ds = n.zeros((thkl.shape[0],1))

        for i in range(thkl.shape[0]):
            ds[i] = 2*tools.sintl(unit_cell,thkl[i,0:3])
    
        #Add ds values to the thkl array    
        thkl = n.concatenate((thkl,ds),1)
    
        # sort rows according to ds, descending
        thkl = thkl[n.argsort(thkl,0)[:,4],:]

        # Output format
        format = "%f "*1 + "%d "*3 +"\n"

        for i in range(thkl.shape[0]):
            out = format %(thkl[i,4],
                           thkl[i,0],
                           thkl[i,1],
                           thkl[i,2]
                           )
            f.write(out)
            
        R_tilt = tools.detect_tilt(param['tilt_x'],param['tilt_y'],param['tilt_z'])

        out = "# xr yr zr xc yc ds eta omega spot3d_id xl yl zl\n" 
        f.write(out)
        format = "%f "*8 + "%i "*1+ "%f "*3+"\n"


        for i in range(nrefl):
            (sc, fc) = detector.detyz2xy([A[i,A_id['dety']],A[i,A_id['detz']]],
                                         param['o11'],
                                         param['o12'],
                                         param['o21'],
                                         param['o22'],
                                         param['dety_size'],
                                         param['detz_size'])
            [xl,yl,zl] = detector.detector2lab(A[i,A_id['dety']],A[i,A_id['detz']],
                                               param['distance'],
                                               param['y_size'],param['z_size'],
                                               param['dety_center'],param['detz_center'],
                                               R_tilt)
            out = format %(A[i,A_id['gv1']]/(2*n.pi),
                           A[i,A_id['gv2']]/(2*n.pi),
                           A[i,A_id['gv3']]/(2*n.pi),
                           sc,#A[i,A_id['detz']],
                           fc,#param['dety_size']-A[i,A_id['dety']],
                           (2*n.sin(.5*A[i,A_id['tth']])/param['wavelength']),
                           A[i,A_id['eta']]*180/n.pi,
                           A[i,A_id['omega']]*180/n.pi,
                           A[i,A_id['spot_id']],
                           xl*1000.,
                           yl*1000.,
                           zl*1000.
                           )
            f.write(out)
    	
        f.close()   


def write_ini(param,hkl):
    """
    Write ini file for grainspotter, see
    http://fable.wiki.sourceforge.net/grainspotter
    
    Henning Osholm Sorensen, RisoeDTU, 2008.
     """
            

    for phase in param['phase_list']:
        out = '! input file for GrainSpotter made by PolyXSim\n'

        if param['no_phases'] > 1:
            filename = '%s/%s_phase_%i.ini' %(param['direc'],param['stem'],phase)
            out = out + 'filespecs %s/%s_phase_%i.gve %s/%s_phase_%i.log\n' %(param['direc'],
                                                                              param['stem'],
                                                                              phase,
                                                                              param['direc'],
                                                                              param['stem'],
                                                                              phase)
        else:
            filename = '%s/%s.ini' %(param['direc'],param['stem'])
            out = out + 'filespecs %s/%s.gve %s/%s.log\n' %(param['direc'],
                                                            param['stem'], 
                                                            param['direc'],
                                                            param['stem'])

        unit_cell = param['unit_cell_phase_%i' %phase]
        
        thkl = hkl[param['phase_list'].index(phase)].copy()
        ds = n.zeros((thkl.shape[0],1))
            
        for i in range(thkl.shape[0]):
            ds[i] = 2*tools.sintl(unit_cell,thkl[i,0:3])
                
        ds.sort()
        ds = ds.round(9)
        ds = n.unique(ds)
        families=  len(ds)
        Nhkls = n.min([families, 8])
        extra_hkls = 4
        if families-Nhkls < extra_hkls:
            tth_max = 2*param['theta_max']
        else:
            ds_max = ds[Nhkls+extra_hkls-1]+ 0.001
            tth_max = 2*n.arcsin(ds_max*param['wavelength']/2.)*180./n.pi

        ds_min = ds[0] - 0.001
        if ds_min < 0.0:
            ds_min = 0.
        tth_min = 2*n.arcsin(ds_min*param['wavelength']/2.)*180./n.pi
                                                            
        out = out + 'spacegroup %i\n' %param['sgno_phase_%i' %phase]
        out = out + 'etarange %f %f\n'%(0.0, 360.0)
        out = out + 'domega %f\n' %param['omega_step']
        out = out + 'omegarange %f %f\n'   %(param['omega_start'],param['omega_end'])
        out = out + 'cuts %i %f %f\n' %(8, 0.6, 0.75)
        out = out + 'eulerstep %f\n' %(5.0)
        out = out + 'uncertainties %f %f %f\n' %(.05, 0.5, 1.0)
        out = out + 'nsigmas %f\n' %(2.0)
        out = out + 'Nhkls_in_indexing %i\n' %(Nhkls)
        out = out + 'tthrange %f %f\n' %(tth_min,tth_max)
        out = out + 'minfracg %f\n'%(0.95)

        f = open(filename,'w')
        f.write(out)

def write_par(param):
    """
    Save the detector parameters

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
    dout = dout + "wedge %f\n" %param['wedge']
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



def write_ref(param,grain,grainno=None):
    """
    write PolyXSim ref file 
    """

    if grainno == None:
        savegrains = range(len(grain))
    else:
        savegrains = grainno

    for grainno in savegrains:
        A = grain[grainno].refs
        setno = 0
        filename = '%s/%s_gr%0.4d_set%0.4d.ref' \
            %(param['direc'],param['stem'],param['grain_list'][grainno],setno)
        f = open(filename,'w')
        format = "%d "*6 + "%f "*14 + "%d "*1 + "\n"
        out = "#"
        A_col = dict([[v,k] for k,v in A_id.items()])
        for col in A_col:
            out = out + ' %s' %A_col[col]
        out = out +"\n"

        f.write(out)
        # Only write reflections to file if some present
        if len(A) > 0:
            ( nrefl, ncol ) = A.shape
            for i in range(nrefl):
                out = format %(A[i,A_id['grain_id']],
                               A[i,A_id['ref_id']],
                               A[i,A_id['spot_id']],   
                               A[i,A_id['h']],
                               A[i,A_id['k']],
                               A[i,A_id['l']],
                               A[i,A_id['tth']]*180/n.pi,
                               A[i,A_id['omega']]*180/n.pi,
                               A[i,A_id['eta']]*180/n.pi,
                               A[i,A_id['dety']],
                               A[i,A_id['detz']],
                               A[i,A_id['detyd']],
                               A[i,A_id['detzd']],
                               A[i,A_id['gv1']],
                               A[i,A_id['gv2']],
                               A[i,A_id['gv3']],
                               A[i,A_id['L']],
                               A[i,A_id['P']],
                               A[i,A_id['F2']],
                               A[i,A_id['Int']],
                               A[i,A_id['overlaps']]
                       )
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

	
