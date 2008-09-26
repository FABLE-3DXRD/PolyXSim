import numpy as n
from xfab import tools,detector
import variables

A_id = variables.refarray().A_id

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
        



def write_ini(param,hkl):
    """
    Write gvector (gve) file, for format see
    http://fable.wiki.sourceforge.net/imaged11+-+file+formats
    
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

def write_gve(param,grain,hkl):
    """
    Write gvector (gve) file, for format see
    http://fable.wiki.sourceforge.net/imaged11+-+file+formats
    
    Henning Osholm Sorensen, RisoeDTU, 2008.
    python translation: Jette Oddershede, Risoe DTU, March 31 2008
    """

    A = grain[0].refs
    for grainno in range(1,param['no_grains']):
        A = n.concatenate((A,grain[grainno].refs))
    		
    nrefl = A.shape[0]

    for phase in param['phase_list']:
        if param['no_phases'] > 1:
            filename = '%s/%s_phase_%i.gve' %(param['direc'],param['stem'],phase)
        else:
            filename = '%s/%s.gve' %(param['direc'],param['stem'])
        f = open(filename,'w')
        lattice = param['sgname_phase_%i' %phase]
        format = "%f "*6 + "%s "*1 +"\n"
        unit_cell = param['unit_cell_phase_%i' %phase]
        out = format %(unit_cell[0],unit_cell[1],unit_cell[2],
                       unit_cell[3],unit_cell[4],unit_cell[5],
                       lattice)
        f.write(out)
        out = "# wavelength = %s\n" %(param['wavelength'])
        f.write(out)
        out = "# wedge = %f\n" %param['wedge']
        f.write(out)
        out = "# ds h k l\n" 
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

        out = "# xr yr zr xc yc ds eta omega\n" 
        f.write(out)
        format = "%f "*8 + "%i"*1+"\n"


        for i in range(nrefl):
            (sc, fc) = detector.detyz2xy([A[i,A_id['dety']],A[i,A_id['detz']]],
                                         param['o11'],
                                         param['o12'],
                                         param['o21'],
                                         param['o22'],
                                         param['dety_size'],
                                         param['detz_size'])
            out = format %(A[i,A_id['gv1']]/(2*n.pi),
                           A[i,A_id['gv2']]/(2*n.pi),
                           A[i,A_id['gv3']]/(2*n.pi),
                           sc,#A[i,A_id['detz']],
                           fc,#param['dety_size']-A[i,A_id['dety']],
                           (2*n.sin(.5*A[i,A_id['tth']])/param['wavelength']),
                           A[i,A_id['eta']]*180/n.pi,
                           A[i,A_id['omega']]*180/n.pi,
                           A[i,A_id['spot_id']]
                           )
            f.write(out)
    	
        f.close()   


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
                       A[i,A_id['Int']]**2,
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
                       A[i,A_id['dety']]-param['dety_size'],
                       A[i,A_id['detz']],
                       0,
                       0,
                       A[i,A_id['spot_id']]
                      )
        f.write(out)

    	
    f.close()   
            
