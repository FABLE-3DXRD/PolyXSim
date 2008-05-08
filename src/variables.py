#!/usr/bin/env python


class frameinfo_cont:
    def __init__(self, frameno):
        self.no = frameno
        self.generic = None
        
    def refl_add(self):
        self.refl.append(refl_cont())
        
    def generic_add(self):
        if self.generic == None:
            self.generic = []
        self.generic.append(generic_cont())

class grain_cont:
    def __init__(self, U=None):
        self.U = U
        self.refl = []
        self.generic = None
        
    def refl_add(self):
        self.refl.append(refl_cont())
        
    def generic_add(self):
        if self.generic == None:
            self.generic = []
        self.generic.append(generic_cont())

class refl_cont:
    def __init__(self,hkl=None, omega=None, tth=None, eta=None, dety=None, detz=None, gt = None):
        self.hkl = hkl
        self.omega = omega
        self.tth = tth
        self.eta = eta
        self.dety = dety
        self.detz = detz
        self.gt = gt

class generic_cont:
    def __init__(self):
        pass

class shoebox_cont:
    def __init__(self):
        pass
        
class refarray:
    def __init__(self):
        self.A_id = { 'grain_id'       :0, 
                      'ref_id'         :1, 
                      'spot_id'        :2, 
                      'h'              :3, 
                      'k'              :4, 
                      'l'              :5, 
                      'tth'            :6, 
                      'omega'          :7, 
                      'eta'            :8, 
                      'dety'           :9, 
                      'detz'           :10,
                      'detyd'          :11,
                      'detzd'          :12,
                      'gv1'            :13,
                      'gv2'            :14,
                      'gv3'            :15,
                      'L'              :16,
                      'P'              :17,
                      'F2'    	       :18,
                      'Int'            :19,
                      'overlaps'       :20
                      }


if __name__=='__main__':
    grain = []
    grain.append(grain_cont())
    grain[0].refl_add()
    print grain[0].refl[0].hkl
#     mycont = container()
#     mycont.grainno = 2
#     mycont.grain = []
#     print mycont.name
#     mycont.grain.append({'reflno': 20})
#     mycont.grain.append({'reflno': 30})
#     mycont.grain.append({'reflno': 40})
#     mycont.grain.append({'reflno': 50})
#     mycont.grain[1]['reflno'] = 100
#     mycont.grain[1][1] = {'tth': 12.3 , 'omega':-0.1}
    
#     print'var ', mycont.grain
#     print mycont.grain[1][1]['tth'] 


#     grain = []
#     grain.append(container())
#     grain[0].nref = []
#     grain[0].nref.append(container())
#     grain[0].nref[0].tth = 12.3
    
#     print grain[0].nref[0].tth
#     print grain[0].nref



# This used to save grain data in find_refl earlier

# First make constainer for grain(X)
 #self.grain.append(variables.grain_cont(U))
# add new reflection             
 #self.grain[grainno].refl_add()
# add specs for reflection
 #self.grain[grainno].refl[nrefl].hkl = hkl
 #self.grain[grainno].refl[nrefl].eta = eta
 #self.grain[grainno].refl[nrefl].omega = omega
 #self.grain[grainno].refl[nrefl].tth = 2*theta
 #self.grain[grainno].refl[nrefl].dety = dety
 #self.grain[grainno].refl[nrefl].detz = detz
 #self.grain[grainno].refl[nrefl].Gt = Gt
 #self.grain[grainno].refl[nrefl].ybox = self.param['ybox']
 #self.grain[grainno].refl[nrefl].zbox = self.param['zbox']
 #self.grain[grainno].refl[nrefl].wbox = self.param['wbox']
 #self.grain[grainno].refl[nrefl].framectr =\
 #                               floor(((omega*180/N.pi-self.param['omega_start'])/self.param['omega_step'])+1)
 #self.grain[grainno].refl[nrefl].AOI  = [round(dety-self.yboxhalf), round(detz-self.zboxhalf),\
 #                                        round(dety+self.yboxhalf), round(detz+self.zboxhalf)]
 #self.grain[grainno].refl[nrefl].frames = range(\
 #    int(self.grain[grainno].refl[nrefl].framectr-(self.param['wbox']-1)/2),\
 #    int(self.grain[grainno].refl[nrefl].framectr+(self.param['wbox']-1)/2+1))
 #Determine if box is on the border
 #if (self.grain[grainno].refl[nrefl].frames[0] < 1) or\
 #       (self.grain[grainno].refl[nrefl].frames[-1] >self.nframes):
 #    self.grain[grainno].refl[nrefl].border = 'yes'
 #else:
 #    self.grain[grainno].refl[nrefl].border = 'no'
