#necessary imports

import numpy as np
import matplotlib.pyplot as plt
from ReadFile import Read
from CenterOfMass import CenterOfMass as COM
import astropy.units as u
import astropy.constants as const

class MassProfile:
    #mass profile and rotation curves based on the mass profile
    
    def __init__(self, galaxy, snap):
            
        #add a string of the filenumber to the value "000"
        ilbl = '000' + str(snap)
        #remove all but the last 3 digits
        ilbl = ilbl[-3:]
        self.filename = "%s_"%(galaxy) + ilbl + ".txt"

        # read data in the given file using Read
        _, _, self.data = Read(self.filename)

        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
        self.mass = self.data['m']

        self.gname = galaxy
            
    def MassEnclosed(self, particle_type, radius_array):
        
        myCOM = COM(self.filename, 2)
        p_COM = myCOM.COM_P(0.1)
        
        idxs = np.where(self.data['type'] == particle_type)
        
        mass_enclosed_array = np.zeros(len(radius_array))
          
        x1 = self.x[idxs] - p_COM[0]
        y1 = self.y[idxs] - p_COM[1]
        z1 = self.z[idxs] - p_COM[2]
        r1 = np.sqrt(x1**2 + y1**2 + z1**2)
        mass1 = self.mass[idxs]
        
        radius_array = radius_array*u.kpc
        
        for i in range(len(radius_array)):

            idxs = r1 < radius_array[i]
            mass_enclosed = np.sum(mass1[idxs])
            mass_enclosed_array[i] = mass_enclosed
            
        mass_enclosed_array = np.round(mass_enclosed_array, 2)*1e10*u.Msun
        
        return mass_enclosed_array
    
    def MassEnclosedTotal(self, radius_array):
        
        mtot_enclosed_halo = self.MassEnclosed(1, radius_array)
        mtot_enclosed_disk = self.MassEnclosed(2, radius_array)
        
        if(self.gname == 'M33'):
            mtot_enclosed_bulge = np.zeros(len(radius_array))*u.Msun
        else:
            mtot_enclosed_bulge = self.MassEnclosed(3, radius_array)
            
        mtot_enclosed_array = mtot_enclosed_halo + mtot_enclosed_disk + mtot_enclosed_bulge
        
        return mtot_enclosed_array
        
    def HernquistMass(self, r, a, Mhalo):
        
        M_Hernquist = Mhalo*(r**2)/((a + r)**2)
        
        return M_Hernquist*u.Msun
            
    def CircularVelocity(self, particle_type, radius_array):
        
        mass_enclosed_array = self.MassEnclosed(particle_type, radius_array) #in Msun
        
        vcirc_array = np.sqrt(const.G*mass_enclosed_array/(radius_array*u.kpc))
        vcirc_array = np.round(vcirc_array.to(u.km/u.s), 2)
        
        return vcirc_array
    
    def CircularVelocityTotal(self, radius_array):
        
        mass_enclosed_array = self.MassEnclosedTotal(radius_array) #in Msun
        
        vcirc_array = np.sqrt(const.G*mass_enclosed_array/(radius_array*u.kpc))
        vcirc_array = np.round(vcirc_array.to(u.km/u.s), 2)
        
        return vcirc_array
    
    def HernquistVCirc(self, r, a, Mhalo):
        
        hern_vcirc = np.sqrt(const.G * self.HernquistMass(r, a, Mhalo)/(r*u.kpc))
        hern_vcirc = np.round(hern_vcirc.to(u.km/u.s), 2)
        
        return hern_vcirc
    
#end