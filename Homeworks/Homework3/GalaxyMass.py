#to calculate the mass of a galaxy

#necessary imports
import numpy as np
from ReadFile import Read

def ComponentMass(filename, particle_type):
    '''
    Returns the total mass of any desired galaxy component.
    Inputs:
    filename -> path to the snapshot file, str
    particle_type -> particle type, int
    Returns:
    Mass -> float, units of 1e12 Msun
    '''
    
    #reading the file
    _, _, data = Read(filename)
    
    #extract the data corresponding to the particle type
    index = np.where(data['type'] == particle_type)
    data_ptype = data[index]
    
    #extract the masses
    masses = data_ptype['m']
    
    #total mass
    tot_mass = np.round(np.sum(masses)/1e2, 3) #required units
    
    return tot_mass

#end