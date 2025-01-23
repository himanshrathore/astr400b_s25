#Reading particle properties of a snaphot

#necessary imports
import numpy as np
import astropy.units as u
from ReadFile import Read

def ParticleInfo(filename, particle_type, particle_number):
    '''
    Obtain information of any given particle from a snapshot.
    Inputs:
    :filename -> string, name of the snapshot file
    :particle_type -> int, particle type out of 1, 2 and 3. 1 is dm, 2 is disk and 3 is bulge.
    :particle_number -> int, the particle number. Note that particle numbers start from 0.
    Returns:
    3-tuple, which contains the following:
    :distance -> magnitude of the position vector of the particle, astropy qty, kpc
    :speed -> magnitude of the velocity vector of the particle, astropy qty, km/s
    :mass -> particle mass in solar masses, astropy qty
    '''
    
    #read the file
    _, _, data = Read(filename)
    
    #extract the data corresponding to the particle type
    index = np.where(data['type'] == particle_type)
    data_ptype = data[index]
  
    #extract the particular particle
    data_particle = data_ptype[particle_number]
    
    #calculate the distance
    x, y, z = data_particle['x'], data_particle['y'], data_particle['z']
    distance = np.round(np.sqrt(x**2 + y**2 + z**2), 3)*u.kpc #kpc
    
    #calculate the speed
    vx, vy, vz = data_particle['vx'], data_particle['vy'], data_particle['vz']
    speed = np.round(np.sqrt(vx**2 + vy**2 + vz**2), 3)*u.km/u.s #km/s
    
    #calculate the mass
    m = data_particle['m']
    mass = np.round(m*1e10, 3)*u.Msun #solar masses
    
    return distance, speed, mass

#end
    