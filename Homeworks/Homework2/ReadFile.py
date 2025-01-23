# Read file function for reading MW-M31 simulation snapshots.

#necessary imports
import numpy as np
import astropy.units as u

def Read(filename): 
    '''
    Function to read simulation snapshots.
    Inputs:
    :filename -> string, path of the snapshot file
    Returns:
    3-tuple containing:
    :time -> snapshot time in Myr, astropy quantity
    :n_tot -> total number of particles, int
    :data -> the datatable of particles. Columns: "type, m, x, y, z, vx, vy, vz"
    '''
    
    #open the file
    file = open(filename, 'r')
    
    #Read the first line and store the time in units of Myr.
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Myr #time of the snapshot in Myr
    
    #Read the second line and store the total number of particles.
    #Next time you run readline(), it will read the second line.
    line2 = file.readline()
    label, value = line2.split()
    n_tot = int(value) #total number of particles
    
    file.close() #close file
    
    #store the remainder of the file
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    
    return time, n_tot, data

#end