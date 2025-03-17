#HW7, orbit integration

import numpy as np
G = 4.50e-6 #kpc^3/M_sun/Gyr

class M33AnalyticOrbit:
	def __init__(self, filename):
		
		self.filename = filename
		self.r = np.array([-98.56, -119.99, -127.76])
		self.v = np.array([-29, 174, 93])
		self.rdisk = 5
		self.Mdisk = 0.120e12
		self.rbulge = 1
		self.Mbulge = 0.019e12
		self.rhalo = 62
		self.Mhalo = 1.921e12
		
	def HernquistAccel(self, M, r_a, r):
		
		r_mag = np.linalg.norm(r)
		x, y, z = r[0], r[1], r[2]
		
		accel_x = -(G*M/(r_mag*((r_a + r_mag)**2)))*x
		accel_y = -(G*M/(r_mag*((r_a + r_mag)**2)))*y
		accel_z = -(G*M/(r_mag*((r_a + r_mag)**2)))*z
		
		accel = np.array([accel_x, accel_y, accel_z])
		
		return accel	
	
	def MiyamotoNagaiAccel(self, M, r_d, r):
		
		x, y, z = r[0], r[1], r[2]
		R = np.sqrt(x**2 + y**2)
		z_d = self.rdisk/5
		B = self.rdisk + np.sqrt(z**2 + z_d**2)
		
		accel = (-G*M/((R**2 + B**2)**1.5))*r*np.array([1, 1, B/np.sqrt(z**2 + z_d**2)])
		
		return accel
		
	def M31Accel(self, r):
		
		accel = self.HernquistAccel(self.Mhalo, self.rhalo, r) + self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, r) + self.HernquistAccel(self.Mbulge, self.rbulge, r)
		return accel
		
	def LeapFrog(self, dt, r, v):
		
		rhalf = r + v*dt/2
		vnew = v + self.M31Accel(rhalf)*dt
		rnew = rhalf + vnew*dt/2
		
		return rnew, vnew
	
	def OrbitIntegrator(self, t0, dt, tmax):
		
		ntimes = int((tmax - t0)/dt) + 1
		orbit = np.zeros((ntimes, 7))
		rnew = self.r
		vnew = self.v
		t = t0
		
		orbit[0, 0] = t0
		orbit[0, 1] = rnew[0]
		orbit[0, 2] = rnew[1]
		orbit[0, 3] = rnew[2]
		orbit[0, 4] = vnew[0]
		orbit[0, 5] = vnew[1]
		orbit[0, 6] = vnew[2]
		
		i = 1
		t = t0 + dt
		while(t <= tmax):
			print(i)
			rnew, vnew = self.LeapFrog(dt, rnew, vnew)
			orbit[i, 0] = t
			orbit[i, 1] = rnew[0]
			orbit[i, 2] = rnew[1]
			orbit[i, 3] = rnew[2]
			orbit[i, 4] = vnew[0]
			orbit[i, 5] = vnew[1]
			orbit[i, 6] = vnew[2]
			t = t + dt
			i = i + 1
			
		np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#',
			   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
			          .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))


#Running the code

M33_orbit = M33AnalyticOrbit('analytic_orbit.txt')
M33_orbit.OrbitIntegrator(0, 0.1, 10)

#end
