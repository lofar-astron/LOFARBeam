#!/usr/bin/env python3

import lobes
from scipy.io import loadmat
import numpy as np
import h5py

simdata = loadmat('LBA_CS302_fine')

E1 = simdata['Vdiff_pol1'].copy() # X? dipole, shape: N_ant x [Etheta, Ephi] x Ntheta x Nphi x Nfreq
E2 = simdata['Vdiff_pol2'].copy() # Y? dipole, shape: N_ant x [Etheta, Ephi] x Ntheta x Nphi x Nfreq


E1 = E1[:,:,1:91,:-1,:]
E2 = E2[:,:,1:91,:-1,:]


#with h5py.File('LBA_CS302.hdf5', 'r') as f:
    #nms = f["nms"][:]


theta = simdata['Theta']/180.0*np.pi
phi = simdata['Phi']/180.0*np.pi

theta = theta[0, 1:91]
phi=phi[0, 0:-1]

theta_mesh, phi_mesh = np.meshgrid(theta, phi)

theta_mesh_flat = theta_mesh.flatten()
phi_mesh_flat = phi_mesh.flatten()

lbm = lobes.LobesBeamModel("LBA_CS302.hdf5")

r = lbm.eval(theta_mesh_flat, phi_mesh_flat)

Ntheta = len(theta)
Nphi = len(phi)
Nfreq = 4
Nant = 96

rr = r.reshape((Nphi, Ntheta, 2, Nfreq, Nant)).transpose((4,2,1,0,3))


