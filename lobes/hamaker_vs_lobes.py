import matplotlib
matplotlib.use('tkagg')
from matplotlib import pyplot as plt
from scipy.io import loadmat
import lobes
import numpy as np

# create a vectorized element_response_lba function from the scalar function in lobes module
# results are stacked in one array
element_response_lba = lambda *args : np.stack(np.frompyfunc(lobes.element_response_lba, 3, 1)(*args))

if 'simdata' not in globals():
    simdata = loadmat('LBA_CS302_fine')

if 'coords' not in globals():
    coords = loadmat('CS302_coords')

P = coords['LBA']['P'][0,0][:,0]
Q = coords['LBA']['Q'][0,0][:,0]
positions = np.array([P, Q])

plt.figure()
for i, (p, q) in enumerate(zip(P,Q)):
    plt.plot(p,q,'x')
    plt.text(p,q, str(i))

theta_idx = -2

freq_idx = 0

theta_deg = simdata['Theta'][0, theta_idx]
phi_deg = simdata['Phi'][0]

E1 = simdata['Vdiff_pol1'].copy()
E2 = simdata['Vdiff_pol2'].copy()

E_norm = np.sqrt(np.mean(abs(E1[:,:,0,:,:]**2 + abs(E2[:,:,0,:,:])**2), axis=(0,1,2)))

E1 /= E_norm[np.newaxis, np.newaxis, np.newaxis, np.newaxis, :]
E2 /= E_norm[np.newaxis, np.newaxis, np.newaxis, np.newaxis, :]

for el in [0,1]:
    plt.figure()
    plt.plot(phi_deg, abs(E1[el,:,theta_idx,:,freq_idx].T))
    plt.plot(phi_deg, abs(E2[el,:,theta_idx,:,freq_idx].T))

    plt.legend(['X dipole, Etheta', 'X dipole, Ephi','Y dipole, Etheta', 'Y dipole, Ephi'])
    plt.title("element = {element}, theta = {theta:f}".format(element=el, theta=theta_deg))
    plt.xlabel('phi(deg)')
    plt.ylabel('E (V)')

theta = theta_deg / 180.0 * np.pi
#phi = np.mod((simdata['Phi'][0, :] - 45) / 180.0 * np.pi, 2*np.pi)

phi_deg = np.linspace(0,360, 361)
phi = (phi_deg - 45) / 180.0 * np.pi

freq = simdata['Freq'][0,freq_idx]

r = element_response_lba(freq, theta, phi)

plt.figure()
plt.plot(phi_deg, abs(r[:,0,0]))
plt.plot(phi_deg, abs(r[:,1,0]))
plt.plot(phi_deg, abs(r[:,0,1]))
plt.plot(phi_deg, abs(r[:,1,1]))
plt.title("Hamaker, theta = {theta:f}".format(element=el, theta=theta_deg))
plt.xlabel('phi(deg)')
plt.ylabel('E (V)')

plt.show()


