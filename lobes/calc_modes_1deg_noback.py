import numpy as np
from scipy.constants import speed_of_light
from scipy.io import loadmat
from F4far_new import F4far_new
import lobes
import h5py


import matplotlib
matplotlib.use('tkagg')
from matplotlib import pyplot as plt

coords = loadmat('CS302_coords')

simdata = loadmat('LBA_CS302_fine')

#%
#% Compensate phase for element position
#%

P = coords['LBA']['P'][0,0][:,0]
Q = coords['LBA']['Q'][0,0][:,0]
positions = np.array([P, Q])

V_pol1 = simdata['Vdiff_pol1']

theta = simdata['Theta']/180.0*np.pi
phi = simdata['Phi']/180.0*np.pi

Freq = simdata['Freq'][0,:]

K = 2.0*np.pi*Freq/speed_of_light

sintheta = np.sin(theta)
cosphi = np.cos(phi)
sinphi = np.sin(phi)


x = sintheta.T*cosphi
y = sintheta.T*sinphi

dirs = np.array([x,y])

distance_diff = np.tensordot(positions, dirs, ((0,),(0,)))

phase_diff = np.multiply.outer(distance_diff, K)

phasor = np.exp( -1j * phase_diff)

V_pol1 *= phasor[:,np.newaxis,:,:,:]

El = 0;
FreqN = 2

nr_elements = V_pol1.shape[0]
nr_frequencies = V_pol1.shape[-1]

Etheta = V_pol1[:,0,1:91,:-1,:].transpose((0,3,1,2))
Ephi = V_pol1[:,1,1:91,:-1,:].transpose((0,3,1,2))

Etheta = Etheta.reshape((nr_elements, nr_frequencies,-1))
Ephi = Ephi.reshape((nr_elements, nr_frequencies,-1))

E = np.concatenate((Etheta, Ephi), axis=2)

theta = theta[0, 1:91]
phi=phi[0, 0:-1]

theta_mesh, phi_mesh = np.meshgrid(theta, phi)

beta = 2.0 * np.pi * Freq[FreqN] / speed_of_light

Nmax = 21

theta_mesh_flat = theta_mesh.flatten()
phi_mesh_flat = phi_mesh.flatten()

F = np.zeros((E.shape[2], 2*Nmax*(Nmax+2)), dtype=np.complex128)

jn = 0
nms = []
for n in range(1, Nmax+1):
    print(n)
    for m in range(-n, n+1):
        for s in [1,2]:
            #q2, q3 = F4far_new(s,m,n,theta_mesh_flat, phi_mesh_flat,beta);
            q2, q3 = lobes.F4far_new(s,m,n,theta_mesh_flat, phi_mesh_flat,beta);
            F[:,jn] = np.concatenate((q2, q3))
            nms.append((n,m,s))
            jn += 1

print("Inverting matrix")
Finv = np.linalg.pinv(F)

print("Compute product")
p = np.tensordot(Finv, E, axes = ([1],[2]))

idx = np.argsort(np.mean(abs(p)**2,axis=(1,2)))[::-1]

plt.figure()
plt.semilogy(abs(p[idx,:,:]).reshape((p.shape[0], -1)))

plt.show()

with h5py.File('LBA_CS302.hdf5', 'w') as f:
    f.create_dataset("coefficients", data=p.T)
    f.create_dataset("nms", data=nms)
    f.create_dataset("frequencies", data=Freq.astype(np.double))

#disp(' ')
#Q=pinv(F)*E;
#Erecon=F*Q;
#E=reshape(E,length(Theta),length(Phi),2);
#Erecon=reshape(Erecon,length(Theta),length(Phi),2);
#%
#h1=figure;
#plot(Theta,20.*log10(abs(E(:,Phi==45,1))),'b')
#grid on
#hold on
#plot(Theta,20.*log10(abs(Erecon(:,Phi==45,1))),'r')
#xlabel('Theta (degrees)')
#ylabel('dB')
#title({['EEP element ' num2str(El) ', \phi=45 deg., freq. ' num2str(Freq(FreqN)./1e6) ' MHz'],'\theta-component (amplitude)'})
#h2=figure;
#plot(Theta,angle(E(:,Phi==45,1))./pi.*180,'b')
#grid on
#hold on
#plot(Theta,angle(Erecon(:,Phi==45,1))./pi.*180,'r')
#xlabel('Theta (degrees)')
#ylabel('degrees')
#title({['EEP element ' num2str(El) ', \phi=45 deg., freq. ' num2str(Freq(FreqN)./1e6) ' MHz'],'\theta-component (phase)'})
#h3=figure;
#plot(Theta,20.*log10(abs(E(:,Phi==135,2))),'b')
#grid on
#hold on
#plot(Theta,20.*log10(abs(Erecon(:,Phi==135,2))),'r')
#xlabel('Theta (degrees)')
#ylabel('dB')
#title({['EEP element 1, \phi=135 deg., freq. ' num2str(Freq(FreqN)./1e6) ' MHz'],'\phi-component (amplitude)'})
#h4=figure;
#plot(Theta,angle(E(:,Phi==135,2))./pi.*180,'b')
#grid on
#hold on
#plot(Theta,angle(Erecon(:,Phi==135,2))./pi.*180,'r')
#xlabel('Theta (degrees)')
#ylabel('degrees')
#title({['EEP element 1, \phi=135 deg., freq. ' num2str(Freq(FreqN)./1e6) ' MHz'],'\phi-component (phase)'})
#%
#Diff=Erecon./E;
#h5=figure;
#plot(Theta,20.*log10(abs(Diff(:,Phi==45,1))),'b')
#grid on
#hold on
#plot(Theta,20.*log10(abs(Diff(:,Phi==135,2))),'r')
#xlabel('Theta (degrees)')
#ylabel('dB')
#legend('\theta-comp., \phi=45^\circ','\phi-comp., \phi=135^\circ')
#title({'Difference between original and','reconstructed EEP (amplitude)', ...
    #['Element ' num2str(El) ', freq. ' num2str(Freq(FreqN)./1e6) ' MHz']})
#h6=figure;
#plot(Theta,angle(Diff(:,Phi==45,1))./pi.*180,'b')
#grid on
#hold on
#plot(Theta,angle(Diff(:,Phi==135,2))./pi.*180,'r')
#xlabel('Theta (degrees)')
#ylabel('degrees')
#legend('\theta-comp., \phi=45^\circ','\phi-comp., \phi=135^\circ')
#title({'Difference between original and','reconstructed EEP (phase)', ...
    #['Element ' num2str(El) ', freq. ' num2str(Freq(FreqN)./1e6) ' MHz']})
#%
#figure(h1)
#legend('original','reconstructed')
#figure(h2)
#legend('original','reconstructed')
#figure(h3)
#legend('original','reconstructed')
#figure(h4)
#legend('original','reconstructed')
#figure
#stem(abs(Q),'b')
#grid on
#xlabel('Index number')
#title('Absolute value of calculated coefficients')
#%
#disp(['Least square error = ' num2str(norm(E(:)-F*pinv(F)*E(:),2))])
#disp(['Condition number F = ' num2str(cond(F))])
#disp(['Size F = ' num2str(size(F))])
#disp(['Rank F = ' num2str(rank(F))])
#disp(' ')
