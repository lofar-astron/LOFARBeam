# Original Matlab code by Michel Arts
# Ported to python by Sebastiaan van der Tol
# Interface change: theta, phi are now in radians
# 2019-10-16

import numpy as np
from scipy.special import lpmn
import lobes

def F4far_new(s,m,n,theta,phi,beta):
    if s==1:
        q1, q2, q3 = F41(m,n,theta,phi,beta)
        return q2,q3
    elif s==2:
        q1, q2, q3 = F42(m,n,theta,phi,beta)
        return q2,q3
    raise RuntimeError("Error arguments to F4far_new: s should be either 1 or 2")

def legendre(n, x):
    results = np.zeros((n+1, len(x)))
    for i in range(len(x)):
        assoc_legendre, assoc_legendre_derivative = lpmn(n,n, x[i])
        results[:, i] = assoc_legendre[:, -1]
    return results

def P(m,n,x):
    q=legendre(n,x)
    if n>0:
        q=q[abs(m),:]

    q1 = lobes.P(m, n, x)

    if m<0:
        q = (-1.**-m) * np.math.factorial(n+m)/np.math.factorial(n-m) * q

    return q

def Pacc(m,n,z):
    q = (-(n+m) * (n-m+1) * np.sqrt(1-z*z) * P(m-1, n, z) - m * z * P(m, n, z))/(z*z-1)
    q1 = lobes.Pacc(m, n, z)

    return q;

def F41(m,n,theta,phi,beta):
    if m != 0:
        C = beta * np.sqrt(60) * 1.0/np.sqrt(n*(n+1)) * (-m/abs(m))**m
    else:
        C = beta * np.sqrt(60) * 1.0/np.sqrt(n * (n+1))
    q1 = np.zeros((1,len(theta)))
    q2 = C * (-1j)**(-n-1)/beta *1j *m / (np.sin(theta)) * np.sqrt((2.*n+1)/2 * np.math.factorial(n-abs(m)) / np.math.factorial(n+abs(m))) * P(abs(m),n,np.cos(theta))*np.exp(1j*m*phi)
    q3 = C * (-1j)**(-n-1)/beta * np.sqrt((2.*n+1) / 2.0 * np.math.factorial(n-abs(m)) / np.math.factorial(n+abs(m))) * Pacc(abs(m),n,np.cos(theta))*np.sin(theta)*np.exp(1j*m*phi)
    return (q1, q2, q3)


def F42(m,n,theta,phi,beta):

    if m != 0:
        C = beta * np.sqrt(60) * 1.0/np.sqrt(n*(n+1))*(-m/abs(m))**m
    else:
        C = beta * np.sqrt(60) * 1.0/np.sqrt(n*(n+1))

    q1 = np.zeros((1,len(theta)))
    q2 = -C * (-1j)**(-n) / beta * np.sqrt((2.*n+1)/2.0 * np.math.factorial(n-abs(m)) / np.math.factorial(n+abs(m))) * Pacc(abs(m),n,np.cos(theta)) * np.sin(theta) * np.exp(1j*m*phi)
    q3 = C * (-1j)**(-n) / beta * 1j * m/np.sin(theta) * np.sqrt((2.*n+1) / 2.0 * np.math.factorial(n-abs(m)) / np.math.factorial(n+abs(m))) * P(abs(m),n,np.cos(theta)) * np.exp(1j*m*phi)

    return (q1, q2, q3)


