import numpy as np
import math
import cmath
import matplotlib.pyplot as plt

def gabor_wave_space_2d(n,sigma,zeta,eta,theta,a,m):
    # generate one gabor wavelet with specified scale and rotation in space
    x = np.arange(-n/2,n/2)
    y = x
    X, Y = np.meshgrid(x,y)
    psi = np.zeros((n,n))

    X_new = a**(-m) * (X * np.cos(theta) + Y * np.sin(theta)); #rotate x and y
    Y_new = a**(-m) * (- X * np.sin(theta) + Y * np.cos(theta));
    #calculate gabor wavelet
    psi = a**(-2 * m) / (2 * math.pi * sigma**2 * zeta)* \
        np.exp((- X_new**2 - (Y_new/zeta)**2)/(2 * sigma**2)) * np.exp(1j * eta * X_new);
    return psi
    
def gabor_wave_family_space_2d(n,K,Q,S,sigma,zeta,eta,a):
    #generate a family of gabor wavelets with specified scales and rotations in space
    psi = np.zeros((n,n,K,Q*S),dtype=complex)
    
    for i in range(K):
        for k in range(S):
            for l in range(Q):
                psi[:, :, i,k*Q+l] = gabor_wave_space_2d(n,sigma,zeta,eta,i*math.pi/K,a,k+l/Q)
    return psi
    
def gabor_wave_freq_2d(n, sigma, zeta, eta, a,j,theta):
    # generate one gabor wavelet with specified scale and rotation in frequency
    pi = math.pi
    omega1 = np.linspace(-3*pi, 3*pi-(2*pi)/n, 3*n)
    omega2 = omega1

    omega1,omega2 = np.meshgrid(omega1,omega2)
    omega1_new = omega1 * np.cos(theta) + omega2 * np.sin(theta)
    omega2_new = - omega1 * np.sin(theta) + omega2 * np.cos(theta)

    psi_hat = np.exp(- 1/2 * sigma**2 * (a**j * omega1_new - eta)**2) * \
              np.exp(- 1/2 * sigma**2 * zeta**2 * (a**j * omega2_new)**2)

    psi_add = np.zeros((n,n))
    for i in range(3):
        for k in range(3):
            psi_add = psi_add + psi_hat[i * n :(i + 1) * n ,k * n:(k + 1) * n ]

    return psi_add
    
def gabor_wave_freq_family_2d(n,K,S,Q,sigma, zeta, eta,a):
    #generate a family of gabor wavelets with specified scales and rotations in space
    psi_hat = np.zeros((n,n,K,S*Q), dtype = complex)
    for i in range(K):
        for j in range(S):
            for l in range(Q):
                psi_hat[:,:,i,j*Q+l] = gabor_wave_freq_2d(n,sigma,zeta,eta,a,j+l/Q,i*math.pi/K)
    return psi_hat
