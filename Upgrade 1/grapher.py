#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 17:34:55 2020

@author: Calum
"""
# =============================================================================
# Package Importing
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
# =============================================================================
# Data Importing
# =============================================================================
lam1A, R1A, T1A = np.loadtxt("RT1A.txt",
                       unpack=True,
                       skiprows=1,
                       delimiter='       ')
lam2A, R2A, T2A = np.loadtxt("RT2A.txt",
                         unpack=True,
                         skiprows=1,
                         delimiter='       ')
lam1B, R1B, T1B = np.loadtxt("RT1B.txt",
                       unpack=True,
                       skiprows=1,
                       delimiter='       ')
lam2B, R2B, T2B = np.loadtxt("RT2B.txt",
                         unpack=True,
                         skiprows=1,
                         delimiter='       ')
# =============================================================================
# Plotting
# =============================================================================
# Initialise
fig1 = plt.figure(figsize=(15,10))
gs = fig1.add_gridspec(2, 2)

# Bottom Left Plot
ax1 = fig1.add_subplot(gs[1,0])
ax1.plot(lam1A, R1A, c='k', label='Version 1')
ax1.plot(lam2A, R2A, c='g', label='Version 2')
plt.xlim(400,525)
plt.ylim(0,1)
plt.xlabel("Wavelength, λ($nm$)",
           fontsize=10)
plt.ylabel("ℛ",
           fontsize=14)
plt.legend(loc='best')

# Bottom Right Plot
ax2 = fig1.add_subplot(gs[1,1])
ax2.plot(lam1B, R1B, c='k', label='Version 1')
ax2.plot(lam2B, R2B, c='g', label='Version 2')
plt.xlim(400,525)
plt.ylim(0,1)
plt.xlabel("Wavelength, λ($nm$)",
           fontsize=10)
plt.ylabel("ℛ",
           fontsize=14)
plt.legend(loc='best')

# Top Plot
ax3 = fig1.add_subplot(gs[0,:])
ax3.plot(lam2A, R2A, c='g', label='TMM')
ax3.plot(lam2A+18, R2A, c='k', label='Ng 2014')
plt.xlim(min(lam2A)+18,max(lam2A)-18)
plt.ylim(0,1)
plt.xlabel("Wavelength, λ($nm$)",
           fontsize=10)
plt.ylabel("ℛ",
           fontsize=14)
plt.legend(loc='best')

# Savefile
plt.savefig("figure.png",
            bbox_inches='tight',
            dpi=300)