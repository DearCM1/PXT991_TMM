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
R, w = np.loadtxt("R_i.txt",
                  unpack=True,
                  skiprows=1,
                  delimiter='    ')
T, w = np.loadtxt("T_i.txt",
                  unpack=True,
                  skiprows=1,
                  delimiter='  ')
# =============================================================================
# Plotting
# =============================================================================
# Initialise
fig1 = plt.figure(figsize=(11, 4))

# Left Plot
ax1 = fig1.add_subplot(121)
ax1.plot(w, R, c='g')
plt.xlim(0.5, 1.5)
plt.ylim(0,1)
plt.xlabel("$\omega$/$\omega_{o}$", fontsize=12)
plt.ylabel("R", family='cursive', fontsize=14)

# Right Plot
ax2 = fig1.add_subplot(122)
ax2.semilogy(w, T, c='c',)
plt.xlim(0.5, 1.5)
plt.ylim(min(T),1)
plt.xlabel("$\omega$/$\omega_{o}$", fontsize=12)
plt.ylabel("T",family='cursive',fontsize=14)

plt.savefig("figure.png", bbox_inches='tight', dpi=300)