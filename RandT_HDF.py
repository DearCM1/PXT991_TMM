#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 14:19:57 2020

@author: Calum
"""
# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# FUNCTIONS
# =============================================================================
# reflection coefficient
def rij(ni, nj):
    return (ni - nj)/(ni + nj)

# transmission coefficient
def tij(ni, nj):
    return (2 * ni)/(ni + nj)

# reflectivity
def R(n1, n2, n3, lam0, theta2, h):
    beta = (2 * np.pi / lam0) * (n2 * h) * np.cos(theta2)
    r12 = rij(n1, n2)
    r23 = rij(n2, n3)
    numerator = (r12**2) + (r23**2) + (2 * r12 * r23 * np.cos(2 * beta))
    denominator = 1 + ((r12**2)*(r23**2)) + (2  * r12 * r23 * np.cos(2 * beta))
    return numerator / denominator

# =============================================================================
# PARAMETERS
# =============================================================================
n1 = 1.
n2 = np.array([1.0, 1.2, 1.4, 1.5, 1.7, 2.0, 3.0])
n3 = 1.5
lam = 1500e-9
h = lam / n2
theta2 = 0

# =============================================================================
# PLOTTING
# =============================================================================
# plot initialise
fig = plt.figure()
ax = plt.subplot(121)
colour = np.array([['green', 'solid'],
                   ['black', 'dashed'],
                   ['black', 'solid'],
                   ['black', (0, (5, 5))],
                   ['black', 'dotted'],
                   ['black', 'dashdot'],
                   ['black', 'solid']], dtype=object)

# plot execute
x = np.linspace(0, 1, 1000)
for i in range(len(n2)):
    H = np.linspace(0, h[i], 1000)
    r = np.array([])
    for j in range(len(H)):
        r = np.append(r, R(n1, n2[i], n3, lam, theta2, H[j]))
    ax.plot(x, r, label="$n_2$ = {:.1f}".format(n2[i]),
             c = colour[i,0],
             linestyle = colour[i,1])
    
# ranging adjustments
plt.xlim(min(x), max(x))
plt.ylim(0, 0.55)

# plot formatting
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.xlabel("Optical Thickness, [$h/(λ_{0}/n_{2})$]")
plt.ylabel("ℛ",fontsize=14)

# export figure
plt.savefig("fig1.png", bbox_inches='tight', dpi=300)