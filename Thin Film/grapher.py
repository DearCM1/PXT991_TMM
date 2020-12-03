#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 22:05:51 2020

@author: Calum
"""
import numpy as np
import matplotlib.pyplot as plt

fig1 = plt.figure(figsize=(22.5,10))
gs = fig1.add_gridspec(4, 3)
plt.tight_layout()
# =============================================================================
# 
# =============================================================================
H, TMM, Theory = np.loadtxt("R_i1.txt",
                            unpack=True,
                            delimiter='       ',
                            skiprows=1)
fig1.add_subplot(gs[0:1,0])
plt.scatter(H[0:-1:20], TMM[0:-1:20], 
            c='k', label='TMM Evaluation',
            marker='x')
plt.plot(H, Theory, c='g', label='Theory Prediction')
plt.xlim(0,1)
plt.ylim(0.03, 0.05)
plt.xlabel("Optical Thickness, [$h / (λ_{0} / n_{1})$]",
       fontsize=12)
plt.ylabel("ℛ",
       fontsize=14)
plt.title("n = 1.00")
plt.legend(loc='best')
# =============================================================================
# 
# =============================================================================
H, TMM, Theory = np.loadtxt("R_i4.txt",
                            unpack=True,
                            delimiter='       ',
                            skiprows=1)
fig1.add_subplot(gs[1:2,0])
plt.scatter(H[0:-1:20], TMM[0:-1:20], 
            c='k', label='TMM Evaluation',
            marker='x')
plt.plot(H, Theory, c='g', label='Theory Prediction')
plt.xlim(0,1)
plt.ylim(0.03, 0.05)
plt.xlabel("Optical Thickness, [$h / (λ_{0} / n_{1})$]",
       fontsize=12)
plt.ylabel("ℛ",
       fontsize=14)
plt.title("n = 1.50")
plt.legend(loc='best')
# =============================================================================
# 
# =============================================================================
H, TMM, Theory = np.loadtxt("R_i2.txt",
                            unpack=True,
                            delimiter='       ',
                            skiprows=1)
fig1.add_subplot(gs[0:2,1])
plt.scatter(H[0:-1:20], TMM[0:-1:20], 
            c='k', label='TMM Evaluation',
            marker='x')
plt.plot(H, Theory, c='g', label='Theory Prediction')
plt.xlim(0,1)
plt.ylim(0,0.042)
plt.xlabel("Optical Thickness, [$h / (λ_{0} / n_{1})$]",
       fontsize=12)
plt.ylabel("ℛ",
       fontsize=14)
plt.title("n = 1.20")
plt.legend(loc='best')
# =============================================================================
# 
# =============================================================================
H, TMM, Theory = np.loadtxt("R_i3.txt",
                            unpack=True,
                            delimiter='       ',
                            skiprows=1)
fig1.add_subplot(gs[0:2,2])
plt.scatter(H[0:-1:20], TMM[0:-1:20], 
            c='k', label='TMM Evaluation',
            marker='x')
plt.plot(H, Theory, c='g', label='Theory Prediction')
plt.xlim(0,1)
plt.ylim(0.016,0.042)
plt.xlabel("Optical Thickness, [$h / (λ_{0} / n_{1})$]",
       fontsize=12)
plt.ylabel("ℛ",
       fontsize=14)
plt.title("n = 1.40")
plt.legend(loc='best')
# =============================================================================
# 
# =============================================================================
H, TMM, Theory = np.loadtxt("R_i5.txt",
                            unpack=True,
                            delimiter='       ',
                            skiprows=1)
fig1.add_subplot(gs[2:4,0])
plt.scatter(H[0:-1:20], TMM[0:-1:20], 
            c='k', label='TMM Evaluation',
            marker='x')
plt.plot(H, Theory, c='g', label='Theory Prediction')
plt.xlim(0,1)
plt.ylim(0.038,0.11)
plt.xlabel("Optical Thickness, [$h / (λ_{0} / n_{1})$]",
       fontsize=12)
plt.ylabel("ℛ",
       fontsize=14)
plt.title("n = 1.70")
plt.legend(loc='best')
# =============================================================================
# 
# =============================================================================
H, TMM, Theory = np.loadtxt("R_i6.txt",
                            unpack=True,
                            delimiter='       ',
                            skiprows=1)
fig1.add_subplot(gs[2:4,1])
plt.scatter(H[0:-1:20], TMM[0:-1:20], 
            c='k', label='TMM Evaluation',
            marker='x')
plt.plot(H, Theory, c='g', label='Theory Prediction')
plt.xlim(0,1)
plt.ylim(0.025,0.225)
plt.xlabel("Optical Thickness, [$h / (λ_{0} / n_{1})$]",
       fontsize=12)
plt.ylabel("ℛ",
       fontsize=14)
plt.title("n = 2.00")
plt.legend(loc='best')
# =============================================================================
# 
# =============================================================================
H, TMM, Theory = np.loadtxt("R_i7.txt",
                            unpack=True,
                            delimiter='       ',
                            skiprows=1)
fig1.add_subplot(gs[2:4,2])
plt.scatter(H[0:-1:20], TMM[0:-1:20], 
            c='k', label='TMM Evaluation',
            marker='x')
plt.plot(H, Theory, c='g', label='Theory Prediction')
plt.xlim(0,1)
plt.ylim(0.03,0.53)
plt.xlabel("Optical Thickness, [$h / (λ_{0} / n_{1})$]",
       fontsize=12)
plt.ylabel("ℛ",
       fontsize=14)
plt.title("n = 3.00")
plt.legend(loc='best')

plt.tight_layout()
plt.savefig("figure.png",
            bbox_inches='tight',
            dpi=300)