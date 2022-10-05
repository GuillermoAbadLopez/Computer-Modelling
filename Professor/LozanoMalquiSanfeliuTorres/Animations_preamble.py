# -*- coding: utf-8 -*-
"""
Created on Sat Jan  1 14:50:37 2022

@author: misae
"""

##############################################################################
# LIBRARIES
##############################################################################

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from matplotlib.animation import PillowWriter
from cycler import cycler

##############################################################################
# GRAPHS CONFIGURATIONS
##############################################################################

def graph_style():
    for param in ['figure.facecolor', 'axes.facecolor', 'savefig.facecolor']:
        #plt.rcParams[param] = '#212946'  # bluish dark grey
        plt.rcParams[param] = '#000000'
        
    for param in ['text.color', 'axes.labelcolor', 'xtick.color', 'ytick.color', 'axes.edgecolor']:
        plt.rcParams[param] = '0.9'  # very light grey
    colors = [
        '#08F7FE',  # teal/cyan
        '#FE53BB',  # pink
        '#F5D300',  # yellow
        '#00ff41', # matrix green
    ]
    plt.rcParams['axes.prop_cycle'] = cycler(color=colors)
    
    return None


plt.rcParams.update({"text.usetex": True})
pre = r"""  \usepackage{amsmath}
            \usepackage{amssymb}
            \usepackage[decimalsymbol=comma]{siunitx}"""
plt.rc('text.latex', preamble=pre)

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True

##############################################################################
# CONSTANTS
##############################################################################

hbar = 1.0
m_e = 1.0
a_0 = 1.0
q_e = 1.0
Eh = 1.0

Nx = 400
Ny = 300
M = 250
Lx = 20.0
Ly = 15.0
iterations = 100
a1 = 24
a2 = 25

i0 = 110

##############################################################################
# GRAPHING FUNCTIONS
##############################################################################

def get_positions(name):
    f = open(name, 'r')
    lines = f.readlines()
    f.close()
    
    spl_x = lines[0].split('\t')
    spl_y = lines[1].split('\t')
    
    xxs1 = np.zeros((Ny, Nx))
    yys1 = np.zeros((Ny, Nx))

    for j in range(Ny):
        for i in range(Nx):
            xxs1[j][i] = float(spl_x[i + Nx*j])
            yys1[j][i] = float(spl_y[i + Nx*j])
            
    xs = np.linspace(0.0, Lx, Nx + 1)
    ys = np.linspace(0, Ly, Ny + 1)
    xxs2, yys2 = np.meshgrid(xs, ys)
            
    return xxs1, yys1, xxs2, yys2


def get_V(name):
    f = open(name, 'r')
    lines = f.readlines()
    f.close()
    
    spl = lines[0].split('\t')
    Vs = np.zeros((Ny, Nx))
    
    for j in range(Ny):
        for i in range(Nx):
            Vs[j][i] = float(spl[i + Nx*j])
            
    return Vs


def get_psi(name):
    f = open(name, 'r')
    lines = f.readlines()
    f.close()
    
    spl_psi = lines[0].split('\t')
    psi_re = np.zeros((Ny, Nx))
    psi_im = np.zeros((Ny, Nx))
    psi_norm = np.zeros((Ny, Nx))
    
    for j in range(Ny):
        for i in range(Nx):
            psi_re[j][i] = float(spl_psi[2*(i + Nx*j)])
            psi_im[j][i] = float(spl_psi[2*(i + Nx*j) + 1])
            psi_norm[j][i] = np.sqrt(psi_re[j][i]**2 + psi_im[j][i]**2)
        
    return psi_re, psi_im, psi_norm


def get_results(name):
    f = open(name, 'r')
    lines = f.readlines()
    f.close()
    M = len(lines)
    print('Nombre d\'iteracions:', M)
    
    psi_re = np.zeros((M, Ny, Nx))
    psi_im = np.zeros((M, Ny, Nx))
    psi_norm = np.zeros((M, Ny, Nx))
    psi_deg = np.zeros((M, Ny, Nx))
    
    for n in range(M):
        spl_psi = lines[n].split('\t')
        for j in range(Ny):
            for i in range(Nx):
                psi_re[n][j][i] = float(spl_psi[2*(i + Nx*j)])
                psi_im[n][j][i] = float(spl_psi[2*(i + Nx*j) + 1])
                psi_norm[n][j][i] = np.sqrt(psi_re[n][j][i]**2 + psi_im[n][j][i]**2)
                psi_deg[n][j][i] = np.arctan2(psi_im[n][j][i], psi_re[n][j][i])

    return psi_re, psi_im, psi_norm, psi_deg, M

##############################################################################
# RANGE DETECTOR FUNCTIONS
##############################################################################

def get_t0(name):
    f = open(name, 'r')
    lines = f.readlines()
    f.close()
    
    psi_re = np.zeros((Ny, Nx))
    psi_im = np.zeros((Ny, Nx))
    psi_prob = np.zeros((Ny, Nx))
    psi_deg = np.zeros((Ny, Nx))
    
    spl_psi = lines[0].split('\t')
    for j in range(Ny):
        for i in range(Nx):
            psi_re[j][i] = float(spl_psi[2*(i + Nx*j)])
            psi_im[j][i] = float(spl_psi[2*(i + Nx*j) + 1])
            psi_prob[j][i] = psi_re[j][i]**2 + psi_im[j][i]**2
            psi_deg[j][i] = np.arctan2(psi_im[j][i], psi_re[j][i])
            
    return psi_re, psi_im, psi_prob, psi_deg

def get_diffraction(psi_prob):
    return [psi_prob[j][Nx - 10] for j in range(Ny)]

def generate_locations(Z):
    d1s = ['-d1', '-d2', '-d3', '-d4', '-d5', '-d6', '-d7', '-d8', '-d9']
    d2s = ['-00', '-50']
    locations = []
    
    for d1 in d1s:
        for d2 in d2s:
            locations.append('Z' + str(Z) + '/results' + d1 + d2 + '.txt')
    locations.append('Z' + str(Z) + '/results-d10-00.txt')
            
    return locations

def generate_destinations(Z):
    d1s = ['-d1', '-d2', '-d3', '-d4', '-d5', '-d6', '-d7', '-d8', '-d9']
    d2s = ['-00', '-50']
    destinations = []
    
    for d1 in d1s:
        for d2 in d2s:
            destinations.append('Z' + str(Z) + '/diffraction' + d1 + d2 + '.png')
    destinations.append('Z' + str(Z) + '/diffraction-d10-00.png')
            
    return destinations