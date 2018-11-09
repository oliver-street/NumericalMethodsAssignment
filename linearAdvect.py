#!/usr/bin/python3

# Outer code for setting up the linear advection problem on a uniform
# grid and calling the function to perform the linear advection and plot.



### The matplotlib package contains plotting functions              ###
import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from initialConditions import *
from advectionSchemes import *
from diagnostics import *

### The main code is inside a function to avoid global variables    ###
def main():
    "Advect the initial conditions using various advection schemes and"
    "compare results"
    "Comment out any schemes which you do not wish to impliment and"
    "adjust the file name appropriately"

    # Parameters
    xmin = 0
    xmax = 1
    nx = 150
    nt = 40
    c = 0.2

    # Derived parameters
    dx = (xmax - xmin)/nx

    # spatial points for plotting and for defining initial conditions
    x = np.arange(xmin, xmax, dx)

    # Initial conditions
    phiOld = squareWave(x, 0.2, 0.4)
    # Exact solution is the initial condition shifted around the domain
    phiAnalytic = squareWave((x - c*nt*dx)%(xmax - xmin), 0.2, 0.4)

    # Advect the profile using finite difference for all the time steps
    phiCTCS = CTCS(phiOld.copy(), c, nt)
    phiFTCS = FTCS(phiOld.copy(), c, nt)
    phiFTBS = FTBS(phiOld.copy(), c, nt)
    phiFTFS = FTFS(phiOld.copy(), c, nt)
    phiCTFS = CTFS(phiOld.copy(), c, nt)
    phiCTBS = CTBS(phiOld.copy(), c, nt)
    phiBTCS = BTCS(phiOld.copy(), c, nt)
    phiLaxWendroff = LaxWendroff(phiOld.copy(), c, nt)

    # Calculate and print out error norms
    print("CTCS l2 error norm = ", l2ErrorNorm(phiCTCS, phiAnalytic))
    print("CTCS linf error norm = ", lInfErrorNorm(phiCTCS, phiAnalytic))
    #print("FTBS l2 error norm = ", l2ErrorNorm(phiFTBS, phiAnalytic))
    #print("FTBS linf error norm = ", lInfErrorNorm(phiFTBS, phiAnalytic))
    #print("BTCS l2 error norm = ", l2ErrorNorm(phiBTCS, phiAnalytic))
    #print("BTCS linf error norm = ", lInfErrorNorm(phiBTCS, phiAnalytic))
    print("Lax-Wendroff l2 error norm = ", l2ErrorNorm(phiLaxWendroff, phiAnalytic))
    print("Lax-Wendroff linf error norm = ", lInfErrorNorm(phiLaxWendroff, phiAnalytic))

    # Plot the solutions
    font = {'size'   : 20}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black',
             linestyle='--', linewidth=2)
    plt.plot(x, phiCTCS, label='CTCS', color='red')
    #plt.plot(x, phiFTBS, label='FTBS', color='green')
    #plt.plot(x, phiBTCS, label='BTCS', color='blue')
    plt.plot(x, phiLaxWendroff, label='Lax-Wendroff', color='magenta')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2,1.2])
    plt.legend(bbox_to_anchor=(1.15 , 1.1))
    plt.xlabel('$x$')
    input('press return to save file and continue')
    plt.savefig('plots/LaxWendroff_FTBS_CTCS.pdf')
    plt.show('plots/LaxWendroff_CTCS.pdf')

### Run the function main defined in this file                      ###
main()
