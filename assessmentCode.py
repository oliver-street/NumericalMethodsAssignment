import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from initialConditions import *
from advectionSchemes import *
from diagnostics import *


def main():
    # Parameters
    xmin = 0
    xmax = 1
    nx = 40
    nt = 30
    c = 0.9

    # Derived parameters
    dx = (xmax - xmin)/nx

    # spatial points for plotting and for defining initial conditions
    x = np.arange(xmin, xmax, dx)

    t = np.arange(nt+1)

    # Initial conditions
    phiOld = squareWave(x, 0.2, 0.4)
    # Exact solution is the initial condition shifted around the domain
    phiAnalytic = squareWave((x - c*nt*dx)%(xmax - xmin), 0.2, 0.4)


    # Plot the solutions
    font = {'size'   : 20}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(t, tvFTBS(phiOld,c,nt), ".-",label='TV FTBS', color='black')
    plt.plot(t, tvLaxWendroff(phiOld,c,nt), ".-", label='TV LW', color='red')
    plt.plot(t, tvCTCS(phiOld,c,nt), ".-", label='TV CTCS', color='blue')
    #plt.plot(x, phiAnalytic, label='Analytic', color='black',
            # linestyle='--', linewidth=2)
    plt.axhline(2, linestyle=':', color='black')
    plt.legend(bbox_to_anchor=(1.15 , 1.1))
    plt.ylim(0,8)
    plt.xlabel('$t$')
    plt.ylabel('$TV$')
    input('press return to continue')
    plt.savefig('plots/TV.pdf')
    plt.show('plots/TV.pdf')

### Run the function main defined in this file                      ###
main()
