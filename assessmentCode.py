import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from initialConditions import *
from advectionSchemes import *
from diagnostics import *

def mainLinearAdvect():
    "Advect the initial conditions using various advection schemes and"
    "compare results"

    # Parameters
    xmin = 0
    xmax = 1
    nx = 100
    nt = 100
    c = 0.4

    # Derived parameters
    dx = (xmax - xmin)/nx

    # spatial points for plotting and for defining initial conditions
    x = np.arange(xmin, xmax, dx)

    # Initial conditions
    phiOld = cosBell(x, 0.2, 0.5)

    # Exact solution is the initial condition shifted around the domain
    phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0.2, 0.5)

    # Advect the profile using finite difference for all the time steps
    phiCTCS = CTCS(phiOld.copy(), c, nt)
    phiFTCS = FTCS(phiOld.copy(), c, nt)
    phiFTBS = FTBS(phiOld.copy(), c, nt)
    phiBTCS = BTCS(phiOld.copy(), c, nt)
    phiLaxWendroff = LaxWendroff(phiOld.copy(), c, nt)

    # Calculate and print out error norms
    print("CTCS l2 error norm = ", l2ErrorNorm(phiCTCS, phiAnalytic))
    print("CTCS linf error norm = ", lInfErrorNorm(phiCTCS, phiAnalytic))
    print("FTBS l2 error norm = ", l2ErrorNorm(phiFTBS, phiAnalytic))
    print("FTBS linf error norm = ", lInfErrorNorm(phiFTBS, phiAnalytic))
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
    plt.plot(x, phiFTBS, label='FTBS', color='green')
    plt.plot(x, phiLaxWendroff, label='Lax-Wendroff', color='magenta')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2,1.2])
    plt.legend(bbox_to_anchor=(1.15 , 1.1))
    plt.xlabel('$x$')
    input('press return to save file and continue')
    plt.savefig('plots/LaxWendroff_FTBS_CTCS_cosBell.pdf')
    plt.show('plots/LaxWendroff_FTBS_CTCS_cosBell.pdf')


    # New Parameters
    xmin = 0
    xmax = 1
    nx = 200
    nt = 100
    c = 0.1

    #  New derived parameters
    dx = (xmax - xmin)/nx

    # spatial points for plotting and for defining initial conditions
    x = np.arange(xmin, xmax, dx)

    # New Initial conditions
    phiOld = squareWave(x, 0.2, 0.5)

    # Exact solution is the New initial condition shifted around the domain
    phiAnalytic = squareWave((x - c*nt*dx)%(xmax - xmin), 0.2, 0.5)

    # Advect the profile using finite difference for all the time steps
    phiCTCS = CTCS(phiOld.copy(), c, nt)
    phiFTCS = FTCS(phiOld.copy(), c, nt)
    phiFTBS = FTBS(phiOld.copy(), c, nt)
    phiBTCS = BTCS(phiOld.copy(), c, nt)
    phiLaxWendroff = LaxWendroff(phiOld.copy(), c, nt)

    # Calculate and print out error norms
    print("CTCS l2 error norm = ", l2ErrorNorm(phiCTCS, phiAnalytic))
    print("CTCS linf error norm = ", lInfErrorNorm(phiCTCS, phiAnalytic))
    print("FTBS l2 error norm = ", l2ErrorNorm(phiFTBS, phiAnalytic))
    print("FTBS linf error norm = ", lInfErrorNorm(phiFTBS, phiAnalytic))
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
    plt.plot(x, phiFTBS, label='FTBS', color='green')
    plt.plot(x, phiLaxWendroff, label='Lax-Wendroff', color='magenta')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2,1.2])
    plt.legend(bbox_to_anchor=(1.15 , 1.1))
    plt.xlabel('$x$')
    input('press return to save file and continue')
    plt.savefig('plots/LaxWendroff_FTBS_CTCS_Squarewave.pdf')
    plt.show('plots/LaxWendroff_FTBS_CTCS_SquareWave.pdf')

### Run the function main defined in this file                      ###
mainLinearAdvect()



def mainTV():
    # Parameters
    xmin = 0
    xmax = 1
    nx = 40
    nt = 50
    c = 0.05

    # Derived parameters
    dx = (xmax - xmin)/nx

    # spatial points for plotting and for defining initial conditions
    x = np.arange(xmin, xmax, dx)

    t = np.arange(nt+1)

    # Initial conditions
    phiOld = squareWave(x, 0.2, 0.4)


    # Plot the solutions
    font = {'size'   : 20}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(t, tvFTBS(phiOld,c,nt), ".-",label='TV FTBS', color='black')
    plt.plot(t, tvLaxWendroff(phiOld,c,nt), ".-", label='TV LW', color='red')
    plt.plot(t, tvCTCS(phiOld,c,nt), ".-", label='TV CTCS', color='blue')
    plt.axhline(2, linestyle=':', color='black')
    plt.legend(bbox_to_anchor=(1.15 , 1.1))
    plt.ylim(0,8)
    plt.xlabel('$t$')
    plt.ylabel('$TV$')
    plt.savefig('plots/TV_CourantNumber05.pdf')
    plt.show('plots/TV_CourantNumber05.pdf')

    c = 0.2

    # Plot the solutions
    font = {'size'   : 20}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(t, tvFTBS(phiOld,c,nt), ".-",label='TV FTBS', color='black')
    plt.plot(t, tvLaxWendroff(phiOld,c,nt), ".-", label='TV LW', color='red')
    plt.plot(t, tvCTCS(phiOld,c,nt), ".-", label='TV CTCS', color='blue')
    plt.axhline(2, linestyle=':', color='black')
    plt.legend(bbox_to_anchor=(1.15 , 1.1))
    plt.ylim(0,8)
    plt.xlabel('$t$')
    plt.ylabel('$TV$')
    plt.savefig('plots/TV_CourantNumber2.pdf')
    plt.show('plots/TV_CourantNumber2.pdf')

    c = 0.4

    # Plot the solutions
    font = {'size'   : 20}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(t, tvFTBS(phiOld,c,nt), ".-",label='TV FTBS', color='black')
    plt.plot(t, tvLaxWendroff(phiOld,c,nt), ".-", label='TV LW', color='red')
    plt.plot(t, tvCTCS(phiOld,c,nt), ".-", label='TV CTCS', color='blue')
    plt.axhline(2, linestyle=':', color='black')
    plt.legend(bbox_to_anchor=(1.15 , 1.1))
    plt.ylim(0,8)
    plt.xlabel('$t$')
    plt.ylabel('$TV$')
    plt.savefig('plots/TV_CourantNumber4.pdf')
    plt.show('plots/TV_CourantNumber4.pdf')

    c = 0.9

    # Plot the solutions
    font = {'size'   : 20}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(t, tvFTBS(phiOld,c,nt), ".-",label='TV FTBS', color='black')
    plt.plot(t, tvLaxWendroff(phiOld,c,nt), ".-", label='TV LW', color='red')
    plt.plot(t, tvCTCS(phiOld,c,nt), ".-", label='TV CTCS', color='blue')
    plt.axhline(2, linestyle=':', color='black')
    plt.legend(bbox_to_anchor=(1.15 , 1.1))
    plt.ylim(0,8)
    plt.xlabel('$t$')
    plt.ylabel('$TV$')
    plt.savefig('plots/TV_CourantNumber9.pdf')
    plt.show('plots/TV_CourantNumber9.pdf')
### Run the function main defined in this file                      ###
mainTV()

def mainTime():
    xmin = 0
    xmax = 1
    nx = 150
    nt = 50
    c = 0.2

    # Derived parameters
    dx = (xmax - xmin)/nx

    # spatial points for plotting and for defining initial conditions
    x = np.arange(xmin, xmax, dx)

    # Initial conditions
    phiOld = squareWave(x, 0.2, 0.4)
    # Exact solution is the initial condition shifted around the domain
    phiAnalytic = squareWave((x - c*nt*dx)%(xmax - xmin), 0.2, 0.4)

    def wrapper(func, *args, **kwargs):
        def wrapped():
            return func(*args, **kwargs)
        return wrapped
    wrapped = wrapper(BTCS, phiOld, c, nt)
    import timeit
    print('BTCS times at', timeit.timeit(wrapped, number=1000))
    wrapped2 = wrapper(FTBS, phiOld, c, nt)
    print('FTBS times at', timeit.timeit(wrapped2, number=1000))

mainTime()


def mainl2():
    l2FTBS = np.zeros(16)
    l2LW = np.zeros(16)
    l2CTCS = np.zeros(16)
    t = np.asarray([1/n for n in np.arange(50,210,10)])
    # Parameters
    i = 0
    for nx in np.arange(50, 210, 10):
        xmin = 0
        xmax = 1
        u = 0.7
        # Derived parameters
        dx = (xmax - xmin)/nx
        dt = dx
        c = (u*dt)/dx
        nt = nx
        # spatial points for plotting and for defining initial conditions
        x = np.arange(xmin, xmax, dx)

        # Initial conditions
        phiOld = sineWave(x)
        # Exact solution is the initial condition shifted around the domain
        phiAnalytic = sineWave((x - c*nt*dx)%(xmax - xmin))

        l2FTBS[i] = l2ErrorNorm(FTBS(phiOld,c,nt), phiAnalytic)
        l2LW[i] = l2ErrorNorm(LaxWendroff(phiOld,c,nt), phiAnalytic)
        l2CTCS[i] = l2ErrorNorm(CTCS(phiOld,c,nt), phiAnalytic)
        i += 1


    # Plot the solutions
    font = {'size'   : 20}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(t, l2FTBS, ".-",label='l2 FTBS', color='black')
    plt.plot(t, l2LW, ".-",label='l2 LW', color='red')
    plt.plot(t, l2CTCS, ".-",label='l2 CTCS', color='blue')
    plt.plot(t, 58*t**2, ":",label='$y \sim 2x$', color='0.3')
    plt.plot(t, 10*t, "--",label='$y \sim x$', color='0.6')
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.xlabel('$dx$')
    plt.ylabel('$l2 Error$')
    input('press return to continue')
    plt.savefig('plots/l2Error.pdf')
    plt.show('plots/l2Error.pdf')

### Run the function main defined in this file                      ###
mainl2()
