# Numerical schemes for simulating linear advection for outer code
# linearAdvect.py


# The numpy package for numerical functions and pi
import numpy as np

def FTCS(phiOld, c, nt):
    "Linear advection of profile in phiOld using FTCS, Courant number c"
    "for nt time-steps"

    nx = len(phiOld)

    # new time-step array for phi
    phi = phiOld.copy()

    # FTCS for each time-step
    for it in range(nt):
        # Loop through all space using remainder after division (%)
        # to cope with periodic boundary conditions
        for j in range(nx):
            phi[j] = phiOld[j] - 0.5*c*\
                     (phiOld[(j+1)%nx] - phiOld[(j-1)%nx])

        # update arrays for next time-step
        phiOld = phi.copy()

    return phi


def FTBS(phiOld, c, nt):
    "Linear advection of profile in phiOld using FTBS, Courant number c"
    "for nt time-steps"

    nx = len(phiOld)

    # new time-step array for phi
    phi = phiOld.copy()


    # FTCS for each time-step
    for it in range(nt):
        # Loop through all space using remainder after division (%)
        # to cope with periodic boundary conditions


        for j in range(nx):
            phi[j] = phiOld[j] - c*\
                     (phiOld[(j)] - phiOld[(j-1)%nx])

        # update arrays for next time-step
        phiOld = phi.copy()

    return phi


def CTCS(phiOld, c, nt):
    "Linear advection of profile in phiOld using CTCS, Courant number c"
    "for nt time-steps"

    nx = len(phiOld)

    if nt == 0:
        phi = phiOld

    else:
    # value of phi generated  at time t_1
        phiOld2 = FTCS(phiOld, c, 1)

        phi = phiOld2.copy()

        for it in range(nt-1):

            for j in range(nx):
                phi[j] = phiOld[j] - c*\
                        (phiOld2[(j+1)%nx] - phiOld2[(j-1)%nx])

            phiOld = phiOld2.copy()
            phiOld2 = phi.copy()


    return phi




def BTCS(phiOld, c, nt):
    "Linear advection of profile in phiOld using BTCS, Courant number c"
    "for nt time-steps"


    nx = len(phiOld)

    #create empty array
    A = np.zeros(shape=(nx,nx))

    #fill diagonal with elemtents 1
    np.fill_diagonal(A, 1)

    #fill 'corner' elements as +-c/2
    A[0,nx-1] = -c/2
    A[nx-1,0] = c/2

    #fill elements above/below the diagonal with +-c/2
    for j in range(nx-1):
        A[j,j+1] = c/2
        A[j+1,j] = -c/2


    for it in range(nt):
        phi = np.linalg.solve(A,phiOld)
        phiOld = phi.copy()

    return phi


def LaxWendroff(phiOld, c, nt):
    "Linear advection of profile in phiOld using Lax-Wendroff"
    "Courant number c, for nt time-steps"

    nx = len(phiOld)

    phi = phiOld.copy()

    for it in range(nt):

        for j in range(nx):

            phi[j] = (1-c**2)*phiOld[j] - (c/2)*(1-c)*phiOld[(j+1)%nx] + \
                    (c/2)*(1+c)*phiOld[(j-1)%nx]

        phiOld = phi.copy()

    return phi


def tvFTBS(phiOld, c, nt):
    "calculation of total variation (at each time step) for the"
    "linear advection of profile in phiOld using FTBS, Courant number c"
    "for nt time-steps"

    nx = len(phiOld)

    # new time-step array for phi
    phi = phiOld.copy()
    #initialise TV as an array
    TV_FTBS = np.zeros(nt+1)
    for j in range(nx):
        TV_FTBS[0] += abs(phiOld[j] - phiOld[(j-1)%nx])

    # FTCS for each time-step
    for it in range(nt):
        # Loop through all space using remainder after division (%)
        # to cope with periodic boundary conditions

        #initialise total variation for time step nt
        TV = 0

        for j in range(nx):
            phi[j] = phiOld[j] - c*\
                     (phiOld[(j)] - phiOld[(j-1)%nx])
            TV += abs(phi[j] - phi[(j-1)%nx])

        #calculate total variation for this time step
        TV_FTBS[it+1] = TV
        # update arrays for next time-step
        phiOld = phi.copy()

    return TV_FTBS


def tvCTCS(phiOld, c, nt):
    "calculation of total variation (at each time step) for the"
    "linear advection of profile in phiOld using CTCS, Courant number c"
    "for nt time-steps"

    nx = len(phiOld)
    TV_CTCS = np.zeros(nt+1)
    for j in range(nx):
        TV_CTCS[0] += abs(phiOld[j] - phiOld[(j-1)%nx])
    if nt == 0:
        phi = phiOld

    else:
        # value of phi generated  at time t_1
        phiOld2 = FTCS(phiOld, c, 1)
        for j in range(nx):
            TV_CTCS[1] += abs(phiOld2[j] - phiOld2[(j-1)%nx])

        phi = phiOld2.copy()

        for it in range(nt-1):

            #initialise TV
            TV = 0

            for j in range(nx):
                phi[j] = phiOld[j] - c*\
                        (phiOld2[(j+1)%nx] - phiOld2[(j-1)%nx])
                TV += abs(phi[j] - phi[(j-1)%nx])

            TV_CTCS[it+2] = TV
            phiOld = phiOld2.copy()
            phiOld2 = phi.copy()

    return TV_CTCS


def tvLaxWendroff(phiOld, c, nt):
    "calculation of total variation (at each time step) for the"
    "linear advection of profile in phiOld using Lax-Wendroff"
    "Courant number c, for nt time-steps"

    nx = len(phiOld)
    TV_LW = np.zeros(nt+1)
    for j in range(nx):
        TV_LW[0] += abs(phiOld[j] - phiOld[(j-1)%nx])

    phi = phiOld.copy()

    for it in range(nt):
        TV = 0
        for j in range(nx):

            phi[j] = (1-c**2)*phiOld[j] - (c/2)*(1-c)*phiOld[(j+1)%nx] + \
                    (c/2)*(1+c)*phiOld[(j-1)%nx]
            TV += abs(phi[j] - phi[(j-1)%nx])
        TV_LW[it+1] = TV
        phiOld = phi.copy()

    return TV_LW
