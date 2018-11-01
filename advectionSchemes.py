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


def FTFS(phiOld, c, nt):
    "Linear advection of profile in phiOld using FTFS, Courant number c"
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
                     (phiOld[(j+1)%nx] - phiOld[(j)%nx])

        # update arrays for next time-step
        phiOld = phi.copy()

    return phi


def CTFS(phiOld, c, nt):
    "Linear advection of profile in phiOld using CTFS, Courant number c"
    "for nt time-steps"

    nx = len(phiOld)

    # value of phi generated  at time t_1
    phiOld2 = FTFS(phiOld, c, 1)

    phi = phiOld2.copy()

    for it in range(nt-1):

        for j in range(nx):
            phi[j] = phiOld[j] - 2*c*\
                     (phiOld2[(j+1)%nx] - phiOld2[(j)%nx])

        phiOld = phiOld2.copy()
        phiOld2 = phi.copy()

    return phi


def CTBS(phiOld, c, nt):
    "Linear advection of profile in phiOld using CTBS, Courant number c"
    "for nt time-steps"

    nx = len(phiOld)

    # value of phi generated  at time t_1
    phiOld2 = FTBS(phiOld, c, 1)

    phi = phiOld2.copy()

    for it in range(nt-1):

        for j in range(nx):
            phi[j] = phiOld[j] - 2*c*\
                     (phiOld2[(j)%nx] - phiOld2[(j-1)%nx])

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
