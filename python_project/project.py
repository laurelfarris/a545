import math
import numpy as np
from scipy.interpolate import interp1d

# Getting error here... :(
#from astropy import constants as const


''' Constants '''
G = 6.67e-8
M_sun = 1.99e33
R_sun = 6.96e10
M_star = 2.02 * M_sun
R_star = 1.711 * R_sun


'''
(1)  Model inputs. All values in cgs units.
'''
T_eff = 5800.
g = G*M_sun/R_sun

''' Mass fractions of hydrogen (X), helium (Y), and "metals" (Z) '''
X = 0.70
Y = 0.28
Z = 0.02

''' tolerance level '''
epsilon = 0.001

''' Rosseland mean opacity table (given X, Y, and X) '''
opacity_table = np.loadtxt("opacity.txt")
T_6 = opacity_table[1:, 0]
R = opacity_table[0, 1:]
chi = opacity_table[1:, 1:]

'''
(2)  Boundary condition (pressure at top layer)
'''
const = (10**1.6)/(10**(8./3.))
P_top = const * g**(2./3.)
P_guess = np.zeros(2)


'''
(3)  Set up optical depth (tau) and temperature (T) arrays
'''
tau = np.logspace(-3, 2, 50)
#T = np.zeros(tau.size)

''' Hopf function q = q(tau). Using single value for now '''
#q_tau = ...
q_tau = 2./3.
T = ( (3./4.) * T_eff**4 * (tau + q_tau) )**(1./4.)


'''
(4)  Solve top layer using optical depth (tau), pressure (P), density (rho)
'''

''' initialize arrays '''
P = np.zeros(tau.size)
n_e = np.zeros(tau.size)
rho = np.zeros(tau.size)
opacity = np.zeros(tau.size)

''' Solve detailed balance of i=1 layer  '''
P[0] = P_top
n_e[0] = P_top
rho[0] = P_top


'''
(5)  Start looping through the layers
'''

for i in range(1, tau.size):

    ''' 6  Estimate pressure in layer, based on opacity of previous layer '''
    P_guess[0] = P[i-1] + g * ( tau[i] - tau[i-1] ) / ( opacity[i-1] )

    ''' 7  Solve detailed balancing using Pressure guess from (6) '''
    #opacity[i] = interpolation from opacity table

    ''' 8  Refine pressure guess '''
    P_guess[1] = P[i-1] + 2. * g * ( tau[i] - tau[i-1] ) / ( opacity[i-1] + opacity[i] )

    ''' 9  Check for convergence '''
    if ( (np.absolute( P_guess[1] - P_guess[0])) / P_guess[0] <= epsilon ):
        P[i] = P_guess[1]
        ''' solve for detailed balancing, then go to step 5 '''
    else:
        P_guess[0] = P_guess[1]
        ''' go to step 7 '''
