#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculator for Vibrations in Course Content Part B of MMAN2300

@author: Steven Watts
"""

import numpy as np
import math
from sympy import symbols, cos, sin, N, exp, solve, simplify
from sig import sig


# number of significant figures to use for answers
accuracy = 4
# acceleration due to gravity
g=9.81 #m/s^2


# defining symbols, do not change
t = symbols('t')
omega = symbols('w')

### FUNCTIONS

### FREE VIBRATION ###
def Omega_N(k, m):
    """
    Parameters
    ----------
    k : Input: k* or k_t*.
    m : Input: m* or I*.

    Returns
    -------
    natural frequency of oscillator - omega_n

    """
    return np.sqrt(k / m)


def Damp_Ratio(c, m, omega_n):
    """
    Parameters
    ----------
    c : Input c* or c_t*
    m : Input m* or I*
    omega_n : the natural frequency

    Returns
    -------
    damping ratio of oscillator

    """
    
    return c / (2*m*omega_n)


def Free_Vibration(x0, xdot0, omega_n, dampRatio):
    """
    Solution for free vibration equation of motion
    Determines if system is undamped, underdamped, critically damped or overdamped
    Finds the constants of the wave function

    Parameters
    ----------
    x0 : TYPE
        DESCRIPTION.
    xdot0 : TYPE
        DESCRIPTION.
    omega_n : TYPE
        DESCRIPTION.
    dampRatio : TYPE
        DESCRIPTION.

    Returns
    -------
    If underdamped :
    A : amplitude of wave function
    phi : phase factor of wave function
    
    If overdamped or critically damped :
        Result is two wave functions
    A1 : amplitude of first wave function
    A2 : amplitude of second wave function
    
    eqn: The Equation of Motion using SymPy

    """
    
   # underdamped
    # NOTE: the case of undamped is inclusive
    if dampRatio < 1:
        
        omega_d = omega_n*np.sqrt(1 - dampRatio**2)
        
        
        A = np.sqrt(x0**2 + (dampRatio*omega_n*x0 + xdot0)**2 / omega_d**2 )
        phi = -math.atan((dampRatio*omega_n*x0 + xdot0) / (omega_d*x0) )
        if x0 < 0:
            phi += np.pi
        
        A = sig(A, accuracy)
        phi = sig(phi, accuracy)
        
        eqn = A * exp(sig(-dampRatio*omega_n,accuracy)*t) * cos(sig(omega_d,accuracy)*t + phi )
        
        
        return A, phi, eqn
    
    # critically damped
    elif sig(dampRatio, 3) == 1:
        A1 = sig(x0, accuracy)
        A2 = sig(xdot0 + omega_n*x0, accuracy)
        
        eqn = (A1 + A2*t)*exp(-sig(dampRatio,accuracy)*sig(omega_n, accuracy)*t)
        return A1, A2, eqn
    
    # overdamped
    elif sig(dampRatio, accuracy) > 1:
        #denominator - same for both A1 and A2
        A_den = 2*omega_n*np.sqrt(dampRatio**2 - 1)
        
        #numerator for A1
        A1_num = (xdot0 + omega_n*x0*(dampRatio + np.sqrt(dampRatio**2 - 1)))
        A1 = sig(A1_num / A_den, accuracy)
        
        #numerator for A2
        A2_num = (xdot0 + omega_n*x0*(dampRatio - np.sqrt(dampRatio**2 - 1)))
        A2 = sig(A2_num / A_den, accuracy)
        
        exponent1 = (-dampRatio + np.sqrt(dampRatio**2 - 1) ) * omega_n
        exponent1 = sig(exponent1, accuracy)
        exponent2 = (-dampRatio - np.sqrt(dampRatio**2 - 1) ) * omega_n
        exponent2 = sig(exponent2, accuracy)
        
        eqn = A1*exp(exponent1*t) + A2*exp(exponent2*t)
        
        return A1, A2, eqn










### FORCED VIBRATION ###
def Forced_Vibration(F0, phi0, k, omegaN, dampRatio, omegaSub):
    """
    Solution for forced vibration motion
    Works for either translation or rotation

    Parameters
    ----------
    F : Input: F_0*, M_0*, inital force/moment
    phi0 : initial phase
    
    k : Input: k* or k_t* 
    omega_n : the natural frequency
    damp_ratio : damping ratio

    omega : variable frequency

    Returns
    -------
    Steady state amplitude

    """
    omegaRatio = omegaSub / omegaN
    
    X = (F0 / k) / np.sqrt((1 - omegaRatio**2 )**2 + (2*dampRatio * omegaRatio)**2 )
    X = sig(X, accuracy)
    phi = phi0 - math.atan((2*dampRatio*omegaRatio) / (1 - omegaRatio**2 ) )
    phi = sig(phi, accuracy)
    
    eqn = X*cos(omega*t + phi)
    return X, phi, eqn


### SYSTEM IDENTIFICATION
def OMEGA_D(N, t1, t2):
    """
    Parameters
    ----------
    N : number of periods
    t1 : time at first peak
    t2 : time at last peak ( t_(N+1) )
    
    Returns
    -------
    Damped Natural Frequency, omega_d

    """
    return (2*np.pi*N) / (t2 - t1)

def Log_Decrement(N, x1, x2):
    """

    Parameters
    ----------
    N : number of periods
    x1 : displacement at first peak
    x2 : displacement at last peak ( x_(N+1) )

    Returns
    -------
    The Logarithmic Decrement, delta. 

    """
    
    return 1/N * np.log(x1 / x2)


# other
def Force_Amplitude(X, k, dampRatio, omegaN, omegaSub):
    
    omegaRatio = omegaSub / omegaN
    
    F0 = X * k * np.sqrt( (1 - omegaRatio**2)**2 + (2*dampRatio*omegaRatio)**2 )
    
    return F0

### 2 DEGREES OF FREEDOM ###

def OMEGA_2DOF(K_STAR, M_STAR):
    
    matrix = K_STAR - omega**2 * M_STAR
    
    det = simplify(matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0])
    
    omega_matrix = solve(det, omega)
    
    omega1 = sig( abs(omega_matrix[0]), accuracy)
    omega2 = sig( abs(omega_matrix[1]), accuracy)
    
    return matrix, det, omega_matrix, omega1, omega2
    

def Mode_Shape(K_STAR, M_STAR):
    
    matrix = OMEGA_2DOF(K_STAR, M_STAR)[0]
    omega1 = OMEGA_2DOF(K_STAR, M_STAR)[-2]
    omega2 = OMEGA_2DOF(K_STAR, M_STAR)[-1]
    
    Ki = symbols('X')
    
    eqn0 = matrix[0][0] + matrix[0][1]*Ki
    
    eqn1 = eqn0.subs(omega, omega1)
    
    eqn2 = eqn0.subs(omega, omega2)
    
    mode1 = sig(solve(eqn1, Ki)[0], accuracy)
    mode2 = sig(solve(eqn2, Ki)[0], accuracy)
    
    return eqn1, eqn2, mode1, mode2


def Forced_Amplitude_2DOF(K_STAR, M_STAR, omegaSub, F0):
    
    matrix = K_STAR - omegaSub**2 * M_STAR
    
    matrixInv = np.linalg.inv(matrix)
    
    X1 = matrixInv[0][1]*F0
    X2 = matrixInv[1][1]*F0
    
    return X1, X2
    


### END OF FUNCTIONS


m = 100 #kg
I = 58.33 #kg
L = 2
c = 20
k = 100
F_0 = 10 #N
omegaSub = 2 #rad/s

N = 19
x1 = 8.42879
x2 = 1.40519
t1 = 1.56667
t2 = 115.867

delta = Log_Decrement(N, x1, x2)
dampRatioIdent = delta / np.sqrt(4*np.pi**2 + delta**2)


m_star = 292 # kg
c_star = 98 # Ns/m
k_star = 27700 #N/m #kg*m/s^2

x0 = 0.000001 # m
xdot0 = 2 # m/s
time = 0.3 #s

omegaN = Omega_N(k_star, m_star)
dampRatio = Damp_Ratio(c_star, m_star, omegaN)
dampRatio = 0



print("Omega_N: ", sig(omegaN, accuracy), "rad/s")
print("Damping Ratio: ", sig(100*dampRatio, accuracy),"%")



initialCond = Free_Vibration(x0, xdot0, omegaN, dampRatio)

if dampRatio < 1:
    
    A = initialCond[0]
    phi = initialCond[1]
    if sig(dampRatio, accuracy) == 0:
        print("System is Undamped")
    else:
        print("System is Underdamped")
    print("Free Vibration Amplitude: ", A, "m")
    print("Free Vibration Phase Factor: ", phi, "rad")
else:
    A1 = initialCond[0]
    A2 = initialCond[1]
    if sig(dampRatio, accuracy) == 1:
        print("System is critically damped")
    else:
        print("System is overdamped")
    print("Free Vibration Amplitude 1: ", A1, "m")
    print("Free Vibration Amplitude 2: ", A2, "rad")
    
eqn = initialCond[2].subs(t, time)



print("Free Vibration Equation: ", initialCond[2])
print("Position at t =",time,"s: ",sig(eqn,accuracy), "m")



M_0 = F_0*(L - L/4)
phi0 = 0

steadyState = Forced_Vibration(M_0, phi0, k_star, omegaN, dampRatio, omegaSub)



### Q10 of quiz 7 ###
L = 1.2
m = 12
I0 = 5.76
c = 10
k = 1010
g = 9.81

k_star = L*(-m*g*L/2 + k*L)
c_star = c*L**2
I_star = I0

theta = 0.01 #vibration amplitude
omegaSub = 22

omegaN = Omega_N(k_star, I_star)
dampRatio = Damp_Ratio(c_star, I_star, omegaN)




momentAmp = Force_Amplitude(theta, k_star, dampRatio, omegaN, omegaSub)

forceAmp = 2*momentAmp / L





## Q1 of qiz 9 ##

r = 1.2
I0 = 2.4
m2 = 2
c = 1
k_t = 240
k = 300

#matrix elements
K_STAR = np.zeros([2,2])
M_STAR = np.zeros_like(K_STAR)
C_STAR = np.zeros_like(K_STAR)

K_STAR[0][0] = k_t + k*r**2
K_STAR[0][1] = -k*r
K_STAR[1][0] = -k*r
K_STAR[1][1] = k


C_STAR[0][0] = c*r**2

M_STAR[0][0] = I0
M_STAR[1][1] = m2


## finding natural freq

temp = OMEGA_2DOF(K_STAR, M_STAR)

omega1 = temp[-2]
omega2 = temp[-1]


mode1 = Mode_Shape(K_STAR, M_STAR)[2]
mode2 = Mode_Shape(K_STAR, M_STAR)[3]
    
    


## Q6 of Quiz 9 ###
r = 1.1
I1 = 2.2
m2 = 2.5
k = 300

#matrix elements
K_STAR = np.zeros([2,2])
M_STAR = np.zeros_like(K_STAR)
C_STAR = np.zeros_like(K_STAR)

M_STAR[0][0] = I1
M_STAR[1][1] = m2

K_STAR[0][0] = 2*k*r**2
K_STAR[0][1] = k*r
K_STAR[1][0] = k*r
K_STAR[1][1] = 2*k

temp = OMEGA_2DOF(K_STAR, M_STAR)


omegaSub = np.sqrt(2*k*r**2 / I1 )
F0 = 28

steadyState = Forced_Amplitude_2DOF(K_STAR, M_STAR, omegaSub, F0)




