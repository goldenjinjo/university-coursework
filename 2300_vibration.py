#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculator for Vibrations in Course Content Part B of MMAN2300

@author: Steven Watts
"""

import numpy as np
import math
from sympy import symbols, cos, sin, N, exp
from sig import sig


# number of significant figures to use for answers
accuracy = 4
# acceleration due to gravity
g=9.81 #m/s^2


# defining symbols, do not change
t = symbols('t')

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
    omega_d = omega_n*np.sqrt(1 - dampRatio**2)
   
     
   
   # underdamped
    # NOTE: the case of undamped is inclusive
    if dampRatio < 1:
        A = np.sqrt(x0**2 + (dampRatio*omega_n*x0 + xdot0)**2 / omega_d**2 )
        phi = -math.atan((dampRatio*omega_n*x0 + xdot0) / (omega_d*x0) )
        if x0 < 0:
            phi += np.pi
        
        eqn = sig(A,accuracy) * exp(sig(-dampRatio*omega_n,accuracy)*t) * cos(sig(omega_d,accuracy)*t + sig(phi,accuracy) )
        
        
        return A, phi, eqn
    
    # critically damped
    elif sig(dampRatio, 3) == 1:
        A1 = x0
        A2 = xdot0 + omega_n*x0
        
        #eqn = 
        return A1, A2
    
    # overdamped
    elif sig(dampRatio, accuracy) > 1:
        #denominator - same for both A1 and A2
        A_den = 2*omega_n*np.sqrt(dampRatio**2 - 1)
        
        #numerator for A1
        A1_num = (xdot0 + omega_n*x0(dampRatio + np.sqrt(dampRatio**2 - 1)))
        A1 = A1_num / A_den
        
        #numerator for A2
        A2_num = (xdot0 + omega_n*x0(dampRatio - np.sqrt(dampRatio**2 - 1)))
        A2 = A2_num / A_den
        
        
        return A1, A2


# symbol testing







### FORCED VIBRATION ###
def Forced_Vibration(F0, phi0, k, omega_n, dampRatio, omega):
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
    omegaRatio = omega / omega_n
    
    X = (F0 / k) / np.sqrt((1 - omegaRatio**2 )**2 + (2*dampRatio * omegaRatio)**2 )
    
    phi = phi0 - math.atan((2*dampRatio*omegaRatio) / (1 - omegaRatio**2 ) )
    
    return X, phi


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

### END OF FUNCTIONS


m = 100 #kg
I = 58.33 #kg
L = 2
c = 20
k = 100
F_0 = 10 #N
omega = 2 #rad/s

N = 5
x1 = 0.1
x2 = 0.4


delta = Log_Decrement(N, x1, x2)
dampRatioIdent = delta / np.sqrt(4*np.pi**2 + delta**2)


m_star = 292 # kg
c_star = 98 # Ns/m
k_star = 27700 #N/m #kg*m/s^2

x0 = 0.000001 # m
xdot0 = 2 # m/s
time = 1.6 #s

omegaN = Omega_N(k_star, m_star)
dampRatio = Damp_Ratio(c_star, m_star, omegaN)

initialCond = Free_Vibration(x0, xdot0, omegaN, dampRatio)
A = initialCond[0]
phi = initialCond[1]
eqn = initialCond[2].subs(t, time)


print("Omega_N: ", sig(omegaN, accuracy), "rad/s")
print("Damping Ratio: ", sig(100*dampRatio, accuracy),"%")
print("Free Vibration Amplitude: ", sig(A, accuracy), "m")
print("Free Vibration Phase Factor: ", sig(phi, accuracy), "rad")
print("Free Vibration Equation: ", initialCond[2])
print("Position at t =",time,"s: ",sig(eqn,accuracy), "m")



M_0 = F_0*(L - L/4)
phi0 = 0

steadyState = Forced_Vibration(M_0, phi0, k_star, omegaN, dampRatio, omega)

    