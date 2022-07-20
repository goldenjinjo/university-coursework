#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 15:01:01 2022

@author: jinjo
"""

import numpy as np
import math

g=9.81
m = 100 #kg
I = 58.33 #kg
L = 2
c = 20
k = 100
F_0 = 10 #N
omega = 2 #rad/s


I_star = I

c_star = c*(L/4)**2

k_star = k*(L/2)**2 + L/4*m*g


def omega_n(k, m):
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


def damp_ratio(c, m, omega_n):
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


def forced_amplitude(F, k, omega_n, damp_ratio, omega):
    """
    Equation to calculate steady state amplitude of forced vibration motion
    Works for either translation or rotation

    Parameters
    ----------
    F : Input: F_0*, M_0*, inital force/moment
    k : Input: k* or k_t* 
    omega_n : the natural frequency
    damp_ratio : damping ratio

    omega : variable frequency

    Returns
    -------
    Steady state amplitude

    """
    omega_ratio = omega / omega_n
    
    X = (F / k) / np.sqrt((1 - omega_ratio**2 )**2 + (2*damp_ratio * omega_ratio)**2 )
    
    return X


omegaN = omega_n(k_star, I_star)
dampRatio = damp_ratio(c_star, I_star, omegaN)


M_0 = F_0*(L - L/4)

steadyState = forced_amplitude(M_0, k_star, omegaN, dampRatio, omega)

    