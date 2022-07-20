# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import sig as s
import math
import numpy as np

mu_N = 5.050783699e-27 #nuclear magneton

mu_F = 2.626868*mu_N #flouride

mu_p = 1.4106067873e-26 # J / T
mu_e = 9.284764e-24 # J / T

g = 2 # g factor
h = 6.62607015e-34 # Js
k_B = 1.380649e-23

f = 14e6 #Hz
f_e = 50e6 #Hz - electron

def B_0(f, mu):
    return (h * f) / (g * mu) * 10**4

def f_0(B, mu):
    return (g * B * mu) / h


B = 0.3
T = 25 + 273.15
energy_ratio = math.exp(-g*mu_p*B / (k_B * T))

energy_ratio = s.sig(energy_ratio, 6)

B_0_p = B_0(f, mu_p)
B_0_e = B_0(f_e, mu_e)
B_0_F = B_0(f, mu_F)

f_0_F = f_0(B_0_p, mu_F)