#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 21:41:24 2022

@author: Brice
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters of the problem ---------------------------------------------------

a = 0.05
b = 0.01

### FV Method -----------------------------------------------------------------

I = 30 # Number of cells
h1 = 1/I # "Volume" of a cell

# Construction of the FV matrix:
diag1 = np.ones(I-1)
A1 = 2*np.eye(I)-np.diag(diag1,-1)-np.diag(diag1,1)
A1[0][0] += 1
A1[-1][-1] += 1
A1 /= h1**2
A1 += np.eye(I)

# Right hand side (f ≡ 1):
B1 = np.ones(I)
B1[0] += 2*a/h1/h1
B1[-1] += 2*b/h1/h1

# Computation of the solution:
X1 = np.linspace(0, 1, I+1)[:-1]+h1/2 # Points of P
U1 = np.linalg.solve(A1, B1) # Solution on each cell

### FD Method -----------------------------------------------------------------

N = 31 # Number of points
X2 = np.linspace(0, 1, N)
h2 = X2[1]-X2[0] # Space step

# Construction of the FD matrix:
diag2 = np.ones(N-3)
A2 = 2*np.eye(N-2)-np.diag(diag2,-1)-np.diag(diag2,1)
A2 /= h2**2
A2 += np.eye(N-2)

# Right hand side (f ≡ 1):
B2 = np.ones(N-2)
B2[0] += a/h2/h2
B2[-1] += b/h2/h2

# Computation of the solution:
U2 = np.linalg.solve(A2, B2) # Solution at the points of the grid

### If we want to take into account the nonlinearity --------------------------

# The nonlinear system can be rewritten as F(U) := AU-B = 0.
# Let's apply Newton's method to solve it:
    
F = lambda U : A2@U + np.sin(U) - B2
jacF = lambda U : A2 + np.diag(np.cos(U))
    
Uold = U2
Unew = Uold - np.linalg.inv(jacF(Uold))@F(Uold)
    
while (np.linalg.norm(Uold-Unew) > 1e-10):
    Uold = Unew
    Unew = Uold - np.linalg.inv(jacF(Uold))@F(Uold)

### Plot of the solutions -----------------------------------------------------

plt.step(X1, U1, where='mid', label="FV solution")
plt.plot(X2[1:-1], U2, '.', label="FD solution")
plt.plot(X2[1:-1], Unew, label="FD w/ nonlinearity")
plt.legend()
plt.show()