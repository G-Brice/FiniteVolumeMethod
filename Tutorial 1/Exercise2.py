#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 17:09:14 2022

@author: Brice
"""

import numpy as np
import matplotlib.pyplot as plt

### FV Method -----------------------------------------------------------------

I = 20 # Number of cells
h1 = 1/I # "Volume" of a cell

# Construction of the FV matrix:
diag1 = np.ones(I-1)
A1 = 2*np.eye(I)-np.diag(diag1,-1)-np.diag(diag1,1)
A1[0][0] += 1
A1[-1][-1] += 1
A1 /= h1**2

# Computation of the solution :
X1 = np.linspace(0, 1, I+1)[:-1]+h1/2 # Points of P
U1 = np.linalg.solve(A1, np.ones(I)) # Solution on each cell

### FD Method -----------------------------------------------------------------

N = 21 # Number of points
X2 = np.linspace(0, 1, N)
h2 = X2[1]-X2[0] # Space step

# Construction of the FD matrix :
diag2 = np.ones(N-3)
A2 = 2*np.eye(N-2)-np.diag(diag2,-1)-np.diag(diag2,1)
A2 /= h2**2

# Computation of the solution:
U2 = np.linalg.solve(A2, np.ones(N-2)) # Solution at the points of the grid

### Plot of the solutions -----------------------------------------------------

plt.plot(X2, 0.5*X2*(1-X2), label="Exact solution")
plt.step(X1, U1, where='mid', label="FV solution")
plt.plot(X2[1:-1], U2, 'o', label="FD solution")
plt.legend()
plt.show()