#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 21:41:24 2022

@author: Brice
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters of the problem ---------------------------------------------------

a = 0.5
b = 0.4

### FV Method -----------------------------------------------------------------

I = 20 # Number of cells
h1 = 1/I # "Volume" of a cell

# Construction of the FV matrix:
diag1 = np.ones(I-1)
A1 = 2*np.eye(I)-np.diag(diag1,-1)-np.diag(diag1,1)
A1[0][0] += 1
A1[-1][-1] += 1-4/(2+h1)
A1 /= h1**2
A1 += 2*np.eye(I)

# Right hand side:
B1 = h1*(np.arange(1, I+1)-1/2)
B1[0] += 2/h1/h1

# Computation of the solution:
X1 = np.linspace(0, 1, I+1)[:-1]+h1/2 # Points of P
U1 = np.linalg.solve(A1, B1) # Solution on each cell

### FD Method -----------------------------------------------------------------

N = 21 # Number of points
X2 = np.linspace(0, 1, N)
h2 = X2[1]-X2[0] # Space step

# Construction of the FD matrix:
diag2 = np.ones(N-3)
A2 = 2*np.eye(N-2)-np.diag(diag2,-1)-np.diag(diag2,1)
A2[-1][-1] -= 1/(h2+1)
A2 /= h2**2
A2 += 2*np.eye(N-2)

# Right hand side:
B2 = h2*np.arange(1, N-1)
B2[0] += 1/h2/h2

# Computation of the solution:
U2 = np.linalg.solve(A2, B2) # Solution at the points of the grid

### Plot of the solutions -----------------------------------------------------

plt.step(X1, U1, where='mid', label="FV solution")
plt.plot(X2[1:-1], U2, 'o', label="FD solution")
plt.legend()
plt.show()