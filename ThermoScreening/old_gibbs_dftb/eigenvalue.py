import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import os


# Inertia tensor initialization
I = [[11.0,-11.0,-7.0],[7.0,-7.0,-5.0],[3.0,-3.0,-1.0]]

### calculate the eigenvecotrs and eigenvalues ###

eigenvalues = np.linalg.eig(I) # in unit u*A^2
# eigenvalue
evl = np.abs(eigenvalues[0])
# eigenvector
evv = eigenvalues[1]

print("Eigenwerte und Eigenvektoren in u °A²: {}".format(eigenvalues))
print("\n")
print("Eigenwert berechnet in u °A²:{}".format(evl)) 
print("Eigenvektor in u °A²: {}".format(evv))
print("\n")
