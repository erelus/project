import numpy as np

# Matrix H0

we = 2860.3691
wx = 1.35066959

def Dirac_delta(i,j):
  return int(i==j)

H0 = np.zeros([64,64])

def sept_to_dec(i,j):
  return i*(8) + j

for i in range(8):
  for n in range(8):
    for j in range(8):
      for m in range(8):
        H0[sept_to_dec(i,n),sept_to_dec(j,m)] = (we*(j+0.5) - wx*(j+0.5)**2)*Dirac_delta(i,j)*Dirac_delta(n,m) + (we*(m+0.5) - wx*(m+0.5)**2)*Dirac_delta(i,j)*Dirac_delta(n,m)
        #H0[sept_to_dec(i,n),sept_to_dec(j,m)] = str(i)+str(n)+str(j)+str(m), i did this just to check whether the matrix is correct

print('H0:')
print(H0)

import math

# Matrix H1

g = 0.03094417965
nu = 0.9771277
f = -987.05

H1 = np.zeros([64,64])

for i in range(8):
  for n in range(8):
    for j in range(8):
      for m in range(8):
        H1[sept_to_dec(i,n),sept_to_dec(j,m)] = ((f/(2*nu*we)) - (g*nu*we)/2)*(math.sqrt(m*j)*Dirac_delta(i,j-1)*Dirac_delta(n,m-1))
        + ((f/(2*nu*we)) + (g*nu*we)/2)*(math.sqrt(m*(j+1))*Dirac_delta(i,j+1)*Dirac_delta(n,m-1))
        + ((f/(2*nu*we)) + (g*nu*we)/2)*(math.sqrt((m+1)*j)*Dirac_delta(i,j-1)*Dirac_delta(n,m+1))
        + ((f/(2*nu*we)) - (g*nu*we)/2)*(math.sqrt((m+1)*(j+1))*Dirac_delta(i,j+1)*Dirac_delta(n,m+1))

print('H1:')
print(H1)

np.sum(1*(H1!=0))

H = H0+H1
print('H = H0 + H1')
print('H:')
print(H)

# eigenvalus and eigenvectors
eigval, eigvect = np.linalg.eig(H)

# Print the eigenvalues
print("Eigenvalues:")
print(eigval)

#print("Eigenvectors:")
#print(eigvect)

normalised_eigvect = []
for i in range(64):
  normalised_eigvect.append(eigvect[:,i]/np.linalg.norm(eigvect[:,i]))

eig_vect = np.array(normalised_eigvect).T

# Print the normalized Eigenvectors
print("Eigenvectors:")
print(eig_vect)




