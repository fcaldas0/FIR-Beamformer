"""
@author: Fernanda Caldas
"""
import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

c = 340
M = 11
N = 64
fs = 16000
d = 0.035
thetaBW = math.pi/6

#Weights Matrix
Mmin = M%2 + 2
pmax = np.int((M-Mmin)/2)
W = np.ones((N,M))
f_k = np.arange(0, fs, fs/N)
f_m_shape = np.int((M-Mmin)/2)+1 
f_m = np.zeros((f_m_shape))
find = np.zeros((f_m_shape))
w = np.zeros(N)

for m in range(0, f_m_shape):
    f_m[m] = c/(np.sin(thetaBW/2)*(M-2*m)*d)
    a = np.where(f_k >= f_m[m])
    find[m] = a[0][0]

for i in range(1, f_m_shape):
    W[np.int(find[i])-1:N, i-1] = 0
    W[np.int(find[i])-1:N, -i] = 0
    
W_esp = np.copy(W)
    
for i in range(np.int(N/2)):
    W_esp[-(i+1),:] = W_esp[i,:]
    
for k in range(0, N):
    w[k] = np.sum(np.absolute(W_esp[k][:]))
    W_esp[k] = W_esp[k][:]/w[k]

#Array manifold vector
aux = fs/N
Ts = 1/fs
v_s = np.zeros((N*M, 91), dtype=complex)

for i in range(0, N*M):
    for j in range(0, 91):
        theta = (j*2 - 90)*math.pi/180
        m = i%M
        k = i/M
        tau = np.sin(theta)*m*d/c
        arg = -2*math.pi*tau*k/(N*Ts)
        v_s[i][j] = cmath.exp(-complex(0,2*math.pi*tau*k/(N*Ts)))
        
#Beampattern
D = np.zeros((N, 91))

for k in range(0, N):
    for theta in range(0, 91):
        D[-(k+1)][theta] = 20*np.log10(np.absolute(np.inner(np.conjugate(W_esp[k,:]), v_s[k*M:k*M+M,theta])))    

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

D = np.zeros((N, 91))

for k in range(0, N):
    #print(k)
    for theta in range(0, 91):
        D[-(k+1)][theta] = 20*np.log10(np.absolute(np.inner(np.conjugate(W_esp[k,:]), v_s[k*M:k*M+M,theta])))

plt.figure(figsize=(20,20))
ax = plt.gca()
im = ax.imshow(D, cmap='jet')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
#im.set_clim(vmin=-3, vmax=0)