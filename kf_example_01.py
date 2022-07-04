#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  2 12:21:25 2022

@author: 1517suj
"""

import numpy as np
from numpy import dot
import matplotlib.pyplot as plt
from numpy.linalg import inv


# sampling time
T = 1

A = np.array([[1, T, 0, 0],
              [0, 1, 0, 0],
              [0, 0, 1, T],
              [0, 0, 0, 1]])

F = np.array([[T**2/2, 0],
              [T,     0],
              [0,     T**2/2],
              [0,     T]])

C = np.array([[1, 0, 0, 0],
              [0, 0, 1, 0]])

G = np.diag([1, 1])

V = np.diag([0.85, 0.3])

W = np.diag([0.28, 0.65])


x = []
y = []

# initial x value
x_k = np.array([[0],
                [5], 
                [0], 
                [7]])

kmax = 200

# simulate the system
for k in range(kmax):
    x_k = dot(A,x_k) + dot(dot(F,V), np.random.randn(2,1)) #x(k+1) = Ax(k) + Fv(k)
    y_k = dot(C,x_k) + dot(dot(G,W), np.random.randn(2,1)) #y(k) = Cx(k) + Fw(k)
    x.append(x_k)
    y.append(y_k)

#get states
x1 = [x[0] for x in x]
x2 = [x[2] for x in x]
#get measurements
y1 = [y[0] for y in y]
y2 = [y[1] for y in y]

#plot
# plot1 = plt.figure()
# plt.plot(x1, x2)
# plt.xlabel("x position")
# plt.ylabel("y piosition")
# plt.title("Traj. of Discrete-time NCV Model")
# plt.grid(True)
# plt.show()

# plot2 = plt.figure()
# plt.plot(y1, y2)
# plt.xlabel("x position")
# plt.ylabel("y piosition")
# plt.title("Traj. of Discrete-time NCV Model")
# plt.grid(True)
# plt.show()


###############################################################################
###############################################################################
# setup Kalman filter
# initial state esimtate --> setup arbitrarily 
x_hat = np.array([[0.5],
                  [20], 
                  [7], 
                  [4]])

P = 10*np.eye(4) #covariance matrix

xhat = []

# V  = 0.5
# W = 0.2

Am = np.array([[1, T, 0, 0],
              [0, 0.87, 0, 0],
              [0, 0, 1, T],
              [0, 0, 0, 0.87]])


V  = 1
W = 1

# Intermittent measurements
y_intermittent = list(y)
y_intermittent[50:80] = [None]*30
y_intermittent[140:170] = [None]*30

for k in range(kmax):
    
    if y_intermittent[k] is None:
        
        # K = dot(dot(Am,dot(P,C.T)), inv(dot(C, dot(P,C.T)) + dot(G, dot(W, G.T)))) #K = APC' / (CPC' + GWG')
        
        x_hat = dot(Am, x_hat)  #x(k+1) = Ax(k)
        
        P = dot(Am, dot(P, Am.T)) + dot(F, dot(V, F.T)) # P = APA' - K*CPA' + FVF'
        
        xhat.append(x_hat)
    
    else:
        
        K = dot(dot(Am,dot(P,C.T)), inv(dot(C, dot(P,C.T)) + dot(G, dot(W, G.T)))) #K = APC' / (CPC' + GWG')
        
        x_hat = dot(Am, x_hat) + dot(K, (y_intermittent[k] - dot(C, x_hat))) #x(k+1) = Ax(k) + K*(y - C*x_hat)
        
        P = dot(Am, dot(P, Am.T)) - dot(K, dot(C, dot(P,Am.T))) + dot(F, dot(V, F.T)) # P = APA' - K*CPA' + FVF'
        
        xhat.append(x_hat)
        
        
    
    
    # K = dot(dot(Am,dot(P,C.T)), inv(dot(C, dot(P,C.T)) + dot(G, dot(W, G.T)))) #K = APC' / (CPC' + GWG')
    
    # x_hat = dot(Am, x_hat) + dot(K, (y[k] - dot(C, x_hat))) #x(k+1) = Ax(k) + K*(y - C*x_hat)
    
    # P = dot(Am, dot(P, Am.T)) - dot(K, dot(C, dot(P,Am.T))) + dot(F, dot(V, F.T)) # P = APA' - K*CPA' + FVF'
    
    # xhat.append(x_hat)
    
    
#get state estimates
xhat1 = [xhat[0] for xhat in xhat]
xhat2 = [xhat[2] for xhat in xhat]


#plot
plot3 = plt.figure()
plt.plot(x1, x2, linewidth=3, label = 'raw traj')
plt.plot(xhat1, xhat2, linewidth=1.5, label = 'estimated traj')

plt.legend()
plt.xlabel("xhat1")
plt.ylabel("xhat2")
plt.title("Traj. of Discrete-time NCV Model")
plt.grid(True)
plt.show()

