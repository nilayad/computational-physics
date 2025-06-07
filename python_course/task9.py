# task 9
#

import matplotlib.pyplot as plt
def dyn(t):
    
    y = h0 - 0.5*g*(t**2)
    v = v0 - g*t
    return y,h
h0 = float(input('starting height h0'))
v0 = float(input('starting velocity v0'))
g =9.81
h =[]
vel =[]
for t in range(30):
    yi,vi = dyn(t)
    h.append(yi)
    vel.append(vi)
#print(h,vel)
plt.plot(h)
plt.scatter(h,range(30))