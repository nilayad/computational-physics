import sys
import math
p = list(range(200))


import random
n =[]
for i in range(100):
    
    u =  float(random.gauss(0.005,200))
    n.append(u)
    
    
def fak(n):
    #n = float(input('enter the decimal number'))
    lim = sys.float_info.max
    try:
        factn = math.gamma(n+1)
        return factn
    except OverflowError as x:
        a = x
        return a
        
           
    
if __name__ == "__main__" :

fac=[]
for g in n :
    t =[g,fak(g)]
    fac.append(t)
    
print('for value g the factorial is ', fac)

    