import random
r = []
for i in range(100):
    u = [ random.gauss(1,10), random.gauss(-2,-6), random.gauss(2,6)]
    r.append(u)

r.append([1,-4,4])
verror=[]
Aerror = []
roots =[]
for i in r:
    try:
        #print(root(*i))
        t = (r.index(i), root(*i))
        roots.append(t)
    except ValueError as x :
        #print ("disc cant be less than zero::.", x)
        t=(r.index(i))
        verror.append(t)
    except ArithmeticError as x:
        #print("A cannot be zero::", x)
        t=(r.index(i),x)
        Aerror.append(t)
        
print('roots found for :', roots )
print('___________________________________')
print('valueerror for index in n:',verror)
print('___________________________________')
print('arithmatic error for index in n :', Aerror)