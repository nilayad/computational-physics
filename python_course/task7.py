# task 7
#a prime number tester

n = int(input('enter the decimal number'))

c=0
for b in range(2,n):
    for i in range(2, b) :
        if ((b%i) ==0):
           # print (b, 'is not a prime no')
            break
        
    else:
        print (b)#, 'is a prime no.')
        c +=1
print('there are ',c,'no. of prime numbers between 1 and', n)
#______________________________________#
#extension to count no. of primes between 2 and the input number n
'''
n = int(input('enter the decimal number'))
i=1
b = 2
c=0
while b<n:
    
    while i< b :
    
        if ((b%i) ==0):
            #print (n, 'is not a prime no')
            i +=1
            b +=1
            continue
        else:
            print (b, 'is a prime no.')
            c +=1
            b +=1
            break
print('total prime nos between 2 and ',n,' are ',c)'''

#  c:-
#Sieve of the Erasthone

n = int(input('enter the decimal number'))
a = list(range(2,n+1))

for i in a:
          
    for m in a :
        if m>i :
       
            if (m%i ==0):
                a.remove(m)
                
#a.insert(0,2)
print ('the list of prime numbers between 2 and', n,' are', a)

#
#________________________________________#


#alternate and compect way to find  prime nos upto 100::

[ x for x in range(2,101) if all(x % y != 0 for y in range(2, x)) ]