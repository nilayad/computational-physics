# task 5

#fibonacci number
fib_n = []
fib_n.append(0)
fib_n.append(1)
i = 2
while (i<12):
    fib_n_i = fib_n[i-1] +fib_n[i-2]
    fib_n.append(fib_n_i)
    i +=1
print (fib_n)
    
    # task5b
    
    #def func fibn
    
def fibn(n):
    fib_n = []
    fib_n.append(0)
    fib_n.append(1)
    i = 2
    while (i<=n):
        fib_n_i = fib_n[i-1] +fib_n[i-2]
        fib_n.append(fib_n_i)
        i +=1
    fibno_n = fib_n[n]
    return fibno_n

# sol to task b
    
    
a = []
b = 0
i=1
a.append('NA')

while (i<10) :
    
    b = 0
    fn_1 = fibn(i+1)
    f1_n = fibn(i-1)
    f_n = fibn(i)
    b = (fn_1)*(f1_n)- (f_n**2)
    i +=1
    print (b)
    print ('LHS= ',b, '& RHS =', (-1)**(i+1))
    a.append(b)
print (a)