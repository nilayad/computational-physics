# task 2 
# a:- squares cubes and sums

min = 1
max = 150
sum_s= 0
sum_c = 0
i = min
print ('squares') 
while (i**2 < 150) :
    print(i**2)
    sum_s +=i**2
    i+=1
print ('sum of all squares is =', sum_s )    

print ('cubes') 
i = min
while (i**3 < 150) :
    print(i**3)
    sum_c +=i**3
    i+=1
print ('sum of all cubes is =', sum_c ) 
#___________________________________________________#
# b))identity proof

def sum_i3(i):
    sum_c = 0
    a=0
    while(a<=i):
        sum_c += a**3
        a +=1
    return sum_c
    
def sum_i(i):
    sum_n = 0
    a=0
    while(a<=i):
        sum_n +=a
        a+=1
        
    return sum_n

s_i= 0
s_i3 = 0
for i in range(200):
    s_i = sum_i(i)
    
    s_i3 = sum_i3(i)
    if (s_i**2 != s_i3):
        print ('Error')
    if(i%20 ==0):
        print ('for', i, 'LHS =',s_i**2,'RHS=', s_i3 )
    
         
    
print( 'identity is proved for upto 200' )
    #print(s_i,s_i3)
    
    
# sol proper

lower = 1
upper= 20

snum = 0
sqnum=0
sknum = 0

for num in range(lower, upper+1):
    snum += num
    sqnum += num**2
    sknum+= num**3
    print(num,num**2,num**3,snum,snum**2,sknum,sqnum)