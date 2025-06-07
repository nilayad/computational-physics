##
#
#


f = open ("numbers.dat")
z = f.readlines()
nums = []
for i in z:
    nums.append(float(i))
print (len(z))
print(min(z), max(z))
#
#alternate
minv =nums[0]
maxv = nums[0]
for x in nums:
    if x<minv:
        minv = x
    if x>maxv:
        maxv = x
print(minv,maxv)
#
a = open("numbers.dat").readlines()
b = [float(t) for t in a]

#
# clossest number

n = float(input('enter the decimal number'))
c= [ abs(n-x) for x in b]
#print(c)
l=c.index(min(c))
print('the number in numbers.dat clossest to', n, 'is', b[l] , 'at index no.', l, 'in sorted list b' )
