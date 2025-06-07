n = int(input('enter the decimal number'))
binrep = ''
a = n
while n > 0 :
    remainder = n%2
    binrep = str(remainder) + binrep
    n = n//2
print ('the binary rep for ', a, 'is' ,binrep) 