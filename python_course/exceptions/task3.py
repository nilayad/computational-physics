a = 10000
while b < 100000:
    b = a*9
    
    if (len(set(str(a)))==5):
        if (len(set(str(b)))==5):
            print(a,b)
            
    a+=1