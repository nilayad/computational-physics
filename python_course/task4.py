# task4

sum_a = 0.
sum_tot = 1.0
for i in range(10):
  sum_a += 0.1;
if ( sum_a != sum_tot ):
    print (" Python kann nicht rechnen: 1 - 1 = %g" %  (sum_a - sum_tot))
#

# error corrrection

    eps = 1.
    a = True
    while (a):
        eps /= 2. 
        onePlusEps = 1.0 + eps
        print(onePlusEps)
        if ( onePlusEps == 1.0 ) :
            a = False