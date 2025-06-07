# print Fahrenheit->Celsius table 
lower =  0 
upper =  300 
step =  20

fahr = lower
while  ( fahr < upper ):  
    celsius =  5./ 9  *  ( fahr -  32 )  # conversion 
    # print(f"Fahrenheit {fahr:6.1f} = Celsius {celsius:8.3f}") # formatted output 
    print ( "Fahrenheit " , fahr , " = Celsius " , celsius ) 
    fahr += step   # increment

# end of loop


# print Fahrenheit<--Celsius table 
lower =  0 
upper =  100 
step =  5

#fahr = lower
celsius = lower
while  ( celsius < upper ):  
    fahr =  (9/5  * celsius)+32  # conversion 
    # print(f"= Celsius {celsius:8.3f}=Fahrenheit {fahr:6.1f} ") # formatted output 
    print (  " Celsius " , celsius ,"= Fahrenheit " , fahr  ) 
    celsius += step   # increment
# end of loop