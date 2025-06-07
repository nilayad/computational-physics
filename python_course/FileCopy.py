fil1 = str(input('name of the file to be copied'))
fil2 = str(input('name of the new copy file to be created '  ) )
         
fck = open(fil1).readlines()
# print(fck)
copiedtext=fck
fwk = open(fil2 , "w")
for i in range(len(copiedtext)):
    fwk.write(copiedtext[i])
fwk.close()