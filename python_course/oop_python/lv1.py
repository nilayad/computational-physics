import math
class lorvec (object):

    def __init__ (self, E , x=0,y=0,z=0 ):
        self.e = E
        self.x = x
        self.y = y
        self.z = z
        

    def SkalProd( self, tv ):
        return( self.x*tv.x + self.y*tv.y + self.z*tv.z)
    def addlv(self,lv):
        return lorvec( self.e+ lv.e, self.x+lv.x, self.y+lv.y, self.z+lv.z)
    
    def Length( self ):
        return math.sqrt( self.x**2 + self.y**2 + self.z**2 )
    def Angle(self, tv):
        return (math.acos(self.SkalProd(tv)/(self.Length()*tv.Length())))*(180/(math.pi))

        
    def mass(self):
        return math.sqrt(((self.e)**2)-(self.Length())**2)
    
    def minv(self,lv):
        #return math.sqrt(((self.e + lv.e)**2)-((self.addlv(lv).Length())**2))
        return math.sqrt(((self.e + lv.e)**2)-((self.x+lv.x)**2 +(self.y+lv.y)**2 +(self.z+lv.z)**2 ))
    def __str__( self ) :
        return "( %6.4f, %6.3f, %6.3f, %6.3f )" % (self.e, self.x, self.y, self.z)
    
    
#
#
if __name__ == "__main__" :
a = lorvec(45.0002,0.,0.,45.0)
b = lorvec(45.0002,31.8198,0.,31.8198)

print ("Angle = " , a.Angle(b))
print ("Mass a = " , a.mass())
print ("Mass b = " , b.mass())
print ("Mass a+b = " , b.minv(a))
print("mass c = a+b", b.addlv(a))
c = b.addlv(a)
print("Mass c = " , c.mass())