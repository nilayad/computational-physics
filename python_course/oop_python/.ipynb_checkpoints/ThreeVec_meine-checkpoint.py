import math
import numpy as np

class ThreeVec(object):
    def __init__( self, x=0., y=0., z=0.):
        self.x = x
        self.y = y
        self.z= z
    def add(self, tv):
        newvec = ThreeVec( self.x+tv.x, self.y+tv.y, self.z+tv.z)
        return newvec
    def __sub__(self, tv):
        newvec = ThreeVec( self.x-tv.x, self.y-tv.y, self.z-tv.z)
        return newvec
    
    def Subtract(self, tv):
        subvec = ThreeVec( self.x-tv.x, self.y-tv.y, self.z-tv.z)
        return subvec
    def Length( self ): # 
        return math.sqrt( self.x**2 + self.y**2 + self.z**2 )
    def Scale( self, a ): # Skalierung mit float
        scalevec = ThreeVec( (self.x)*a, (self.y)*a, (self.z)*a)
        return scalevec
    

    def Angle( self, tv ): 
        angle = np.arccos((self.ScalProd(tv))/((self.Length())*(tv.Length())))
        return np.degrees(angle)
    

    def ScalProd( self, tv ): 
        dotval = np.dot([self.x, self.y, self.z], [tv.x,tv.y,tv.z])
        return dotval
 
    def CrossProd( self, tv ): 
        cpvec = np.cross([self.x, self.y, self.z], [tv.x, tv.y, tv.z])
        return ThreeVec(*cpvec)