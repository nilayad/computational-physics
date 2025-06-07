from ThreeVec import ThreeVec
import math
class lvt(ThreeVec):
    
    def __init__( self,t=0., x=0., y=0., z=0.):
        ThreeVec.__init__( self, x, y, z )  # initialize  ThreeVec
        self.t = t
    
    def 