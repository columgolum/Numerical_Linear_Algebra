from juliacall import Main as jl

class vec():
    def __init__(self,*comp) : # comp stands for components
        self.__comp = jl.convert(jl.Vector[jl.Number],list(comp))
    def len(self): 
        return jl.size(self.__comp)[0]
    def get_vec(self):
        return self.__comp
    def __repr__(self):
        # Leveraging Julia's own function 
        return jl.sprint(jl.show, jl.MIME("text/plain"), self.__comp)
    def __add__(p1,p2) :
        return vec(*jl.add(p1.get_vec(),p2.get_vec()))
    def __sub__(p1,p2) :  # Gives p1-p2
        return vec(*jl.sub(p1.get_vec(),p2.get_vec())) 
    def scale(v,a = 1) :
        return vec(*jl.scale(v.get_vec(),a))
    def mag(v) :
       return jl.mag(v.get_vec())
    
    def dot(u,v) :
        return jl.dot(u.get_vec(),v.get_vec())
    def cosine(u,v):
        w=  vec.dot(u,v)/(u.mag()*v.mag())
        return   w
    def normalise(self) :
        w = self.mag()
        z = self.scale(1/w)
        return z
