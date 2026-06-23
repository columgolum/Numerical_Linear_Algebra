from juliacall import Main as jl
# NLA_atsh.vector_engine ki jagah do dots (..) ka use karein parent folder par jaane ke liye
from ..vector_engine import vec


class mat ():
    def __init__(self, *rv): # rv = row vectors
        # If Julia gives exactly one matrix/column vector, then it converts it to a matrix
        if len(rv) == 1 and "juliacall" in str(type(rv[0])):
            self.__rv = jl.convert(jl.Matrix[jl.Number], rv[0])
            
        # if multiple items are pass by julia then it horizontally concatenates them and then transposes the result to get a matrix
        else:
            self.__rv = jl.convert(jl.Matrix[jl.Number], jl.transpose(jl.hcat(*rv)))
        
    def iter(self):
        return self.__rv  
        
    def __repr__(self):
        # Leveraging Julia's own function 
        return jl.sprint(jl.show, jl.MIME("text/plain"), self.__rv)
        
    def dim(self) :
        return jl.size(self.iter())
        
    def transpose(self):
       return mat(jl.transpose_(self.iter()))
    
    def __mul__(u,v): 
        return mat(jl.matmult(u.iter(),v.iter()))

    ''' in julia matrix multiplication
    first column is build first ( using the 3rd loop { k is changing fastest})
    then using 2nd loo the column is updated again and again ( cuz j changes)
    then i changes and we move to the next column ( then above 2 steps are repeated)'''
    def __add__(u,v):
        return mat((jl.add(u.iter(),v.iter())))
    def __sub__(u,v):
        return mat((jl.sub(u.iter(),v.iter())))
    def trace(u):
        return jl.trace(u.iter())
    def ref(u):
       output,rank,nullity = jl.ref(u.iter())
       return mat(output),rank,nullity
    def rref(u):
        output,rank,nullity = jl.rref(u.iter())
        return mat(output),rank,nullity
    def solveLSE(self):
        result = jl.solveLSE(self.iter())
        if len(result) == 1:
            return "No solution exists"
        elif len(result) == 2:
            
            status, sol = result
            return status, mat(sol)
        else:
            status, part_sol, basis_vecs = result
            return status, mat(part_sol), [mat(v) for v in basis_vecs]

    def inverse(u):
        return mat(jl.inverse(u.iter()))
    def det(u):
        return jl.det(u.iter())
    
