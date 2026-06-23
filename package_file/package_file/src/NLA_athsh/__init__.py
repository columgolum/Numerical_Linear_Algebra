import os
from juliacall import Main as jl

# Automatically track the root directory of the package
base_dir = os.path.dirname(os.path.abspath(__file__))


# Single-Graph Compilation Engine: Loads everything into a unified global memory context
jl.include(os.path.join(base_dir, "vector_engine", "backend2.jl"))
jl.include(os.path.join(base_dir, "matrix_engine", "backend.jl"))
jl.include(os.path.join(base_dir, "matrix_engine", "RDIS.jl"))

print("Julia computation pipelines successfully inter-linked!")
