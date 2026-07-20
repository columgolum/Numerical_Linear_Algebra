import os
from juliacall import Main as jl

base_dir = os.path.dirname(os.path.abspath(__file__))

# Define paths to your engine and decomposition files
engine_path = os.path.join(base_dir, "matrix_engine")
decomp_path = os.path.join(base_dir, "decompositions")

# Load all files into the same global Julia memory context
jl.include(os.path.join(base_dir, "vector_engine", "backend2.jl"))
jl.include(os.path.join(engine_path, "backend.jl"))
jl.include(os.path.join(engine_path, "RDIS.jl"))
jl.include(os.path.join(decomp_path, "LU_QR.jl"))

print("Julia computation pipelines successfully inter-linked!")
