using Pkg
Pkg.activate(".")
Pkg.instantiate()

using IJulia
using PyPlot # Seems to install matplotlib only when being used explicitly for the first time

notebook(; dir=pwd())
