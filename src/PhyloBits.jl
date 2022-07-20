module PhyloBits
__precompile__(true)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/
export hello_PhyloBits, add_one_PhyloBits

print("\nPhyloBits: loading PhyloBits.jl.\n")

using Distributed 	# for workers, spawnat :any, etc.
using Hwloc					# for Hwloc.num_physical_cores(), Hwloc.num_virtual_cores()
using SpecialFunctions
using StatsBase
using StaticArrays	# for SVector
using DataFrames
using Printf


print("For multithreading purposes, Threads.nthreads() = ")
print(Threads.nthreads())
print("\n")
print("For multiprocessor purposes\n:")
print("Hwloc.num_physical_cores() = ")
print(Hwloc.num_physical_cores())
print("\n")
print("Hwloc.num_virtual_cores() = ")
print(Hwloc.num_virtual_cores())
print("\n")
print("Number of workers turned on in this session? Distributed.nworkers() = ")
print(Distributed.nworkers())
print("\n")
print("List of workers turned on in this session? Distributed.workers() = ")
print(Distributed.workers())





# List each PhyloBits code file here
# NOTE: LOAD THE DEPENDENCY .jl FILES *FIRST*, or you get "not recognized" errors

include("PNtypes.jl")					# Types: Tree, Node, Edge etc.
include("PNmanipulateNet.jl")	# 
include("PNauxiliary.jl")			# Helper functions
include("PNreadwrite.jl")			# Reading and writing trees
include("PNdescriptive.jl")		# show() commands for trees etc.
include("TrUtils.jl")					# basic utility functions, R-like functions
include("TreeTable.jl")				# for prt() tree tables (tree dataframes, trdfs), bd_liks(), etc.





"""
# Local nstallation:
using Pkg
Pkg.add(PackageSpec(path="/GitHub/PhyloBits.jl"))

# GitHub Installation
using Pkg
Pkg.add(Pkg.PackageSpec(url="https://github.com/nmatzke/PhyloBits.jl"))

# Add to another package:
cd("/GitHub/BioGeoJulia.jl")
]
activate .
add https://github.com/BioJulia/BioJuliaRegistry.git

# Loading:
Pkg.instantiate()
using PhyloBits

using PhyloBits.TrUtils			# basic utility functions 
using PhyloBits.TreeTable			# for prt() tree tables (DFs), bd_liks(), etc.
"""


#######################################################
# Put functions here
#######################################################
"""
    hello(who::String)

Return "PhyloBits says, hi `who`".
"""
hello_PhyloBits(who::String) = "hello_PhyloBits() says '$who'"

"""
    add_one_PhyloBits(x::Number)

Return `x + 1`.
"""
add_one_PhyloBits(x::Number) = x + 1

end # Ends the module command
