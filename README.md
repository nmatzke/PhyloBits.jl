# PhyloBits
 Phylogenetic tree tables, PhyloNetworks tree reader without dependencies

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaLang.github.io/PhyloBits.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaLang.github.io/PhyloBits.jl/dev)

Linux and macOS: [![Build Status](https://travis-ci.org/JuliaLang/PhyBEARS.jl.svg?branch=master)](https://travis-ci.org/JuliaLang/PhyloBits.jl)

Windows: [![Build Status](https://ci.appveyor.com/api/projects/status/github/JuliaLang/PhyloBits.jl?branch=master&svg=true)](https://ci.appveyor.com/project/tkelman/example-jl/branch/master)

[![Coverage Status](https://coveralls.io/repos/JuliaLang/PhyBEARS.jl/badge.svg?branch=master)](https://coveralls.io/r/JuliaLang/PhyloBits.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaLang/PhyloBits.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaLang/PhyloBits.jl?branch=master)



# Local add instructions:
Pkg.add(PackageSpec(path="/GitHub/PhyloBits.jl"))
using PhyloBits

# Get the list of installed packages:
x = Pkg.installed()
list_of_installed_packages = collect(x);
println.(list_of_installed_packages)	# Messy

# Or:
x["PhyloBits"]
# v"0.0.1"

