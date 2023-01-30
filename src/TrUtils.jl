module TrUtils
__precompile__(true)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/

print("PhyloBits: loading TrUtils.jl dependencies...")

using CSV						# for CSV.read(file, DataFrame; delim="\t")
using DataFrames
#using Plots  						# for Plots.savefig
using StatsBase					# for countmap(a; alg=:auto)
using Distributed 	# for workers, spawnat :any, etc.
using Hwloc					# for Hwloc.num_physical_cores(), Hwloc.num_virtual_cores()

#using RCall

print("...done.\n")

export hello_world_TrUtils, offdiag, make_diag_TF, make_offdiag_TF, convert_df_datatypes!, nthreads_procs, getwd, Rgetwd, setwd, getfn, readtable, recursive_find, include_jls, source, get_a_most_common_value, indexed_Dict_to_DF, convert_is_js_to_single_index, pair_of_indices_to_single_index_column_first, dim, Rdim, seq, Rchoose, Rcbind, Rrbind, Rpaste, Rpaste0, paste, paste0, type, class, Rclass, odds, evens, slashslash, addslash, df_to_Rdata, Reval, Rdput, julian_dput, Rnames, rnames, rn, Rtypes, rtypes, compare_dfs, get_max_df_diffs_byCol, vector_of_vectors_to_df, vvdf, vfft, ont, Rnrow, Rncol, Rsize, Rorder, headLR, flat2, rowSums, single_element_array_to_scalar, headf, moref, get_alphabets, LETTERS, letters, GREEKLETTERS, greekletters, greekletters2, scr2str, lagrange_to_tip

# cutting as it requires the loading of Plots (slow)
# saveopen, 

# Basic checks during Julia startup
function hello_world_TrUtils()
	display("PhyloBits.TrUtils.hello_world_TrUtils() says hello on 2022-06-13_abc")
end


"""
Get the off-diagonal elements of a matrix. 
i.e., the opposite of diag()

Source of function:
https://discourse.julialang.org/t/off-diagonal-elements-of-matrix/41169/7
"""
function offdiag(A::AbstractMatrix)
    [A[ι] for ι in CartesianIndices(A) if ι[1] ≠ ι[2]]
end 


"""
# Create a BitArray (a Matrix of Booleans) with
# true on the diagonal
"""
function make_diag_TF(n::Int64)
	diagTF = BitArray(undef, n, n)
	for i in 1:n
		diagTF[i,i] = true
	end
	return(diagTF)
end

"""
# Create a BitArray (a Matrix of Booleans) with
# true on the diagonal
"""
function make_diag_TF(Qmat::Matrix)
	if (length(dim(Qmat)) != 2)
		errtxt = "STOP ERROR in make_diag_TF(): the input must have 2 dimensions. Instead, dim() gives:"
		print("\n\n")
		print(errtxt)
		print(dim(Qmat))
		print("\n\n")
		error(errtxt)
	end
	if (dim(Qmat)[1] != dim(Qmat)[2])
		errtxt = "STOP ERROR in make_diag_TF(): the Matrix must be square. Instead, dim() gives:"
		print("\n\n")
		print(errtxt)
		print(dim(Qmat))
		print("\n\n")
		error(errtxt)
	end

	n = dim(Qmat)[1]
	diagTF = BitArray(undef, n, n)
	for i in 1:n
		diagTF[i,i] = true
	end
	return(diagTF)
end

"""
# Create a BitArray (a Matrix of Booleans) with
# true on the off-diagonal
"""
function make_offdiag_TF(Qmat::Matrix)
	if (length(dim(Qmat)) != 2)
		errtxt = "STOP ERROR in make_offdiag_TF(): the input must have 2 dimensions. Instead, dim() gives:"
		print("\n\n")
		print(errtxt)
		print(dim(Qmat))
		print("\n\n")
		error(errtxt)
	end
	if (dim(Qmat)[1] != dim(Qmat)[2])
		errtxt = "STOP ERROR in make_offdiag_TF(): the Matrix must be square. Instead, dim() gives:"
		print("\n\n")
		print(errtxt)
		print(dim(Qmat))
		print("\n\n")
		error(errtxt)
	end

	n = dim(Qmat)[1]
	offdiagTF = BitArray(undef, n, n)
	offdiagTF .= true
	for i in 1:n
		offdiagTF[i,i] = false
	end
	return(offdiagTF)
end


"""
# Create a BitArray (a Matrix of Booleans) with
# true on the off-diagonal
"""
function make_offdiag_TF(n::Int64)
	offdiagTF = BitArray(undef, n, n)
	offdiagTF .= true
	for i in 1:n
		offdiagTF[i,i] = false
	end
	return(offdiagTF)
end





"""
dftxt = "Any[\"y\" 2 2 2 1 1.0 1.0 0.2222222222222222 0.2222222222222222 0.0; \"y\" 3 3 3 1 1.0 1.0 0.2222222222222222 0.2222222222222222 0.0; \"v\" 4 3 2 2 2.0 0.16666666666666666 0.037037037037037035 0.037037037037037035 0.0; \"s\" 4 4 2 2 2.0 0.16666666666666666 0.037037037037037035 0.037037037037037035 0.0; \"s\" 4 4 3 2 2.0 0.16666666666666666 0.037037037037037035 0.037037037037037035 0.0; \"v\" 4 2 3 2 2.0 0.16666666666666666 0.037037037037037035 0.037037037037037035 0.0; \"s\" 4 2 4 2 2.0 0.16666666666666666 0.037037037037037035 0.037037037037037035 0.0; \"s\" 4 3 4 2 2.0 0.16666666666666666 0.037037037037037035 0.037037037037037035 0.0]"

dfnames_txt = "[\"event\", \"i\", \"j\", \"k\", \"pair\", \"wt\", \"prob\", \"rate\", \"val\", \"rates_t\"]"
dfnames = Reval(dfnames_txt)

datatypes = Reval("DataType[String, Int64, Int64, Int64, Int64, Float64, Float64, Float64, Float64, Float64]")

ctable = DataFrame(Reval(dftxt), dfnames)
convert_df_datatypes!(ctable, datatypes)
ctable

"""
function convert_df_datatypes!(df, datatypes)
	for i in 1:Rncol(df)
		df[!,i] = convert.(datatypes[i], df[:, i])
	end
end # END function convert_df_datatypes(df, datatypes)




# R-like utilities, and other short functions

# Handy aliases
#function next()
#	continue
#end


# Reload PhyBEARS
reload_cmds = """
function re()
	# Remove and re-install
	Pkg.rm("PhyBEARS")
	Pkg.add(PackageSpec(path="/GitHub/PhyBEARS.jl"))
	@eval using PhyBEARS
	@eval using PhyBEARS.MaxentInterp
	@eval using PhyBEARS.BGExample
	#@eval using PhyBEARS.TrUtils
	#@eval using PhyBEARS.TreeTable
	@eval using PhyBEARS.StateSpace
	@eval using PhyBEARS.Parsers
	@eval using PhyBEARS.SSEs
	@eval using PhyBEARS.TreePass

	# Refresh Revise's look
	atreplinit() do repl
	try
	@eval using Revise
	@async Revise.wait_steal_repl_backend()
	catch
	end
	end
end
"""


# Number of threads and processors

function nthreads_procs()
	print("For multithreading purposes, Threads.nthreads() = ")
	print(Threads.nthreads())
	print("\n")
	print("For multiprocessor purposes:\n")
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
	print("\n")
end # END function nthreads_procs()




# getwd
function getwd()
	pwd()
end

function Rgetwd()
	pwd()
end

# setwd
function setwd(path=expanduser("~"))
	cd(path)
end

"""
# Extract the filename from a path
fn = "/GitHub/PhyBEARS.jl/notes/res_table_df_Ci_eq_i_GflowMap_v9.txt"
getfn(fn)
"""
function getfn(fn)
	words = split(fn, "/")
	lastword = words[length(words)]
	return lastword
end



"""
# Read a tab-delimited text file
fn = "/GitHub/PhyBEARS.jl/notes/area_of_areas_NZ_Oz_v1.txt"
using DataFrames
using CSV
df = readtable(fn)
"""
function readtable(fn; delim="\t", header=1)
	df = CSV.read(fn, DataFrame; delim=delim, header=header)
	return df
end



# Find all code *.jl files in a package

"""
package_path = "/GitHub/PhyBEARS.jl"
recursive_find(package_path)
include_jls(package_path)
"""

function recursive_find(package_path)
	# Look in the source directory, "src"
	srcpath = joinpath(package_path, "src")
	srcfiles = readdir(srcpath)
	jl_files = []  # empty array
	for fn in srcfiles
		TF = endswith(fn, ".jl")
		if (TF == true)
			tmpfn = joinpath(srcpath, fn)
			push!(jl_files, tmpfn)
		end
	end
	
	# Look in the tests directory, "test"
# 	srcpath = joinpath(package_path, "test")
# 	srcfiles = readdir(srcpath)
# 	for fn in srcfiles
# 		TF = endswith(fn, ".jl")
# 		if (TF == true)
# 			tmpfn = joinpath(srcpath, fn)
# 			push!(jl_files, tmpfn)
# 		end
# 	end
	
	return jl_files
end

"""
package_path = "/GitHub/PhyBEARS.jl"
recursive_find(package_path)
include_jls(package_path)
"""
function include_jls(package_path)
	srcfiles = recursive_find(package_path)
	for fn in srcfiles
		include(fn)
	end
end

# source
function source(str)
	include(str)
end


#######################################################
# Gets (a) most common value in a vector
# Note: if there is a tie, the first one in the
# list is returned.
#######################################################
function get_a_most_common_value(vals)
	vals_count_dict = StatsBase.countmap(vals; alg = :auto)
	nonzeroTF = vals_count_dict.slots .!= 0
	maxTF = vals_count_dict.vals[nonzeroTF] .== maximum(vals_count_dict.vals[nonzeroTF])
	most_common_val = vals_count_dict.keys[nonzeroTF][maxTF][1]
	return most_common_val
end


#######################################################
# Take a dictionary of results (e.g., timing results, tres)
# that are indexed by integers. Convert to a DataFrame.
#
# This function assumes the values are themselves 
# DataFrame rows(!).
#######################################################
function indexed_Dict_to_DF(tres)
	tdf = tres[1]
	for i in 2:tres.count
		tdf = vcat(tdf, tres[i])
	end
	return tdf
end # END function indexed_Dict_to_DF(tres)




# Just convert is, js to a single-number index for a matrix
"""
include("/GitHub/PhyBEARS.jl/notes/convert_IsJs.jl")

A=[1 2 3; 4 5 6; 7 8 9]
# See how the single-number indices count down the columns
A[1:9]
n = dim(A)[1]
is= [[1, 1, 1, 1, 1, 1, 1, 1],
 [2, 2, 2, 2, 2, 2, 2, 2],
 [3, 3, 3, 3, 3, 3, 3, 3]]

js= [[1, 1, 1, 1, 1, 1, 1, 1],
 [2, 2, 2, 2, 2, 2, 2, 2],
 [3, 3, 3, 3, 3, 3, 3, 3]]

indices = convert_is_js_to_single_index(is, js, n)
"""

function convert_is_js_to_single_index(is, js, n)
	if length(is) != length(js)
		txt = "STOP ERROR in convert_is_js_to_single_index(is, js, n). The lengths of 'is' and 'js' must be identical. Check these inputs and re-run."
		throw(txt)
	end
	
	
	tmpvec = Any[]
	indices = collect(repeat([tmpvec], length(is)))
	for tmpi in 1:length(indices)
		# Check for Int64[], i.e. a null range with no possibility
		# of being ancestral in a Qmat or Cmat/Carray
		if length(is[tmpi]) == 0
			indices[tmpi] = Int64[]
			continue  # equivalent of next()
		end
		
		#print(tmpi)
		if maximum(is[tmpi])[1] > n
			txt = "STOP ERROR in convert_is_js_to_single_index(is, js, n). The maximum of 'is' cannot be larger than 'n'.  Check these inputs and re-run."
			throw(txt)
		end
	
		if maximum(js[tmpi])[1] > n
			txt = "STOP ERROR in convert_is_js_to_single_index(is, js, n). The maximum of 'js' cannot be larger than 'n'.  Check these inputs and re-run."
			throw(txt)
		end

		#indices[i] .= ((js[i].-1) * n) + is[i]
		single_num_indices = is[tmpi] .+ (js[tmpi] .-1) .* n
		indices[tmpi] = single_num_indices
		# Adding the [1]'s mean that this error is avoided:
		# ERROR: MethodError: no method matching -(::Vector{Int64}, ::Int64)
		# ...which comes from mistmatched types
	end # END for i in 1:length(indices)
	
	return indices
end # END function convert_is_js_to_single_index(is, js, n)


"""
# Convert a list of tuples (representing pairs of i,j indices in a matrix)
# into a list of single-number indices into the matrix. Note that Julia 
# counts down the columns.
#
# This function works fine in the example, but when broadcasting it crashes with:
# # ERROR: MethodError: no method matching -(::Vector{Int64}, ::Int64)
#	# For element-wise subtraction, use broadcasting with dot syntax: array .- scalar
#	# Closest candidates are:
#
# See convert_is_js_to_single_index(is, js, n) for a better option.
# 
A=[1 2 3; 4 5 6; 7 8 9]
# See how the single-number indices count down the columns
A[1:9]
n = dim(A)[1]
is = [2,2,2]
js = [1,2,3]
ij_tuples = tuple.(is,js)
pair_of_indices_to_single_index_column_first.(ij_tuples, n)
"""
function pair_of_indices_to_single_index_column_first(ij_tuple, n)
	ind = 0
	#ind = ((ij_tuple[2]-1) * n) + ij_tuple[1]   # This line attracts: 
	# ERROR: MethodError: no method matching -(::Vector{Int64}, ::Int64)
	# For element-wise subtraction, use broadcasting with dot syntax: array .- scalar
	# Closest candidates are:

# This should work:
	ind = ((ij_tuple[2][1]-1) * n) + ij_tuple[1][1] 


	return ind
end # END function pair_of_indices_to_single_index_column_first(ij_tuple, n)
	



# dimensions
"""
A = ones(3,3)
dim(A)
Rdim(A)
dim(A)[1]
dim(A)[2]
"""
function dim(A)
	size(A)
end

# dimensions
function Rdim(A)
	size(A)
end


# seq
function seq(from, to, by=1)
	return(collect(from:by:to))
end


# choose
# n choose k
function Rchoose(n,k)
	return(binomial(n,k))
end


# cbind()
function Rcbind(A...)
	hcat(A...)
end

# rbind
# c()
# concatenate
function Rrbind(A...)
	vcat(A...)
end

# paste
function Rpaste(array_of_strings; delim)
	newtxt = join(array_of_strings, delim)
	return(newtxt)
end

# paste0
function Rpaste0(array_of_strings; delim="")
	newtxt = join(array_of_strings, delim)
	return(newtxt)
end


# paste
function paste(array_of_strings; delim)
	newtxt = join(array_of_strings, delim)
	return(newtxt)
end

# paste0
function paste0(array_of_strings; delim="")
	newtxt = join(array_of_strings, delim)
	return(newtxt)
end

"""
# type
# Shortcut for julia's typeof()

obj = construct_Res()
Rnames(obj)
Rtypes(obj)
Rcbind(Rnames(obj), Rtypes(obj))
ont(obj)

type(obj)
"""
function type(obj)
	typeof(obj)
end

# class
# Returns a plain-test version of the type/class
function class(obj)
	string(typeof(obj))
end

# Rclass
# Returns a plain-test version of the type/class
function Rclass(obj)
	string(typeof(obj))
end

"""
# Returns the odd numbers
x = 1:5
odds(x)
evens(x)
"""
function odds(x=1:5)
	filter(x->isodd(x), x)
end

"""
# Returns the even numbers
x = 1:5
odds(x)
evens(x)
"""
function evens(x=1:5)
	filter(x->iseven(x), x)
end



# Convert any multiple slashes to single slashes
function slashslash(txt)
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	return(txt)
end

# Add a slash to the end of a string, if it is not there
# (Unless the txt string is "", then return "")
function addslash(txt)
	if (txt == "")
		return("")
	end
	if (endswith(txt, "/") == false)
		txt = join([txt, "/"])
	end
	return(txt)
end


# Save a julia DataFrame to an R data.frame
# as an Rdata file that can be easily loaded.
# Source: https://stackoverflow.com/questions/28084403/saving-julia-dataframe-to-read-in-r-using-hdf5/57903489#57903489

#######################################################
# NOTE: YOU MAY HAVE TO PASTE THE FUNCTION DEFINITION
# OF df_to_Rdata(), ALONG WITH THE "using" STATEMENTS,
# INTO YOUR MAIN ENVIRONMENT FOR IT TO WORK!!!
#######################################################

#using DataFrames
#using RCall
function df_to_Rdata(df; fn="dfjulia.RData", path=expanduser("~"))
	# Create the output fn
	pathfn = slashslash(join([addslash(path), fn], ""))
		
	# R environment in a session started from Julia
	g = globalEnv  # This requires "using RCall" previously, or you get:
								 # ERROR: UndefVarError: globalEnv not defined
	reval(rparse("dfls <- NULL"))

	# add columns one at a time converting Julia vectors to R-types via   RCall.sexp
	#  https://github.com/JuliaStats/RCall.jl/blob/master/src/sexp.jl
	for cnm in DataFrames._names(df)
		g[:colcnm] = sexp(convert(Array, df[!,cnm]))
		reval(rparse("dfls\$$cnm <- colcnm"))
	end
	reval(rparse("df <- data.frame(dfls)"))
	
	# Make and run the command to save the .Rdata file
	# (the .Rdata file will load to object df in R)
	txt = join(["save(file='", pathfn, "', df)"], "")
	reval(rparse(txt))
	return(pathfn)
end


# Input an object from a string representation from repr()
# dget, load, eval()
# eval
"""
tmpmatrix = [3 1; 3 2; 5 3; 5 4; 7 5; 7 6]
tmpstr = repr(tmpmatrix)
tmpstr2 = eval(Meta.parse(tmpstr))
tmpstr2
"""
function Reval(tmpstr)
	eval(Meta.parse(tmpstr))
end

# Output an object to a string representation with repr()
# dput, dump, str()
"""
tmpstr = "[[1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]"
tmpstr2 = tmpstr
states_list = Reval(tmpstr)
@test Rdput(states_list) == tmpstr2
"""
function Rdput(item)
	tmpstr = repr(item)
end

"""
Output a DataFrame to a text representation (to screen) that can be 
read back in.
Source: https://stackoverflow.com/questions/61723347/how-to-provide-reproducible-sample-data-in-julia

# Example:
df = DataFrame(AbstractVector[["y", "j", "j", "y"], [1, 1, 2, 2], [1, 1, 2, 2], [1, 2, 1, 2], [1, 2, 2, 1], [0.9166666666666666, 0.5, 0.5, 0.9166666666666666], [0.6470588235294118, 0.35294117647058826, 0.35294117647058826, 0.6470588235294118], [0.14379084952941176, 0.07843137247058823, 0.07843137247058823, 0.14379084952941176], [0.14379084952941176, 0.07843137247058823, 0.07843137247058823, 0.14379084952941176], [0.0, 0.0, 0.0, 0.0]], DataFrames.Index(Dict(:event => 1, :wt => 6, :val => 9, :j => 3, :k => 4, :rates_t => 10, :prob => 7, :i => 2, :rate => 8, :pair => 5), [:event, :i, :j, :k, :pair, :wt, :prob, :rate, :val, :rates_t]); copycols=true)

# Output to screen:
julian_dput(df)

# DataFrame(AbstractVector[["y", "j", "j", "y"], [1, 1, 2, 2], [1, 1, 2, 2], [1, 2, 1, 2], [1, 2, 2, 1], [0.9166666666666666, 0.5, 0.5, 0.9166666666666666], [0.6470588235294118, 0.35294117647058826, 0.35294117647058826, 0.6470588235294118], [0.14379084952941176, 0.07843137247058823, 0.07843137247058823, 0.14379084952941176], [0.14379084952941176, 0.07843137247058823, 0.07843137247058823, 0.14379084952941176], [0.0, 0.0, 0.0, 0.0]], DataFrames.Index(Dict(:event => 1, :wt => 6, :val => 9, :j => 3, :k => 4, :rates_t => 10, :prob => 7, :i => 2, :rate => 8, :pair => 5), [:event, :i, :j, :k, :pair, :wt, :prob, :rate, :val, :rates_t]), nothing, nothing, true)

# Cut ", nothing, nothing, true", and it works!
df2 = DataFrame(AbstractVector[["y", "j", "j", "y"], [1, 1, 2, 2], [1, 1, 2, 2], [1, 2, 1, 2], [1, 2, 2, 1], [0.9166666666666666, 0.5, 0.5, 0.9166666666666666], [0.6470588235294118, 0.35294117647058826, 0.35294117647058826, 0.6470588235294118], [0.14379084952941176, 0.07843137247058823, 0.07843137247058823, 0.14379084952941176], [0.14379084952941176, 0.07843137247058823, 0.07843137247058823, 0.14379084952941176], [0.0, 0.0, 0.0, 0.0]], DataFrames.Index(Dict(:event => 1, :wt => 6, :val => 9, :j => 3, :k => 4, :rates_t => 10, :prob => 7, :i => 2, :rate => 8, :pair => 5), [:event, :i, :j, :k, :pair, :wt, :prob, :rate, :val, :rates_t]))
"""
function julian_dput(df)
	invoke(show, Tuple{typeof(stdout), Any}, stdout, df);
end

# fields / "names()" of an object
# https://stackoverflow.com/questions/41687418/how-to-get-fields-of-a-julia-object
"""
obj = construct_Res()
Rnames(obj)
Rtypes(obj)
Rcbind(Rnames(obj), Rtypes(obj))
ont(obj)
"""
function Rnames(obj)
	flat2(fieldnames(typeof(obj)))
end

# fields / "names" of an object
# (shortcut for Rnames)
# https://stackoverflow.com/questions/41687418/how-to-get-fields-of-a-julia-object
"""
obj = construct_Res()
rnames(obj)
rtypes(obj)
Rcbind(Rnames(obj), rtypes(obj))
ont(obj)
"""
function rnames(obj)
	Rnames(obj)
end


# fields / "names" of an object
# (shortcut for Rnames)
# https://stackoverflow.com/questions/41687418/how-to-get-fields-of-a-julia-object
"""
obj = construct_Res()
rn(obj)
rtypes(obj)
Rcbind(Rnames(obj), rtypes(obj))
ont(obj)
"""
function rn(obj)
	Rnames(obj)
end




"""
obj = construct_Res()
Rnames(obj)
Rtypes(obj)
Rcbind(Rnames(obj), Rtypes(obj))
ont(obj)
"""
function Rtypes(obj)
	tmpnames = Rnames(obj)
	types = collect(repeat([DataType],length(tmpnames))) # empty array
	for i in 1:length(tmpnames)
		types[i] = typeof(getfield(obj, tmpnames[i]))		
	end
	return types
end

"""
obj = construct_Res()
Rnames(obj)
Rtypes(obj)
Rcbind(Rnames(obj), rtypes(obj))
ont(obj)
"""
function rtypes(obj)
	types = Rtypes(obj)
	return types
end



"""
Compare 2 DataFrames, elementwise, for numerical differences


df1 = bgb_ancstates_AT_nodes_df
df2 = vfft(res.anc_estimates_at_each_nodeIndex_branchTop[R_order])
tol = 1e-4
"""
function compare_dfs(df1, df2, tol=1e-4)
	difs = eachcol(df1) .- eachcol(df2)
	TFs = zeros(dim(df1))
	for col in 1:length(difs)
		abs_diffs = abs.(difs[col])
		abs_diffs[isnan.(abs_diffs)] .= 0.0
		TFs[:,col] .= abs_diffs .< tol
	end
	TFs
	
	return(TFs)
end

function get_max_df_diffs_byCol(df1, df2, tol=1e-4)
	max_diffs = repeat([0.0], Rncol(df1))
	difs = eachcol(df1) .- eachcol(df2)
	for col in 1:length(difs)
		abs_diffs = abs.(difs[col])
		abs_diffs[isnan.(abs_diffs)] .= 0.0
		max_diffs[col] = maximum(abs_diffs)
	end
	return(max_diffs)
end




"""
# Convert a vector of vectors (e.g. ancestral states) to DataFrame
# Source: https://stackoverflow.com/questions/72957438/how-to-convert-a-vector-of-vectors-into-a-dataframe-in-julia-without-for-loop
"""
function vector_of_vectors_to_df(vec_of_vecs)
	df = DataFrame( mapreduce(permutedims, vcat, vec_of_vecs), :auto)
	return(df)
end

"""
# Convert a vector of vectors (e.g. ancestral states) to DataFrame
# Source: https://stackoverflow.com/questions/72957438/how-to-convert-a-vector-of-vectors-into-a-dataframe-in-julia-without-for-loop
"""
function vvdf(vec_of_vecs)
	df = vector_of_vectors_to_df(vec_of_vecs)
	return(df)
end




"""
dfft: conditional format table (BioGeoBEARS "cft"), converting
a Julia vector of vectors to a DataFrame, then rounding
"""
function vfft(vv; digits=4)
	tdf = vvdf(vv)
	df = round.(tdf; digits=digits)
	return(df)
end




"""
# Rename the column names in a DataFrame
"""
function rename_df(df, newnames)
	symbols=Array{Symbol,1}(newnames)
	rename!(df, symbols)
	return(df)
end



# ont = object names and types
"""
obj = construct_Res()
Rnames(obj)
Rtypes(obj)
Rcbind(Rnames(obj), Rtypes(obj))
ont(obj)
"""
function ont(obj)
	Rcbind(Rnames(obj), Rtypes(obj))
end



# Send the just-done plot to PDF, and open
#function saveopen(fn)
#	print("Saving ", """'""", fn, """'""")
#	Plots.savefig(fn)
#	cmdtxt = join(["`", "open ", fn, "`"], "")
#	print("""Running """, cmdtxt, "\n")
#	run(`open $fn`)
#end


function Rnrow(obj)
	return size(obj)[1]
end

function Rncol(obj)
	return size(obj)[2]
end

function Rsize(obj)
	return size(obj)
end

function Rorder(obj)
	return sortperm(obj)
end

# Print the left and rightmost columns of a table
function headLR(df, num_startcols=4, num_endcols=4)
	ncols = Rncol(df)
	startcols = collect(1:num_startcols)
	endcols = collect((ncols-(num_endcols-1)):ncols)
	colnums = flat2(collect([startcols, endcols]))
	colnums
	print(df[:,colnums])
end


# Flattens an array of arrays into a vector
# Similar to R's unlist()
function flat2(arr)
    rst = Any[]
    grep(v) = for x in v
        if isa(x, Array) grep(x) else push!(rst, x) end
    end
    grep(arr)
    rst
end

"""
Version of R's rowSums
Designed for 2-dimensional matrices. Technically, it will return 
a vector of the 2nd dimension of the Matrix/Array

Qmat = [0.0 0.0 0.0 0.0;
e_val 0.0 0.0 d_val;
e_val 0.0 0.0 d_val;
0.0 e_val e_val 0.0]


Qmat[make_diag_TF(Qmat)] .= flat2(sum(Qmat; dims=2))

u0 = [0, 0.125, 0.75, 0.125]

"""
function rowSums(Qmat)
	return flat2(sum(Qmat; dims=2))
end





# Convert a single-element array to scalar
# Julia often produces single-element arrays. 
# To convert to scalar, just take item [1]
# https://stackoverflow.com/questions/39079428/1-element-array-to-scalar-in-julia
function single_element_array_to_scalar(tmparray)
	if length(tmparray) != 1
		txt = ["STOP ERROR in single_element_array_to_scalar().\nThe input 'tmparray' has to be of length 1, but length(tmparray)=", string(length(tmparray)), ".\nPrinting input tmparray...\n"]
		errortxt = join(txt, "")
		println(errortxt)
		print(tmparray)
		error(errortxt)
	end
	
	# If check passed, go ahead.
	tmpscalar = tmparray[1]
	return tmpscalar
end



# Print the file to screen, with line numbers
function headf(fn; numlines=5)
	open(fn, "r") do f
		for (i,ln) in enumerate(eachline(f))
			if i > numlines
				break
			end
			println("$i $ln")
		end
	end
end

# Print the file to screen, with line numbers
function moref(fn)
	open(fn, "r") do f
	 for (i,ln) in enumerate(eachline(f))
		 println("$i $ln")
	 end
	end
end




"""
# Gives a tuple containing the letters, then Greek symbols available in Julia
# (taken from: https://en.wikipedia.org/wiki/Greek_alphabet )
#
# Example:

alphabets = get_alphabets()
rn(alphabets)
hcat(alphabets.greek_letters_in_english, alphabets.greek_alphabet, alphabets.greek_alphabet2, alphabets.GREEK_ALPHABET)
"""

function get_alphabets(return_val="all")
	# https://en.wikipedia.org/wiki/Greek_alphabet

	alphabet = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"]

	ALPHABET = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]
	

	greek_alphabet = ["α", "β", "γ", "δ", "ε", "ζ", "η", "θ", "ι", "κ", "λ", "μ", "ν", "ξ", "ο", "π", "ρ", "σ", "τ", "υ", "φ", "χ", "ψ", "ω"]
	greek_alphabet2 = ["α", "β", "γ", "δ", "ε", "ζ", "η", "θ", "ι", "κ", "λ", "μ", "ν", "ξ", "ο", "π", "ρ", "ς", "τ", "υ", "φ", "χ", "ψ", "ω"]
	GREEK_ALPHABET = ["Α", "Β", "Γ", "Δ", "Ε", "Ζ", "Η", "Θ", "Ι", "Κ", "Λ", "Μ", "Ν", "Ξ", "Ο", "Π", "Ρ", "Σ", "Τ", "Υ", "Φ", "Χ", "Ψ", "Ω"]

	additional_greek = ["ς", "Γ", "Δ", "Θ", "Λ", "Ξ", "Π", "Σ", "Φ", "Ψ", "Ω"]

	greek_letters_in_english = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega"]

	corresponding_english_letters = ["a", "b", "g", "d", "e",  "z", "h", "th", "i", "k", "l", "m", "n", "ks", "o", "p",  "r", "s", "t", "u", "v", "x", "y", "o"]

	alphabets = (alphabet=alphabet, ALPHABET=ALPHABET, greek_alphabet=greek_alphabet, greek_alphabet2=greek_alphabet2, GREEK_ALPHABET=GREEK_ALPHABET, additional_greek=additional_greek, greek_letters_in_english=greek_letters_in_english, corresponding_english_letters=corresponding_english_letters)

	
	if (return_val == "all")
		return alphabets
	end

	if (return_val == "greek_df")
		greek_df = DataFrame(greek_alphabet=greek_alphabet, greek_alphabet2=greek_alphabet2, GREEK_ALPHABET=GREEK_ALPHABET, additional_greek=additional_greek, greek_letters_in_english=greek_letters_in_english, corresponding_english_letters=corresponding_english_letters)
		return greek_df
	end
	return alphabets
end

"""
# After the 26 standard letters, the function uses the opposite case, then Greek symbols available in Julia
# (taken from: https://en.wikipedia.org/wiki/Greek_alphabet )
"""
function LETTERS(i=1:26; output_char=false)
	alphabet = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"]

	ALPHABET = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]

	greek_alphabet = ["α", "β", "γ", "δ", "ε", "ζ", "η", "θ", "ι", "κ", "λ", "μ", "ν", "ξ", "ο", "π", "ρ", "σ", "τ", "υ", "φ", "χ", "ψ", "ω"]
	additional_unique_greek = ["ς", "Γ", "Δ", "Θ", "Λ", "Ξ", "Π", "Σ", "Φ", "Ψ", "Ω"]
	alphabet = vcat(ALPHABET, alphabet, greek_alphabet, additional_unique_greek)
	k = 1
	tmp = join.(collect(Iterators.product(ntuple(_ -> alphabet, k)...))[:])
	letters_vec = tmp[i]
	
	if output_char == true
		letters_vec = collect(Iterators.flatten(vec.(collect.(area_names))))
	end
	
	return(letters_vec)
end

function letters(i=1:26; output_char=false)
	alphabet = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"]

	ALPHABET = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]
	greek_alphabet = ["α", "β", "γ", "δ", "ε", "ζ", "η", "θ", "ι", "κ", "λ", "μ", "ν", "ξ", "ο", "π", "ρ", "σ", "τ", "υ", "φ", "χ", "ψ", "ω"]
	additional_unique_greek = ["ς", "Γ", "Δ", "Θ", "Λ", "Ξ", "Π", "Σ", "Φ", "Ψ", "Ω"]
	alphabet = vcat(alphabet, ALPHABET, greek_alphabet, additional_unique_greek)
	k = 1
	tmp = join.(collect(Iterators.product(ntuple(_ -> alphabet, k)...))[:])
	letters_vec = tmp[i]

	if output_char == true
		letters_vec = collect(Iterators.flatten(vec.(collect.(area_names))))
	end

	return(letters_vec)
end



function GREEKLETTERS(i=1:26; output_char=false)
	GREEK_ALPHABET = ["Α", "Β", "Γ", "Δ", "Ε", "Ζ", "Η", "Θ", "Ι", "Κ", "Λ", "Μ", "Ν", "Ξ", "Ο", "Π", "Ρ", "Σ", "Τ", "Υ", "Φ", "Χ", "Ψ", "Ω"]
	alphabet = GREEK_ALPHABET
	k = 1
	tmp = join.(collect(Iterators.product(ntuple(_ -> alphabet, k)...))[:])
	letters_vec = tmp[i]
	
	if output_char == true
		letters_vec = collect(Iterators.flatten(vec.(collect.(area_names))))
	end
	
	return(letters_vec)
end

function greekletters(i=1:26; output_char=false)
	greek_alphabet = ["α", "β", "γ", "δ", "ε", "ζ", "η", "θ", "ι", "κ", "λ", "μ", "ν", "ξ", "ο", "π", "ρ", "σ", "τ", "υ", "φ", "χ", "ψ", "ω"]
	alphabet = greek_alphabet
	k = 1
	tmp = join.(collect(Iterators.product(ntuple(_ -> alphabet, k)...))[:])
	letters_vec = tmp[i]

	if output_char == true
		letters_vec = collect(Iterators.flatten(vec.(collect.(area_names))))
	end

	return(letters_vec)
end


function greekletters2(i=1:26; output_char=false)
	greek_alphabet2 = ["α", "β", "γ", "δ", "ε", "ζ", "η", "θ", "ι", "κ", "λ", "μ", "ν", "ξ", "ο", "π", "ρ", "ς", "τ", "υ", "φ", "χ", "ψ", "ω"]
	alphabet = greek_alphabet2
	k = 1
	tmp = join.(collect(Iterators.product(ntuple(_ -> alphabet, k)...))[:])
	letters_vec = tmp[i]

	if output_char == true
		letters_vec = collect(Iterators.flatten(vec.(collect.(area_names))))
	end

	return(letters_vec)
end











function scr2str(obj)
	io = IOBuffer()
	show(io, "text/plain", obj)
	str = String(take!(io))
	return str
end

function lagrange_to_tip(inputs, geog_df)
	dfnames = names(geog_df)
	area_column_nums = 2:length(dfnames)
	areas_txt_list = dfnames[area_column_nums]
	numareas = length(inputs.setup.areas_list)
	areas_list = collect(1:numareas)
	maxareas = numareas
	include_null_range = true
	states_list = areas_list_to_states_list(areas_list, maxareas, include_null_range)
	statenums = collect(1:length(states_list))
	observed_statenums = collect(repeat([0], nrow(geog_df)))
	trdf_nodenums = collect(1:nrow(trdf))

	for i in 1:nrow(geog_df)
		tmprow = geog_df[i,area_column_nums]
		tmpnums = parse.(Int, flat2(tmprow))
		range_as_areanums = inputs.setup.areas_list[tmpnums .== 1]
		# Compare observed range_as_areanums to full states_list
		TF = [range_as_areanums] .== inputs.setup.states_list
		if (sum(TF) != 1)
			txt = paste0(["STOP ERROR: An observed range in your geography file, from tipname '", geog_df[i,:tipnames], "', is not found in the list of states in inputs.setup.states_list. Printing range_as_areanums, then inputs.setup.states_list"])
			print(txt)
			print("\nrange_as_areanums (this is the observed range that was not found):\n")
			print(range_as_areanums)
			print("\nstates_list:\n")
			print(inputs.setup.states_list)
			error(txt)
		end
	
		# Yay, a single match to the states_list was found!
		# Convert to a state index number
		observed_statenums[i] = statenums[TF][1]
	end

	# Go through the geography file tips, match each one to the correct node of trdf,
	# then update the tip likelihoods at that node.

	for i in 1:nrow(geog_df)
		spname = geog_df[i,:tipnames]
		TF = spname .== inputs.trdf[!,:nodeName]
		nodeNum = trdf_nodenums[TF][1]
	
		# Input likelihoods of 1 for the observed state, 0 otherwise
		inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum] .= 0.0         # zero out
		inputs.res.normlikes_at_each_nodeIndex_branchTop[nodeNum] .= 0.0     # zero out
		inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum][inputs.setup.observed_statenums[i]] = 1.0
		inputs.res.normlikes_at_each_nodeIndex_branchTop[nodeNum][inputs.setup.observed_statenums[i]] = 1.0
		inputs.res.sumLikes_at_node_at_branchTop[nodeNum] = 1.0
	end
end # END function lagrange_to_tip(inputs, geog_df)




end # end of module