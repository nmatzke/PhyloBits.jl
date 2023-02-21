#######################################################
# TreeTable.jl
# 
# Contains prt(), "print tree table", which writes
# a Julia PhyloNetworks tree object to a DataFrames.jl
# table, much like the prt() function of BioGeoBEARS.
# 
# prt() also includes the R node numbers, so that the
# tree table (usually trtable or trdf) can be reordered
# to match the needs of some function needing an 
# ape-like tree structure.
#
# Various utility functions that feed into prt() are
# also included (tree traversals to trace node 
# ancestors etc.
#
# Also includes bd_liks(), which produces the lnL
# of a tree under a constant-rate birth-death process.
# This is derived from ape's birthdeath() function.
# However, all of the parts of the calculation are 
# saved separately, as various other programs include
# or exclude various constant terms (e.g. diversitree's
# ClaSSE does not include the logfactorial(nb_node) 
# term). 
# 
# bd_liks_trdf does the same, where the input is a 
# trtable/trdf instead of the actual tree.
#######################################################

module TreeTable
__precompile__(true)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/
print("PhyloBits: loading TreeTable.jl dependencies...")

using PhyloBits.PNtypes # for e.g. HybridNetwork
using PhyloBits.TrUtils # for e.g. flat2
using DataFrames			# for e.g. DataFrame()
#using PhyloBits		# for e.g. readTopology()
using SpecialFunctions			# for e.g. logfactorial()

print("...done.\n")

export get_nodenumbers_above_node, get_postorder_nodenumbers_above_node, initialize_edgematrix, get_pruningwise_postorder_edgematrix, get_LR_uppass_edgematrix, get_LR_downpass_edgematrix, get_LR_uppass_nodeIndexes, get_LR_downpass_nodeIndexes, get_Rnodenums, get_nodeIndex_PNnumber, get_nodeIndex_from_PNnumber, get_nonrootnodes, get_nonrootnodes_trdf, nodetimes, get_fossils, get_hooks, prt, get_taxa_descending_from_each_node, isTip_TF, get_NodeIndexes_from_edge, get_NodeIndex_df_by_tree_edges, get_node_heights, get_node_ages, get_root_age, get_tree_height, branching_times, bd_liks, bd_liks_trdf, get_treelength, ML_yule_birthRate, ML_yule_birthRate_wRoot, get_num_speciation_nodes, get_num_tips_from_speciation





# Get the 2 nodeIndexes descending from a node; iterates for uppass
# (assumes a binary, PhyloNetwork, rooted tree)
# Reverse for downpass

"""
using DataFrames
using PhyloBits

# For Nick's editing (ignore)
include("/GitHub/BioGeoJulia.jl/notes/TreeTable.jl")

#######################################################
# Typical bifurcating (binary) tree
#######################################################
great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr

nodeIndex_array = collect(repeat([0], tr.numNodes))

iterNum = 1

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

TreeTable.get_nodenumbers_above_node(tr, tr.root, nodeIndex_array, iterNum, indexNum_table=indexNum_table )


#######################################################
# Tree with a 2-degree node inside a branch
#######################################################
include("/GitHub/BioGeoJulia.jl/notes/TreeTable.jl")

great_ape_newick_string = "((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
tr

nodeIndex_array = collect(repeat([0], tr.numNodes))

iterNum = 1

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

TreeTable.get_nodenumbers_above_node(tr, tr.root, nodeIndex_array, iterNum, indexNum_table=indexNum_table )
"""

function get_nodenumbers_above_node(tr, rootnodenum, nodeIndex_array, iterNum; indexNum_table)
  if (tr.node[rootnodenum].leaf != true)
  	# * A typical internal node will be attached to 3 edges
  	#   (left descendant, right descendant, ancestor edge)
  	# * A root node will be attached to 2 edges (left, right descendant edges
  	# * A degree-2 (mid-branch) node will have 1 descendant edge

  	if (length(tr.node[rootnodenum].edge) > 3)
  		txt = join(["STOP ERROR in get_nodenumbers_above_node(): tr.node[rootnodenum=", string(rootnodenum), "] has more than 3 edges. Probably this is a multifurcating node, which is not allowed."])
  		error(txt)
  	end

  	
  	# Is the current node a typical internal node?
  	typicalTF = length(tr.node[rootnodenum].edge) == 3
  	
  	# Is the current node the root?
  	root_PNnumber = tr.node[tr.root].number  # PhyloNetworks node number of root
  	# rootnodenum = current node being examined
  	current_PNnumber = tr.node[rootnodenum].number
  	rootTF = root_PNnumber == current_PNnumber
		
		# If typical or root, proceed
		typical_or_root_TF = typicalTF || rootTF
  	if (typical_or_root_TF == true)
			# Left descendant edge
			one_edge = tr.node[rootnodenum].edge[1]
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			left_dec_PNnumber = anc_decPNnumbers[2]
			left_dec_nodeIndex = get_nodeIndex_from_PNnumber(left_dec_PNnumber, indexNum_table=indexNum_table)
		
			one_edge = tr.node[rootnodenum].edge[2]
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			right_dec_PNnumber = anc_decPNnumbers[2]
			right_dec_nodeIndex = get_nodeIndex_from_PNnumber(right_dec_PNnumber, indexNum_table=indexNum_table)
		
			# Then, iterate through left and right clades
			#println(rootnodenum)
			nodeIndex_array[iterNum] = rootnodenum
			#print(nodeIndex_array)
			iterNum = iterNum + 1
			(nodeIndex_array, iterNum) = get_nodenumbers_above_node(tr, right_dec_nodeIndex, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
			(nodeIndex_array, iterNum) = get_nodenumbers_above_node(tr, left_dec_nodeIndex, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
			return (nodeIndex_array, iterNum)
		else # It must be a nontypical node (i.e. a 2-degree, inside branch)
			# Check number of edges
			num_edges_attached_to_this_node = length(tr.node[rootnodenum].edge)
			if (num_edges_attached_to_this_node != 2)
				txt = ["ERROR in get_nodenumbers_above_node(): node ", string(rootnodenum), " has ", string(num_edges_attached_to_this_node), " edges attached. It should have exactly 2. Exiting with error."]
				error(join(txt, ""))
			end # END if (num_edges_attached_to_this_node != 2)
			
			one_edge = tr.node[rootnodenum].edge[1]

			# Ancestor node PNnumber for the descendant edge
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			dec_PNnumber = anc_decPNnumbers[2]
			dec_nodeIndex = get_nodeIndex_from_PNnumber(dec_PNnumber, indexNum_table=indexNum_table)
			
			# Ancestor edge not needed
			
			# Then, iterate through single descendant clade
			nodeIndex_array[iterNum] = rootnodenum
			iterNum = iterNum + 1
			(nodeIndex_array, iterNum) = get_nodenumbers_above_node(tr, dec_nodeIndex, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
			return (nodeIndex_array, iterNum)
		end # END if (typical_or_root_TF == true)
  else
  	# Leaf node:
  	#println(rootnodenum)
  	nodeIndex_array[iterNum] = rootnodenum
  	#print(nodeIndex_array)
  	iterNum = iterNum + 1
  	return (nodeIndex_array, iterNum)
  end

	# Error check
	txt = ["ERROR in get_nodenumbers_above_node(): You shouldn't be able to get to the last line of this function."]
	error(join(txt, ""))

end # END get_nodenumbers_above_node


"""
using DataFrames
using PhyloBits

# For Nick's editing (ignore)
include("/GitHub/BioGeoJulia.jl/notes/TreeTable.jl")

#######################################################
# Typical bifurcating (binary) tree
#######################################################
great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr

nodeIndex_array = collect(repeat([0], tr.numNodes))

iterNum = 1

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

TreeTable.get_postorder_nodenumbers_above_node(tr, tr.root, nodeIndex_array, iterNum, indexNum_table=indexNum_table)


#######################################################
# Tree with a 2-degree node inside a branch
#######################################################
include("/GitHub/BioGeoJulia.jl/notes/TreeTable.jl")

great_ape_newick_string = "((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
tr

nodeIndex_array = collect(repeat([0], tr.numNodes))

iterNum = 1

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

TreeTable.get_postorder_nodenumbers_above_node(tr, tr.root, nodeIndex_array, iterNum, indexNum_table=indexNum_table )
"""

function get_postorder_nodenumbers_above_node(tr, rootnodenum, nodeIndex_array, iterNum; indexNum_table)
  if (tr.node[rootnodenum].leaf != true)
  	# * A typical internal node will be attached to 3 edges
  	#   (left descendant, right descendant, ancestor edge)
  	# * A root node will be attached to 2 edges (left, right descendant edges
  	# * A degree-2 (mid-branch) node will have 1 descendant edge

  	if (length(tr.node[rootnodenum].edge) > 3)
  		txt = join(["STOP ERROR in get_postorder_nodenumbers_above_node(): tr.node[rootnodenum=", string(rootnodenum), "] has more than 3 edges. Probably this is a multifurcating node, which is not allowed."])
  		error(txt)
  	end
  	
  	# Is the current node a typical internal node?
  	typicalTF = length(tr.node[rootnodenum].edge) == 3
  	
  	# Is the current node the root?
  	root_PNnumber = tr.node[tr.root].number  # PhyloNetworks node number of root
  	# rootnodenum = current node being examined
  	current_PNnumber = tr.node[rootnodenum].number
  	rootTF = root_PNnumber == current_PNnumber
		
		# If typical or root, proceed
		typical_or_root_TF = typicalTF || rootTF
  	if (typical_or_root_TF == true)
			# Left descendant edge
			one_edge = tr.node[rootnodenum].edge[1]
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			left_dec_PNnumber = anc_decPNnumbers[2]
			left_dec_nodeIndex = get_nodeIndex_from_PNnumber(left_dec_PNnumber, indexNum_table=indexNum_table)
		
			one_edge = tr.node[rootnodenum].edge[2]
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			right_dec_PNnumber = anc_decPNnumbers[2]
			right_dec_nodeIndex = get_nodeIndex_from_PNnumber(right_dec_PNnumber, indexNum_table=indexNum_table)
		
			# Then, iterate through left and right clades
			#println(rootnodenum)
			nodeIndex_array[iterNum] = rootnodenum
			iterNum = iterNum + 1
			(nodeIndex_array, iterNum) = get_postorder_nodenumbers_above_node(tr, right_dec_nodeIndex, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
			(nodeIndex_array, iterNum) = get_postorder_nodenumbers_above_node(tr, left_dec_nodeIndex, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
			#print(nodeIndex_array)

			return (nodeIndex_array, iterNum)
		else # It must be a nontypical node (i.e. a 2-degree, inside branch)
			# Check number of edges
			num_edges_attached_to_this_node = length(tr.node[rootnodenum].edge)
			if (num_edges_attached_to_this_node != 2)
				txt = ["ERROR in get_postorder_nodenumbers_above_node(): node ", string(rootnodenum), " has ", string(num_edges_attached_to_this_node), " edges attached. It should have exactly 2. Exiting with error."]
				error(join(txt, ""))
			end # END if (num_edges_attached_to_this_node != 2)
			
			one_edge = tr.node[rootnodenum].edge[1]

			# Ancestor node PNnumber for the descendant edge
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			dec_PNnumber = anc_decPNnumbers[2]
			dec_nodeIndex = get_nodeIndex_from_PNnumber(dec_PNnumber, indexNum_table=indexNum_table)
			
			# Ancestor edge not needed
			
			# Then, iterate through single descendant clade
			nodeIndex_array[iterNum] = rootnodenum
			iterNum = iterNum + 1
			(nodeIndex_array, iterNum) = get_postorder_nodenumbers_above_node(tr, dec_nodeIndex, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
			return (nodeIndex_array, iterNum)
		end # END if (typical_or_root_TF == true)
  else
  	# Leaf node:
  	#println(rootnodenum)
  	nodeIndex_array[iterNum] = rootnodenum
  	iterNum = iterNum + 1
  	#print(nodeIndex_array)
  	return (nodeIndex_array, iterNum)
  end

	# Error check
	txt = ["ERROR in get_postorder_nodenumbers_above_node(): You shouldn't be able to get to the last line of this function."]
	error(join(txt, ""))

end # END get_postorder_nodenumbers_above_node


function initialize_edgematrix(tr)
	#ancNodeIndex = collect(repeat([0], 2*(tr.numNodes-tr.numTaxa)))
  #decNodeIndex = collect(repeat([0], 2*(tr.numNodes-tr.numTaxa)))
  ancNodeIndex = collect(repeat([0], length(tr.edge)))
	decNodeIndex = collect(repeat([0], length(tr.edge)))
  edgematrix = hcat(ancNodeIndex,decNodeIndex)
  return(edgematrix)
end


# Returns an edgematrix, with
# 1st column gives ancestral nodeIndex for the edge
# 2nd column gives descendant nodeIndex for the edge
# The edges are in order:
#   pair of edges descending from the root
#   pairs following up the right branch above root
#   pairs following up the left branch above root
#
# Iterate backwards for a postorder traversal in "pruningwise" order 
# (as in APE phylo.reorder "pruningwise")
#

"""
using DataFrames
using PhyloBits
#using PhyloPlots


great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,orangutan:12);"
great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,(orang1:1,orang2:1):11);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root





rootnodenum = tr.root
#tipNodeIndex_array = collect(repeat([0], length(tr.numTaxa)))
NodeIndex_array = collect(repeat([0], 2*(tr.numNodes-tr.numTaxa)));
iterNum = 0;
indexNum_table = get_nodeIndex_PNnumber(tr);
uppass_edgematrix = initialize_edgematrix(tr)
res = get_pruningwise_postorder_edgematrix(tr, rootnodenum);
uppass_edgematrix = res[1]

get_LR_uppass_edgematrix(tr)
get_LR_downpass_edgematrix(tr)

"""

function get_pruningwise_postorder_edgematrix(tr, rootnodenum, iterNum=1; edgematrix=initialize_edgematrix(tr), indexNum_table=get_nodeIndex_PNnumber(tr))
  if (tr.node[rootnodenum].leaf != true)
  	# * A typical internal node will be attached to 3 edges
  	#   (left descendant, right descendant, ancestor edge)
  	# * A root node will be attached to 2 edges (left, right descendant edges
  	# * A degree-2 (mid-branch) node will have 1 descendant edge
  	
  	if (length(tr.node[rootnodenum].edge) > 3)
  		txt = join(["STOP ERROR in get_pruningwise_postorder_edgematrix(): tr.node[rootnodenum=", string(rootnodenum), "] has more than 3 edges. Probably this is a multifurcating node, which is not allowed."])
  		error(txt)
  	end
  	
  	
  	# Is the current node a typical internal node?
  	typicalTF = length(tr.node[rootnodenum].edge) == 3
  	
  	# Is the current node the root?
  	root_PNnumber = tr.node[tr.root].number  # PhyloNetworks node number of root
  	# rootnodenum = current node being examined
  	current_PNnumber = tr.node[rootnodenum].number
  	rootTF = root_PNnumber == current_PNnumber
		
		# If typical or root, proceed
		typical_or_root_TF = typicalTF || rootTF
  	if (typical_or_root_TF == true)
			# Left descendant edge
			one_edge = tr.node[rootnodenum].edge[1]
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			left_anc_PNnumber = anc_decPNnumbers[1]
			left_dec_PNnumber = anc_decPNnumbers[2]
			left_anc_nodeIndex = get_nodeIndex_from_PNnumber(left_anc_PNnumber, indexNum_table=indexNum_table)
			left_dec_nodeIndex = get_nodeIndex_from_PNnumber(left_dec_PNnumber, indexNum_table=indexNum_table)
		
			# Right descendant edge
			one_edge = tr.node[rootnodenum].edge[2]
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			right_anc_PNnumber = anc_decPNnumbers[1]
			right_dec_PNnumber = anc_decPNnumbers[2]
			right_anc_nodeIndex = get_nodeIndex_from_PNnumber(right_anc_PNnumber, indexNum_table=indexNum_table)
			right_dec_nodeIndex = get_nodeIndex_from_PNnumber(right_dec_PNnumber, indexNum_table=indexNum_table)
		
			# Then, iterate through left and right clades
			#println(rootnodenum)
			edgematrix[iterNum,1] = right_anc_nodeIndex
			edgematrix[iterNum,2] = right_dec_nodeIndex
			iterNum = iterNum + 1
			edgematrix[iterNum,1] = left_anc_nodeIndex
			edgematrix[iterNum,2] = left_dec_nodeIndex
			iterNum = iterNum + 1
		
			(edgematrix, iterNum) = get_pruningwise_postorder_edgematrix(tr, right_dec_nodeIndex, iterNum, edgematrix=edgematrix, indexNum_table=indexNum_table)
			(edgematrix, iterNum) = get_pruningwise_postorder_edgematrix(tr, left_dec_nodeIndex, iterNum, edgematrix=edgematrix, indexNum_table=indexNum_table)
			#print(nodeIndex_array)
			return (edgematrix, iterNum)
		else
			# Single descendant edge
			one_edge = tr.node[rootnodenum].edge[1]
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			anc_PNnumber = anc_decPNnumbers[1]
			dec_PNnumber = anc_decPNnumbers[2]
			anc_nodeIndex = get_nodeIndex_from_PNnumber(anc_PNnumber, indexNum_table=indexNum_table)
			dec_nodeIndex = get_nodeIndex_from_PNnumber(dec_PNnumber, indexNum_table=indexNum_table)

			# Then, iterate through the single descendant clade
			#println(rootnodenum)
			edgematrix[iterNum,1] = anc_nodeIndex
			edgematrix[iterNum,2] = dec_nodeIndex
			iterNum = iterNum + 1
			
			(edgematrix, iterNum) = get_pruningwise_postorder_edgematrix(tr, dec_nodeIndex, iterNum, edgematrix=edgematrix, indexNum_table=indexNum_table)
			return (edgematrix, iterNum)
		end # END if (typical_or_root_TF == true)
  else
  	#println(rootnodenum)
  	#iterNum = iterNum + 1
  	#nodeIndex_array[iterNum] = rootnodenum
  	#print(nodeIndex_array)
  	return (edgematrix, iterNum)
  end
	
	# Shouldn't get here
	error("get_pruningwise_postorder_edgematrix(): Shouldn't get here!")
end




# Get the node indexes for an uppass from the root to the tips
# Reverse for downpass (postorder)
# Uses get_nodenumbers_above_node() recursively to iterate up tree

"""
using DataFrames
using PhyloBits
#using PhyloPlots


great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

uppass_edgematrix = get_LR_uppass_edgematrix(tr)
uppass_edgematrix
"""

function get_LR_uppass_edgematrix(tr)
	rootnodenum = tr.root
	iterNum = 1
	edgematrix = initialize_edgematrix(tr)
	indexNum_table = get_nodeIndex_PNnumber(tr)

	res = get_pruningwise_postorder_edgematrix(tr, rootnodenum, iterNum, edgematrix=edgematrix, indexNum_table=indexNum_table)
	uppass_edgematrix = res[1]
	return(uppass_edgematrix)
end


# Get the node indexes for an downpass from the root to the tips
# Reverse of get_LR_uppass_nodeIndexes()
# 

"""
using DataFrames
using PhyloBits
#using PhyloPlots


great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

uppass_edgematrix = get_LR_uppass_edgematrix(tr)
uppass_edgematrix

downpass_edgematrix = get_LR_downpass_edgematrix(tr)
downpass_edgematrix
"""
function get_LR_downpass_edgematrix(tr)
	rootnodenum = tr.root
	iterNum = 1
	edgematrix = initialize_edgematrix(tr)
	indexNum_table = get_nodeIndex_PNnumber(tr)

	res = get_pruningwise_postorder_edgematrix(tr, rootnodenum, iterNum, edgematrix=edgematrix, indexNum_table=indexNum_table)
	uppass_edgematrix = res[1]
	
	# Reverse
	numrows = size(uppass_edgematrix)[1]
	reverse_rownums = seq(numrows,1,-1)
	downpass_edgematrix = uppass_edgematrix[reverse_rownums,:]
	return(downpass_edgematrix)
end



# Preorder traversal (I think)
"""
using DataFrames
using PhyloBits
#using PhyloPlots


great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
#great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,orangutan:12);"
#great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,(orang1:1,orang2:1):11);"
tr = readTopology(great_ape_newick_string)

# Preorder traversal
uppass_nodeIndexes = get_LR_uppass_nodeIndexes(tr)

# Reverse preorder traversal (not the same as postorder!)
downpass_nodeIndexes = get_LR_downpass_nodeIndexes(tr)
"""
function get_LR_uppass_nodeIndexes(tr)
	rootnodenum = tr.root
	iterNum = 1
	nodeIndex_array = collect(repeat([0], tr.numNodes))
	indexNum_table = get_nodeIndex_PNnumber(tr)

	res = get_nodenumbers_above_node(tr, rootnodenum, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
	uppass_nodeIndexes = res[1]
	return(uppass_nodeIndexes)
end

# Reverse-Preorder traversal (I think)
"""
using DataFrames
using PhyloBits
#using PhyloPlots


great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
#great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,orangutan:12);"
#great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,(orang1:1,orang2:1):11);"
tr = readTopology(great_ape_newick_string)

# Preorder traversal
uppass_nodeIndexes = get_LR_uppass_nodeIndexes(tr)

# Reverse preorder traversal (not the same as postorder!)
downpass_nodeIndexes = get_LR_downpass_nodeIndexes(tr)
"""
function get_LR_downpass_nodeIndexes(tr)
	rootnodenum = tr.root
	iterNum = 1
	nodeIndex_array = collect(repeat([0], tr.numNodes))
	indexNum_table = get_nodeIndex_PNnumber(tr)

	res = get_nodenumbers_above_node(tr, rootnodenum, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
	downpass_nodeIndexes = reverse(res[1])
	return(downpass_nodeIndexes)
end


# Get Rnodenums
# They seem to be just 
# sort(tip_PNnumbers) + sort(abs(internal_PNnumbers))

"""
using DataFrames
using PhyloBits
#using PhyloPlots


great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root
trdf = prt(tr, rootnodenum)
df

# Get the edges, get the nodeIndexes corresponding to these
edge = tr.edge
edge

one_edge = edge[1]

edge_df = get_NodeIndex_df_by_tree_edges(tr, indexNum_table=indexNum_table)
edge_df

uppass_edgematrix = get_LR_uppass_edgematrix(tr)
uppass_edgematrix

downpass_edgematrix = get_LR_downpass_edgematrix(tr)
downpass_edgematrix


# Get the R node numbers, append to indexNum_table
Rnodenums_in_indexNum_table_order = get_Rnodenums(tr, indexNum_table)
indexNum_table2 = hcat(indexNum_table, Rnodenums_in_indexNum_table_order)
indexNum_table2

indexNum_table3 = indexNum_table2[sortperm(indexNum_table2[:,3]),:]
indexNum_table3

"""


function get_Rnodenums(tr, indexNum_table)
	numnodes = length(tr.node)
	Rnodenums = collect(1:numnodes)
	tipsTF = indexNum_table[:,2] .> 0
	internalTF = tipsTF .== false
	
	numtips = sum(tipsTF)
	numinternal = sum(internalTF)
	
	
	PNnumbers_in_Rnodenums_order = collect(repeat([0], numnodes))
	PNnumbers_in_Rnodenums_order[1:numtips] = sort(indexNum_table[:,2][tipsTF])
	
	PNnumbers_in_Rnodenums_order[(numtips+1):(numtips+numinternal)] = reverse(sort(indexNum_table[:,2][internalTF]))
	
	
	tmpmat = hcat(Rnodenums, PNnumbers_in_Rnodenums_order)
	
	# Match to indexNum_table
	indices = collect(1:numnodes)
	Rnodenums_in_indexNum_table_order = collect(repeat([0], numnodes))
	for i in 1:numnodes
		current_PNnumber = indexNum_table[i,2] 
		TF = tmpmat[:,2] .== current_PNnumber
		rownum_in_tmpmat = indices[TF][1]
		Rnodenums_in_indexNum_table_order[i] = tmpmat[rownum_in_tmpmat,1]
	end
	return(Rnodenums_in_indexNum_table_order)
end





"""
using DataFrames
using PhyloBits
#using PhyloPlots

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

"""
function get_nodeIndex_PNnumber(tr)
	# Get a numNode x 2 table
	# Index, then ".number"
	numnodes = length(tr.node)
	indexNum_table = Array{Int}(undef, numnodes, 2)
	for i in 1:numnodes
		indexNum_table[i,1] = i
		indexNum_table[i,2] = tr.node[i].number
	end
	return(indexNum_table)
end



# Go from a PhyloNetwork Node Number (PNnumber) to the node index 
# (i.e., the index of that node in the list of nodes)
"""
using DataFrames
using PhyloBits


great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

PNnumber = 3;
get_nodeIndex_from_PNnumber(PNnumber; indexNum_table)

PNnumber = -2;
get_nodeIndex_from_PNnumber(PNnumber; indexNum_table)

PNnumber = -4;
get_nodeIndex_from_PNnumber(PNnumber; indexNum_table)


using DataFrames
using PhyloBits

#######################################################
# Typical bifurcating (binary) tree
#######################################################
great_ape_newick_string = "((human:1.0,chimp:0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
tr

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

PNnumber = 3;
get_nodeIndex_from_PNnumber(PNnumber; indexNum_table)

PNnumber = -2;
get_nodeIndex_from_PNnumber(PNnumber; indexNum_table)

PNnumber = -3;
get_nodeIndex_from_PNnumber(PNnumber; indexNum_table)

# Node doesn't exist:
# PNnumber = -4;
# get_nodeIndex_from_PNnumber(PNnumber; indexNum_table)



#######################################################
# Tree with a 2-degree node inside a branch
#######################################################
great_ape_newick_string = "((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
tr

nodeIndex_array = collect(repeat([0], tr.numNodes))

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

PNnumber = 3;
get_nodeIndex_from_PNnumber(PNnumber; indexNum_table)

PNnumber = -2;
get_nodeIndex_from_PNnumber(PNnumber; indexNum_table)

PNnumber = -4;
get_nodeIndex_from_PNnumber(PNnumber; indexNum_table)

"""
function get_nodeIndex_from_PNnumber(PNnumber; indexNum_table)
	TF01 = indexNum_table[:,2] .== PNnumber
	# Adding the [1] returns a scalar
	nodeIndex = indexNum_table[:,1][TF01][1]
	return nodeIndex
end


"""
# Tree with a direct-ancestor node (degree-2)
great_ape_newick_string = "((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
nonroot_nodes = get_nonrootnodes(tr, include_direct_ancestors=false)
nonroot_nodes = get_nonrootnodes(tr, include_direct_ancestors=true)
trdf = prt(tr)
nonroot_nodes = get_nonrootnodes_trdf2(trdf, include_direct_ancestors=false)
nonroot_nodes = get_nonrootnodes_trdf2(trdf, include_direct_ancestors=true)
"""
# Return a Vector of node numbers for a tree, without the root node
function get_nonrootnodes(tr; include_direct_ancestors=false)
	nonroot_nodeIndexes = deleteat!(collect(1:tr.numNodes), tr.root)
	if include_direct_ancestors == true
		return nonroot_nodeIndexes
	else
		TF = collect(repeat([true], length(nonroot_nodeIndexes)))
		# Remove nonroot nodes with only 2 edges, as they are direct ancestors
		for i in 1:length(nonroot_nodeIndexes)
			TF[i] = length(tr.node[i].edge) != 2  # degree-2 nodes
		end
		nonroot_nodes = nonroot_nodeIndexes[TF]
		return nonroot_nodes
	end
end # END function get_nonrootnodes(tr)

"""
# Tree with a direct-ancestor node (degree-2)
great_ape_newick_string = "((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
nonroot_nodes = get_nonrootnodes(tr, include_direct_ancestors=false)
nonroot_nodes = get_nonrootnodes(tr, include_direct_ancestors=true)
trdf = prt(tr)
nonroot_nodes = get_nonrootnodes_trdf2(trdf, include_direct_ancestors=false)
nonroot_nodes = get_nonrootnodes_trdf2(trdf, include_direct_ancestors=true)
"""
function get_nonrootnodes_trdf(trdf; include_direct_ancestors=false)
	TF1 = trdf[!,:nodeType] .== "tip"
	TF2 = trdf[!,:nodeType] .== "intern"
	TF3 = trdf[!,:nodeType] .== "direct" # degree-2 node
	TF4 = trdf[!,:nodeType] .== "root"	 # degree-2 node, but root
	
	# If include_direct_ancestors = false
	# Remove nonroot nodes with only 2 edges, as they are direct ancestors
	if include_direct_ancestors == false
		TF = (TF1 .+ TF2 .+ TF4) .== 1
	else
		# If include_direct_ancestors = true
		TF = (TF1 .+ TF2 .+ TF3 .+ TF4) .== 1
	end

	num_tip_int_root_nodes = sum(TF)
	nodenums_tip_int_root = trdf[!,:nodeIndex][TF]
	root_nodenum = trdf[!,:nodeIndex][TF4]
	nonroot_nodes = nodenums_tip_int_root[nodenums_tip_int_root .!= root_nodenum]
	return(nonroot_nodes)
end # END function get_nonrootnodes_trdf(trdf)

"""
# Get the unique node times, in ascending order from 0.0 (the tips).
using DataFrames
using PhyloBits

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
trdf = prt(tr)
nodetimes(trdf)
nodetimes(tr)

"""
function nodetimes(trdf::DataFrame; tipDateEq0=1.0e-7)
	node_ages = trdf.node_age
	# Remove effective zeros; sort to unique
	uniq_node_ages = sort(unique(node_ages[node_ages .> tipDateEq0]))
	return(uniq_node_ages)
end


function nodetimes(tr::HybridNetwork; tipDateEq0=1.0e-7)
	trdf = prt(tr)
	node_ages = trdf.node_age
	# Remove effective zeros; sort to unique
	uniq_node_ages = sort(unique(node_ages[node_ages .> tipDateEq0]))
	return(uniq_node_ages)
end


function get_fossils(tr::HybridNetwork; fossils_older_than=1.0e-5)
	trdf = prt(tr, rootnodenum=tr.root, get_taxa_by_node=false);
	fossils_TF = get_fossils(trdf; fossils_older_than=fossils_older_than)
	return fossils_TF
end

# Get NODES (TF for all nodes, including tips) which are fossils (tips, including hooktips, not hooknodes)
function get_fossils(trdf::DataFrame; fossils_older_than=1.0e-5)
	TF1 = trdf.nodeType .== "tip"
	TF2 = trdf.node_age .>= fossils_older_than
	fossils_TF = (TF1 .+ TF2) .== 2
	return fossils_TF
end

# Get NODES (TF for all nodes, including tips) which are fossils (tips, including hooktips, not hooknodes)
function get_hooks(tr::HybridNetwork; hooks_below=1.0e-6)
	trdf = prt(tr, rootnodenum=tr.root, get_taxa_by_node=false);
	hooknode_TF = get_hooks(trdf; hooks_below=hooks_below)
	return hooknode_TF
end

function get_hooks(trdf::DataFrame; hooks_below=1.0e-6)
	hooknode_TF = repeat([false], nrow(trdf))
	for nodeIndex in 1:nrow(trdf)
		# Only internal nodes
		if (trdf.nodeType[nodeIndex] == "intern") || (trdf.nodeType[nodeIndex] == "root")
			lnode = trdf.leftNodeIndex[nodeIndex]
			rnode = trdf.rightNodeIndex[nodeIndex]
			brlen_above_Left_corner = trdf.brlen[lnode]
			brlen_above_Right_corner = trdf.brlen[rnode]
			
			# Is this a hooknode?
			if (brlen_above_Left_corner <= hooks_below) || (brlen_above_Right_corner <= hooks_below)
				hooknode_TF[nodeIndex] = true
			end
			
			# Is this a hooktip?
			if (brlen_above_Left_corner <= hooks_below)
				hooknode_TF[lnode] = true
			end
			if (brlen_above_Right_corner <= hooks_below)
				hooknode_TF[rnode] = true
			end
			
		end # END if (trdf.nodeType[nodeIndex] == "intern") |
	end # END for nodeIndex in 1:nrow(trdf)
	return hooknode_TF
end

"""
using DataFrames
using PhyloBits


great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root
trdf = prt(tr, rootnodenum)
trdf


using DataFrames
using PhyloBits

#######################################################
# Typical bifurcating (binary) tree
#######################################################
great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr

nodeIndex_array = collect(repeat([0], tr.numNodes))

iterNum = 1

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

TreeTable.get_postorder_nodenumbers_above_node(tr, tr.root, nodeIndex_array, iterNum, indexNum_table=indexNum_table)


#######################################################
# Typical tree with a fossil
#######################################################
great_ape_newick_string = "((((human:2.5,Lucy:2.5):3.5,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr

nodeIndex_array = collect(repeat([0], tr.numNodes))

iterNum = 1

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

TreeTable.get_postorder_nodenumbers_above_node(tr, tr.root, nodeIndex_array, iterNum, indexNum_table=indexNum_table)


#######################################################
# Tree with a 2-degree node inside a branch
#######################################################
great_ape_newick_string = "((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);";
tr = readTopology(great_ape_newick_string);
trdf = prt(tr);
trdf

nodeIndex_array = collect(repeat([0], tr.numNodes))

iterNum = 1

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

TreeTable.get_postorder_nodenumbers_above_node(tr, tr.root, nodeIndex_array, iterNum, indexNum_table=indexNum_table )

"""
# Return a DataFrame with the edge numbers
function prt(tr::HybridNetwork; rootnodenum::Int64=tr.root, get_taxa_by_node::Bool=true,  hooks_below::Float64=1.0e-6, fossils_older_than::Float64=1.0e-5)
	#using DataFrames
	numnodes = length(tr.node)
	# number of digits for internal node numbers
	numdigits = length(digits(numnodes))
	if (numdigits < 2)
		numdigits = 2
	end
	
	# Initialize the dataframe
	trdf = DataFrames.DataFrame(nodeIndex=collect(1:numnodes), PNnumber=collect(repeat([0],numnodes)))
	
	# Fill in the PNnumber node numbers
	indexNum_table = get_nodeIndex_PNnumber(tr)
	trdf[!, :PNnumber] = indexNum_table[:,2]
	trdf
	
	# Add the R node numbers
	Rnodenums = get_Rnodenums(tr, indexNum_table)
	trdf[!,:Rnodenums] = Rnodenums
	trdf
	
	
	# Add the node ages
	node_age = get_node_ages(tr)
	trdf[!,:node_age] = node_age
	
	# Add branch lengths
	brlen = collect(repeat([0.0], numnodes))
	ancNodeIndex = collect(repeat([0], numnodes))
	leftNodeIndex = collect(repeat([0], numnodes))
	rightNodeIndex = collect(repeat([0], numnodes))
	nodeName = collect(repeat([""], numnodes))
	nodeType = collect(repeat([""], numnodes))
	
	edge_df = get_NodeIndex_df_by_tree_edges(tr, indexNum_table=indexNum_table)
	edge_df_rownums = collect(1:Rnrow(edge_df))
	for i in 1:Rnrow(trdf)
		nodeIndex = trdf[i,:nodeIndex]
		TF = edge_df[:,:edge_decNodeIndex] .== nodeIndex
		if (sum(TF) > 0)
			edge_df_rownum = edge_df_rownums[TF][1]
			brlen[i] = edge_df[:,:edge_length][edge_df_rownum]
			ancNodeIndex[i] = edge_df[:,:edge_ancNodeIndex][edge_df_rownum]
		else
			brlen[i] = -999.0
			ancNodeIndex[i] = -999
		end
		
		# Tip node
		if (tr.node[nodeIndex].leaf == true)
			# Tip node
			leftNodeIndex[i] = -999
			rightNodeIndex[i] = -999
			nodeType[i] = "tip"
			nodeName[i] = tr.node[nodeIndex].name
		end
		
		# Left and right descendant nodeIndexes
		TF_edgerows_descend_from_decNodeIndex = edge_df[:,:edge_ancNodeIndex] .== nodeIndex
		
		# Direct-ancestor node
		if ((sum(TF_edgerows_descend_from_decNodeIndex) == 1) && (tr.node[nodeIndex].leaf == false))
			leftNodeIndex[i] = edge_df[TF_edgerows_descend_from_decNodeIndex,:edge_decNodeIndex][1]
			rightNodeIndex[i] = -99999 # No right descendant, as it's a direct ancestor
			nodeType[i] = "direct"  # da = direct ancestor node

			# Node names
			tmp_nodeName = tr.node[nodeIndex].name
			if (tmp_nodeName == "")
				internal_nodeIndex_as_string = lpad(string(nodeIndex), numdigits, "0")
				nodeName[i] = join(["da", internal_nodeIndex_as_string], "")
			else
				nodeName[i] = tmp_nodeName
			end
		end # END if (sum(TF_edgerows_descend_from_decNodeIndex) == 1)
		
		# Bifurcating node
		if (sum(TF_edgerows_descend_from_decNodeIndex) == 2)
			# Internal node
			leftNodeIndex[i] = edge_df[TF_edgerows_descend_from_decNodeIndex,:edge_decNodeIndex][1]
			rightNodeIndex[i] = edge_df[TF_edgerows_descend_from_decNodeIndex,:edge_decNodeIndex][2]
			nodeType[i] = "intern"
			if (nodeIndex == tr.root)
				nodeType[i] = "root"
			end # END if (nodeIndex == tr.root)
			
			# Node names
			tmp_nodeName = tr.node[nodeIndex].name
			if (tmp_nodeName == "")
				internal_nodeIndex_as_string = lpad(string(nodeIndex), numdigits, "0")
				nodeName[i] = join(["in", internal_nodeIndex_as_string], "")
			else
				nodeName[i] = tmp_nodeName
			end
		end # END if (sum(TF_edgerows_descend_from_decNodeIndex) == 2)
	end # END for i in 1:Rnrow(trdf)
	
	# Regimes: the field "reg" holds the regimes
	reg = collect(repeat([1], numnodes))
	
	# Add the fields
	trdf[!,:brlen] = brlen
	trdf[!,:ancNodeIndex] = ancNodeIndex
	trdf[!,:leftNodeIndex] = leftNodeIndex
	trdf[!,:rightNodeIndex] = rightNodeIndex
	trdf[!,:nodeName] = nodeName
	trdf[!,:nodeType] = nodeType
	
	# Identify hook nodes and tips
	hook = hooknode_TF = get_hooks(trdf; hooks_below=hooks_below)
	trdf[!,:hook] = hook
	
	# Identify fossil tips (including hooktips)
	TF1 = trdf.nodeType .== "tip"
	TF2 = trdf.node_age .>= fossils_older_than
	foss = (TF1 .+ TF2) .== 2
	trdf[!,:foss] = foss

	trdf[!,:reg] = reg
	
	if get_taxa_by_node == true
		downpass_edgematrix = get_LR_downpass_edgematrix(tr)
		taxa = get_taxa_descending_from_each_node(tr, trdf, downpass_edgematrix=downpass_edgematrix)
		trdf[!,:taxa] = taxa
	end
	
	return(trdf)
end # END function prt(tr, rootnodenum=tr.root, get_taxa_by_node=true)



# Get tipnames descending from each node
# The list is comma-delimited and alphabetical
"""
using DataFrames
using PhyloBits
#using PhyloPlots

include("/GitHub/BioGeoJulia.jl/notes/TreeTable.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root

# Get trdf WITHOUT "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, false)

# Not all columns print now that trdf is big, let's print just the left and right
# Get trdf WITH "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, true);
headLR(trdf)

downpass_edgematrix = get_LR_downpass_edgematrix(tr)
taxa = get_taxa_descending_from_each_node(tr, trdf, downpass_edgematrix=get_LR_downpass_edgematrix(tr))



#######################################################
# Tree with a 2-degree node inside a branch
#######################################################
#include("/GitHub/BioGeoJulia.jl/notes/TreeTable.jl")

great_ape_newick_string = "((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root

# Get trdf WITHOUT "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, false)

# Not all columns print now that trdf is big, let's print just the left and right
# Get trdf WITH "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, true);
headLR(trdf)

downpass_edgematrix = get_LR_downpass_edgematrix(tr)
taxa = get_taxa_descending_from_each_node(tr, trdf, downpass_edgematrix=get_LR_downpass_edgematrix(tr))

"""

function get_taxa_descending_from_each_node(tr, trdf; downpass_edgematrix=get_LR_downpass_edgematrix(tr))
	# Get sizes
	numnodes = length(tr.node)
	node_indices_in_downpass_edgematrix_tmp = flat2(downpass_edgematrix)
	node_indices_in_downpass_edgematrix = node_indices_in_downpass_edgematrix_tmp[node_indices_in_downpass_edgematrix_tmp .!== 0]
	numIndices = length(sort(unique(node_indices_in_downpass_edgematrix)))
	
	# Remove 0s from edge matrix, as they constitute "fake" right branches
	
	
	
	# Error check
	if (numnodes != numIndices)
		txt = ["ERROR in get_taxa_descending_from_each_node(): the number of nodes in tr (", string(numnodes), ") does not match the number of rows in downpass_edgematrix (", string(numIndices), ")."]
		error(join(txt, ""))
	end
	
	taxa = collect(repeat([""], numnodes))
	
	# Assumes a binary tree, and a nodeIndexes df in downpass order, from
	# downpass_edgematrix = get_LR_downpass_edgematrix(tr)
	
	# Step through the downpass_edgematrix, in pairs
	edgematrix_rows_to_visit = collect(1:2:Rnrow(downpass_edgematrix))
	
	numrows_in_downpass_edgematrix = Rnrow(downpass_edgematrix)
	i = 1
	while i <= numrows_in_downpass_edgematrix
		j = i+1
		nodeIndex_ancL = downpass_edgematrix[i,1]
		nodeIndex_ancR = downpass_edgematrix[j,1]
		
		# If direct ancestor node
		if (nodeIndex_ancL != nodeIndex_ancR)
			nodeIndex_left = downpass_edgematrix[i,2]
			nodeIndex_right = -99999
			nodeIndex_anc = nodeIndex_ancL
			
			# Get the taxa
			tmp_taxa = [""]
			if (tr.node[nodeIndex_left].leaf == true)
				taxa[nodeIndex_left] = tr.node[nodeIndex_left].name
				tmp_taxa[1] = tr.node[nodeIndex_left].name
			else
				tmp_taxa[1] = taxa[nodeIndex_left]
			end

			# Put the taxa in the ancestral node
			ancNodeIndex_of_Left = nodeIndex_ancL
			ancNodeIndex = ancNodeIndex_of_Left
			taxa_unordered = flat2(split.(tmp_taxa, ","))
			taxa_ordered = sort(taxa_unordered)
			taxa_at_ancNode = join(taxa_ordered, ",")
			taxa[ancNodeIndex] = taxa_at_ancNode
		
			# Advance i by just 1 (because this was a node with a single descendant)
			i = i+1
		end # END if (nodeIndex_ancL != nodeIndex_ancR) # (end direct ancestor node)

		# Bifurcating node, 2 adjacent edges
		if (nodeIndex_ancL == nodeIndex_ancR)
			nodeIndex_anc = nodeIndex_ancL
			nodeIndex_left = downpass_edgematrix[i,2]
			nodeIndex_right = downpass_edgematrix[j,2]

			tmp_taxa = collect(repeat([""], 2))
			if (tr.node[nodeIndex_left].leaf == true)
				taxa[nodeIndex_left] = tr.node[nodeIndex_left].name
				tmp_taxa[1] = tr.node[nodeIndex_left].name
			else
				tmp_taxa[1] = taxa[nodeIndex_left]
			end
			if (tr.node[nodeIndex_right].leaf == true)
				taxa[nodeIndex_right] = tr.node[nodeIndex_right].name
				tmp_taxa[2] = tr.node[nodeIndex_right].name
			else
				tmp_taxa[2] = taxa[nodeIndex_right]
			end
		
			# Join the two lists and put in the ancestral node
			ancNodeIndex_of_Left = nodeIndex_ancL
			ancNodeIndex_of_Right = nodeIndex_ancR
			if (ancNodeIndex_of_Left != ancNodeIndex_of_Right)
				txt = ["ERROR in get_taxa_descending_from_each_node(): ancNodeIndex_of_Left must match ancNodeIndex_of_Right in trdf, but doesn't.\nAncestor of nodeIndex_left ", string(nodeIndex_left), " is ", string(ancNodeIndex_of_Left), ", but\nancestor of nodeIndex_right ", string(nodeIndex_right), " is ", string(ancNodeIndex_of_Right), "\n"]
				msg = join(txt, "")
				println(msg)
				error(msg)
			end
			# No error, so ancNodeIndex was found:
			ancNodeIndex = ancNodeIndex_of_Left
			taxa_unordered = flat2(split.(tmp_taxa, ","))
			taxa_ordered = sort(taxa_unordered)
			taxa_at_ancNode = join(taxa_ordered, ",")
			taxa[ancNodeIndex] = taxa_at_ancNode
			
			# Advance i by 2 (because this was a node with 2 descendants)
			i = i+2
		end # END if (nodeIndex_ancL == nodeIndex_ancR) # (bifurcating node)
	end # END while
	
	
	return(taxa)
end # END function get_taxa_descending_from_each_node









# Is the node a tip?
"""
#######################################################
# Standard bifurcating tree
#######################################################
great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root

# Get trdf WITHOUT "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, false)

# Not all columns print now that trdf is big, let's print just the left and right
# Get trdf WITH "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, true);
headLR(trdf)

downpass_edgematrix = get_LR_downpass_edgematrix(tr)
taxa = get_taxa_descending_from_each_node(tr, trdf, downpass_edgematrix=get_LR_downpass_edgematrix(tr))

# Make a trdf, a DataFrame holding all of the key tree information
cumulative_height_at_each_node = get_node_heights(tr)
node_age = get_node_ages(tr)

nodeIndex = 1
isTip_TF(tr, nodeIndex)
nodeIndex = 5
isTip_TF(tr, nodeIndex)

#######################################################
# Tree with a 2-degree node inside a branch
#######################################################

great_ape_newick_string = "((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root

# Get trdf WITHOUT "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, false)

# Not all columns print now that trdf is big, let's print just the left and right
# Get trdf WITH "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, true);
headLR(trdf)

downpass_edgematrix = get_LR_downpass_edgematrix(tr)
taxa = get_taxa_descending_from_each_node(tr, trdf, downpass_edgematrix=get_LR_downpass_edgematrix(tr))

# Make a trdf, a DataFrame holding all of the key tree information
cumulative_height_at_each_node = get_node_heights(tr)
node_age = get_node_ages(tr)

nodeIndex = 1
isTip_TF(tr, nodeIndex)
nodeIndex = 3
isTip_TF(tr, nodeIndex)
"""

function isTip_TF(tr, nodeIndex)
	# Returns true if it's a tip (a "leaf"), false if not
	return(tr.node[nodeIndex].leaf)
end

# Get the nodeIndex ancestral to 2 nodes



# Elaborate function to MAKE SURE we are getting the ancestor & descendant PNnumbers correct

"""
using DataFrames
using PhyloBits
#using PhyloPlots


great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root
trdf = prt(tr, rootnodenum)
df

# Get the edges, get the nodeIndexes corresponding to these
edge = tr.edge
edge

one_edge = edge[1]

# Every edge has only 2 nodes. 
for i in 1:length(edge)
	one_edge = edge[i]
	anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
	print(anc_decPNnumbers)
end
"""

function get_NodeIndexes_from_edge(one_edge)
	# Make sure there are only 2 nodes for each edge.
	numnodec_on_this_edge = length(one_edge.node)
	if ( numnodec_on_this_edge != 2 )
		array_of_strings = ["\nError in get_NodeIndexes_from_edge(): this function assumes each edge has only 2 nodes, ancestor & descendant.\nInstead, this edge gives length(edge.node)=", string(numnodec_on_this_edge), ".\nPrinting edge to screen...\n"]
		txt = join(array_of_strings, "")
		println(txt)
		print(edge)
		error(txt)
	end
	
	# Initialize
	tmp_nodeIndices = collect(repeat([0], numnodec_on_this_edge))
	tmp_nodePNnumbers = collect(repeat([0], numnodec_on_this_edge))
	ancPNnumber = 0
	decPNnumber = 0
	
	# Record the PN (PhyloNetworks) node numbers
	tmp_nodePNnumbers[1] = one_edge.node[1].number
	tmp_nodePNnumbers[2] = one_edge.node[2].number
	tmp_nodeIndices_on_edge = collect(1:2)
	
	# Positive node numbers are tips; negative node numbers are internal

	# Declare error if both node numbers are the same
	PNnumber_eq_0_TF = tmp_nodePNnumbers .== 0
	if (sum(PNnumber_eq_0_TF) > 0)
		error("Error in get_NodeIndexes_from_edge(): both node numbers attached to this edge are the same. This should be impossible.")
	end

	
	# Declare error if both node numbers are positive (both are tips)
	positiveTF = tmp_nodePNnumbers .> 0
	if (sum(positiveTF) == 2)
		error("Error in get_NodeIndexes_from_edge(): both node numbers attached to this edge are tips (PhyloNetworks node number > 0). This should be impossible. It may mean that you read in a Newick string that already had node labels. You can remove these in R with tr$node.label=NULL.")
	end
	
	# If one PNnumber is positive and one negative, then you have a descendant tip, and ancestor internal node
	if (sum(positiveTF) == 1)
		anc_nodeIndex_on_edge = tmp_nodeIndices_on_edge[positiveTF .== false]
		ancPNnumber = tmp_nodePNnumbers[anc_nodeIndex_on_edge] # internal nodenum
		dec_nodeIndex_on_edge = tmp_nodeIndices_on_edge[positiveTF .== true]
		decPNnumber = tmp_nodePNnumbers[dec_nodeIndex_on_edge]
	end
	
	# If both are negative, then you have 2 internal nodes. The one closest to 0 is the root
	if (sum(positiveTF) == 0)
		min_absval = -1 * minimum(abs.(tmp_nodePNnumbers))
		matches_min_TF = tmp_nodePNnumbers .== min_absval
		anc_nodeIndex_on_edge = tmp_nodeIndices_on_edge[matches_min_TF .== true]
		dec_nodeIndex_on_edge = tmp_nodeIndices_on_edge[matches_min_TF .== false]
		ancPNnumber = tmp_nodePNnumbers[anc_nodeIndex_on_edge] # ancestral internal nodenum
		decPNnumber = tmp_nodePNnumbers[dec_nodeIndex_on_edge] # ancestral internal nodenum
	end
	
	anc_decPNnumbers = [ancPNnumber[1], decPNnumber[1]]
	return(anc_decPNnumbers)
end


"""
using DataFrames
using PhyloBits
#using PhyloPlots


great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

#######################################################
# Standard bifurcating tree
#######################################################
great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root

# Get trdf WITHOUT "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, false)

# Not all columns print now that trdf is big, let's print just the left and right
# Get trdf WITH "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, true);
headLR(trdf)

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

edge_df = get_NodeIndex_df_by_tree_edges(tr, indexNum_table=indexNum_table)

#######################################################
# Tree with a 2-degree node inside a branch
#######################################################

great_ape_newick_string = "((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root

# Get trdf WITHOUT "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, false)

# Not all columns print now that trdf is big, let's print just the left and right
# Get trdf WITH "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, true);
headLR(trdf)

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

edge_df = get_NodeIndex_df_by_tree_edges(tr, indexNum_table=indexNum_table)
"""
function get_NodeIndex_df_by_tree_edges(tr; indexNum_table)
	# Get the ancestral and descendant node numbers for each edge
	num_edges = length(tr.edge)
	edgeIndex = collect(1:num_edges)
	edge_ancNodeIndex = collect(repeat([0], num_edges))
	edge_decNodeIndex = collect(repeat([0], num_edges))
	edge_ancPNnumber = collect(repeat([0], num_edges))
	edge_decPNnumber = collect(repeat([0], num_edges))
	edge_length = collect(repeat([0.0], num_edges))

	# Initialize the dataframe
	edge_df = DataFrames.DataFrame(edgeIndex=edgeIndex, edge_ancNodeIndex=edge_ancNodeIndex, edge_decNodeIndex=edge_decNodeIndex, edge_ancPNnumber=edge_ancPNnumber, edge_decPNnumber=edge_decPNnumber, edge_length=edge_length)
	edge_df

	# Every edge has only 2 nodes. 
	for i in 1:length(tr.edge)
		one_edge = tr.edge[i]
		anc_dec_PNnumbers = get_NodeIndexes_from_edge(one_edge)
	
		# Use a named keyword argument, indexNum_table, to provide the translation from PNnumber to node index
		anc_dec_nodeIndices = get_nodeIndex_from_PNnumber.(anc_dec_PNnumbers, indexNum_table=indexNum_table)
		edge_df[i, :edge_ancNodeIndex] = anc_dec_nodeIndices[1]
		edge_df[i, :edge_decNodeIndex] = anc_dec_nodeIndices[2]
		edge_df[i, :edge_ancPNnumber] = anc_dec_PNnumbers[1]
		edge_df[i, :edge_decPNnumber] = anc_dec_PNnumbers[2]
		edge_df[i, :edge_length] = one_edge.length
	end

	return(edge_df)
end


"""
using DataFrames
using PhyloBits
#using PhyloPlots


# Make a trdf, a DataFrame holding all of the key tree information
cumulative_height_at_each_node = get_node_heights(tr)
node_age = get_node_ages(tr)
"""
function get_node_heights(tr)
	indexNum_table = get_nodeIndex_PNnumber(tr)
	uppass_nodeIndexes = get_LR_uppass_nodeIndexes(tr)
	edge_df = get_NodeIndex_df_by_tree_edges(tr, indexNum_table=indexNum_table)
	cumulative_height_at_each_node = collect(repeat([0.0], length(tr.node)))
	
	#print("\nuppass_nodeIndexes:\n")
	#print(uppass_nodeIndexes)
	#print("\n")
	#print("\nedge_df:\n")
	#print(edge_df)
	
	# Iterate up through the nodes from the root
	for i in 1:length(uppass_nodeIndexes)
		#print("\n")
		#print("\n")
		#print(i)
		
		nodeIndex = uppass_nodeIndexes[i]
		if (nodeIndex == tr.root)
			cumulative_height_at_each_node[nodeIndex] = 0.0
		else
			#print("\nedge_df[:,:edge_decNodeIndex]:\n")
			#print(edge_df[:,:edge_decNodeIndex])
			#print("\nnodeIndex:\n")
			#print(nodeIndex)
			ancTF = edge_df[:,:edge_decNodeIndex] .== nodeIndex
			#print("\nancTF:\n")
			#print(ancTF)
			anc_nodeIndex = edge_df[:,:edge_ancNodeIndex][ancTF][1]
			# Get previous age
			previous_height_above_root = cumulative_height_at_each_node[anc_nodeIndex]
			length_of_this_branch = edge_df[:,:edge_length][ancTF][1]
			current_height = length_of_this_branch + previous_height_above_root
			cumulative_height_at_each_node[nodeIndex] = current_height
		end
	end
	return(cumulative_height_at_each_node)
end # END function get_node_heights(tr)




"""
using DataFrames
using PhyloBits
#using PhyloPlots

#######################################################
# Standard bifurcating tree
#######################################################
great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root

# Get trdf WITHOUT "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, false)

# Not all columns print now that trdf is big, let's print just the left and right
# Get trdf WITH "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, true);
headLR(trdf)

downpass_edgematrix = get_LR_downpass_edgematrix(tr)
taxa = get_taxa_descending_from_each_node(tr, trdf, downpass_edgematrix=get_LR_downpass_edgematrix(tr))

# Make a trdf, a DataFrame holding all of the key tree information
cumulative_height_at_each_node = get_node_heights(tr)
node_age = get_node_ages(tr)

#######################################################
# Tree with a 2-degree node inside a branch
#######################################################

great_ape_newick_string = "((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root

# Get trdf WITHOUT "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, false)

# Not all columns print now that trdf is big, let's print just the left and right
# Get trdf WITH "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, true);
headLR(trdf)

downpass_edgematrix = get_LR_downpass_edgematrix(tr)
taxa = get_taxa_descending_from_each_node(tr, trdf, downpass_edgematrix=get_LR_downpass_edgematrix(tr))

# Make a trdf, a DataFrame holding all of the key tree information
cumulative_height_at_each_node = get_node_heights(tr)
node_age = get_node_ages(tr)

"""
function get_node_ages(tr)
	cumulative_height_at_each_node = get_node_heights(tr)
	tree_height = maximum(cumulative_height_at_each_node)
	node_age = tree_height .- cumulative_height_at_each_node
	return node_age
end

"""
# Returns the age of the root, i.e. the height of the highest tip
# above the root.

# Ultrametric tree
great_ape_newick_string = "((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
tr

get_root_age(tr)


# Non-ultrametric tree
great_ape_newick_string = "((human:0.75,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
tr

get_root_age(tr)

"""
function get_root_age(tr)
	cumulative_height_at_each_node = get_node_heights(tr)
	root_age = maximum(cumulative_height_at_each_node)
	return root_age
end


"""
# Returns the age of the root, i.e. the height of the highest tip
# above the root.

# Ultrametric tree
great_ape_newick_string = "((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
tr

get_root_age(tr)


# Non-ultrametric tree
great_ape_newick_string = "((human:0.75,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
tr

get_root_age(tr)

"""
function get_tree_height(tr)
	root_age = get_root_age(tr)
	return root_age
end


















# Get the branching times, of the non-tip nodes
# (used in the likelihood calculation of a birthdeath process)
# Modeled on ape's branching.times()
function branching_times(tr)
"""
tr = readTopology("((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);")
split_times = branching_times(tr)

tr = readTopology("((chimp:1,human:1):1,gorilla:2);")
split_times = branching_times(tr)
"""
	trdf2 = prt(tr);
	# Sort the trtable by Rnodenums
	trdf3 = sort(trdf2, :Rnodenums)
	
	TF = trdf3[!,:nodeType] .!= "tip"
	split_times = trdf3[!,:node_age][TF]
	prepend!(split_times, NaN) # The NaN is added so that the number of nodes equals the number of tips(?)
	return(split_times)
end

#######################################################
# Likelihood equation in the birthdeath function
# (derived by pulling apart the birthdeath() function from ape)
# This version stores all of the pieces of the calculation, for comparison
#######################################################
function bd_liks(tr, birthRate=1.0, deathRate=0.0)
	runtxt="""
	tr = readTopology("((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);")
	birthRate=0.3288164   # ML birthRate for Psychotria tree
	deathRate=0.0					# ML deathRate for Psychotria tree
	bd = bd_liks(tr, birthRate, deathRate)
	Rnames(bd)
	bd.lnL
	bd.deviance
	bd.lnl_topology
	bd.lnl_numBirths
	bd.lnl_Births_above_root
	bd.lnl_numtips_wOneMinusDeathRate
	bd.lnl_branching_times

	tr = readTopology("((chimp:1,human:1):1,gorilla:2);")
	birthRate=1.0
	deathRate=0.0
	bd = bd_liks(tr, birthRate, deathRate)
	Rnames(bd)
	bd.lnL
	bd.deviance
	bd.lnl_topology
	bd.lnl_numBirths
	bd.lnl_Births_above_root
	bd.lnl_numtips_wOneMinusDeathRate
	bd.lnl_branching_times
	
	# Error cases when deathRate > birthRate:
	bd = bd_liks(tr, birthRate, 0.999)
	bd = bd_liks(tr, birthRate, 1.0)
	bd = bd_liks(tr, birthRate, 1.1)
	
	"""
	numTips = tr.numTaxa
	nb_node = tr.numNodes - tr.numTaxa
	nb_node_minusRoot = tr.numNodes - tr.numTaxa - 1
	split_times = branching_times(tr)

	a = deathRate / birthRate	# relative death rate
	r = birthRate - deathRate	# net diversification rate

	# Error trap -- this birthdeath function is used in Maximum Likelihood, 
	# thus it ignores the possibility of deathRate > birthRate
	if (r <= 0 || a >= 1)
		lnL = -1e+100
		dev =  -2 * lnL  # deviance, if desired, is -2 * the log-likelihood
		lnl_topology = logfactorial(nb_node)
		lnl_numBirths = NaN # nb_node_minusRoot * log(r)
		lnl_Births_above_root = r * sum(split_times[3:numTips])
		lnl_numtips_wOneMinusDeathRate = NaN # numTips * log(1 - a)
		lnl_branching_times = NaN #-2 * sum(log.(exp.(r .* split_times[2:numTips]) .- a))
	
		# results tuple
		bd = (lnL=lnL, deviance=dev, lnl_topology=lnl_topology, lnl_numBirths=lnl_numBirths, lnl_Births_above_root=lnl_Births_above_root, lnl_numtips_wOneMinusDeathRate=lnl_numtips_wOneMinusDeathRate, lnl_branching_times=lnl_branching_times)

		return(bd)
	end

	
	lnl_topology = logfactorial(nb_node)
	lnl_numBirths = nb_node_minusRoot * log(r)
	lnl_Births_above_root = r * sum(split_times[3:numTips])
	
	lnl_numtips_wOneMinusDeathRate = numTips * log(1 - a)
	# Interpretation: more tips are less likely, if relativeDeathRate is >0
	# If relativeDeathRate = 1, a=0, and lnl=-Inf... 
	#    CLEARLY WRONG EXCEPT IN A MAXIMUM LIKELIHOOD CONTEXT!!!
	# If relativeDeathRate = 0, a=0, and lnl=0, i.e. any number of tips is equiprobable
	
	lnl_branching_times = -2 * sum(log.(exp.(r .* split_times[2:numTips]) .- a))
	# For each observed branching,
	# netDiversificationRate * timeOfEvent <- take exponent of that ; this means recorded events are less likely in the past
	# <- subtract "a", a constant (relativeDeathRate)
	#
	# This would be a straight likelihood as:
	# 1/
	# (exp(r*branching_time)-a)^2
	#
	# Sum the logs of these
	#
	# If deathRate = 0
	# lnl_branching_times = -2 * sum(log(exp(birthRate*sptimes[2:N]) - 0))
	# lnl_branching_times = -2 * sum(log( exp(birthRate*sptimes[2:N]) )
	# lnl_branching_times = -2 * sum( birthRate*sptimes[2:N] )
	#
	# Note: sum(X) = 9 = total branchlength of tr
	# In BD:
	# -2*sum(sptimes[2:N]) = -12
	# sum(sptimes[3:N]) = 3
	# So, lnl_branching_times + lnl_Births_above_root = yule's -lambda * X
	lnL = lnl_topology + lnl_numBirths + lnl_Births_above_root + lnl_numtips_wOneMinusDeathRate + lnl_branching_times
	dev =  -2 * lnL  # deviance, if desired, is -2 * the log-likelihood
	
	# results tuple
	bd = (lnL=lnL, deviance=dev, lnl_topology=lnl_topology, lnl_numBirths=lnl_numBirths, lnl_Births_above_root=lnl_Births_above_root, lnl_numtips_wOneMinusDeathRate=lnl_numtips_wOneMinusDeathRate, lnl_branching_times=lnl_branching_times)
	
	extract="""
	bd = (lnL=lnL, deviance=dev, lnl_topology=lnl_topology, lnl_numBirths=lnl_numBirths, lnl_Births_above_root=lnl_Births_above_root, lnl_numtips_wOneMinusDeathRate=lnl_numtips_wOneMinusDeathRate, lnl_branching_times=lnl_branching_times)
	
	lnL = bd.lnL
	dev = bd.deviance
	lnl_topology = bd.lnl_topology
	lnl_numBirths = bd.lnl_numBirths
	lnl_Births_above_root = bd.lnl_Births_above_root
	lnl_numtips_wOneMinusDeathRate = bd.lnl_numtips_wOneMinusDeathRate
	lnl_branching_times = bd.lnl_branching_times
	"""
	
	return(bd)
end # END function bd_liks(tr, birthRate=1.0, deathRate=0.0)



# Function that does bd_liks, but on a trdf (tree table from prt(tr), tree data.frame)
function bd_liks_trdf(trdf, birthRate=1.0, deathRate=0.0)
	runtxt="""
	tr = readTopology("((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);")
	trdf = prt(tr)
	birthRate=0.3288164   # ML birthRate for Psychotria tree
	deathRate=0.0					# ML deathRate for Psychotria tree
	bd = bd_liks_trdf(trdf, birthRate, deathRate)
	Rnames(bd)
	bd.lnL
	bd.deviance
	bd.lnl_topology
	bd.lnl_numBirths
	bd.lnl_Births_above_root
	bd.lnl_numtips_wOneMinusDeathRate
	bd.lnl_branching_times

	tr = readTopology("((chimp:1,human:1):1,gorilla:2);")
	trdf = prt(tr)
	birthRate=1.0
	deathRate=0.0
	bd = bd_liks_trdf(trdf, birthRate, deathRate)
	Rnames(bd)
	bd.lnL
	bd.deviance
	bd.lnl_topology
	bd.lnl_numBirths
	bd.lnl_Births_above_root
	bd.lnl_numtips_wOneMinusDeathRate
	bd.lnl_branching_times
	
	# Error cases when deathRate > birthRate:
	bd = bd_liks_trdf(trdf, birthRate, 0.999)
	bd = bd_liks_trdf(trdf, birthRate, 1.0)
	bd = bd_liks_trdf(trdf, birthRate, 1.1)
	
	"""
	# Get basic tree info
	numTips = sum(trdf.nodeType .== "tip")
	nb_node = sum(trdf.nodeType .== "intern") + sum(trdf.nodeType .== "root")
	nb_node_minusRoot = nb_node - 1

	trdf2 = sort(trdf, :Rnodenums)
	TF = trdf2[!,:nodeType] .!= "tip"
	split_times = trdf2[!,:node_age][TF]
	prepend!(split_times, NaN) # The NaN is added so that the number of nodes equals the number of tips(?)
	
	a = deathRate / birthRate	# relative death rate
	r = birthRate - deathRate	# net diversification rate

	# Error trap -- this birthdeath function is used in Maximum Likelihood, 
	# thus it ignores the possibility of deathRate > birthRate
	if (r <= 0 || a >= 1)
		lnL = -1e+100
		dev =  -2 * lnL  # deviance, if desired, is -2 * the log-likelihood
		lnl_topology = logfactorial(nb_node)
		lnl_numBirths = NaN # nb_node_minusRoot * log(r)
		lnl_Births_above_root = r * sum(split_times[3:numTips])
		lnl_numtips_wOneMinusDeathRate = NaN # numTips * log(1 - a)
		lnl_branching_times = NaN #-2 * sum(log.(exp.(r .* split_times[2:numTips]) .- a))
	
		# results tuple
		bd = (lnL=lnL, deviance=dev, lnl_topology=lnl_topology, lnl_numBirths=lnl_numBirths, lnl_Births_above_root=lnl_Births_above_root, lnl_numtips_wOneMinusDeathRate=lnl_numtips_wOneMinusDeathRate, lnl_branching_times=lnl_branching_times)

		return(bd)
	end

	
	lnl_topology = logfactorial(nb_node)
	lnl_numBirths = nb_node_minusRoot * log(r)
	lnl_Births_above_root = r * sum(split_times[3:numTips])
	
	lnl_numtips_wOneMinusDeathRate = numTips * log(1 - a)
	# Interpretation: more tips are less likely, if relativeDeathRate is >0
	# If relativeDeathRate = 1, a=0, and lnl=-Inf... 
	#    CLEARLY WRONG EXCEPT IN A MAXIMUM LIKELIHOOD CONTEXT!!!
	# If relativeDeathRate = 0, a=0, and lnl=0, i.e. any number of tips is equiprobable
	
	lnl_branching_times = -2 * sum(log.(exp.(r .* split_times[2:numTips]) .- a))
	# For each observed branching,
	# netDiversificationRate * timeOfEvent <- take exponent of that ; this means recorded events are less likely in the past
	# <- subtract "a", a constant (relativeDeathRate)
	#
	# This would be a straight likelihood as:
	# 1/
	# (exp(r*branching_time)-a)^2
	#
	# Sum the logs of these
	#
	# If deathRate = 0
	# lnl_branching_times = -2 * sum(log(exp(birthRate*sptimes[2:N]) - 0))
	# lnl_branching_times = -2 * sum(log( exp(birthRate*sptimes[2:N]) )
	# lnl_branching_times = -2 * sum( birthRate*sptimes[2:N] )
	#
	# Note: sum(X) = 9 = total branchlength of tr
	# In BD:
	# -2*sum(sptimes[2:N]) = -12
	# sum(sptimes[3:N]) = 3
	# So, lnl_branching_times + lnl_Births_above_root = yule's -lambda * X
	lnL = lnl_topology + lnl_numBirths + lnl_Births_above_root + lnl_numtips_wOneMinusDeathRate + lnl_branching_times
	dev =  -2 * lnL  # deviance, if desired, is -2 * the log-likelihood
	
	# results tuple
	bd = (lnL=lnL, deviance=dev, lnl_topology=lnl_topology, lnl_numBirths=lnl_numBirths, lnl_Births_above_root=lnl_Births_above_root, lnl_numtips_wOneMinusDeathRate=lnl_numtips_wOneMinusDeathRate, lnl_branching_times=lnl_branching_times)
	
	extract="""
	bd = (lnL=lnL, deviance=dev, lnl_topology=lnl_topology, lnl_numBirths=lnl_numBirths, lnl_Births_above_root=lnl_Births_above_root, lnl_numtips_wOneMinusDeathRate=lnl_numtips_wOneMinusDeathRate, lnl_branching_times=lnl_branching_times)
	
	lnL = bd.lnL
	dev = bd.deviance
	lnl_topology = bd.lnl_topology
	lnl_numBirths = bd.lnl_numBirths
	lnl_Births_above_root = bd.lnl_Births_above_root
	lnl_numtips_wOneMinusDeathRate = bd.lnl_numtips_wOneMinusDeathRate
	lnl_branching_times = bd.lnl_branching_times
	"""
	
	return(bd)
end # END function bd_liks(tr, birthRate=1.0, deathRate=0.0)



"""
using PhyloBits
tr = readTopology("((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);")
total_branch_length = get_treelength(tr)
"""
function get_treelength(tr)
	treelength = 0.0
	@inbounds for i in 1:length(tr.edge)
		treelength += tr.edge[i].length
	end
	return treelength
end

"""
using PhyloBits
tr = readTopology("((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);")
total_branch_length = get_treelength(tr)
ML_yule_birthRate(tr)
ML_yule_birthRate_wRoot(tr)

# A simple great ape tree, shortest branch lengths of 1.0
tr = readTopology("(((chimp:1,human:1):1,gorilla:2):1,orang:3);")
total_branch_length = get_treelength(tr)
ML_yule_birthRate(tr)
ML_yule_birthRate_wRoot(tr)
"""
function ML_yule_birthRate(tr::HybridNetwork; hooks_below=1.0e-6)
	# You *HAVE* to remove (1) direct ancestor nodes (singletons) and
	#                      (2) hook nodes & tips; neither are direct-ancestor nodes
	trdf = prt(tr, rootnodenum=tr.root, get_taxa_by_node=false; hooks_below=hooks_below);
	
	num_speciation_nodes = get_num_speciation_nodes(trdf; hooks_below=hooks_below)
	num_tips_from_speciation = get_num_tips_from_speciation(trdf; hooks_below=hooks_below)
	
	#num_sp_events = tr.numNodes-tr.numTaxa
	num_sp_events = num_speciation_nodes - num_tips_from_speciation
	total_branch_length = get_treelength(tr)
	lambda = (num_sp_events-1.0) / total_branch_length # -1.0, following Yule birthdeath function
	return lambda
end


"""
using PhyloBits

#################################################################
# Psychotria example tree, from Ree & Smith (2008)
#################################################################
tr = readTopology("((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);")
total_branch_length = get_treelength(tr)
ML_yule_birthRate(tr)
ML_yule_birthRate_wRoot(tr)

trdf = prt(tr)
birthRate = ML_yule_birthRate(tr)
deathRate = 0.0
bd = bd_liks_trdf(trdf, birthRate, deathRate)
bd_lnL_noTopo = bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times

birthRate = ML_yule_birthRate_wRoot(tr)
bd_wRoot = bd_liks_trdf(trdf, birthRate, deathRate)
bd_lnL_noTopo_wRoot = bd_wRoot.lnl_numBirths + bd_wRoot.lnl_Births_above_root + bd_wRoot.lnl_numtips_wOneMinusDeathRate + bd_wRoot.lnl_branching_times

bd.lnL - bd_wRoot.lnL
bd_lnL_noTopo - bd_lnL_noTopo_wRoot

flat2(bd) ./ flat2(bd_wRoot)
difs = flat2(bd) .- flat2(bd_wRoot)
difs[1]
sum(difs[2:7])
1 + difs[4]
difs[2] / 2

#################################################################
# A simple great ape tree, shortest branch lengths of 1.0
#################################################################
tr = readTopology("(((chimp:1,human:1):1,gorilla:2):1,orang:3);")
total_branch_length = get_treelength(tr)
ML_yule_birthRate(tr)
ML_yule_birthRate_wRoot(tr)

trdf = prt(tr)
birthRate = ML_yule_birthRate(tr)
deathRate = 0.0
bd = bd_liks_trdf(trdf, birthRate, deathRate)
bd_lnL_noTopo = bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times

birthRate = ML_yule_birthRate_wRoot(tr)
bd_wRoot = bd_liks_trdf(trdf, birthRate, deathRate)
bd_lnL_noTopo_wRoot = bd_wRoot.lnl_numBirths + bd_wRoot.lnl_Births_above_root + bd_wRoot.lnl_numtips_wOneMinusDeathRate + bd_wRoot.lnl_branching_times

bd.lnL - bd_wRoot.lnL
bd_lnL_noTopo - bd_lnL_noTopo_wRoot

difs = flat2(bd) .- flat2(bd_wRoot)
difs[1]
sum(difs[2:7])
1 + difs[4]
difs[2] / 2

"""
function ML_yule_birthRate_wRoot(tr::HybridNetwork; hooks_below=1.0e-6)
	# You *HAVE* to remove (1) direct ancestor nodes (singletons) and
	#                      (2) hook nodes & tips; neither are direct-ancestor nodes
	trdf = prt(tr, rootnodenum=tr.root, get_taxa_by_node=false; hooks_below=hooks_below);
	
	num_speciation_nodes = get_num_speciation_nodes(trdf; hooks_below=hooks_below)
	num_tips_from_speciation = get_num_tips_from_speciation(trdf; hooks_below=hooks_below)
	
	#num_sp_events = tr.numNodes-tr.numTaxa
	num_sp_events = num_speciation_nodes - num_tips_from_speciation
	total_branch_length = get_treelength(tr)
	lambda = (num_sp_events-0.0) / total_branch_length # -0.0, meaning the root node 
																										 # speciation event is counted.
	return lambda
end

# Get the number of nodes resulting from speciation events (tips+nodes)
# i.e. traditional nb.node in Yule process calcs etc.
# eliminates hooknodes, hooktips, and direct ancestor nodes
function get_num_speciation_nodes(trdf::DataFrame; hooks_below=1.0e-6)
	numNodes = nrow(trdf)
	hooktips_TF = (trdf.hook .+ trdf.nodeType.=="tip") .== 2
	hooknodes_TF = (trdf.hook .+ trdf.nodeType.=="intern" .+ trdf.nodeType.=="root") .== 2
	direct_TF = trdf.nodeType .== "direct"
	num_speciation_nodes = numNodes - sum(hooknodes_TF) - sum(hooktips_TF) - sum(direct_TF)
	return(num_speciation_nodes)
end

# Get the number of nodes resulting from speciation events (tips+nodes)
# i.e. traditional nb.node in Yule process calcs etc.
# eliminates hooknodes, hooktips, and direct ancestor nodes
function get_num_speciation_nodes(tr::HybridNetwork; hooks_below=1.0e-6)
	trdf = prt(tr, rootnodenum=tr.root, get_taxa_by_node=false; hooks_below=hooks_below);
	#numNodes = nrow(trdf)
	#numNodes = tr.numNodes
	num_speciation_nodes = get_num_speciation_nodes(trdf; hooks_below=hooks_below)
	return(num_speciation_nodes)
end

# Get the number of tips resulting from speciation events
# i.e. traditional nb.node in Yule process calcs etc.
# eliminates hooknodes, hooktips, and direct ancestor nodes
function get_num_tips_from_speciation(trdf::DataFrame; hooks_below=1.0e-6)
	numTaxa = sum(trdf.nodeType .== "tip")
	hooktips_TF = (trdf.hook .+ trdf.nodeType.=="tip") .== 2
	num_tips_from_speciation = numTaxa - sum(hooktips_TF)
	return(num_tips_from_speciation)
end

function get_num_tips_from_speciation(tr::HybridNetwork; hooks_below=1.0e-6)
	trdf = prt(tr, rootnodenum=tr.root, get_taxa_by_node=false; hooks_below=hooks_below);
	num_tips_from_speciation = get_num_tips_from_speciation(trdf; hooks_below=hooks_below)
	return(num_tips_from_speciation)
end

end # end of module