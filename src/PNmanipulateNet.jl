module PNmanipulateNet
__precompile__(true)

export isTree, directEdges!, traverseDirectEdges!, getParents, getMajorParentEdge, getChildren


# Source:
# https://github.com/crsl4/PhyloNetworks.jl/blob/08eab8b17a171db74e08b92708913582cd30b526/src/snaq_optimization.jl
function isTree(net::HybridNetwork)
    net.numHybrids == length(net.hybrid) || error("numHybrids does not match to length of net.hybrid")
    net.numHybrids != 0 || return true
    return false
end






# Source: https://github.com/crsl4/PhyloNetworks.jl/blob/904c248720325b2e603918d9441aeacf0d2729af/src/manipulateNet.jl

# Claudia SL & Paul Bastide: November 2015, Cecile: Feb 2016

#################################################
# Direct Edges
#################################################

"""
    directEdges!(net::HybridNetwork; checkMajor=true::Bool)
Updates the edges' attribute `isChild1`, according to the root placement.
Also updates edges' attribute `containRoot`, for other possible root placements
compatible with the direction of existing hybrid edges.
Relies on hybrid nodes having exactly 1 major hybrid parent edge,
but checks for that if checkMajor=true.
Warnings:
1. Assumes that isChild1 is correct on hybrid edges
(to avoid changing the identity of which nodes are hybrids and which are not).
2. Does not check for cycles (to maintain a network's DAG status)
Returns the network. Throws a 'RootMismatch' Exception if the root was found to
conflict with the direction of any hybrid edge.
"""
function directEdges!(net::HybridNetwork; checkMajor=true::Bool)
    if checkMajor # check each node has 2+ hybrid parent edges (if any), and exactly one major.
        for n in net.node
            nparents = 0 # 0 or 2 normally, but could be >2 if polytomy.
            nmajor = 0   # there should be exactly 1 major parent if nparents>0
            for e in n.edge
                if e.hybrid && n == getChild(e)
                    nparents += 1
                    if (e.isMajor) nmajor +=1; end
                end
            end
            (nparents!=1) || error("node $(n.number) has exactly 1 hybrid parent edge")
            (nparents==0 || nmajor == 1) ||
              error("hybrid node $(n.number) has 0 or 2+ major hybrid parents")
            (nparents!=2 || n.hybrid) ||
              @warn "node $(n.number) has 2 parents but its hybrid attribute is false.
It is not used in directEdges!, but might cause an error elsewhere."
            # to fix this: change n.hybrid, net.hybrid, net.numHybrids etc.
            # none of those attributes are used here.
        end
    end
    net.cleaned = false # attributed used by snaq! Will change isChild1 and containRoot
    for e in net.node[net.root].edge
        traverseDirectEdges!(net.node[net.root],e,true)
    end
    net.isRooted = true
    return net
end

# containroot = true until the path goes through a hybrid node, below which
# containroot is turned to false.
function traverseDirectEdges!(node::Node, edge::Edge, containroot::Bool)
    if edge.hybrid && node==getChild(edge)
        throw(RootMismatch(
"direction (isChild1) of hybrid edge $(edge.number) conflicts with the root.
isChild1 and containRoot were updated for a subset of edges in the network only."))
    end
    if node == edge.node[1]
        edge.isChild1 = false
        cn = edge.node[2] # cn = child node
    else
        edge.isChild1 = true
        cn = edge.node[1]
    end
    edge.containRoot = containroot
    if !cn.leaf && (!edge.hybrid || edge.isMajor) # continue down recursion
        if edge.hybrid containroot=false; end # changes containroot locally, intentional.
        nchildren=0
        for e in cn.edge
            if e==edge continue; end
            if (e.hybrid && cn == getChild(e)) continue; end
            traverseDirectEdges!(cn,e,containroot)
            nchildren += 1
        end
        if nchildren==0
            throw(RootMismatch("non-leaf node $(cn.number) had 0 children.
Could be a hybrid whose parents' direction conflicts with the root.
isChild1 and containRoot were updated for a subset of edges in the network only."))
        end
    end
    return nothing
end

#################################################
## Topological sorting
#################################################

"""
    getParents(node)
Get vector of all parent nodes of `n`, based on `isChild1` field (for edges).
To get the parent node of an edge: see [`getParent`](@ref).
To get individual parent edges (rather than all parent *nodes*):
see [`getMajorParentEdge`](@ref) and `getMinorParentEdge`.
"""
@inline function getParents(node::Node)
    parents = Node[]
    for e in node.edge
            if node == getChild(e)
                push!(parents, getParent(e))
            end
    end
    return parents
end

# getParent, getMajorParent, getMinorParent: defined in auxiliary.jl

"""
    getMajorParentEdge(node)
    getMinorParentEdge(node)
return the parent edge of a given node: the major / minor if hybrid.
**warning**: assume isChild1 and isMajor attributes are correct
To get all parent *nodes*: see [`getParents`](@ref).
"""
@inline function getMajorParentEdge(n::Node)
    for ee in n.edge
        if n == ee.node[(ee.isChild1 ? 1 : 2)] && ee.isMajor
            return ee
        end
    end
    error("node $(n.number) has no major parent")
end
@doc (@doc getMajorParentEdge) getMinorParentEdge
@inline function getMinorParentEdge(n::Node)
    for ee in n.edge
        if !ee.isMajor && n == ee.node[(ee.isChild1 ? 1 : 2)]
            return ee
        end
    end
    error("node $(n.number) has no minor parent")
end

"""
    getChildren(node)
return a vector with all children *nodes* of `node`.
**warning**: assume `isChild1` field (for edges) are correct
To get all parent *nodes*: see [`getParents`](@ref).
"""
function getChildren(node::Node)
    children = Node[]
    for e in node.edge
        if node === getParent(e)
            push!(children, getChild(e))
        end
    end
    return children
end


end # END module PNmanipulateNet
