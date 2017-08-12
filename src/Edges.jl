#__precompile__() 
#module Edges

#include("MeshPrimitives.jl")

#using .MeshPrimitives

include("edges/curvededge.jl")
include("edges/lineedge.jl")

export CurvedEdge, compare, line_divide,
       LineEdge, poistion, mag

#end
