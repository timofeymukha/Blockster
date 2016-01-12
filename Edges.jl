module Edges

using ImmutableArrays

include("edges/curvededge.jl")
include("edges/lineedge.jl")

export CurvedEdge, compare, line_divide,
       LineEdge, poistion, mag

end
