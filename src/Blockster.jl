__precompile__() 
module Blockster

using StaticArrays: SVector

include("MeshPrimitives.jl")
include("Read.jl")
include("Edges.jl")
include("Blocks.jl")
include("MeshCreate.jl")
include("Write.jl")

end
