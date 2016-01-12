module blockMesh

using JSON
using ImmutableArrays
using Edges
using Blocks

# Parse the dictionary
dict = JSON.parsefile("dictJSON")

# Get the number of blocks
nBlocks = size(dict["blocks"], 1)
nVertices = size(dict["vertices"], 1)

vertices = Vector{Vector3{Float64}}(nVertices)

for i in 1:nVertices
    vertices[i] = dict["vertices"][i]
end


function read_boundary!(meshDict, 
                        patchNames::Array{ASCIIString, 1},
                        patchBlockFaces::Vector{Vector{Vector{Int64}}})

    nPatches = size(meshDict["boundary"], 1)

    resize!(patchNames, nPatches)
    resize!(patchBlockFaces, nPatches)

    # Grab the patch names
    for i in 1:nPatches
        patchNames[i] = meshDict["boundary"][i]["name"]
        patchBlockFaces[i] = meshDict["boundary"][i]["faces"] + 1
    end

end

function check_patch_vertex_labels(patchNames,
                                   patchBlockFaces,
                                   vertices)
    for patchI in 1:size(patchNames, 1)
        for faceI in size(patchBlockFaces[patchI], 1)
            if !isempty(find(patchBlockFaces[patchI][faceI] .< 1))
                error("""check_patch_vertex_labels() : face vertex label < 0
                       in patch """, patchNames[patchI])
            elseif !isempty(find(patchBlockFaces[patchI][faceI] .>
                                 size(vertices, 1)))
                error("""check_patch_vertex_labels() : face vertex label 
                out of bounds in patch """, patchNames[patchI])
            end
        end
    end
end

patchNames = ASCIIString[]
patchBlockFaces = Vector{Vector{Vector{Int64}}}(0)

read_boundary!(dict, patchNames, patchBlockFaces)
check_patch_vertex_labels(patchNames, patchBlockFaces, vertices)


println("Creating blocks...")

function create_topology()
    
end

blocks = Vector{Block}(nBlocks)

for blockI in 1:nBlocks
    println("    Block number ", blockI)
    blocks[blockI] = Block()
    blocks[blockI].nCells = dict["blocks"][blockI]["number of cells"]
    blocks[blockI].vertexLabels = dict["blocks"][blockI]["vertices"] + 1

    for j in 1:8
        blocks[blockI].vertices[j] = vertices[blocks[blockI].vertexLabels[j]]
    end
    
    # Create the edge points
    println("        Creating edge-points")
    make_block_edges!(blocks[blockI])

    # Create the points
    println("        Creating points")
    create_points!(blocks[blockI])

    # Create the cells
    println("        Creating cells")
    create_cells!(blocks[blockI])
end

end
