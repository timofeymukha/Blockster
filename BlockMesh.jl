module blockMesh

using JSON
using MeshPrimitives
using Edges
using Blocks

# Parse the dictionary
dict = JSON.parsefile("dictJSON")

# Get the number of blocks
nBlocks = size(dict["blocks"], 1)
nVertices = size(dict["vertices"], 1)

vertices = Vector{point}(nVertices)

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



function create_topology()
    
end



"Type defining the mesh"
type Mesh
    "The json file to read the mesh properties from"
    dict

    "The blocks defining the mesh" 
    blocks::Vector{Block}

    "The points of the mesh"
    points::Vector{point}

    "The vertices of the mesh, used to define the blocks"
    vertices::Vector{point}
end

println("Creating blocks...")

blocks = Vector{Block}(nBlocks)

function point_cell_addressing(cells::Array{Int64, 2}, nPoints::Int64)

    # Number of cells each point is included in
    pointCellAddressing = Vector{Vector{Int64}}(nPoints)

    for i in 1:nPoints
        pointCellAddressing[i] = Vector{Int64}(0)
    end

    # For each cells
    for cellI in 1:size(cells, 1)

        # For each point in the cell
        for pointI in 1:size(cells[cellI], 2)

            # Add current cell to the addressing list
            push!(pointCellAddressing[pointI], cellI)
        end

    end

    return pointCellAddressing
end

function patch_face_cells(faces::Vector{face}, patchId::Int64,
                          pointCellAddressing::Vector{Vector{Int64}})

    faceCells = Vector{Int64}(size(faces,1))

    # For each face determine the cell
    for faceI in 1:size(faces, 1)
        currFace = faces[faceI]

        # For each point of the face
        for facePointI in 1:size(currFace, 1)
            
            facePointCells = pointCellAddressing[currFace[facePointI]]

            # For each cell that this point is part of
            for cellI in 1:size(facePointCells, 1)

                #For each face of the cells check if it is the same as faceI

            end

        end



    end

    return facecells
end

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
    println("        Done")
end

pointCellsAddressing = point_cell_addressing(blocks[1].cells, size(blocks[1].points, 1))

end
